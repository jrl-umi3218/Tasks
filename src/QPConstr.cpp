// This file is part of Tasks.
//
// Tasks is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Tasks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Tasks.  If not, see <http://www.gnu.org/licenses/>.

// associated header
#include "QPConstr.h"

// includes
// std
#include <set>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// SCD
#include <SCD/CD/CD_Pair.h>
#include <SCD/S_Object/S_Object.h>

namespace tasks
{

namespace qp
{

MotionConstr::ContactData::ContactData(const rbd::MultiBody& mb, int b,
  std::vector<Eigen::Vector3d> p, int nrGen):
  jac(mb, b),
  body(jac.jointsPath().back()),
  points(p),
  generators(3, nrGen),
  jacTrans(6, jac.dof()),
  generatorsComp(3, nrGen)
{ }


MotionConstr::MotionConstr(const rbd::MultiBody& mb):
	fd_(mb),
	cont_(),
	fullJac_(6, mb.nrDof()),
	AEq_(),
	BEq_(),
	XL_(),
	XU_(),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0)
{
}


void MotionConstr::updateNrVars(const rbd::MultiBody& mb,
	const SolverData& data)
{
	const auto& uniCont = data.unilateralContacts();
	const auto& biCont = data.bilateralContacts();
	cont_.resize(data.nrContacts());

	nrDof_ = data.alphaD();
	nrFor_ = data.lambda();
	nrTor_ = data.torque();

	std::size_t iCont = 0;
	for(const UnilateralContact& c: uniCont)
	{
		cont_[iCont] = ContactData(mb, c.bodyId, c.points,
			static_cast<int>(c.cone.generators.size()));
		for(std::size_t i = 0; i < c.cone.generators.size(); ++i)
		{
			cont_[iCont].generators.col(i) = c.cone.generators[i];
		}

		++iCont;
	}

	for(const BilateralContact& c: biCont)
	{
		cont_[iCont] = ContactData(mb, c.bodyId, c.points, 3);
		cont_[iCont].generators = c.frame;

		++iCont;
	}

	AEq_.resize(nrDof_, nrDof_ + nrFor_ + nrTor_);
	BEq_.resize(nrDof_);

	AEq_.setZero();
	BEq_.setZero();

	XL_.resize(data.nrUnilateralLambda());
	XU_.resize(data.nrUnilateralLambda());

	XL_.fill(0.);
	XU_.fill(std::numeric_limits<double>::infinity());
}


void MotionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	fd_.computeH(mb, mbc);
	fd_.computeC(mb, mbc);

	// H*alphaD - tau - tau_c = -C

	// AEq
	//         nrDof      nrFor            nrTor
	// nrDof [   H      -Sum J_i^t*ni     [0 ... -1]

	AEq_.block(0, 0, nrDof_, nrDof_) = fd_.H();

	int contPos = nrDof_;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		const MatrixXd& jac = cont_[i].jac.jacobian(mb, mbc);
		cont_[i].generatorsComp =
			mbc.bodyPosW[cont_[i].body].rotation().transpose()*cont_[i].generators;

		// for each contact point we compute all the torques
		// due to each generator of the friction cone
		for(std::size_t j = 0; j < cont_[i].points.size(); ++j)
		{
			cont_[i].jac.translateJacobian(jac, mbc,
				cont_[i].points[j], cont_[i].jacTrans);
			cont_[i].jac.fullJacobian(mb, cont_[i].jacTrans, fullJac_);

			AEq_.block(0, contPos, nrDof_, cont_[i].generatorsComp.cols()) =
				-fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*
					cont_[i].generatorsComp;

			contPos += cont_[i].generatorsComp.cols();
		}
	}

	AEq_.block(mb.joint(0).dof(), contPos, nrTor_, nrTor_) =
		-MatrixXd::Identity(nrTor_, nrTor_);

	// BEq = -C
	BEq_ = -fd_.C();
}


int MotionConstr::nrEqLine()
{
	return nrDof_;
}


const Eigen::MatrixXd& MotionConstr::AEq() const
{
	return AEq_;
}


const Eigen::VectorXd& MotionConstr::BEq() const
{
	return BEq_;
}


int MotionConstr::beginVar()
{
	return nrDof_;
}


const Eigen::VectorXd& MotionConstr::Lower() const
{
	return XL_;
}


const Eigen::VectorXd& MotionConstr::Upper() const
{
	return XU_;
}


/**
	*															ContactAccConstr
	*/


std::set<int> bodyIdInContact(const rbd::MultiBody& mb,
	const SolverData& data)
{
	std::set<int> ret;
	auto isValid = [&mb](int bodyId)
	{
		// if fixed base and support body we don't add the contact
		return !(bodyId == mb.body(0).id() &&
			mb.joint(0).type() == rbd::Joint::Fixed);
	};

	for(const UnilateralContact& c: data.unilateralContacts())
	{
		if(isValid(c.bodyId))
		{
			ret.insert(c.bodyId);
		}
	};

	for(const BilateralContact& c: data.bilateralContacts())
	{
		if(isValid(c.bodyId))
		{
			ret.insert(c.bodyId);
		}
	}

	return std::move(ret);
}


ContactAccConstr::ContactAccConstr(const rbd::MultiBody& mb):
	cont_(),
	fullJac_(6, mb.nrDof()),
	alphaVec_(mb.nrDof()),
	AEq_(),
	BEq_(),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0)
{}


void ContactAccConstr::updateNrVars(const rbd::MultiBody& mb,
	const SolverData& data)
{
	cont_.clear();
	nrDof_ = data.alphaD();
	nrFor_ = data.lambda();
	nrTor_ = data.torque();

	std::set<int> bodyIdSet = bodyIdInContact(mb, data);
	for(int bodyId: bodyIdSet)
	{
		cont_.emplace_back(rbd::Jacobian(mb, bodyId));
	}

	AEq_.resize(cont_.size()*6 , nrDof_ + nrFor_ + nrTor_);
	BEq_.resize(cont_.size()*6);

	AEq_.setZero();
	BEq_.setZero();
}


void ContactAccConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	rbd::paramToVector(mbc.alpha, alphaVec_);

	// J_i*alphaD + JD_i*alpha = 0

	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		// AEq = J_i
		const MatrixXd& jac = cont_[i].jac.jacobian(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jac, fullJac_);
		AEq_.block(i*6, 0, 6, mb.nrDof()) = fullJac_;

		// BEq = -JD_i*alpha
		const MatrixXd& jacDot = cont_[i].jac.jacobianDot(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jacDot, fullJac_);
		BEq_.segment(i*6, 6) = -fullJac_*alphaVec_;
	}
}


int ContactAccConstr::nrEqLine()
{
	return AEq_.rows();
}


const Eigen::MatrixXd& ContactAccConstr::AEq() const
{
	return AEq_;
}


const Eigen::VectorXd& ContactAccConstr::BEq() const
{
	return BEq_;
}


/**
	*															ContactSpeedConstr
	*/


ContactSpeedConstr::ContactSpeedConstr(const rbd::MultiBody& mb, double timeStep):
	cont_(),
	fullJac_(6, mb.nrDof()),
	alphaVec_(mb.nrDof()),
	AEq_(),
	BEq_(),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0),
	timeStep_(timeStep)
{}


void ContactSpeedConstr::updateNrVars(const rbd::MultiBody& mb,
	const SolverData& data)
{
	cont_.clear();
	nrDof_ = data.alphaD();
	nrFor_ = data.lambda();
	nrTor_ = data.torque();

	std::set<int> bodyIdSet = bodyIdInContact(mb, data);
	for(int bodyId: bodyIdSet)
	{
		cont_.emplace_back(rbd::Jacobian(mb, bodyId));
	}

	AEq_.resize(cont_.size()*6 , nrDof_ + nrFor_ + nrTor_);
	BEq_.resize(cont_.size()*6);

	AEq_.setZero();
	BEq_.setZero();
}


void ContactSpeedConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	rbd::paramToVector(mbc.alpha, alphaVec_);

	// J_i*alphaD + JD_i*alpha = 0

	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		// AEq = J_i
		const MatrixXd& jac = cont_[i].jac.jacobian(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jac, fullJac_);
		AEq_.block(i*6, 0, 6, mb.nrDof()) = fullJac_;

		// BEq = -JD_i*alpha
		const MatrixXd& jacDot = cont_[i].jac.jacobianDot(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jacDot, fullJac_);
		BEq_.segment(i*6, 6) = -fullJac_*alphaVec_ -
			mbc.bodyVelW[cont_[i].body].vector()/timeStep_;
	}
}


int ContactSpeedConstr::nrEqLine()
{
	return AEq_.rows();
}


const Eigen::MatrixXd& ContactSpeedConstr::AEq() const
{
	return AEq_;
}


const Eigen::VectorXd& ContactSpeedConstr::BEq() const
{
	return BEq_;
}


/**
	*															JointLimitsConstr
	*/


JointLimitsConstr::JointLimitsConstr(const rbd::MultiBody& mb,
	std::vector<std::vector<double>> lBound,
	std::vector<std::vector<double>> uBound,
	double step):
	lower_(),
	upper_(),
	qMin_(mb.nrParams()),
	qMax_(mb.nrParams()),
	qVec_(mb.nrParams()),
	alphaVec_(mb.nrDof()),
	begin_(mb.joint(0).dof()),
	step_(step)
{
	int vars = mb.nrDof() - mb.joint(0).dof();
	qMin_.resize(vars);
	qMax_.resize(vars);
	lower_.resize(vars);
	upper_.resize(vars);

	// remove the joint 0
	lBound[0] = {};
	uBound[0] = {};

	rbd::paramToVector(lBound, qMin_);
	rbd::paramToVector(uBound, qMax_);
}


void JointLimitsConstr::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& /* data */)
{
}


void JointLimitsConstr::update(const rbd::MultiBody& /* mb */, const rbd::MultiBodyConfig& mbc)
{
	double dts = step_*step_*0.5;
	int vars = lower_.rows();

	rbd::paramToVector(mbc.q, qVec_);
	rbd::paramToVector(mbc.alpha, alphaVec_);

	lower_ = qMin_ - qVec_.tail(vars) - alphaVec_.tail(vars)*step_;
	lower_ /= dts;

	upper_ = qMax_ - qVec_.tail(vars) - alphaVec_.tail(vars)*step_;
	upper_ /= dts;
}


int JointLimitsConstr::beginVar()
{
	return begin_;
}


const Eigen::VectorXd& JointLimitsConstr::Lower() const
{
	return lower_;
}


const Eigen::VectorXd& JointLimitsConstr::Upper() const
{
	return upper_;
}



/**
	*															TorqueLimitsConstr
	*/


TorqueLimitsConstr::TorqueLimitsConstr(const rbd::MultiBody& mb,
	std::vector<std::vector<double>> lBound,
	std::vector<std::vector<double>> uBound):
	lower_(),
	upper_(),
	begin_()
{
	int vars = mb.nrDof() - mb.joint(0).dof();
	lower_.resize(vars);
	upper_.resize(vars);

	// remove the joint 0
	lBound[0] = {};
	uBound[0] = {};

	rbd::paramToVector(lBound, lower_);
	rbd::paramToVector(uBound, upper_);
}


void TorqueLimitsConstr::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& data)
{
	begin_ = data.torqueBegin();
}


void TorqueLimitsConstr::update(const rbd::MultiBody& /* mb */,
	const rbd::MultiBodyConfig& /* mbc */)
{
}


int TorqueLimitsConstr::beginVar()
{
	return begin_;
}


const Eigen::VectorXd& TorqueLimitsConstr::Lower() const
{
	return lower_;
}


const Eigen::VectorXd& TorqueLimitsConstr::Upper() const
{
	return upper_;
}

/**
	*													SelfCollisionConstr
	*/


SCD::Matrix4x4 toSCD(const sva::PTransform& t)
{
	SCD::Matrix4x4 m;
	const Eigen::Matrix3d& rot = t.rotation();
	const Eigen::Vector3d& tran = t.translation();

	for(int i = 0; i < 3; ++i)
	{
		for(int j = 0; j < 3; ++j)
		{
			m(i,j) = rot(j,i);
		}
	}

	m(0,3) = tran(0);
	m(1,3) = tran(1);
	m(2,3) = tran(2);

	return m;
}


SelfCollisionConstr::CollData::CollData(const rbd::MultiBody& mb,
	int body1Id, SCD::S_Object* body1, const sva::PTransform& body1T,
	int body2Id, SCD::S_Object* body2, const sva::PTransform& body2T,
	double di, double ds, double damping):
		pair(new SCD::CD_Pair(body1, body2)),
		body1T(body1T),
		body2T(body2T),
		normVecDist(Eigen::Vector3d::Zero()),
		jacB1(rbd::Jacobian(mb, body1Id)),
		jacB2(rbd::Jacobian(mb, body2Id)),
		di(di),
		ds(ds),
		damping(damping),
		body1Id(body1Id),
		body2Id(body2Id),
		body1(mb.bodyIndexById(body1Id)),
		body2(mb.bodyIndexById(body2Id))
{
}



SelfCollisionConstr::SelfCollisionConstr(const rbd::MultiBody& mb, double step):
  dataVec_(),
  step_(step),
  nrVars_(0),
  AInEq_(),
  BInEq_(),
  fullJac_(6, mb.nrDof()),
  fullJacDot_(6, mb.nrDof()),
  alphaVec_(mb.nrDof()),
  calcVec_(mb.nrDof())
{
}


void SelfCollisionConstr::addCollision(const rbd::MultiBody& mb,
	int body1Id, SCD::S_Object* body1, const sva::PTransform& body1T,
	int body2Id, SCD::S_Object* body2, const sva::PTransform& body2T,
	double di, double ds, double damping)
{
	dataVec_.emplace_back(mb, body1Id, body1, body1T,
		body2Id, body2, body2T, di, ds, damping);
}


void SelfCollisionConstr::rmCollision(int body1Id, int body2Id)
{
	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
		[body1Id, body2Id](const CollData& data)
		{
			return data.body1Id == body1Id && data.body2Id == body2Id;
		});

	if(it != dataVec_.end())
	{
		delete it->pair;
		dataVec_.erase(it);
	}
}


void SelfCollisionConstr::reset()
{
	dataVec_.clear();
}


void SelfCollisionConstr::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& data)
{
	nrVars_ = data.nrVars();
}


void SelfCollisionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	if(static_cast<unsigned int>(AInEq_.rows()) != dataVec_.size()
		 || AInEq_.cols() != nrVars_)
	{
		AInEq_.resize(dataVec_.size(), nrVars_);
		BInEq_.resize(dataVec_.size());
		AInEq_.setZero();
		BInEq_.setZero();
	}

	rbd::paramToVector(mbc.alpha, alphaVec_);

	int i = 0;
	for(CollData& d: dataVec_)
	{
		SCD::Point3 pb1Tmp, pb2Tmp;

		d.pair->operator[](0)->setTransformation(toSCD(d.body1T*mbc.bodyPosW[d.body1]));
		d.pair->operator[](1)->setTransformation(toSCD(d.body2T*mbc.bodyPosW[d.body2]));

		double dist = d.pair->getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);

		Vector3d pb1(pb1Tmp[0], pb1Tmp[1], pb1Tmp[2]);
		Vector3d pb2(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]);

		Eigen::Vector3d normVecDist = (pb1 - pb2)/dist;

		pb1 = (sva::PTransform(pb1)*mbc.bodyPosW[d.body1].inv()).translation();
		pb2 = (sva::PTransform(pb2)*mbc.bodyPosW[d.body2].inv()).translation();

		if(dist < d.di)
		{
			double dampers = d.damping*((dist - d.ds)/(d.di - d.ds));

			Vector3d nf = normVecDist;
			Vector3d onf = d.normVecDist;
			Vector3d dnf = (nf - onf)/step_;

			// Compute body1
			d.jacB1.point(pb1);
			const MatrixXd& jac1 = d.jacB1.jacobian(mb, mbc);
			const MatrixXd& jacDot1 = d.jacB1.jacobianDot(mb, mbc);

			d.jacB1.fullJacobian(mb, jac1, fullJac_);
			d.jacB1.fullJacobian(mb, jacDot1, fullJacDot_);

			double jqdn = ((fullJac_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*nf)(0);
			double jqdnd = ((fullJac_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*dnf*step_)(0);
			double jdqdn = ((fullJacDot_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*nf*step_)(0);

			calcVec_ = -fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*nf*step_;

			// Compute body2
			d.jacB2.point(pb2);
			const MatrixXd& jac2 = d.jacB2.jacobian(mb, mbc);
			const MatrixXd& jacDot2 = d.jacB2.jacobianDot(mb, mbc);

			d.jacB2.fullJacobian(mb, jac2, fullJac_);
			d.jacB2.fullJacobian(mb, jacDot2, fullJacDot_);

			jqdn -= ((fullJac_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*nf)(0);
			jqdnd -= ((fullJac_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*dnf*step_)(0);
			jdqdn -= ((fullJacDot_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*nf*step_)(0);

			calcVec_ += fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*nf*step_;

			// distdot + distdotdot*dt > -damp*((d - ds)/(di - ds))
			AInEq_.block(i, 0, 1, mb.nrDof()) = calcVec_.transpose();
			BInEq_(i) = dampers + jqdn + jqdnd + jdqdn;
		}
		else
		{
			AInEq_.block(i, 0, 1, mb.nrDof()).setZero();
			BInEq_(i) = 0.;
		}

		d.normVecDist = normVecDist;
		++i;
	}
}


int SelfCollisionConstr::nrInEqLine()
{
	return dataVec_.size();
}


const Eigen::MatrixXd& SelfCollisionConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& SelfCollisionConstr::BInEq() const
{
	return BInEq_;
}


/**
	*													StaticEnvCollisionConstr
	*/


StaticEnvCollisionConstr::CollData::CollData(const rbd::MultiBody& mb,
	int bodyId, SCD::S_Object* body, const sva::PTransform& bodyT,
	int envId, SCD::S_Object* env,
	double di, double ds, double damping):
		pair(new SCD::CD_Pair(body, env)),
		bodyT(bodyT),
		normVecDist(Eigen::Vector3d::Zero()),
		jacB1(rbd::Jacobian(mb, bodyId)),
		di(di),
		ds(ds),
		damping(damping),
		bodyId(bodyId),
		envId(envId),
		body(mb.bodyIndexById(bodyId))
{
}


StaticEnvCollisionConstr::StaticEnvCollisionConstr(const rbd::MultiBody& mb, double step):
  dataVec_(),
  step_(step),
  nrVars_(0),
  AInEq_(),
  BInEq_(),
  fullJac_(6, mb.nrDof()),
  fullJacDot_(6, mb.nrDof()),
  alphaVec_(mb.nrDof()),
  calcVec_(mb.nrDof())
{
}


void StaticEnvCollisionConstr::addCollision(const rbd::MultiBody& mb,
	int bodyId, SCD::S_Object* body, const sva::PTransform& bodyT,
	int envId, SCD::S_Object* env,
	double di, double ds, double damping)
{
	dataVec_.emplace_back(mb, bodyId, body, bodyT, envId, env,
		di, ds, damping);
}


void StaticEnvCollisionConstr::rmCollision(int bodyId, int envId)
{
	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
		[bodyId, envId](const CollData& data)
		{
			return data.bodyId == bodyId && data.envId == envId;
		});

	if(it != dataVec_.end())
	{
		delete it->pair;
		dataVec_.erase(it);
	}
}


void StaticEnvCollisionConstr::reset()
{
	dataVec_.clear();
}


void StaticEnvCollisionConstr::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& data)
{
	nrVars_ = data.nrVars();
}


void StaticEnvCollisionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	if(static_cast<unsigned int>(AInEq_.rows()) != dataVec_.size()
		 || AInEq_.cols() != nrVars_)
	{
		AInEq_.resize(dataVec_.size(), nrVars_);
		BInEq_.resize(dataVec_.size());
		AInEq_.setZero();
		BInEq_.setZero();
	}

	rbd::paramToVector(mbc.alpha, alphaVec_);

	int i = 0;
	for(CollData& d: dataVec_)
	{
		SCD::Point3 pb1Tmp, pb2Tmp;

		d.pair->operator[](0)->setTransformation(toSCD(d.bodyT*mbc.bodyPosW[d.body]));

		double dist = d.pair->getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);

		Vector3d pb1(pb1Tmp[0], pb1Tmp[1], pb1Tmp[2]);
		Vector3d pb2(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]);

		Eigen::Vector3d normVecDist = (pb1 - pb2)/dist;

		pb1 = (sva::PTransform(pb1)*mbc.bodyPosW[d.body].inv()).translation();

		if(dist < d.di)
		{
			double dampers = d.damping*((dist - d.ds)/(d.di - d.ds));

			Vector3d nf = normVecDist;
			Vector3d onf = d.normVecDist;
			Vector3d dnf = (nf - onf)/step_;

			// Compute body
			d.jacB1.point(pb1);
			const MatrixXd& jac1 = d.jacB1.jacobian(mb, mbc);
			const MatrixXd& jacDot1 = d.jacB1.jacobianDot(mb, mbc);

			d.jacB1.fullJacobian(mb, jac1, fullJac_);
			d.jacB1.fullJacobian(mb, jacDot1, fullJacDot_);

			double jqdn = ((fullJac_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*nf)(0);
			double jqdnd = ((fullJac_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*dnf*step_)(0);
			double jdqdn = ((fullJacDot_.block(3, 0, 3, fullJac_.cols())*alphaVec_).transpose()*nf*step_)(0);

			calcVec_ = -fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*nf*step_;

			// distdot + distdotdot*dt > -damp*((d - ds)/(di - ds))
			AInEq_.block(i, 0, 1, mb.nrDof()) = calcVec_.transpose();
			BInEq_(i) = dampers + jqdn + jqdnd + jdqdn;
		}
		else
		{
			AInEq_.block(i, 0, 1, mb.nrDof()).setZero();
			BInEq_(i) = 0.;
		}

		d.normVecDist = normVecDist;
		++i;
	}
}


int StaticEnvCollisionConstr::nrInEqLine()
{
	return dataVec_.size();
}


const Eigen::MatrixXd& StaticEnvCollisionConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& StaticEnvCollisionConstr::BInEq() const
{
	return BInEq_;
}

} // namespace qp

} // namespace tasks
