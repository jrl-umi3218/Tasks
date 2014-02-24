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
#include <cmath>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// SCD
#include <SCD/CD/CD_Pair.h>
#include <SCD/S_Object/S_Object.h>

// Tasks
#include "utils.h"

namespace tasks
{

namespace qp
{


/**
	*															ContactCommon
	*/


bool ContactConstrCommon::addVirtualContact(int bodyId)
{
	return virtualContacts_.insert(bodyId).second;
}


bool ContactConstrCommon::removeVirtualContact(int bodyId)
{
	return virtualContacts_.erase(bodyId) == 1;
}


void ContactConstrCommon::resetVirtualContacts()
{
	virtualContacts_.clear();
}


std::set<int> ContactConstrCommon::bodyIdInContact(const rbd::MultiBody& mb,
	const SolverData& data)
{
	std::set<int> ret;
	auto isValid = [&mb, this](int bodyId)
	{
		// if is virtualContacts we don't add it
		// if fixed base and support body we don't add the contact
		return (virtualContacts_.find(bodyId) == virtualContacts_.end()) ||
			(!(bodyId == mb.body(0).id() &&
			mb.joint(0).type() == rbd::Joint::Fixed));
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
	A_(),
	ALU_(),
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

	A_.resize(cont_.size()*6, data.nrVars());
	ALU_.resize(cont_.size()*6);

	A_.setZero();
	ALU_.setZero();
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
		A_.block(i*6, 0, 6, mb.nrDof()) = fullJac_;

		// BEq = -JD_i*alpha
		const MatrixXd& jacDot = cont_[i].jac.jacobianDot(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jacDot, fullJac_);
		ALU_.segment(i*6, 6) = -fullJac_*alphaVec_;
	}
}


std::string ContactAccConstr::nameInEq() const
{
	return "ContactAccConstr";
}


std::string ContactAccConstr::descInEq(const rbd::MultiBody& mb, int line)
{
	int contact = line/6;
	int body = cont_[contact].jac.jointsPath().back();
	return std::string("Contact: ") + mb.body(body).name();
}


int ContactAccConstr::maxInEq() const
{
	return int(A_.rows());
}


const Eigen::MatrixXd& ContactAccConstr::AInEq() const
{
	return A_;
}


const Eigen::VectorXd& ContactAccConstr::LowerInEq() const
{
	return ALU_;
}


const Eigen::VectorXd& ContactAccConstr::UpperInEq() const
{
	return ALU_;
}


/**
	*															ContactSpeedConstr
	*/


ContactSpeedConstr::ContactSpeedConstr(const rbd::MultiBody& mb, double timeStep):
	cont_(),
	fullJac_(6, mb.nrDof()),
	alphaVec_(mb.nrDof()),
	A_(),
	ALU_(),
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

	A_.resize(cont_.size()*6, data.nrVars());
	ALU_.resize(cont_.size()*6);

	A_.setZero();
	ALU_.setZero();
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
		A_.block(i*6, 0, 6, mb.nrDof()) = fullJac_;

		// BEq = -JD_i*alpha
		const MatrixXd& jacDot = cont_[i].jac.jacobianDot(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jacDot, fullJac_);
		ALU_.segment(i*6, 6) = -fullJac_*alphaVec_ -
			mbc.bodyVelW[cont_[i].body].vector()/timeStep_;
	}
}


std::string ContactSpeedConstr::nameInEq() const
{
	return "ContactSpeedConstr";
}


std::string ContactSpeedConstr::descInEq(const rbd::MultiBody& mb, int line)
{
	int contact = line/6;
	int body = cont_[contact].jac.jointsPath().back();
	return std::string("Contact: ") + mb.body(body).name();
}


int ContactSpeedConstr::maxInEq() const
{
	return int(A_.rows());
}


const Eigen::MatrixXd& ContactSpeedConstr::AInEq() const
{
	return A_;
}


const Eigen::VectorXd& ContactSpeedConstr::LowerInEq() const
{
	return ALU_;
}


const Eigen::VectorXd& ContactSpeedConstr::UpperInEq() const
{
	return ALU_;
}


/**
	*															JointLimitsConstr
	*/


JointLimitsConstr::JointLimitsConstr(const rbd::MultiBody& mb,
	std::vector<std::vector<double> > lBound,
	std::vector<std::vector<double> > uBound,
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
	int vars = int(lower_.rows());

	rbd::paramToVector(mbc.q, qVec_);
	rbd::paramToVector(mbc.alpha, alphaVec_);

	lower_ = qMin_ - qVec_.tail(vars) - alphaVec_.tail(vars)*step_;
	lower_ /= dts;

	upper_ = qMax_ - qVec_.tail(vars) - alphaVec_.tail(vars)*step_;
	upper_ /= dts;
}


std::string JointLimitsConstr::nameBound() const
{
	return "JointLimitsConstr";
}


std::string JointLimitsConstr::descBound(const rbd::MultiBody& mb, int line)
{
	int jIndex = findJointFromVector(mb, line, false);
	return std::string("Joint: ") + mb.joint(jIndex).name();
}


int JointLimitsConstr::beginVar() const
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
	*												DamperJointLimitsConstr
	*/


DamperJointLimitsConstr::DamperJointLimitsConstr(const rbd::MultiBody& mb,
	const std::vector<std::vector<double> >& lBound,
	const std::vector<std::vector<double> >& uBound,
	std::vector<std::vector<double> > lVel,
	std::vector<std::vector<double> > uVel,
	double interPercent, double securityPercent,
	double damperOffset, double step):
	data_(),
	lower_(),
	upper_(),
	begin_(mb.joint(0).dof()),
	step_(step),
	damperOff_(damperOffset)
{
	int vars = mb.nrDof() - begin_;
	lower_.resize(vars);
	upper_.resize(vars);

	// remove the joint 0
	lVel[0] = {};
	uVel[0] = {};

	rbd::paramToVector(lVel, lower_);
	rbd::paramToVector(uVel, upper_);

	for(int i = 0; i < mb.nrJoints(); ++i)
	{
		if(mb.joint(i).dof() == 1)
		{
			double dist = (uBound[i][0] - lBound[i][0]);
			data_.emplace_back(lBound[i][0], uBound[i][0], lVel[i][0], uVel[i][0],
				dist*interPercent, dist*securityPercent,
				mb.jointPosInDof(i) - begin_, i);
		}
	}
}


void DamperJointLimitsConstr::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& /* data */)
{
}


void DamperJointLimitsConstr::update(const rbd::MultiBody& /* mb */,
	const rbd::MultiBodyConfig& mbc)
{
	for(DampData& d: data_)
	{
		double ld = mbc.q[d.index][0] - d.min;
		double ud = d.max - mbc.q[d.index][0];
		double alpha = mbc.alpha[d.index][0];

		lower_[d.vecPos] = (d.minVel - alpha)/step_;
		upper_[d.vecPos] = (d.maxVel - alpha)/step_;

		if(ld < d.iDist)
		{
			// damper(dist) < alpha
			// dist > 0 -> negative < alpha -> joint angle can decrease
			// dist < 0 -> positive < alpha -> joint angle must increase
			if(d.state != DampData::Low)
			{
				d.damping =
					std::abs(computeDamping(alpha, ld, d.iDist, d.sDist)) + damperOff_;
				d.state = DampData::Low;
			}

			double damper = -computeDamper(ld, d.iDist, d.sDist, d.damping);
			lower_[d.vecPos] = std::max((damper - alpha)/step_, lower_[d.vecPos]);
		}
		else if(ud < d.iDist)
		{
			// alpha < damper(dist)
			// dist > 0 -> alpha < positive -> joint angle can increase
			// dist < 0 -> alpha < negative -> joint angle must decrease
			if(d.state != DampData::Upp)
			{
				d.damping =
					std::abs(computeDamping(alpha, ud, d.iDist, d.sDist)) + damperOff_;
				d.state = DampData::Upp;
			}

			double damper = computeDamper(ud, d.iDist, d.sDist, d.damping);
			upper_[d.vecPos] = std::min((damper - alpha)/step_, upper_[d.vecPos]);
		}
		else
		{
			d.state = DampData::Free;
		}
	}
}


std::string DamperJointLimitsConstr::nameBound() const
{
	return "DamperJointLimitsConstr";
}


std::string DamperJointLimitsConstr::descBound(const rbd::MultiBody& mb, int line)
{
	int jIndex = findJointFromVector(mb, line, false);
	return std::string("Joint: ") + mb.joint(jIndex).name();
}


int DamperJointLimitsConstr::beginVar() const
{
	return begin_;
}


const Eigen::VectorXd& DamperJointLimitsConstr::Lower() const
{
	return lower_;
}


const Eigen::VectorXd& DamperJointLimitsConstr::Upper() const
{
	return upper_;
}


double DamperJointLimitsConstr::computeDamping(double alpha, double dist,
	double iDist, double sDist)
{
	return ((iDist - sDist)/(dist - sDist))*alpha;
}


double DamperJointLimitsConstr::computeDamper(double dist,
	double iDist, double sDist ,double damping)
{
	return damping*((dist - sDist)/(iDist - sDist));
}


/**
	*													SelfCollisionConstr
	*/


SCD::Matrix4x4 toSCD(const sva::PTransformd& t)
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


SelfCollisionConstr::CollData::CollData(const rbd::MultiBody& mb, int collId,
	int body1Id, SCD::S_Object* body1, const sva::PTransformd& body1T,
	int body2Id, SCD::S_Object* body2, const sva::PTransformd& body2T,
	double di, double ds, double damping, double dampOff):
		pair(new SCD::CD_Pair(body1, body2)),
		body1T(body1T),
		body2T(body2T),
		normVecDist(Eigen::Vector3d::Zero()),
		jacB1(rbd::Jacobian(mb, body1Id)),
		jacB2(rbd::Jacobian(mb, body2Id)),
		di(di),
		ds(ds),
		damping(damping),
		collId(collId),
		body1Id(body1Id),
		body2Id(body2Id),
		body1(mb.bodyIndexById(body1Id)),
		body2(mb.bodyIndexById(body2Id)),
		dampingType(damping > 0. ? DampingType::Hard : DampingType::Soft),
		dampingOff(dampOff)
{
}



SelfCollisionConstr::SelfCollisionConstr(const rbd::MultiBody& mb, double step):
  dataVec_(),
  step_(step),
  nrVars_(0),
  nrActivated_(0),
  AInEq_(),
  AL_(),
  AU_(),
  fullJac_(6, mb.nrDof()),
  fullJacDot_(6, mb.nrDof()),
  alphaVec_(mb.nrDof()),
  calcVec_(mb.nrDof())
{
}


void SelfCollisionConstr::addCollision(const rbd::MultiBody& mb, int collId,
	int body1Id, SCD::S_Object* body1, const sva::PTransformd& body1T,
	int body2Id, SCD::S_Object* body2, const sva::PTransformd& body2T,
	double di, double ds, double damping, double dampingOff)
{
	dataVec_.emplace_back(mb, collId, body1Id, body1, body1T,
		body2Id, body2, body2T, di, ds, damping, dampingOff);
}


bool SelfCollisionConstr::rmCollision(int collId)
{
	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
		[collId](const CollData& data)
		{
			return data.collId == collId;
		});

	if(it != dataVec_.end())
	{
		delete it->pair;
		dataVec_.erase(it);
		return true;
	}

	return false;
}


std::size_t SelfCollisionConstr::nrCollisions() const
{
	return dataVec_.size();
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
		AL_.resize(dataVec_.size());
		AU_.resize(dataVec_.size());
		AInEq_.setZero();
		AL_.fill(-std::numeric_limits<double>::infinity());
		AU_.setZero();
	}

	rbd::paramToVector(mbc.alpha, alphaVec_);

	nrActivated_ = 0;
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

		pb1 = (sva::PTransformd(pb1)*mbc.bodyPosW[d.body1].inv()).translation();
		pb2 = (sva::PTransformd(pb2)*mbc.bodyPosW[d.body2].inv()).translation();

		if(dist < d.di)
		{
			if(d.dampingType == CollData::DampingType::Free)
			{
				d.dampingType = CollData::DampingType::Soft;
				Vector3d v1(mbc.bodyPosW[d.body1].rotation().transpose()*
						(sva::PTransformd(pb1)*mbc.bodyVelB[d.body1]).linear());
				Vector3d v2(mbc.bodyPosW[d.body2].rotation().transpose()*
						(sva::PTransformd(pb2)*mbc.bodyVelB[d.body2]).linear());
				double distDot = std::abs((v1 - v2).dot(normVecDist));

				/// @todo find a bette solution.
				// use a value slightly upper ds if dist <= ds
				double fixedDist = dist <= d.ds ? d.ds + (d.di - d.ds)*0.2 : dist;
				d.damping = ((d.di - d.ds)/(fixedDist - d.ds))*distDot + d.dampingOff;
			}

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
			AInEq_.block(nrActivated_, 0, 1, mb.nrDof()) = calcVec_.transpose();
			AU_(nrActivated_) = dampers + jqdn + jqdnd + jdqdn;
			++nrActivated_;
		}
		else
		{
			if(d.dampingType == CollData::DampingType::Soft)
			{
				d.dampingType = CollData::DampingType::Free;
			}
		}

		d.normVecDist = normVecDist;
	}
}


std::string SelfCollisionConstr::nameInEq() const
{
	return "SelfCollisionConstr";
}


std::string SelfCollisionConstr::descInEq(const rbd::MultiBody& mb, int line)
{
	int curLine = 0;
	for(CollData& d: dataVec_)
	{
		double dist = d.pair->getDistance();
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
		if(dist < d.di)
		{
			if(curLine == line)
			{
				std::stringstream ss;
				ss << mb.body(d.body1).name() << " / " << mb.body(d.body2).name() << std::endl;
				ss << "collId: " << d.collId << std::endl;
				ss << "dist: " << dist << std::endl;
				ss << "di: " << d.di << std::endl;
				ss << "ds: " << d.ds << std::endl;
				ss << "damp: " << d.damping + d.dampingOff << std::endl;
				return ss.str();
			}
			++curLine;
		}
	}
	return "";
}


int SelfCollisionConstr::nrInEq() const
{
	return nrActivated_;
}


int SelfCollisionConstr::maxInEq() const
{
	return int(dataVec_.size());
}


const Eigen::MatrixXd& SelfCollisionConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& SelfCollisionConstr::LowerInEq() const
{
	return AL_;
}


const Eigen::VectorXd& SelfCollisionConstr::UpperInEq() const
{
	return AU_;
}


/**
	*													StaticEnvCollisionConstr
	*/


StaticEnvCollisionConstr::CollData::CollData(const rbd::MultiBody& mb, int collId,
	int bodyId, SCD::S_Object* body, const sva::PTransformd& bodyT,
	int envId, SCD::S_Object* env,
	double di, double ds, double damping, double dampOff):
		pair(new SCD::CD_Pair(body, env)),
		bodyT(bodyT),
		normVecDist(Eigen::Vector3d::Zero()),
		jacB1(rbd::Jacobian(mb, bodyId)),
		di(di),
		ds(ds),
		damping(damping),
		collId(collId),
		bodyId(bodyId),
		envId(envId),
		body(mb.bodyIndexById(bodyId)),
		dampingType(damping > 0. ? DampingType::Hard : DampingType::Soft),
		dampingOff(dampOff)
{
}


StaticEnvCollisionConstr::StaticEnvCollisionConstr(const rbd::MultiBody& mb, double step):
  dataVec_(),
  step_(step),
  nrVars_(0),
  nrActivated_(0),
  AInEq_(),
  AL_(),
  AU_(),
  fullJac_(6, mb.nrDof()),
  fullJacDot_(6, mb.nrDof()),
  alphaVec_(mb.nrDof()),
  calcVec_(mb.nrDof())
{
}


void StaticEnvCollisionConstr::addCollision(const rbd::MultiBody& mb, int collId,
	int bodyId, SCD::S_Object* body, const sva::PTransformd& bodyT,
	int envId, SCD::S_Object* env,
	double di, double ds, double damping, double dampingOff)
{
	dataVec_.emplace_back(mb, collId, bodyId, body, bodyT, envId, env,
		di, ds, damping, dampingOff);
}


bool StaticEnvCollisionConstr::rmCollision(int collId)
{
	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
		[collId](const CollData& data)
		{
			return data.collId == collId;
		});

	if(it != dataVec_.end())
	{
		delete it->pair;
		dataVec_.erase(it);
		return true;
	}

	return false;
}


std::size_t StaticEnvCollisionConstr::nrCollisions() const
{
	return dataVec_.size();
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
		AL_.resize(dataVec_.size());
		AU_.resize(dataVec_.size());
		AInEq_.setZero();
		AL_.fill(-std::numeric_limits<double>::infinity());
		AU_.setZero();
	}

	rbd::paramToVector(mbc.alpha, alphaVec_);

	nrActivated_ = 0;
	for(CollData& d: dataVec_)
	{
		SCD::Point3 pb1Tmp, pb2Tmp;

		d.pair->operator[](0)->setTransformation(toSCD(d.bodyT*mbc.bodyPosW[d.body]));

		double dist = d.pair->getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);

		Vector3d pb1(pb1Tmp[0], pb1Tmp[1], pb1Tmp[2]);
		Vector3d pb2(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]);

		Eigen::Vector3d normVecDist = (pb1 - pb2)/dist;

		pb1 = (sva::PTransformd(pb1)*mbc.bodyPosW[d.body].inv()).translation();

		if(dist < d.di)
		{
			if(d.dampingType == CollData::DampingType::Free)
			{
				d.dampingType = CollData::DampingType::Soft;
				Vector3d v1(mbc.bodyPosW[d.body].rotation().transpose()*
						(sva::PTransformd(pb1)*mbc.bodyVelB[d.body]).linear());
				double distDot = std::abs(v1.dot(normVecDist));

				/// @todo find a bette solution.
				// use a value slightly upper ds if dist <= ds
				double fixedDist = dist <= d.ds ? d.ds + (d.di - d.ds)*0.2 : dist;
				d.damping = ((d.di - d.ds)/(fixedDist - d.ds))*distDot + d.dampingOff;
			}

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
			AInEq_.block(nrActivated_, 0, 1, mb.nrDof()) = calcVec_.transpose();
			AU_(nrActivated_) = dampers + jqdn + jqdnd + jdqdn;
			++nrActivated_;
		}
		else
		{
			if(d.dampingType == CollData::DampingType::Soft)
			{
				d.dampingType = CollData::DampingType::Free;
			}
		}

		d.normVecDist = normVecDist;
	}
}


std::string StaticEnvCollisionConstr::nameInEq() const
{
	return "StaticEnvCollisionConstr";
}


std::string StaticEnvCollisionConstr::descInEq(const rbd::MultiBody& mb, int line)
{
	int curLine = 0;
	for(CollData& d: dataVec_)
	{
		double dist = d.pair->getDistance();
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
		if(dist < d.di)
		{
			if(curLine == line)
			{
				std::stringstream ss;
				ss << mb.body(d.body).name() << " / " << d.envId << std::endl;
				ss << "collId: " << d.collId << std::endl;
				ss << "dist: " << dist << std::endl;
				ss << "di: " << d.di << std::endl;
				ss << "ds: " << d.ds << std::endl;
				ss << "damp: " << d.damping + d.dampingOff << std::endl;
				return ss.str();
			}
			++curLine;
		}
	}
	return "";
}


int StaticEnvCollisionConstr::nrInEq() const
{
	return nrActivated_;
}


int StaticEnvCollisionConstr::maxInEq() const
{
	return int(dataVec_.size());
}


const Eigen::MatrixXd& StaticEnvCollisionConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& StaticEnvCollisionConstr::LowerInEq() const
{
	return AL_;
}


const Eigen::VectorXd& StaticEnvCollisionConstr::UpperInEq() const
{
	return AU_;
}

/**
	*													CoMCollisionConstr
	*/


CoMCollisionConstr::CollData::CollData(const rbd::MultiBody& mb,
	int collId, SCD::S_Object* env,
	double di, double ds, double damping, double dampOff):
		comSphere_(ds/5.),
		pair(new SCD::CD_Pair(&comSphere_, env)),
		normVecDist(Eigen::Vector3d::Zero()),
		jacCoM(rbd::CoMJacobian(mb)),
		di(di),
		ds(ds),
		damping(damping),
		collId(collId),
		dampingType(damping > 0. ? DampingType::Hard : DampingType::Soft),
		dampingOff(dampOff)
{
}


CoMCollisionConstr::CoMCollisionConstr(const rbd::MultiBody& mb, double step):
  dataVec_(),
  step_(step),
  nrVars_(0),
  nrActivated_(0),
  AInEq_(),
  AL_(),
  AU_(),
  fullJac_(6, mb.nrDof()),
  fullJacDot_(6, mb.nrDof()),
  alphaVec_(mb.nrDof()),
  calcVec_(mb.nrDof())
{
}


void CoMCollisionConstr::addCollision(const rbd::MultiBody& mb, int collId,
	SCD::S_Object* env,
	double di, double ds, double damping, double dampingOff)
{
	dataVec_.emplace_back(mb, collId, env,
		di, ds, damping, dampingOff);
}


bool CoMCollisionConstr::rmCollision(int collId)
{
	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
		[collId](const CollData& data)
		{
			return data.collId == collId;
		});

	if(it != dataVec_.end())
	{
		delete it->pair;
		dataVec_.erase(it);
		return true;
	}

	return false;
}


std::size_t CoMCollisionConstr::nrCollisions() const
{
	return dataVec_.size();
}


void CoMCollisionConstr::reset()
{
	dataVec_.clear();
}


void CoMCollisionConstr::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& data)
{
	nrVars_ = data.nrVars();
}


void CoMCollisionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	if(static_cast<unsigned int>(AInEq_.rows()) != dataVec_.size()
		 || AInEq_.cols() != nrVars_)
	{
		AInEq_.resize(dataVec_.size(), nrVars_);
		AL_.resize(dataVec_.size());
		AU_.resize(dataVec_.size());
		AInEq_.setZero();
		AL_.fill(-std::numeric_limits<double>::infinity());
		AU_.setZero();
	}

	rbd::paramToVector(mbc.alpha, alphaVec_);
	Eigen::Vector3d com = rbd::computeCoM(mb, mbc);

	nrActivated_ = 0;
	for(CollData& d: dataVec_)
	{
		SCD::Point3 pb1Tmp, pb2Tmp;

		d.pair->operator[](0)->setTransformation(toSCD(sva::PTransformd(com)));

		double dist = d.pair->getClosestPoints(pb1Tmp, pb2Tmp);
		Vector3d pb2(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]);

		dist = -std::copysign((com - pb2).norm(), dist);

		Eigen::Vector3d normVecDist = (com - pb2)/dist;

		if(dist < d.di)
		{
			/*if(d.dampingType == CollData::DampingType::Free)
			{
				d.dampingType = CollData::DampingType::Soft;
				Vector3d v1(mbc.bodyPosW[d.body].rotation().transpose()*
						(sva::PTransformd(pb1)*mbc.bodyVelB[d.body]).linear());
				double distDot = std::abs(v1.dot(normVecDist));

				/// @todo find a bette solution.
				// use a value slightly upper ds if dist <= ds
				double fixedDist = dist <= d.ds ? d.ds + (d.di - d.ds)*0.2 : dist;
				d.damping = ((d.di - d.ds)/(fixedDist - d.ds))*distDot + d.dampingOff;
			}*/

			double dampers = d.damping*((dist - d.ds)/(d.di - d.ds));

			Vector3d nf = normVecDist;
			Vector3d onf = d.normVecDist;
			Vector3d dnf = (nf - onf)/step_;

			// Compute body
			//d.jacB1.point(pb1);
			const MatrixXd& jac1 = d.jacCoM.jacobian(mb, mbc);
			const MatrixXd& jacDot1 = d.jacCoM.jacobianDot(mb, mbc);

			//d.jacB1.fullJacobian(mb, jac1, fullJac_);
			//d.jacB1.fullJacobian(mb, jacDot1, fullJacDot_);

			double jqdn = ((jac1*alphaVec_).transpose()*nf)(0);
			double jqdnd = ((jac1*alphaVec_).transpose()*dnf*step_)(0);
			double jdqdn = ((jacDot1*alphaVec_).transpose()*nf*step_)(0);

			calcVec_ = -jac1.transpose()*nf*step_;

			// distdot + distdotdot*dt > -damp*((d - ds)/(di - ds))
			AInEq_.block(nrActivated_, 0, 1, mb.nrDof()) = calcVec_.transpose();
			AU_(nrActivated_) = dampers + jqdn + jqdnd + jdqdn;
			++nrActivated_;
		}
		else
		{
			if(d.dampingType == CollData::DampingType::Soft)
			{
				d.dampingType = CollData::DampingType::Free;
			}
		}

		d.normVecDist = normVecDist;
	}
}


std::string CoMCollisionConstr::nameInEq() const
{
	return "CoMCollisionConstr";
}


std::string CoMCollisionConstr::descInEq(const rbd::MultiBody& /*mb*/, int line)
{
	int curLine = 0;
	for(CollData& d: dataVec_)
	{
		double dist = d.pair->getDistance();
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
		if(dist < d.di)
		{
			if(curLine == line)
			{
				std::stringstream ss;
				ss << "collId: " << d.collId << std::endl;
				ss << "dist: " << dist << std::endl;
				ss << "di: " << d.di << std::endl;
				ss << "ds: " << d.ds << std::endl;
				ss << "damp: " << d.damping + d.dampingOff << std::endl;
				return ss.str();
			}
			++curLine;
		}
	}
	return "";
}


int CoMCollisionConstr::nrInEq() const
{
	return nrActivated_;
}


int CoMCollisionConstr::maxInEq() const
{
	return int(dataVec_.size());
}


const Eigen::MatrixXd& CoMCollisionConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& CoMCollisionConstr::LowerInEq() const
{
	return AL_;
}


const Eigen::VectorXd& CoMCollisionConstr::UpperInEq() const
{
	return AU_;
}


/**
	*													GripperTorqueConstr
	*/


GripperTorqueConstr::GripperData::GripperData(int bId, double tl,
  const Eigen::Vector3d& o, const Eigen::Vector3d& a):
  bodyId(bId),
  torqueLimit(tl),
  origin(o),
  axis(a)
{}


GripperTorqueConstr::GripperTorqueConstr():
  dataVec_(),
  AInEq_(),
  AL_(),
  AU_()
{}


void GripperTorqueConstr::addGripper(int bodyId, double torqueLimit,
	const Eigen::Vector3d& origin, const Eigen::Vector3d& axis)
{
	dataVec_.emplace_back(bodyId, torqueLimit, origin, axis);
}


bool GripperTorqueConstr::rmGripper(int bodyId)
{
	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
		[bodyId](const GripperData& data)
		{
			return data.bodyId == bodyId;
		});

	if(it != dataVec_.end())
	{
		dataVec_.erase(it);
		return true;
	}

	return false;
}


void GripperTorqueConstr::reset()
{
	dataVec_.clear();
}


void GripperTorqueConstr::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& data)
{
	using namespace Eigen;
	AInEq_.setZero(dataVec_.size(), data.nrVars());
	AL_.resize(dataVec_.size());
	AL_.fill(-std::numeric_limits<double>::infinity());
	AU_.setZero(dataVec_.size());

	int line = 0;
	for(const GripperData& gd: dataVec_)
	{
		int begin = data.bilateralBegin();
		for(const BilateralContact& bc: data.bilateralContacts())
		{
			int curLambda = 0;
			// compute the number of lambda needed by the current bilateral
			for(std::size_t i = 0; i < bc.points.size(); ++i)
			{
				curLambda += bc.nrLambda(static_cast<int>(i));
			}

			if(bc.bodyId == gd.bodyId)
			{
				int col = begin;
				// Torque applied on the gripper motor
				// Sum_i^nrF  T_i·( p_i^T_o x f_i)
				for(std::size_t i = 0; i < bc.cones.size(); ++i)
				{
					Vector3d T_o_p = bc.points[i] - gd.origin;
					for(std::size_t j = 0; j < bc.cones[i].generators.size(); ++j)
					{
						// we use abs because the contact force cannot apply
						// negative torque on the gripper
						AInEq_(line, col) = std::abs(
							gd.axis.transpose()*(T_o_p.cross(bc.cones[i].generators[j])));
						++col;
					}
				}
				AU_(line) = gd.torqueLimit;
				break;
			}

			begin += curLambda;
			// if the bodyId is not found the AInEq_ and BInEq_ line stay at zero
		}

		++line;
	}
}


void GripperTorqueConstr::update(const rbd::MultiBody& /* mb */,
	const rbd::MultiBodyConfig& /* mbc */)
{}


std::string GripperTorqueConstr::nameInEq() const
{
	return "GripperTorqueConstr";
}


std::string GripperTorqueConstr::descInEq(const rbd::MultiBody& mb, int line)
{
	std::stringstream ss;
	int bodyIndex = mb.bodyIndexById(dataVec_[line].bodyId);
	ss << mb.body(bodyIndex).name() << std::endl;
	ss << "limits: " << dataVec_[line].torqueLimit << std::endl;
	return ss.str();
}


int GripperTorqueConstr::maxInEq() const
{
	return static_cast<int>(dataVec_.size());
}


const Eigen::MatrixXd& GripperTorqueConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& GripperTorqueConstr::LowerInEq() const
{
	return AL_;
}


const Eigen::VectorXd& GripperTorqueConstr::UpperInEq() const
{
	return AU_;
}

} // namespace qp

} // namespace tasks
