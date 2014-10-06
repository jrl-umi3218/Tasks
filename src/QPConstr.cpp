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

// sch
#include <sch/CD/CD_Pair.h>
#include <sch/S_Object/S_Object.h>

// Tasks
#include "Bounds.h"
#include "utils.h"

namespace tasks
{

namespace qp
{


/**
	*															ContactCommon
	*/


bool ContactConstrCommon::ContactCommon::operator==(const ContactCommon& cc) const
{
	return cId == cc.cId;
}


bool ContactConstrCommon::ContactCommon::operator<(const ContactCommon& cc) const
{
	return cId < cc.cId;
}


bool ContactConstrCommon::addVirtualContact(const ContactId& cId)
{
	return virtualContacts_.insert(cId).second;
}


bool ContactConstrCommon::removeVirtualContact(const ContactId& cId)
{
	return virtualContacts_.erase(cId) == 1;
}


void ContactConstrCommon::resetVirtualContacts()
{
	virtualContacts_.clear();
}


bool ContactConstrCommon::addDofContact(const ContactId& cId,
	const Eigen::MatrixXd& dof)
{
	return dofContacts_.insert({cId, dof}).second;
}


bool ContactConstrCommon::removeDofContact(const ContactId& cId)
{
	return dofContacts_.erase(cId) == 1;
}


void ContactConstrCommon::resetDofContacts()
{
	dofContacts_.clear();
}


std::set<ContactConstrCommon::ContactCommon>
ContactConstrCommon::contactCommonInContact(
	const std::vector<rbd::MultiBody>& mbs, const SolverData& data)
{
	std::set<ContactCommon> ret;
	auto isValid = [&mbs, this](const ContactId& contactId)
	{
		// if is virtualContacts we don't add it
		return (virtualContacts_.find(contactId) == virtualContacts_.end());
	};

	for(const BilateralContact& c: data.allContacts())
	{
		if(isValid(c.contactId))
		{
			ret.insert({c.contactId, c.X_b1_b2});
		}
	}

	return std::move(ret);
}


/**
	*															ContactAccConstr
	*/


ContactAccConstr::ContactAccConstr():
	cont_(),
	fullJac_(),
	A_(),
	b_(),
	nrEq_(0)
{}


void ContactAccConstr::updateDofContacts()
{
	for(ContactData& c: cont_)
	{
		auto it = dofContacts_.find(c.contactId);
		if(it != dofContacts_.end())
		{
			c.dof = it->second;
		}
		else
		{
			c.dof.setIdentity(6, 6);
		}
	}
	updateNrEq();
}


void ContactAccConstr::updateNrVars(const std::vector<rbd::MultiBody>& mbs,
	const SolverData& data)
{
	cont_.clear();
	fullJac_.resize(mbs.size());

	for(std::size_t i = 0; i < mbs.size(); ++i)
	{
		fullJac_[i].resize(6, mbs[i].nrDof());
	}

	std::set<ContactCommon> contactCSet = contactCommonInContact(mbs, data);
	for(const ContactCommon& cC: contactCSet)
	{
		Eigen::MatrixXd dof(Eigen::MatrixXd::Identity(6, 6));
		auto it = dofContacts_.find(cC.cId);
		if(it != dofContacts_.end())
		{
			dof = it->second;
		}
		std::vector<ContactSideData> contacts;
		auto addContact = [&mbs, &data, &contacts](int rIndex, int bId,
			double sign, const Eigen::Vector3d& point)
		{
			if(mbs[rIndex].nrDof() > 0)
			{
				contacts.emplace_back(rIndex, data.alphaDBegin(rIndex), sign,
															rbd::Jacobian(mbs[rIndex], bId, point));
			}
		};
		addContact(cC.cId.r1Index, cC.cId.r1BodyId, 1., cC.X_b1_b2.translation());
		addContact(cC.cId.r2Index, cC.cId.r2BodyId, -1., Eigen::Vector3d::Zero());

		cont_.emplace_back(std::move(contacts), dof, cC.cId);
	}
	updateNrEq();

	A_.setZero(cont_.size()*6, data.nrVars());
	b_.setZero(cont_.size()*6);
}


void ContactAccConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	using namespace Eigen;

	A_.setZero();
	b_.setZero();
	// J_i*alphaD + JD_i*alpha = 0

	int index = 0;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		ContactData& cd = cont_[i];
		int rows = int(cd.dof.rows());

		for(std::size_t j = 0; j < cd.contacts.size(); ++j)
		{
			ContactSideData& csd = cd.contacts[j];
			const rbd::MultiBody& mb = mbs[csd.robotIndex];
			const rbd::MultiBodyConfig& mbc = mbcs[csd.robotIndex];
			Eigen::MatrixXd& fullJac = fullJac_[csd.robotIndex];

			// AEq = J_i
			const MatrixXd& jacMat = csd.jac.jacobian(mb, mbc);
			csd.jac.fullJacobian(mb, jacMat, fullJac);
			/// TODO don't apply dof on full jac
			A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() +=
				csd.sign*cd.dof*fullJac;

			// BEq = -JD_i*alpha
			Vector6d normalAcc = csd.jac.normalAcceleration(
				mb, mbc, data.normalAccB(csd.robotIndex)).vector();
			b_.segment(index, rows).noalias() -= csd.sign*cd.dof*normalAcc;
		}
		index += rows;
	}
}


std::string ContactAccConstr::nameEq() const
{
	return "ContactAccConstr";
}


std::string ContactAccConstr::descEq(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	std::ostringstream oss;
	int contact = line/6;
	for(const ContactSideData& csd: cont_[contact].contacts)
	{
		int body = csd.jac.jointsPath().back();
		oss << "Contact: " << mbs[csd.robotIndex].body(body).name() << std::endl;
	}
	return oss.str();
}


int ContactAccConstr::nrEq() const
{
	return nrEq_;
}


int ContactAccConstr::maxEq() const
{
	return int(A_.rows());
}


const Eigen::MatrixXd& ContactAccConstr::AEq() const
{
	return A_;
}


const Eigen::VectorXd& ContactAccConstr::bEq() const
{
	return b_;
}


void ContactAccConstr::updateNrEq()
{
	nrEq_ = 0;
	for(const ContactData& c: cont_)
	{
		nrEq_ += int(c.dof.rows());
	}
}


/**
	*															ContactSpeedConstr
	*/


ContactSpeedConstr::ContactSpeedConstr(double timeStep):
	cont_(),
	fullJac_(),
	A_(),
	b_(),
	nrEq_(0),
	timeStep_(timeStep)
{}


void ContactSpeedConstr::updateDofContacts()
{
	for(ContactData& c: cont_)
	{
		auto it = dofContacts_.find(c.contactId);
		if(it != dofContacts_.end())
		{
			c.dof = it->second;
		}
		else
		{
			c.dof.setIdentity(6, 6);
		}
	}
	updateNrEq();
}


void ContactSpeedConstr::updateNrVars(const std::vector<rbd::MultiBody>& mbs,
	const SolverData& data)
{
	cont_.clear();
	fullJac_.resize(mbs.size());

	for(std::size_t i = 0; i < mbs.size(); ++i)
	{
		fullJac_[i].resize(6, mbs[i].nrDof());
	}

	std::set<ContactCommon> contactCSet = contactCommonInContact(mbs, data);
	for(const ContactCommon& cC: contactCSet)
	{
		Eigen::MatrixXd dof(Eigen::MatrixXd::Identity(6, 6));
		auto it = dofContacts_.find(cC.cId);
		if(it != dofContacts_.end())
		{
			dof = it->second;
		}
		std::vector<ContactSideData> contacts;
		auto addContact = [&mbs, &data, &contacts](int rIndex, int bId,
			double sign, const Eigen::Vector3d& point)
		{
			if(mbs[rIndex].nrDof() > 0)
			{
				contacts.emplace_back(rIndex, data.alphaDBegin(rIndex), sign,
															rbd::Jacobian(mbs[rIndex], bId, point));
			}
		};
		addContact(cC.cId.r1Index, cC.cId.r1BodyId, 1., cC.X_b1_b2.translation());
		addContact(cC.cId.r2Index, cC.cId.r2BodyId, -1., Eigen::Vector3d::Zero());

		cont_.emplace_back(std::move(contacts), dof, cC.cId);
	}
	updateNrEq();

	A_.setZero(cont_.size()*6, data.nrVars());
	b_.setZero(cont_.size()*6);
}


void ContactSpeedConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	using namespace Eigen;

	A_.setZero();
	b_.setZero();
	// J_i*alphaD + JD_i*alpha = 0

	int index = 0;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		ContactData& cd = cont_[i];
		int rows = int(cd.dof.rows());

		for(std::size_t j = 0; j < cd.contacts.size(); ++j)
		{
			ContactSideData& csd = cd.contacts[j];
			const rbd::MultiBody& mb = mbs[csd.robotIndex];
			const rbd::MultiBodyConfig& mbc = mbcs[csd.robotIndex];
			Eigen::MatrixXd& fullJac = fullJac_[csd.robotIndex];

			// AEq = J_i
			const MatrixXd& jacMat = csd.jac.jacobian(mb, mbc);
			csd.jac.fullJacobian(mb, jacMat, fullJac);
			/// TODO don't apply dof on full jac
			A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() +=
				csd.sign*cd.dof*fullJac;

			// BEq = -JD_i*alpha
			Vector6d normalAcc = csd.jac.normalAcceleration(
				mb, mbc, data.normalAccB(csd.robotIndex)).vector();
			b_.segment(index, rows).noalias() -= csd.sign*cd.dof*(normalAcc +
				mbc.bodyVelW[csd.bodyIndex].vector()/timeStep_);
		}
		index += rows;
	}
}


int ContactSpeedConstr::nrEq() const
{
	return nrEq_;
}


std::string ContactSpeedConstr::nameEq() const
{
	return "ContactSpeedConstr";
}


std::string ContactSpeedConstr::descEq(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	std::ostringstream oss;
	int contact = line/6;
	for(const ContactSideData& csd: cont_[contact].contacts)
	{
		int body = csd.jac.jointsPath().back();
		oss << "Contact: " << mbs[csd.robotIndex].body(body).name() << std::endl;
	}
	return oss.str();
}


int ContactSpeedConstr::maxEq() const
{
	return int(A_.rows());
}


const Eigen::MatrixXd& ContactSpeedConstr::AEq() const
{
	return A_;
}


const Eigen::VectorXd& ContactSpeedConstr::bEq() const
{
	return b_;
}


void ContactSpeedConstr::updateNrEq()
{
	nrEq_ = 0;
	for(const ContactData& c: cont_)
	{
		nrEq_ += int(c.dof.rows());
	}
}


/**
	*															JointLimitsConstr
	*/


JointLimitsConstr::JointLimitsConstr(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, QBound bound, double step):
	robotIndex_(robotIndex),
	alphaDBegin_(-1),
	alphaDOffset_(mbs[robotIndex].joint(0).dof() > 1 ? mbs[robotIndex].nrDof() : 0),
	step_(step),
	qMin_(),
	qMax_(),
	qVec_(),
	alphaVec_(),
	lower_(),
	upper_()
{
	assert(std::size_t(robotIndex_) < mbs.size() && robotIndex_ >= 0);

	const rbd::MultiBody& mb = mbs[robotIndex_];

	// we don't manage joint with more than 1
	/// @todo remove this dirty hack
	int nrVars = mb.nrDof() - alphaDOffset_;
	qMin_.resize(nrVars);
	qMax_.resize(nrVars);
	qVec_.resize(mb.nrParams());
	alphaVec_.resize(mb.nrDof());

	// if first joint is not managed remove it
	if(alphaDOffset_ != 0)
	{
		bound.lQBound[0] = {};
		bound.uQBound[0] = {};
	}

	rbd::paramToVector(bound.lQBound, qMin_);
	rbd::paramToVector(bound.uQBound, qMax_);

	lower_.setConstant(nrVars, -std::numeric_limits<double>::infinity());
	upper_.setConstant(nrVars, std::numeric_limits<double>::infinity());
}


void JointLimitsConstr::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_) + alphaDOffset_;
}


void JointLimitsConstr::update(const std::vector<rbd::MultiBody>& /* mbs */,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& /* data */)
{
	const rbd::MultiBodyConfig& mbc = mbcs[robotIndex_];

	double dts = step_*step_*0.5;

	int vars = int(qMin_.rows());

	rbd::paramToVector(mbc.q, qVec_);
	rbd::paramToVector(mbc.alpha, alphaVec_);

	lower_.noalias() = qMin_ - qVec_.tail(vars) - alphaVec_.tail(vars)*step_;
	lower_ /= dts;

	upper_.noalias() = qMax_ - qVec_.tail(vars) - alphaVec_.tail(vars)*step_;
	upper_ /= dts;
}


std::string JointLimitsConstr::nameBound() const
{
	return "JointLimitsConstr";
}


std::string JointLimitsConstr::descBound(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	int jIndex = findJointFromVector(mbs[robotIndex_], line, false);
	return std::string("Joint: ") + mbs[robotIndex_].joint(jIndex).name();
}


int JointLimitsConstr::beginVar() const
{
	return alphaDBegin_;
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


DamperJointLimitsConstr::DamperJointLimitsConstr(
	const std::vector<rbd::MultiBody>& mbs, int robotIndex,
	const QBound& qBound, const AlphaBound& aBound,
	double interPercent, double securityPercent,
	double damperOffset, double step):
	robotIndex_(robotIndex),
	alphaDBegin_(-1),
	data_(),
	lower_(mbs[robotIndex].nrDof()),
	upper_(mbs[robotIndex].nrDof()),
	step_(step),
	damperOff_(damperOffset)
{
	assert(std::size_t(robotIndex_) < mbs.size() && robotIndex_ >= 0);

	const rbd::MultiBody& mb = mbs[robotIndex_];

	for(int i = 0; i < mb.nrJoints(); ++i)
	{
		if(mb.joint(i).dof() == 1)
		{
			double dist = (qBound.uQBound[i][0] - qBound.lQBound[i][0]);
			data_.emplace_back(qBound.lQBound[i][0], qBound.uQBound[i][0],
				aBound.lAlphaBound[i][0], aBound.uAlphaBound[i][0],
				dist*interPercent, dist*securityPercent,
				mb.jointPosInDof(i), i);
		}
	}

	rbd::paramToVector(aBound.lAlphaBound, lower_);
	rbd::paramToVector(aBound.uAlphaBound, upper_);
}


void DamperJointLimitsConstr::updateNrVars(
	const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
}


void DamperJointLimitsConstr::update(const std::vector<rbd::MultiBody>& /* mbs */,
	const std::vector<rbd::MultiBodyConfig>& mbcs, const SolverData& /* data */)
{
	const rbd::MultiBodyConfig& mbc = mbcs[robotIndex_];

	for(DampData& d: data_)
	{
		double ld = mbc.q[d.jointIndex][0] - d.min;
		double ud = d.max - mbc.q[d.jointIndex][0];
		double alpha = mbc.alpha[d.jointIndex][0];

		lower_[d.alphaDBegin] = (d.minVel - alpha)/step_;
		upper_[d.alphaDBegin] = (d.maxVel - alpha)/step_;

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
			lower_[d.alphaDBegin] = std::max((damper - alpha)/step_,
				lower_[d.alphaDBegin]);
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
			upper_[d.alphaDBegin] = std::min((damper - alpha)/step_,
				upper_[d.alphaDBegin]);
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


std::string DamperJointLimitsConstr::descBound(
	const std::vector<rbd::MultiBody>& mbs, int line)
{
	int jIndex = findJointFromVector(mbs[robotIndex_], line, false);
	return std::string("Joint: ") + mbs[robotIndex_].joint(jIndex).name();
}


int DamperJointLimitsConstr::beginVar() const
{
	return alphaDBegin_;
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


///**
//	*													SelfCollisionConstr
//	*/


//sch::Matrix4x4 tosch(const sva::PTransformd& t)
//{
//	sch::Matrix4x4 m;
//	const Eigen::Matrix3d& rot = t.rotation();
//	const Eigen::Vector3d& tran = t.translation();

//	for(int i = 0; i < 3; ++i)
//	{
//		for(int j = 0; j < 3; ++j)
//		{
//			m(i,j) = rot(j,i);
//		}
//	}

//	m(0,3) = tran(0);
//	m(1,3) = tran(1);
//	m(2,3) = tran(2);

//	return m;
//}


//SelfCollisionConstr::CollData::CollData(const rbd::MultiBody& mb, int collId,
//	int body1Id, sch::S_Object* body1, const sva::PTransformd& body1T,
//	int body2Id, sch::S_Object* body2, const sva::PTransformd& body2T,
//	double di, double ds, double damping, double dampOff):
//		pair(new sch::CD_Pair(body1, body2)),
//		body1T(body1T),
//		body2T(body2T),
//		normVecDist(Eigen::Vector3d::Zero()),
//		jacB1(rbd::Jacobian(mb, body1Id)),
//		jacB2(rbd::Jacobian(mb, body2Id)),
//		di(di),
//		ds(ds),
//		damping(damping),
//		collId(collId),
//		body1Id(body1Id),
//		body2Id(body2Id),
//		body1(mb.bodyIndexById(body1Id)),
//		body2(mb.bodyIndexById(body2Id)),
//		dampingType(damping > 0. ? DampingType::Hard : DampingType::Free),
//		dampingOff(dampOff)
//{
//}



//SelfCollisionConstr::SelfCollisionConstr(const rbd::MultiBody& mb, double step):
//	dataVec_(),
//	step_(step),
//	nrVars_(0),
//	nrActivated_(0),
//	AInEq_(),
//	bInEq_(),
//	fullJac_(3, mb.nrDof()),
//	calcVec_(mb.nrDof())
//{
//}


//void SelfCollisionConstr::addCollision(const rbd::MultiBody& mb, int collId,
//	int body1Id, sch::S_Object* body1, const sva::PTransformd& body1T,
//	int body2Id, sch::S_Object* body2, const sva::PTransformd& body2T,
//	double di, double ds, double damping, double dampingOff)
//{
//	dataVec_.emplace_back(mb, collId, body1Id, body1, body1T,
//		body2Id, body2, body2T, di, ds, damping, dampingOff);
//}


//bool SelfCollisionConstr::rmCollision(int collId)
//{
//	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
//		[collId](const CollData& data)
//		{
//			return data.collId == collId;
//		});

//	if(it != dataVec_.end())
//	{
//		delete it->pair;
//		dataVec_.erase(it);
//		return true;
//	}

//	return false;
//}


//std::size_t SelfCollisionConstr::nrCollisions() const
//{
//	return dataVec_.size();
//}


//void SelfCollisionConstr::reset()
//{
//	dataVec_.clear();
//}


//void SelfCollisionConstr::updateNrVars(const rbd::MultiBody& /* mb */,
//	const SolverData& data)
//{
//	nrVars_ = data.nrVars();
//}


//void SelfCollisionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
//	const SolverData& data)
//{
//	using namespace Eigen;

//	if(static_cast<unsigned int>(AInEq_.rows()) != dataVec_.size()
//		 || AInEq_.cols() != nrVars_)
//	{
//		AInEq_.setZero(dataVec_.size(), nrVars_);
//		bInEq_.setZero(dataVec_.size());
//	}

//	nrActivated_ = 0;
//	for(CollData& d: dataVec_)
//	{
//		sch::Point3 pb1Tmp, pb2Tmp;

//		d.pair->operator[](0)->setTransformation(tosch(d.body1T*mbc.bodyPosW[d.body1]));
//		d.pair->operator[](1)->setTransformation(tosch(d.body2T*mbc.bodyPosW[d.body2]));

//		double dist = d.pair->getClosestPoints(pb1Tmp, pb2Tmp);
//		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);

//		Vector3d pb1(pb1Tmp[0], pb1Tmp[1], pb1Tmp[2]);
//		Vector3d pb2(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]);

//		Eigen::Vector3d normVecDist = (pb1 - pb2)/dist;

//		pb1 = (sva::PTransformd(pb1)*mbc.bodyPosW[d.body1].inv()).translation();
//		pb2 = (sva::PTransformd(pb2)*mbc.bodyPosW[d.body2].inv()).translation();

//		if(dist < d.di)
//		{
//			if(d.dampingType == CollData::DampingType::Free)
//			{
//				d.dampingType = CollData::DampingType::Soft;
//				Vector3d v1(mbc.bodyPosW[d.body1].rotation().transpose()*
//						(sva::PTransformd(pb1)*mbc.bodyVelB[d.body1]).linear());
//				Vector3d v2(mbc.bodyPosW[d.body2].rotation().transpose()*
//						(sva::PTransformd(pb2)*mbc.bodyVelB[d.body2]).linear());
//				double distDot = std::abs((v1 - v2).dot(normVecDist));

//				/// @todo find a bette solution.
//				// use a value slightly upper ds if dist <= ds
//				double fixedDist = dist <= d.ds ? d.ds + (d.di - d.ds)*0.2 : dist;
//				d.damping = ((d.di - d.ds)/(fixedDist - d.ds))*distDot + d.dampingOff;
//			}

//			double dampers = d.damping*((dist - d.ds)/(d.di - d.ds));

//			Vector3d nf = normVecDist;
//			Vector3d onf = d.normVecDist;
//			Vector3d dnf = (nf - onf)/step_;

//			// Compute body1
//			d.jacB1.point(pb1);
//			const MatrixXd& jac1 = d.jacB1.jacobian(mb, mbc);
//			Eigen::Vector3d p1Speed = d.jacB1.velocity(mb, mbc).linear();
//			Eigen::Vector3d p1NormalAcc = d.jacB1.normalAcceleration(
//				mb, mbc, data.normalAccB()).linear();

//			d.jacB1.fullJacobian(mb, jac1.block(3, 0, 3, jac1.cols()), fullJac_);

//			double jqdn = (p1Speed.transpose()*nf)(0);
//			double jqdnd = (p1Speed.transpose()*dnf*step_)(0);
//			double jdqdn = (p1NormalAcc.transpose()*nf*step_)(0);

//			calcVec_.noalias() = -fullJac_.transpose()*(nf*step_);

//			// Compute body2
//			d.jacB2.point(pb2);
//			const MatrixXd& jac2 = d.jacB2.jacobian(mb, mbc);
//			Eigen::Vector3d p2Speed = d.jacB2.velocity(mb, mbc).linear();
//			Eigen::Vector3d p2NormalAcc = d.jacB2.normalAcceleration(
//				mb, mbc, data.normalAccB()).linear();

//			d.jacB2.fullJacobian(mb, jac2.block(3, 0, 3, jac2.cols()), fullJac_);

//			jqdn -= (p2Speed.transpose()*nf)(0);
//			jqdnd -= (p2Speed.transpose()*dnf*step_)(0);
//			jdqdn -= (p2NormalAcc.transpose()*nf*step_)(0);

//			calcVec_.noalias() += fullJac_.transpose()*(nf*step_);

//			// distdot + distdotdot*dt > -damp*((d - ds)/(di - ds))
//			AInEq_.block(nrActivated_, 0, 1, mb.nrDof()).noalias() = calcVec_.transpose();
//			bInEq_(nrActivated_) = dampers + jqdn + jqdnd + jdqdn;
//			++nrActivated_;
//		}
//		else
//		{
//			if(d.dampingType == CollData::DampingType::Soft)
//			{
//				d.dampingType = CollData::DampingType::Free;
//			}
//		}

//		d.normVecDist = normVecDist;
//	}
//}


//std::string SelfCollisionConstr::nameInEq() const
//{
//	return "SelfCollisionConstr";
//}


//std::string SelfCollisionConstr::descInEq(const rbd::MultiBody& mb, int line)
//{
//	int curLine = 0;
//	for(CollData& d: dataVec_)
//	{
//		double dist = d.pair->getDistance();
//		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
//		if(dist < d.di)
//		{
//			if(curLine == line)
//			{
//				std::stringstream ss;
//				ss << mb.body(d.body1).name() << " / " << mb.body(d.body2).name() << std::endl;
//				ss << "collId: " << d.collId << std::endl;
//				ss << "dist: " << dist << std::endl;
//				ss << "di: " << d.di << std::endl;
//				ss << "ds: " << d.ds << std::endl;
//				ss << "damp: " << d.damping + d.dampingOff << std::endl;
//				return ss.str();
//			}
//			++curLine;
//		}
//	}
//	return "";
//}


//int SelfCollisionConstr::nrInEq() const
//{
//	return nrActivated_;
//}


//int SelfCollisionConstr::maxInEq() const
//{
//	return int(dataVec_.size());
//}


//const Eigen::MatrixXd& SelfCollisionConstr::AInEq() const
//{
//	return AInEq_;
//}


//const Eigen::VectorXd& SelfCollisionConstr::bInEq() const
//{
//	return bInEq_;
//}


///**
//	*													StaticEnvCollisionConstr
//	*/


//StaticEnvCollisionConstr::CollData::CollData(const rbd::MultiBody& mb, int collId,
//	int bodyId, sch::S_Object* body, const sva::PTransformd& bodyT,
//	int envId, sch::S_Object* env,
//	double di, double ds, double damping, double dampOff):
//		pair(new sch::CD_Pair(body, env)),
//		bodyT(bodyT),
//		normVecDist(Eigen::Vector3d::Zero()),
//		jacB1(rbd::Jacobian(mb, bodyId)),
//		di(di),
//		ds(ds),
//		damping(damping),
//		collId(collId),
//		bodyId(bodyId),
//		envId(envId),
//		body(mb.bodyIndexById(bodyId)),
//		dampingType(damping > 0. ? DampingType::Hard : DampingType::Free),
//		dampingOff(dampOff)
//{
//}


//StaticEnvCollisionConstr::StaticEnvCollisionConstr(const rbd::MultiBody& mb, double step):
//	dataVec_(),
//	step_(step),
//	nrVars_(0),
//	nrActivated_(0),
//	AInEq_(),
//	bInEq_(),
//	fullJac_(3, mb.nrDof()),
//	calcVec_(mb.nrDof())
//{
//}


//void StaticEnvCollisionConstr::addCollision(const rbd::MultiBody& mb, int collId,
//	int bodyId, sch::S_Object* body, const sva::PTransformd& bodyT,
//	int envId, sch::S_Object* env,
//	double di, double ds, double damping, double dampingOff)
//{
//	dataVec_.emplace_back(mb, collId, bodyId, body, bodyT, envId, env,
//		di, ds, damping, dampingOff);
//}


//bool StaticEnvCollisionConstr::rmCollision(int collId)
//{
//	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
//		[collId](const CollData& data)
//		{
//			return data.collId == collId;
//		});

//	if(it != dataVec_.end())
//	{
//		delete it->pair;
//		dataVec_.erase(it);
//		return true;
//	}

//	return false;
//}


//std::size_t StaticEnvCollisionConstr::nrCollisions() const
//{
//	return dataVec_.size();
//}


//void StaticEnvCollisionConstr::reset()
//{
//	dataVec_.clear();
//}


//void StaticEnvCollisionConstr::updateNrVars(const rbd::MultiBody& /* mb */,
//	const SolverData& data)
//{
//	nrVars_ = data.nrVars();
//}


//void StaticEnvCollisionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
//	const SolverData& data)
//{
//	using namespace Eigen;

//	if(static_cast<unsigned int>(AInEq_.rows()) != dataVec_.size()
//		 || AInEq_.cols() != nrVars_)
//	{
//		AInEq_.setZero(dataVec_.size(), nrVars_);
//		bInEq_.setZero(dataVec_.size());
//	}

//	nrActivated_ = 0;
//	for(CollData& d: dataVec_)
//	{
//		sch::Point3 pb1Tmp, pb2Tmp;

//		d.pair->operator[](0)->setTransformation(tosch(d.bodyT*mbc.bodyPosW[d.body]));

//		double dist = d.pair->getClosestPoints(pb1Tmp, pb2Tmp);
//		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);

//		Vector3d pb1(pb1Tmp[0], pb1Tmp[1], pb1Tmp[2]);
//		Vector3d pb2(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]);

//		Eigen::Vector3d normVecDist = (pb1 - pb2)/dist;

//		pb1 = (sva::PTransformd(pb1)*mbc.bodyPosW[d.body].inv()).translation();

//		if(dist < d.di)
//		{
//			if(d.dampingType == CollData::DampingType::Free)
//			{
//				d.dampingType = CollData::DampingType::Soft;
//				Vector3d v1(mbc.bodyPosW[d.body].rotation().transpose()*
//						(sva::PTransformd(pb1)*mbc.bodyVelB[d.body]).linear());
//				double distDot = std::abs(v1.dot(normVecDist));

//				/// @todo find a bette solution.
//				// use a value slightly upper ds if dist <= ds
//				double fixedDist = dist <= d.ds ? d.ds + (d.di - d.ds)*0.2 : dist;
//				d.damping = ((d.di - d.ds)/(fixedDist - d.ds))*distDot + d.dampingOff;
//			}

//			double dampers = d.damping*((dist - d.ds)/(d.di - d.ds));

//			Vector3d nf = normVecDist;
//			Vector3d onf = d.normVecDist;
//			Vector3d dnf = (nf - onf)/step_;

//			// Compute body
//			d.jacB1.point(pb1);
//			const MatrixXd& jac1 = d.jacB1.jacobian(mb, mbc);

//			d.jacB1.fullJacobian(mb, jac1.block(3, 0, 3, jac1.cols()), fullJac_);
//			Eigen::Vector3d p1Speed = d.jacB1.velocity(mb, mbc).linear();
//			Eigen::Vector3d p1NormalAcc = d.jacB1.normalAcceleration(
//				mb, mbc, data.normalAccB()).linear();

//			double jqdn = (p1Speed.transpose()*nf)(0);
//			double jqdnd = (p1Speed.transpose()*dnf*step_)(0);
//			double jdqdn = (p1NormalAcc.transpose()*nf*step_)(0);

//			calcVec_.noalias() = -fullJac_.transpose()*(nf*step_);

//			// distdot + distdotdot*dt > -damp*((d - ds)/(di - ds))
//			AInEq_.block(nrActivated_, 0, 1, mb.nrDof()).noalias() = calcVec_.transpose();
//			bInEq_(nrActivated_) = dampers + jqdn + jqdnd + jdqdn;
//			++nrActivated_;
//		}
//		else
//		{
//			if(d.dampingType == CollData::DampingType::Soft)
//			{
//				d.dampingType = CollData::DampingType::Free;
//			}
//		}

//		d.normVecDist = normVecDist;
//	}
//}


//std::string StaticEnvCollisionConstr::nameInEq() const
//{
//	return "StaticEnvCollisionConstr";
//}


//std::string StaticEnvCollisionConstr::descInEq(const rbd::MultiBody& mb, int line)
//{
//	int curLine = 0;
//	for(CollData& d: dataVec_)
//	{
//		double dist = d.pair->getDistance();
//		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
//		if(dist < d.di)
//		{
//			if(curLine == line)
//			{
//				std::stringstream ss;
//				ss << mb.body(d.body).name() << " / " << d.envId << std::endl;
//				ss << "collId: " << d.collId << std::endl;
//				ss << "dist: " << dist << std::endl;
//				ss << "di: " << d.di << std::endl;
//				ss << "ds: " << d.ds << std::endl;
//				ss << "damp: " << d.damping + d.dampingOff << std::endl;
//				return ss.str();
//			}
//			++curLine;
//		}
//	}
//	return "";
//}


//int StaticEnvCollisionConstr::nrInEq() const
//{
//	return nrActivated_;
//}


//int StaticEnvCollisionConstr::maxInEq() const
//{
//	return int(dataVec_.size());
//}


//const Eigen::MatrixXd& StaticEnvCollisionConstr::AInEq() const
//{
//	return AInEq_;
//}


//const Eigen::VectorXd& StaticEnvCollisionConstr::bInEq() const
//{
//	return bInEq_;
//}


///**
//	*													CoMCollisionConstr
//	*/


//CoMCollisionConstr::CollData::CollData(const rbd::MultiBody& mb,
//	int collId, sch::S_Object* env,
//	double di, double ds, double damping, double dampOff):
//		comSphere_(ds/5.),
//		pair(new sch::CD_Pair(&comSphere_, env)),
//		normVecDist(Eigen::Vector3d::Zero()),
//		jacCoM(rbd::CoMJacobian(mb)),
//		di(di),
//		ds(ds),
//		damping(damping),
//		collId(collId),
//		dampingType(damping > 0. ? DampingType::Hard : DampingType::Free),
//		dampingOff(dampOff)
//{
//}


//CoMCollisionConstr::CoMCollisionConstr(const rbd::MultiBody& mb, double step):
//	dataVec_(),
//	step_(step),
//	nrVars_(0),
//	nrActivated_(0),
//	AInEq_(),
//	bInEq_(),
//	calcVec_(mb.nrDof())
//{
//}


//void CoMCollisionConstr::addCollision(const rbd::MultiBody& mb, int collId,
//	sch::S_Object* env,
//	double di, double ds, double damping, double dampingOff)
//{
//	dataVec_.emplace_back(mb, collId, env,
//		di, ds, damping, dampingOff);
//}


//bool CoMCollisionConstr::rmCollision(int collId)
//{
//	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
//		[collId](const CollData& data)
//		{
//			return data.collId == collId;
//		});

//	if(it != dataVec_.end())
//	{
//		delete it->pair;
//		dataVec_.erase(it);
//		return true;
//	}

//	return false;
//}


//std::size_t CoMCollisionConstr::nrCollisions() const
//{
//	return dataVec_.size();
//}


//void CoMCollisionConstr::reset()
//{
//	dataVec_.clear();
//}


//void CoMCollisionConstr::updateNrVars(const rbd::MultiBody& /* mb */,
//	const SolverData& data)
//{
//	nrVars_ = data.nrVars();
//}


//void CoMCollisionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
//	const SolverData& data)
//{
//	using namespace Eigen;

//	if(static_cast<unsigned int>(AInEq_.rows()) != dataVec_.size()
//		 || AInEq_.cols() != nrVars_)
//	{
//		AInEq_.setZero(dataVec_.size(), nrVars_);
//		bInEq_.setZero(dataVec_.size());
//	}

//	Eigen::Vector3d com = rbd::computeCoM(mb, mbc);

//	nrActivated_ = 0;
//	for(CollData& d: dataVec_)
//	{
//		sch::Point3 pb1Tmp, pb2Tmp;

//		d.pair->operator[](0)->setTransformation(tosch(sva::PTransformd(com)));

//		double dist = d.pair->getClosestPoints(pb1Tmp, pb2Tmp);
//		Vector3d pb2(pb2Tmp[0], pb2Tmp[1], pb2Tmp[2]);

//		dist = -std::copysign((com - pb2).norm(), dist);

//		Eigen::Vector3d normVecDist = (com - pb2)/dist;

//		if(dist < d.di)
//		{
//			// Compute CoM jacobian, speed and normalAcc
//			const MatrixXd& jac1 = d.jacCoM.jacobian(mb, mbc);
//			Eigen::Vector3d comSpeed = d.jacCoM.velocity(mb, mbc);
//			Eigen::Vector3d comNormalAcc = d.jacCoM.normalAcceleration(
//				mb, mbc, data.normalAccB());

//			if(d.dampingType == CollData::DampingType::Free)
//			{
//				d.dampingType = CollData::DampingType::Soft;
//				double distDot = std::abs(comSpeed.dot(normVecDist));

//				/// @todo find a bette solution.
//				// use a value slightly upper ds if dist <= ds
//				double fixedDist = dist <= d.ds ? d.ds + (d.di - d.ds)*0.2 : dist;
//				d.damping = ((d.di - d.ds)/(fixedDist - d.ds))*distDot + d.dampingOff;
//			}

//			double dampers = d.damping*((dist - d.ds)/(d.di - d.ds));

//			Vector3d nf = normVecDist;
//			Vector3d onf = d.normVecDist;
//			Vector3d dnf = (nf - onf)/step_;

//			double jqdn = (comSpeed.transpose()*nf)(0);
//			double jqdnd = (comSpeed.transpose()*dnf*step_)(0);
//			double jdqdn = (comNormalAcc.transpose()*nf*step_)(0);

//			calcVec_.noalias() = -jac1.transpose()*(nf*step_);

//			// distdot + distdotdot*dt > -damp*((d - ds)/(di - ds))
//			AInEq_.block(nrActivated_, 0, 1, mb.nrDof()) = calcVec_.transpose();
//			bInEq_(nrActivated_) = dampers + jqdn + jqdnd + jdqdn;
//			++nrActivated_;
//		}
//		else
//		{
//			if(d.dampingType == CollData::DampingType::Soft)
//			{
//				d.dampingType = CollData::DampingType::Free;
//			}
//		}

//		d.normVecDist = normVecDist;
//	}
//}


//std::string CoMCollisionConstr::nameInEq() const
//{
//	return "CoMCollisionConstr";
//}


//std::string CoMCollisionConstr::descInEq(const rbd::MultiBody& /*mb*/, int line)
//{
//	int curLine = 0;
//	for(CollData& d: dataVec_)
//	{
//		double dist = d.pair->getDistance();
//		dist = dist >= 0 ? -std::sqrt(dist) : std::sqrt(-dist);
//		if(dist < d.di)
//		{
//			if(curLine == line)
//			{
//				std::stringstream ss;
//				ss << "collId: " << d.collId << std::endl;
//				ss << "dist: " << dist << std::endl;
//				ss << "di: " << d.di << std::endl;
//				ss << "ds: " << d.ds << std::endl;
//				ss << "damp: " << d.damping + d.dampingOff << std::endl;
//				return ss.str();
//			}
//			++curLine;
//		}
//	}
//	return "";
//}


//int CoMCollisionConstr::nrInEq() const
//{
//	return nrActivated_;
//}


//int CoMCollisionConstr::maxInEq() const
//{
//	return int(dataVec_.size());
//}


//const Eigen::MatrixXd& CoMCollisionConstr::AInEq() const
//{
//	return AInEq_;
//}


//const Eigen::VectorXd& CoMCollisionConstr::bInEq() const
//{
//	return bInEq_;
//}


///**
//	*													CoMIncPlaneConstr
//	*/


//CoMIncPlaneConstr::PlaneData::PlaneData(
//	int planeId, const Eigen::Vector3d& normal, double offset,
//	double di, double ds, double damping, double dampOff):
//		normal(normal),
//		offset(offset),
//		dist(0.),
//		di(di),
//		ds(ds),
//		damping(damping),
//		planeId(planeId),
//		dampingType(damping > 0. ? DampingType::Hard : DampingType::Free),
//		dampingOff(dampOff)
//{
//}


//CoMIncPlaneConstr::CoMIncPlaneConstr(const rbd::MultiBody& mb, double step):
//	dataVec_(),
//	step_(step),
//	nrVars_(0),
//	nrActivated_(0),
//	activated_(0),
//	jacCoM_(mb),
//	AInEq_(),
//	bInEq_(),
//	calcVec_(mb.nrDof())
//{
//}


//void CoMIncPlaneConstr::addPlane(int planeId,
//	const Eigen::Vector3d& normal, double offset,
//	double di, double ds, double damping, double dampingOff)
//{
//	dataVec_.emplace_back(planeId, normal, offset,
//		di, ds, damping, dampingOff);
//	activated_.reserve(dataVec_.size());
//}


//bool CoMIncPlaneConstr::rmPlane(int planeId)
//{
//	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
//		[planeId](const PlaneData& data)
//		{
//			return data.planeId == planeId;
//		});

//	if(it != dataVec_.end())
//	{
//		dataVec_.erase(it);
//		return true;
//	}
//	// no need to resize activated, only the max size is important

//	return false;
//}


//std::size_t CoMIncPlaneConstr::nrPlanes() const
//{
//	return dataVec_.size();
//}


//void CoMIncPlaneConstr::reset()
//{
//	dataVec_.clear();
//}


//void CoMIncPlaneConstr::updateNrVars(const rbd::MultiBody& /* mb */,
//	const SolverData& data)
//{
//	nrVars_ = data.nrVars();
//}


//void CoMIncPlaneConstr::update(const rbd::MultiBody& mb,
//	const rbd::MultiBodyConfig& mbc, const SolverData& data)
//{
//	using namespace Eigen;

//	if(static_cast<unsigned int>(AInEq_.rows()) != dataVec_.size()
//		 || AInEq_.cols() != nrVars_)
//	{
//		AInEq_.setZero(dataVec_.size(), nrVars_);
//		bInEq_.setZero(dataVec_.size());
//	}

//	Eigen::Vector3d com = rbd::computeCoM(mb, mbc);

//	for(std::size_t i = 0; i < dataVec_.size(); ++i)
//	{
//		PlaneData& d = dataVec_[i];
//		d.dist = d.normal.dot(com) + d.offset;
//		if(d.dist <= d.di)
//		{
//			// don't allocate since we set capacity to dataVec_ size
//			activated_.push_back(i);
//		}
//		else
//		{
//			if(d.dampingType == PlaneData::DampingType::Soft)
//			{
//				d.dampingType = PlaneData::DampingType::Free;
//			}
//		}
//	}

//	nrActivated_ = 0;
//	if(!activated_.empty())
//	{
//		const MatrixXd& jacComMat = jacCoM_.jacobian(mb, mbc);
//		Eigen::Vector3d comSpeed = jacCoM_.velocity(mb, mbc);
//		Eigen::Vector3d comNormalAcc = jacCoM_.normalAcceleration(
//			mb, mbc, data.normalAccB());

//		for(std::size_t i: activated_)
//		{
//			PlaneData& d = dataVec_[i];
//			double distDot = d.normal.dot(comSpeed);

//			if(d.dampingType == PlaneData::DampingType::Free)
//			{
//				d.dampingType = PlaneData::DampingType::Soft;

//				/// @todo find a bette solution.
//				// use a value slightly upper ds if dist <= ds
//				double fixedDist = d.dist <= d.ds ? d.ds + (d.di - d.ds)*0.2 : d.dist;
//				d.damping = -((d.di - d.ds)/(fixedDist - d.ds))*distDot + d.dampingOff;
//			}

//			double dampers = d.damping*((d.dist - d.ds)/(d.di - d.ds));

//			// -dt*normal^T*J_com
//			AInEq_.block(nrActivated_, 0, 1, mb.nrDof()).noalias() =
//				-(step_*d.normal.transpose())*jacComMat;

//			// dampers + ddot + dt*normal^T*J*qdot
//			bInEq_(nrActivated_) = dampers + distDot + step_*(d.normal.dot(comNormalAcc));
//			++nrActivated_;
//		}
//	}

//	activated_.clear(); // don't free the vector, just say there is 0 elements
//}


//std::string CoMIncPlaneConstr::nameInEq() const
//{
//	return "CoMIncPlaneConstr";
//}


//std::string CoMIncPlaneConstr::descInEq(const rbd::MultiBody& /*mb*/, int line)
//{
//	int curLine = 0;
//	for(PlaneData& d: dataVec_)
//	{
//		if(d.dist < d.di)
//		{
//			if(curLine == line)
//			{
//				std::stringstream ss;
//				ss << "planeId: " << d.planeId << std::endl;
//				ss << "dist: " << d.dist << std::endl;
//				ss << "di: " << d.di << std::endl;
//				ss << "ds: " << d.ds << std::endl;
//				ss << "damp: " << d.damping << std::endl;
//				return ss.str();
//			}
//			++curLine;
//		}
//	}
//	return "";
//}


//int CoMIncPlaneConstr::nrInEq() const
//{
//	return nrActivated_;
//}


//int CoMIncPlaneConstr::maxInEq() const
//{
//	return int(dataVec_.size());
//}


//const Eigen::MatrixXd& CoMIncPlaneConstr::AInEq() const
//{
//	return AInEq_;
//}


//const Eigen::VectorXd& CoMIncPlaneConstr::bInEq() const
//{
//	return bInEq_;
//}


///**
//	*													GripperTorqueConstr
//	*/


//GripperTorqueConstr::GripperData::GripperData(int bId, double tl,
//	const Eigen::Vector3d& o, const Eigen::Vector3d& a):
//	bodyId(bId),
//	torqueLimit(tl),
//	origin(o),
//	axis(a)
//{}


//GripperTorqueConstr::GripperTorqueConstr():
//	dataVec_(),
//	AInEq_(),
//	bInEq_()
//{}


//void GripperTorqueConstr::addGripper(int bodyId, double torqueLimit,
//	const Eigen::Vector3d& origin, const Eigen::Vector3d& axis)
//{
//	dataVec_.emplace_back(bodyId, torqueLimit, origin, axis);
//}


//bool GripperTorqueConstr::rmGripper(int bodyId)
//{
//	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
//		[bodyId](const GripperData& data)
//		{
//			return data.bodyId == bodyId;
//		});

//	if(it != dataVec_.end())
//	{
//		dataVec_.erase(it);
//		return true;
//	}

//	return false;
//}


//void GripperTorqueConstr::reset()
//{
//	dataVec_.clear();
//}


//void GripperTorqueConstr::updateNrVars(const rbd::MultiBody& /* mb */,
//	const SolverData& data)
//{
//	using namespace Eigen;
//	AInEq_.setZero(dataVec_.size(), data.nrVars());
//	bInEq_.setZero(dataVec_.size());

//	int line = 0;
//	for(const GripperData& gd: dataVec_)
//	{
//		int begin = data.bilateralBegin();
//		for(const BilateralContact& bc: data.bilateralContacts())
//		{
//			int curLambda = 0;
//			// compute the number of lambda needed by the current bilateral
//			for(std::size_t i = 0; i < bc.points.size(); ++i)
//			{
//				curLambda += bc.nrLambda(static_cast<int>(i));
//			}

//			if(bc.bodyId == gd.bodyId)
//			{
//				int col = begin;
//				// Torque applied on the gripper motor
//				// Sum_i^nrF  T_i·( p_i^T_o x f_i)
//				for(std::size_t i = 0; i < bc.cones.size(); ++i)
//				{
//					Vector3d T_o_p = bc.points[i] - gd.origin;
//					for(std::size_t j = 0; j < bc.cones[i].generators.size(); ++j)
//					{
//						// we use abs because the contact force cannot apply
//						// negative torque on the gripper
//						AInEq_(line, col) = std::abs(
//							gd.axis.transpose()*(T_o_p.cross(bc.cones[i].generators[j])));
//						++col;
//					}
//				}
//				bInEq_(line) = gd.torqueLimit;
//				break;
//			}

//			begin += curLambda;
//			// if the bodyId is not found the AInEq_ and BInEq_ line stay at zero
//		}

//		++line;
//	}
//}


//void GripperTorqueConstr::update(const rbd::MultiBody& /* mb */,
//	const rbd::MultiBodyConfig& /* mbc */, const SolverData& /* data */)
//{}


//std::string GripperTorqueConstr::nameInEq() const
//{
//	return "GripperTorqueConstr";
//}


//std::string GripperTorqueConstr::descInEq(const rbd::MultiBody& mb, int line)
//{
//	std::stringstream ss;
//	int bodyIndex = mb.bodyIndexById(dataVec_[line].bodyId);
//	ss << mb.body(bodyIndex).name() << std::endl;
//	ss << "limits: " << dataVec_[line].torqueLimit << std::endl;
//	return ss.str();
//}


//int GripperTorqueConstr::maxInEq() const
//{
//	return static_cast<int>(dataVec_.size());
//}


//const Eigen::MatrixXd& GripperTorqueConstr::AInEq() const
//{
//	return AInEq_;
//}


//const Eigen::VectorXd& GripperTorqueConstr::bInEq() const
//{
//	return bInEq_;
//}


///**
//	*															ConstantSpeedConstr
//	*/


//ConstantSpeedConstr::ConstantSpeedConstr(const rbd::MultiBody& mb, double timeStep):
//	cont_(),
//	fullJac_(6, mb.nrDof()),
//	A_(),
//	b_(),
//	nrVars_(0),
//	timeStep_(timeStep)
//{}


//void ConstantSpeedConstr::addConstantSpeed(const rbd::MultiBody& mb, int bodyId,
//																				const Eigen::Vector3d& bodyPoint,
//																				const Eigen::MatrixXd& dof,
//																				const Eigen::VectorXd& speed)
//{
//	rbd::Jacobian jac(mb, bodyId, bodyPoint);
//	cont_.push_back({jac, dof, speed, bodyId});
//	updateNrEq();
//}


//bool ConstantSpeedConstr::removeConstantSpeed(int bodyId)
//{
//	auto it = std::find_if(cont_.begin(), cont_.end(),
//		[bodyId](const ConstantSpeedData& data)
//		{
//			return data.bodyId == bodyId;
//		});

//	if(it != cont_.end())
//	{
//		cont_.erase(it);
//		updateNrEq();
//		return true;
//	}

//	return false;
//}


//void ConstantSpeedConstr::resetConstantSpeed()
//{
//	cont_.clear();
//	updateNrEq();
//}


//std::size_t ConstantSpeedConstr::nrConstantSpeed() const
//{
//	return cont_.size();
//}


//void ConstantSpeedConstr::updateNrVars(const rbd::MultiBody& /* mb */,
//	const SolverData& data)
//{
//	nrVars_ = data.nrVars();
//	updateNrEq();
//}


//void ConstantSpeedConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
//	const SolverData& data)
//{
//	using namespace Eigen;

//	// TargetSpeed = V_k + A_{k+1}*dt
//	// TargetSpeed - V_k = J_k*alphaD_{k+1} + JD_k*alpha_k
//	// (TargetSpeed - V_k)/dt - JD_k*alpha_k = J_k*alphaD_{k+1}

//	int index = 0;
//	for(std::size_t i = 0; i < cont_.size(); ++i)
//	{
//		int rows = int(cont_[i].dof.rows());

//		// AEq
//		const MatrixXd& jac = cont_[i].jac.bodyJacobian(mb, mbc);
//		cont_[i].jac.fullJacobian(mb, jac, fullJac_);
//		A_.block(index, 0, rows, mb.nrDof()).noalias() = cont_[i].dof*fullJac_;

//		// BEq
//		Vector6d speed = cont_[i].jac.bodyVelocity(mb, mbc).vector();
//		Vector6d normalAcc = cont_[i].jac.bodyNormalAcceleration(
//			mb, mbc, data.normalAccB()).vector();
//		b_.segment(index, rows).noalias() = cont_[i].dof*(-normalAcc -
//			(speed/timeStep_));
//		b_.segment(index, rows).noalias() += (cont_[i].speed/timeStep_);
//		index += rows;
//	}
//}


//std::string ConstantSpeedConstr::nameEq() const
//{
//	return "ConstantSpeedConstr";
//}


//std::string ConstantSpeedConstr::descEq(const rbd::MultiBody& mb, int line)
//{
//	int curRow = 0;
//	for(const ConstantSpeedData& c: cont_)
//	{
//		curRow += int(c.dof.rows());
//		if(line < curRow)
//		{
//			return std::string("Body: ") + mb.body(c.body).name();
//		}
//	}
//	return std::string("");
//}


//int ConstantSpeedConstr::maxEq() const
//{
//	return int(A_.rows());
//}


//const Eigen::MatrixXd& ConstantSpeedConstr::AEq() const
//{
//	return A_;
//}


//const Eigen::VectorXd& ConstantSpeedConstr::bEq() const
//{
//	return b_;
//}


//void ConstantSpeedConstr::updateNrEq()
//{
//	int nrEq = 0;
//	for(const ConstantSpeedData& c: cont_)
//	{
//		nrEq += int(c.dof.rows());
//	}

//	A_.setZero(nrEq, nrVars_);
//	b_.setZero(nrEq);
//}



} // namespace qp

} // namespace tasks
