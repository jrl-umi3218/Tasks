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


/**
	*													CollisionConstr
	*/


sch::Matrix4x4 tosch(const sva::PTransformd& t)
{
	sch::Matrix4x4 m;
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



CollisionConstr::BodyCollData::BodyCollData(const rbd::MultiBody& mb,
	int rI, int bId, sch::S_Object* h, const sva::PTransformd& X):
	hull(h),
	jac(mb, bId),
	X_op_o(X),
	rIndex(rI),
	bIndex(mb.bodyIndexById(bId)),
	bodyId(bId)
{}



CollisionConstr::CollData::CollData(
		std::vector<BodyCollData> bcds, int collId,
		sch::S_Object* body1, sch::S_Object* body2,
		double di, double ds, double damp, double dampOff):
		pair(new sch::CD_Pair(body1, body2)),
		normVecDist(Eigen::Vector3d::Zero()),
		di(di),
		ds(ds),
		damping(damp),
		bodies(std::move(bcds)),
		dampingType(damping > 0. ? DampingType::Hard : DampingType::Free),
		dampingOff(dampOff),
		collId(collId)
{
}



CollisionConstr::CollisionConstr(const std::vector<rbd::MultiBody>& mbs, double step):
	dataVec_(),
	step_(step),
	nrActivated_(0),
	AInEq_(),
	bInEq_(),
	fullJac_(mbs.size())
{
	for(std::size_t i = 0; i < mbs.size(); ++i)
	{
		fullJac_[i].resize(3, mbs[i].nrDof());
	}
}


void CollisionConstr::addCollision(const std::vector<rbd::MultiBody>& mbs, int collId,
	int r1Index, int r1BodyId,
	sch::S_Object* body1, const sva::PTransformd& X_op1_o1,
	int r2Index, int r2BodyId,
	sch::S_Object* body2, const sva::PTransformd& X_op2_o2,
	double di, double ds, double damping, double dampingOff)
{
	const rbd::MultiBody mb1 = mbs[r1Index];
	const rbd::MultiBody mb2 = mbs[r2Index];
	std::vector<BodyCollData> bodies;
	if(mb1.nrDof() > 0)
	{
		bodies.emplace_back(mb1, r1Index, r1BodyId, body1, X_op1_o1);
	}
	if(mb2.nrDof() > 0)
	{
		bodies.emplace_back(mb2, r2Index, r2BodyId, body2, X_op2_o2);
	}

	dataVec_.emplace_back(std::move(bodies), collId, body1, body2,
		di, ds, damping, dampingOff);
}


bool CollisionConstr::rmCollision(int collId)
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


std::size_t CollisionConstr::nrCollisions() const
{
	return dataVec_.size();
}


void CollisionConstr::reset()
{
	dataVec_.clear();
}


void CollisionConstr::updateNrCollisions()
{
	AInEq_.setZero(dataVec_.size(), nrVars_);
	bInEq_.setZero(dataVec_.size());
}


void CollisionConstr::updateNrVars(const std::vector<rbd::MultiBody>& /* mb */,
	const SolverData& data)
{
	nrVars_ = data.nrVars();
	updateNrCollisions();
}


void CollisionConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	using namespace Eigen;

	Vector3d nearestPoint[2];

	nrActivated_ = 0;
	for(CollData& d: dataVec_)
	{
		// update moving hull position
		for(BodyCollData& bcd: d.bodies)
		{
			const rbd::MultiBodyConfig& mbc = mbcs[bcd.rIndex];
			bcd.hull->setTransformation(tosch(bcd.X_op_o*mbc.bodyPosW[bcd.bIndex]));
		}

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = d.pair->getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);

		nearestPoint[0] << pb1Tmp[0], pb1Tmp[1], pb1Tmp[2];
		nearestPoint[1] << pb2Tmp[0], pb2Tmp[1], pb2Tmp[2];

		Eigen::Vector3d normVecDist = (nearestPoint[0] - nearestPoint[1])/dist;

		// compute nearestPoint in body coordinate
		for(std::size_t i = 0; i < d.bodies.size(); ++i)
		{
			BodyCollData& bcd = d.bodies[i];
			const rbd::MultiBodyConfig& mbc = mbcs[bcd.rIndex];
			nearestPoint[i] = (sva::PTransformd(nearestPoint[i])*
				mbc.bodyPosW[bcd.bIndex].inv()).translation();

			// change the jacobian end point
			bcd.jac.point(nearestPoint[i]);
		}

		if(dist < d.di)
		{
			// automatic damping computation if needed
			if(d.dampingType == CollData::DampingType::Free)
			{
				d.dampingType = CollData::DampingType::Soft;
				d.damping = computeDamping(mbs, mbcs, d, normVecDist, dist);
			}

			double dampers = d.damping*((dist - d.ds)/(d.di - d.ds));

			Vector3d nf = normVecDist;
			Vector3d onf = d.normVecDist;
			Vector3d dnf = (nf - onf)/step_;

			double sign = 1.;
			bInEq_(nrActivated_) = dampers;
			AInEq_.row(nrActivated_).setZero();
			for(std::size_t i = 0; i < d.bodies.size(); ++i)
			{
				BodyCollData& bcd = d.bodies[i];
				const rbd::MultiBody& mb = mbs[bcd.rIndex];
				const rbd::MultiBodyConfig& mbc = mbcs[bcd.rIndex];

				// Compute body1
				const MatrixXd& jac = bcd.jac.jacobian(mb, mbc);
				Eigen::Vector3d pSpeed = bcd.jac.velocity(mb, mbc).linear();
				Eigen::Vector3d pNormalAcc = bcd.jac.normalAcceleration(
					mb, mbc, data.normalAccB(bcd.rIndex)).linear();

				bcd.jac.fullJacobian(mb, jac.block(3, 0, 3, jac.cols()),
					fullJac_[bcd.rIndex]);

				double jqdn = pSpeed.dot(nf);
				double jqdnd = pSpeed.dot(dnf*step_);
				double jdqdn = pNormalAcc.dot(nf*step_);

				AInEq_.block(nrActivated_, data.alphaDBegin(bcd.rIndex),
					1, mb.nrDof()).noalias() -=
						(nf*step_*sign).transpose()*fullJac_[bcd.rIndex];
				bInEq_(nrActivated_) += sign*(jqdn + jqdnd + jdqdn);
				// little hack
				// the max iteration number is two, so at the second iteration
				// sign will be -1
				sign = -1.;
			}
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


std::string CollisionConstr::nameInEq() const
{
	return "SelfCollisionConstr";
}


std::string CollisionConstr::descInEq(const std::vector<rbd::MultiBody>& mbs,
	int line)
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
				for(const BodyCollData& bcd: d.bodies)
				{
					const rbd::MultiBody& mb = mbs[bcd.rIndex];
					ss << "robot: " << bcd.rIndex << std::endl;
					ss << "body: " << mb.body(bcd.bIndex).name() << std::endl;
				}
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


int CollisionConstr::nrInEq() const
{
	return nrActivated_;
}


int CollisionConstr::maxInEq() const
{
	return int(dataVec_.size());
}


const Eigen::MatrixXd& CollisionConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& CollisionConstr::bInEq() const
{
	return bInEq_;
}


double CollisionConstr::computeDamping(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs, const CollData& cd,
	const Eigen::Vector3d& normVecDist, double dist) const
{
	Eigen::Vector3d diffVel(Eigen::Vector3d::Zero());
	double sign = 1.;
	for(std::size_t i = 0; i < cd.bodies.size(); ++i)
	{
		const BodyCollData& bcd = cd.bodies[i];
		const rbd::MultiBody& mb = mbs[bcd.rIndex];
		const rbd::MultiBodyConfig& mbc = mbcs[bcd.rIndex];

		Eigen::Vector3d velW = bcd.jac.velocity(mb, mbc).linear();

		diffVel += sign*velW;
		// little hack
		// the max iteration number is two, so at the second iteration
		// sign will be -1
		sign = -1;
	}

	double distDot = std::abs((diffVel).dot(normVecDist));

	/// @todo find a bette solution.
	// use a value slightly upper ds if dist <= ds
	double fixedDist = dist <= cd.ds ? cd.ds + (cd.di - cd.ds)*0.2 : dist;
	return ((cd.di - cd.ds)/(fixedDist - cd.ds))*distDot + cd.dampingOff;
}



/**
	*													CoMIncPlaneConstr
	*/


CoMIncPlaneConstr::PlaneData::PlaneData(
	int planeId, const Eigen::Vector3d& normal, double offset,
	double di, double ds, double damping, double dampOff):
		normal(normal),
		offset(offset),
		dist(0.),
		di(di),
		ds(ds),
		damping(damping),
		planeId(planeId),
		dampingType(damping > 0. ? DampingType::Hard : DampingType::Free),
		dampingOff(dampOff)
{
}


CoMIncPlaneConstr::CoMIncPlaneConstr(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, double step):
	robotIndex_(robotIndex),
	alphaDBegin_(-1),
	dataVec_(),
	step_(step),
	nrVars_(0),
	nrActivated_(0),
	activated_(0),
	jacCoM_(mbs[robotIndex]),
	AInEq_(),
	bInEq_()
{
}


void CoMIncPlaneConstr::addPlane(int planeId,
	const Eigen::Vector3d& normal, double offset,
	double di, double ds, double damping, double dampingOff)
{
	dataVec_.emplace_back(planeId, normal, offset,
		di, ds, damping, dampingOff);
	activated_.reserve(dataVec_.size());
}


bool CoMIncPlaneConstr::rmPlane(int planeId)
{
	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
		[planeId](const PlaneData& data)
		{
			return data.planeId == planeId;
		});

	if(it != dataVec_.end())
	{
		dataVec_.erase(it);
		return true;
	}
	// no need to resize activated, only the max size is important

	return false;
}


std::size_t CoMIncPlaneConstr::nrPlanes() const
{
	return dataVec_.size();
}


void CoMIncPlaneConstr::reset()
{
	dataVec_.clear();
}


void CoMIncPlaneConstr::updateNrPlanes()
{
	AInEq_.setZero(dataVec_.size(), nrVars_);
	bInEq_.setZero(dataVec_.size());
}


void CoMIncPlaneConstr::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
	nrVars_ = data.nrVars();
	updateNrPlanes();
}


void CoMIncPlaneConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs, const SolverData& data)
{
	using namespace Eigen;

	const rbd::MultiBody& mb = mbs[robotIndex_];
	const rbd::MultiBodyConfig& mbc = mbcs[robotIndex_];

	Eigen::Vector3d com = rbd::computeCoM(mb, mbc);

	for(std::size_t i = 0; i < dataVec_.size(); ++i)
	{
		PlaneData& d = dataVec_[i];
		d.dist = d.normal.dot(com) + d.offset;
		if(d.dist <= d.di)
		{
			// don't allocate since we set capacity to dataVec_ size
			activated_.push_back(i);
		}
		else
		{
			if(d.dampingType == PlaneData::DampingType::Soft)
			{
				d.dampingType = PlaneData::DampingType::Free;
			}
		}
	}

	nrActivated_ = 0;
	if(!activated_.empty())
	{
		const MatrixXd& jacComMat = jacCoM_.jacobian(mb, mbc);
		Eigen::Vector3d comSpeed = jacCoM_.velocity(mb, mbc);
		Eigen::Vector3d comNormalAcc = jacCoM_.normalAcceleration(
			mb, mbc, data.normalAccB(robotIndex_));

		for(std::size_t i: activated_)
		{
			PlaneData& d = dataVec_[i];
			double distDot = d.normal.dot(comSpeed);

			if(d.dampingType == PlaneData::DampingType::Free)
			{
				d.dampingType = PlaneData::DampingType::Soft;

				/// @todo find a bette solution.
				// use a value slightly upper ds if dist <= ds
				double fixedDist = d.dist <= d.ds ? d.ds + (d.di - d.ds)*0.2 : d.dist;
				d.damping = -((d.di - d.ds)/(fixedDist - d.ds))*distDot + d.dampingOff;
			}

			double dampers = d.damping*((d.dist - d.ds)/(d.di - d.ds));

			// -dt*normal^T*J_com
			AInEq_.block(nrActivated_, alphaDBegin_, 1, mb.nrDof()).noalias() =
				-(step_*d.normal.transpose())*jacComMat;

			// dampers + ddot + dt*normal^T*J*qdot
			bInEq_(nrActivated_) = dampers + distDot + step_*(d.normal.dot(comNormalAcc));
			++nrActivated_;
		}
	}

	activated_.clear(); // don't free the vector, just say there is 0 elements
}


std::string CoMIncPlaneConstr::nameInEq() const
{
	return "CoMIncPlaneConstr";
}


std::string CoMIncPlaneConstr::descInEq(const std::vector<rbd::MultiBody>& /* mbs */,
	int line)
{
	int curLine = 0;
	for(PlaneData& d: dataVec_)
	{
		if(d.dist < d.di)
		{
			if(curLine == line)
			{
				std::stringstream ss;
				ss << "planeId: " << d.planeId << std::endl;
				ss << "dist: " << d.dist << std::endl;
				ss << "di: " << d.di << std::endl;
				ss << "ds: " << d.ds << std::endl;
				ss << "damp: " << d.damping << std::endl;
				return ss.str();
			}
			++curLine;
		}
	}
	return "";
}


int CoMIncPlaneConstr::nrInEq() const
{
	return nrActivated_;
}


int CoMIncPlaneConstr::maxInEq() const
{
	return int(dataVec_.size());
}


const Eigen::MatrixXd& CoMIncPlaneConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& CoMIncPlaneConstr::bInEq() const
{
	return bInEq_;
}


/**
	*													GripperTorqueConstr
	*/


GripperTorqueConstr::GripperData::GripperData(const ContactId& cId, double tl,
	const Eigen::Vector3d& o, const Eigen::Vector3d& a):
	contactId(cId),
	torqueLimit(tl),
	origin(o),
	axis(a)
{}


GripperTorqueConstr::GripperTorqueConstr():
	dataVec_(),
	AInEq_(),
	bInEq_()
{}


void GripperTorqueConstr::addGripper(const ContactId& cId, double torqueLimit,
	const Eigen::Vector3d& origin, const Eigen::Vector3d& axis)
{
	dataVec_.emplace_back(cId, torqueLimit, origin, axis);
}


bool GripperTorqueConstr::rmGripper(const ContactId& contactId)
{
	auto it = std::find_if(dataVec_.begin(), dataVec_.end(),
		[contactId](const GripperData& data)
		{
			return data.contactId == contactId;
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


void GripperTorqueConstr::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	using namespace Eigen;
	AInEq_.setZero(dataVec_.size(), data.nrVars());
	bInEq_.setZero(dataVec_.size());

	int line = 0;
	int nrUni = int(data.unilateralContacts().size());
	for(const GripperData& gd: dataVec_)
	{
		for(std::size_t bi = 0; bi < data.bilateralContacts().size(); ++bi)
		{
			const BilateralContact& bc = data.bilateralContacts()[bi];

			if(bc.contactId == gd.contactId)
			{
				int col = data.lambdaBegin(int(bi) + nrUni);
				// Torque applied on the gripper motor
				// Sum_i^nrF  T_i·( p_i^T_o x f_i)
				for(std::size_t i = 0; i < bc.r1Cones.size(); ++i)
				{
					Vector3d T_o_p = bc.r1Points[i] - gd.origin;
					for(std::size_t j = 0; j < bc.r1Cones[i].generators.size(); ++j)
					{
						// we use abs because the contact force cannot apply
						// negative torque on the gripper
						AInEq_(line, col) = std::abs(
							gd.axis.transpose()*(T_o_p.cross(bc.r1Cones[i].generators[j])));
						++col;
					}
				}
				bInEq_(line) = gd.torqueLimit;
				++line;
				break;
			}
			// if the bodyId is not found the AInEq_ and BInEq_ line stay at zero
		}
	}
}


void GripperTorqueConstr::update(const std::vector<rbd::MultiBody>& /* mbs */,
	const std::vector<rbd::MultiBodyConfig>& /* mbcs */, const SolverData& /* data */)
{}


std::string GripperTorqueConstr::nameInEq() const
{
	return "GripperTorqueConstr";
}


std::string GripperTorqueConstr::descInEq(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	std::stringstream ss;
	const GripperData& gd = dataVec_[line];

	const rbd::MultiBody& mb1 = mbs[gd.contactId.r1Index];
	const rbd::MultiBody& mb2 = mbs[gd.contactId.r1Index];
	int r1BodyIndex = mb1.bodyIndexById(gd.contactId.r1BodyId);
	int r2BodyIndex = mb2.bodyIndexById(gd.contactId.r2BodyId);

	ss << mb1.body(r1BodyIndex).name() << "/" <<
				mb2.body(r2BodyIndex).name() << std::endl;
	ss << "limits: " << gd.torqueLimit << std::endl;
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


const Eigen::VectorXd& GripperTorqueConstr::bInEq() const
{
	return bInEq_;
}


/**
	*															ConstantSpeedConstr
	*/


ConstantSpeedConstr::ConstantSpeedConstr(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, double timeStep):
	robotIndex_(robotIndex),
	cont_(),
	fullJac_(6, mbs[robotIndex_].nrDof()),
	A_(),
	lower_(),
	upper_(),
	nrVars_(0),
	timeStep_(timeStep)
{}


void ConstantSpeedConstr::addConstantSpeed(
	const std::vector<rbd::MultiBody>& mbs, int bodyId,
	const Eigen::Vector3d& bodyPoint, const Eigen::MatrixXd& dof,
	const Eigen::VectorXd& speed)
{
	addConstantSpeed(mbs, bodyId, bodyPoint, dof, speed, speed);
}


void ConstantSpeedConstr::addConstantSpeed(
	const std::vector<rbd::MultiBody>& mbs, int bodyId,
	const Eigen::Vector3d& bodyPoint, const Eigen::MatrixXd& dof,
	const Eigen::VectorXd& lowerSpeed, const Eigen::VectorXd& upperSpeed)
{
	rbd::Jacobian jac(mbs[robotIndex_], bodyId, bodyPoint);
	cont_.push_back({jac, dof, lowerSpeed, upperSpeed, bodyId});
}


bool ConstantSpeedConstr::removeConstantSpeed(int bodyId)
{
	auto it = std::find_if(cont_.begin(), cont_.end(),
		[bodyId](const ConstantSpeedData& data)
		{
			return data.bodyId == bodyId;
		});

	if(it != cont_.end())
	{
		cont_.erase(it);
		return true;
	}

	return false;
}


void ConstantSpeedConstr::resetConstantSpeeds()
{
	cont_.clear();
}


std::size_t ConstantSpeedConstr::nrConstantSpeeds() const
{
	return cont_.size();
}


void ConstantSpeedConstr::updateConstantSpeeds()
{
	updateNrEq();
}


void ConstantSpeedConstr::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
	nrVars_ = data.nrVars();
	updateNrEq();
}


void ConstantSpeedConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	using namespace Eigen;

	const rbd::MultiBody& mb = mbs[robotIndex_];
	const rbd::MultiBodyConfig& mbc = mbcs[robotIndex_];

	// TargetSpeed = V_k + A_{k+1}*dt
	// TargetSpeed - V_k = J_k*alphaD_{k+1} + JD_k*alpha_k
	// (TargetSpeed - V_k)/dt - JD_k*alpha_k = J_k*alphaD_{k+1}

	int index = 0;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		int rows = int(cont_[i].dof.rows());

		// AEq
		const MatrixXd& jac = cont_[i].jac.bodyJacobian(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jac, fullJac_);
		A_.block(index, alphaDBegin_, rows, mb.nrDof()).noalias() =
			cont_[i].dof*fullJac_;

		// BEq
		Vector6d speed = cont_[i].jac.bodyVelocity(mb, mbc).vector();
		Vector6d normalAcc = cont_[i].jac.bodyNormalAcceleration(
			mb, mbc, data.normalAccB(robotIndex_)).vector();

		lower_.segment(index, rows).noalias() = cont_[i].dof*(-normalAcc -
			(speed/timeStep_));
		upper_.segment(index, rows).noalias() = lower_.segment(index, rows);

		lower_.segment(index, rows).noalias() += (cont_[i].lSpeed/timeStep_);
		upper_.segment(index, rows).noalias() += (cont_[i].uSpeed/timeStep_);
		index += rows;
	}
}


std::string ConstantSpeedConstr::nameGenInEq() const
{
	return "ConstantSpeedConstr";
}


std::string ConstantSpeedConstr::descGenInEq(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	int curRow = 0;
	for(const ConstantSpeedData& c: cont_)
	{
		curRow += int(c.dof.rows());
		if(line < curRow)
		{
			return std::string("Body: ") + mbs[robotIndex_].body(c.body).name();
		}
	}
	return std::string("");
}


int ConstantSpeedConstr::maxGenInEq() const
{
	return int(A_.rows());
}


const Eigen::MatrixXd& ConstantSpeedConstr::AGenInEq() const
{
	return A_;
}


const Eigen::VectorXd& ConstantSpeedConstr::LowerGenInEq() const
{
	return lower_;
}


const Eigen::VectorXd& ConstantSpeedConstr::UpperGenInEq() const
{
	return upper_;
}


void ConstantSpeedConstr::updateNrEq()
{
	int nrEq = 0;
	for(const ConstantSpeedData& c: cont_)
	{
		nrEq += int(c.dof.rows());
	}

	A_.setZero(nrEq, nrVars_);
	lower_.setZero(nrEq);
	upper_.setZero(nrEq);
}



} // namespace qp

} // namespace tasks
