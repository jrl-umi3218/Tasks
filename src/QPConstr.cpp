// Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
//
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
#include "Tasks/QPConstr.h"

// includes
// std
#include <cmath>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/VisServo.h>

// sch
#include <sch/CD/CD_Pair.h>
#include <sch/S_Object/S_Object.h>

// Tasks
#include "Tasks/Bounds.h"
#include "utils.h"

namespace tasks
{

namespace qp
{

/**
	*															JointLimitsConstr
	*/


JointLimitsConstr::JointLimitsConstr(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, QBound bound, double step):
	robotIndex_(robotIndex),
	alphaDBegin_(-1),
	alphaDOffset_(mbs[robotIndex].joint(0).dof() > 1 ? mbs[robotIndex].joint(0).dof() : 0),
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
	int rI, const std::string& bName, sch::S_Object* h, const sva::PTransformd& X):
	hull(h),
	jac(mb, bName),
	X_op_o(X),
	rIndex(rI),
	bIndex(mb.bodyIndexByName(bName)),
	bodyName(bName)
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
	totalAlphaD_(-1),
	AInEq_(),
	bInEq_(),
	fullJac_(),
	distJac_()
{
	int maxDof = std::max_element(mbs.begin(), mbs.end(), compareDof)->nrDof();
	fullJac_.resize(1, maxDof);
	distJac_.resize(1, maxDof);
}


void CollisionConstr::addCollision(const std::vector<rbd::MultiBody>& mbs, int collId,
	int r1Index, const std::string& r1BodyName,
	sch::S_Object* body1, const sva::PTransformd& X_op1_o1,
	int r2Index, const std::string& r2BodyName,
	sch::S_Object* body2, const sva::PTransformd& X_op2_o2,
	double di, double ds, double damping, double dampingOff)
{
	const rbd::MultiBody mb1 = mbs[r1Index];
	const rbd::MultiBody mb2 = mbs[r2Index];
	std::vector<BodyCollData> bodies;
	if(mb1.nrDof() > 0)
	{
		bodies.emplace_back(mb1, r1Index, r1BodyName, body1, X_op1_o1);
	}
	if(mb2.nrDof() > 0)
	{
		bodies.emplace_back(mb2, r2Index, r2BodyName, body2, X_op2_o2);
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
	totalAlphaD_ = data.totalAlphaD();
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
			AInEq_.block(nrActivated_, 0, 1, totalAlphaD_).setZero();
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

				distJac_.block(0, 0, 1, bcd.jac.dof()).noalias() =
					(nf*step_*sign).transpose()*jac.block(3, 0, 3, bcd.jac.dof());

				bcd.jac.fullJacobian(mb, distJac_.block(0, 0, 1, bcd.jac.dof()), fullJac_);

				double jqdn = pSpeed.dot(nf);
				double jqdnd = pSpeed.dot(dnf*step_);
				double jdqdn = pNormalAcc.dot(nf*step_);

				AInEq_.block(nrActivated_, data.alphaDBegin(bcd.rIndex),
					1, mb.nrDof()).noalias() -= fullJac_.block(0, 0, 1, mb.nrDof());
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
	double di, double ds, double damping, double dampOff,
	const Eigen::Vector3d& speed,
	const Eigen::Vector3d& normalDot):
		normal(normal),
		normalDot(normalDot),
		offset(offset),
		dist(0.),
		di(di),
		ds(ds),
		damping(damping),
		planeId(planeId),
		dampingType(damping > 0. ? DampingType::Hard : DampingType::Free),
		dampingOff(dampOff),
		speed(speed)
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
		di, ds, damping, dampingOff, Eigen::Vector3d::Zero(),
		Eigen::Vector3d::Zero());
	activated_.reserve(dataVec_.size());
}

void CoMIncPlaneConstr::addPlane(int planeId,
	const Eigen::Vector3d& normal, double offset,
	double di, double ds, double damping,
	const Eigen::Vector3d& speed, const Eigen::Vector3d& normalDot, double dampingOff)
{
	dataVec_.emplace_back(planeId, normal, offset,
		di, ds, damping, dampingOff, speed, normalDot);
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
			double distDot = d.normal.dot(comSpeed-d.speed);
			double distDDot = d.normalDot.dot(comSpeed-d.speed);

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
			bInEq_(nrActivated_) = dampers + distDot + step_*(d.normal.dot(comNormalAcc)+distDDot);
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
				ss << "normal: " << d.normal.transpose();
				ss << "offset: " << d.offset << std::endl;
				ss << "dist: " << d.dist << std::endl;
				ss << "di: " << d.di << std::endl;
				ss << "ds: " << d.ds << std::endl;
				ss << "damp: " << d.damping << std::endl;
				ss << "speed: " << d.speed.transpose() << std::endl;
				ss << "normalDot: " << d.normalDot.transpose() << std::endl;
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


std::string GripperTorqueConstr::descInEq(const std::vector<rbd::MultiBody>& /* mbs */,
	int line)
{
	std::stringstream ss;
	const GripperData& gd = dataVec_[line];

	ss << gd.contactId.r1BodyName << "/" <<
				gd.contactId.r2BodyName << std::endl;
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
	*															BoundedSpeedConstr
	*/


BoundedSpeedConstr::BoundedSpeedConstr(const std::vector<rbd::MultiBody>& mbs,
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


void BoundedSpeedConstr::addBoundedSpeed(
	const std::vector<rbd::MultiBody>& mbs, const std::string& bodyName,
	const Eigen::Vector3d& bodyPoint, const Eigen::MatrixXd& dof,
	const Eigen::VectorXd& speed)
{
	addBoundedSpeed(mbs, bodyName, bodyPoint, dof, speed, speed);
}


void BoundedSpeedConstr::addBoundedSpeed(
	const std::vector<rbd::MultiBody>& mbs, const std::string& bodyName,
	const Eigen::Vector3d& bodyPoint, const Eigen::MatrixXd& dof,
	const Eigen::VectorXd& lowerSpeed, const Eigen::VectorXd& upperSpeed)
{
	rbd::Jacobian jac(mbs[robotIndex_], bodyName, bodyPoint);
	cont_.push_back({jac, dof, lowerSpeed, upperSpeed, bodyName});
}


bool BoundedSpeedConstr::removeBoundedSpeed(const std::string& bodyName)
{
	auto it = std::find_if(cont_.begin(), cont_.end(),
		[bodyName](const BoundedSpeedData& data)
		{
			return data.bodyName == bodyName;
		});

	if(it != cont_.end())
	{
		cont_.erase(it);
		return true;
	}

	return false;
}


void BoundedSpeedConstr::resetBoundedSpeeds()
{
	cont_.clear();
}


std::size_t BoundedSpeedConstr::nrBoundedSpeeds() const
{
	return cont_.size();
}


void BoundedSpeedConstr::updateBoundedSpeeds()
{
	updateNrEq();
}


void BoundedSpeedConstr::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
	nrVars_ = data.nrVars();
	updateNrEq();
}


void BoundedSpeedConstr::update(const std::vector<rbd::MultiBody>& mbs,
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


std::string BoundedSpeedConstr::nameGenInEq() const
{
	return "BoundedSpeedConstr";
}


std::string BoundedSpeedConstr::descGenInEq(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	int curRow = 0;
	for(const BoundedSpeedData& c: cont_)
	{
		curRow += int(c.dof.rows());
		if(line < curRow)
		{
			return std::string("Body: ") + mbs[robotIndex_].body(c.body).name();
		}
	}
	return std::string("");
}


int BoundedSpeedConstr::maxGenInEq() const
{
	return int(A_.rows());
}


const Eigen::MatrixXd& BoundedSpeedConstr::AGenInEq() const
{
	return A_;
}


const Eigen::VectorXd& BoundedSpeedConstr::LowerGenInEq() const
{
	return lower_;
}


const Eigen::VectorXd& BoundedSpeedConstr::UpperGenInEq() const
{
	return upper_;
}


void BoundedSpeedConstr::updateNrEq()
{
	int nrEq = 0;
	for(const BoundedSpeedData& c: cont_)
	{
		nrEq += int(c.dof.rows());
	}

	A_.setZero(nrEq, nrVars_);
	lower_.setZero(nrEq);
	upper_.setZero(nrEq);
}


/**
	*													ImageConstr
	*/


ImageConstr::PointData::PointData(const Eigen::Vector2d& pt, const double d):
	point2d(pt),
	depthEstimate(d)
{}


ImageConstr::RobotPointData::RobotPointData(const std::string& bn, const sva::PTransformd& X, const rbd::Jacobian& j):
	bName(bn),
	X_b_p(X),
	jac(j)
{}


ImageConstr::ImageConstr(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, const std::string& bName, const sva::PTransformd& X_b_gaze,
	double step, double constrDirection):
	dataVec_(),
	dataVecRob_(),
	robotIndex_(robotIndex),
	bodyIndex_(mbs[robotIndex].bodyIndexByName(bName)),
	alphaDBegin_(-1),
	nrVars_(0),
	step_(step),
	accelFactor_(0.5*step*step),
	nrActivated_(0),
	jac_(mbs[robotIndex], bName),
	X_b_gaze_(X_b_gaze),
	L_img_(new Eigen::Matrix<double, 2, 6>(Eigen::Matrix<double, 2, 6>::Zero())),
	surfaceVelocity_(new Eigen::Matrix<double,6,1>(Eigen::Matrix<double,6,1>::Zero())),
	L_Z_dot_(new Eigen::Matrix<double,1,6> (Eigen::Matrix<double,1,6> ::Zero())),
	L_img_dot_(new Eigen::Matrix<double,2,6>(Eigen::Matrix<double,2,6>::Zero())),
	speed_(new Eigen::Vector2d(Eigen::Vector2d::Zero())),
	normalAcc_(new Eigen::Vector2d(Eigen::Vector2d::Zero())),
	jacMat_(2, mbs[robotIndex].nrDof()),
	iDistMin_(new Eigen::Vector2d(Eigen::Vector2d::Zero())),
	iDistMax_(new Eigen::Vector2d(Eigen::Vector2d::Zero())),
	sDistMin_(new Eigen::Vector2d(Eigen::Vector2d::Zero())),
	sDistMax_(new Eigen::Vector2d(Eigen::Vector2d::Zero())),
	damping_(0.),
	dampingOffset_(0.),
	ineqInversion_(1),
	constrDirection_(constrDirection),
	AInEq_(),
	bInEq_()
{}

ImageConstr::ImageConstr(const ImageConstr& rhs):
	dataVec_(rhs.dataVec_),
	dataVecRob_(rhs.dataVecRob_),
	robotIndex_(rhs.robotIndex_),
	bodyIndex_(rhs.bodyIndex_),
	alphaDBegin_(rhs.alphaDBegin_),
	nrVars_(rhs.nrVars_),
	step_(rhs.step_),
	accelFactor_(rhs.accelFactor_),
	nrActivated_(rhs.nrActivated_),
	jac_(rhs.jac_),
	X_b_gaze_(rhs.X_b_gaze_),
	L_img_(new Eigen::Matrix<double, 2, 6>(*rhs.L_img_)),
	surfaceVelocity_(new Eigen::Matrix<double,6,1>(*rhs.surfaceVelocity_)),
	L_Z_dot_(new Eigen::Matrix<double,1,6>(*rhs.L_Z_dot_)),
	L_img_dot_(new Eigen::Matrix<double,2,6>(*rhs.L_img_dot_)),
	speed_(new Eigen::Vector2d(*rhs.speed_)),
	normalAcc_(new Eigen::Vector2d(*rhs.normalAcc_)),
	jacMat_(rhs.jacMat_),
	iDistMin_(new Eigen::Vector2d(*rhs.iDistMin_)),
	iDistMax_(new Eigen::Vector2d(*rhs.iDistMax_)),
	sDistMin_(new Eigen::Vector2d(*rhs.sDistMin_)),
	sDistMax_(new Eigen::Vector2d(*rhs.sDistMax_)),
	damping_(rhs.damping_),
	dampingOffset_(rhs.dampingOffset_),
	ineqInversion_(rhs.ineqInversion_),
	constrDirection_(rhs.constrDirection_),
	AInEq_(rhs.AInEq_),
	bInEq_(rhs.bInEq_)
{
}

ImageConstr& ImageConstr::operator=(const ImageConstr& rhs)
{
	if(&rhs != this)
	{
		dataVec_ = rhs.dataVec_;
		dataVecRob_ = rhs.dataVecRob_;
		robotIndex_ = rhs.robotIndex_;
		bodyIndex_ = rhs.bodyIndex_;
		alphaDBegin_ = rhs.alphaDBegin_;
		nrVars_ = rhs.nrVars_;
		step_ = rhs.step_;
		accelFactor_ = rhs.accelFactor_;
		nrActivated_ = rhs.nrActivated_;
		jac_ = rhs.jac_;
		X_b_gaze_ = rhs.X_b_gaze_;
		*L_img_ = *rhs.L_img_;
		*surfaceVelocity_ = *rhs.surfaceVelocity_;
		*L_Z_dot_ = *rhs.L_Z_dot_;
		*L_img_dot_ = *rhs.L_img_dot_;
		*speed_ = *rhs.speed_;
		*normalAcc_ = *rhs.normalAcc_;
		jacMat_ = rhs.jacMat_;
		*iDistMin_ = *rhs.iDistMin_;
		*iDistMax_ = *rhs.iDistMax_;
		*sDistMin_ = *rhs.sDistMin_;
		*sDistMax_ = *rhs.sDistMax_;
		damping_ = rhs.damping_;
		dampingOffset_ = rhs.dampingOffset_;
		ineqInversion_ = rhs.ineqInversion_;
		constrDirection_ = rhs.constrDirection_;
		AInEq_ = rhs.AInEq_;
		bInEq_ = rhs.bInEq_;
	}
	return *this;
}

int ImageConstr::addPoint(const Eigen::Vector2d& point2d, const double depthEstimate)
{
	dataVec_.emplace_back(point2d, depthEstimate);
	return int(dataVec_.size())-1;
}


int ImageConstr::addPoint(const Eigen::Vector3d& point3d)
{
	Eigen::Vector2d point2d(point3d[0]/point3d[2], point3d[1]/point3d[2]);
	int id = addPoint(point2d, point3d[2]);
	return id;
}


void ImageConstr::addPoint(const std::vector<rbd::MultiBody>& mbs,
	const std::string& bName, const sva::PTransformd& X_b_p)
{
	dataVecRob_.emplace_back(bName, X_b_p, rbd::Jacobian(mbs[robotIndex_], bName));
}


void ImageConstr::reset()
{
	dataVec_.clear();
	dataVecRob_.clear();
}


void ImageConstr::updatePoint(const int pointId, const Eigen::Vector2d& point2d)
{
	dataVec_[pointId].point2d = point2d;
}


void ImageConstr::updatePoint(const int pointId, const Eigen::Vector2d& point2d,
	const double depthEstimate)
{
	dataVec_[pointId] = PointData(point2d, depthEstimate);
}


void ImageConstr::updatePoint(const int pointId, const Eigen::Vector3d& point3d)
{
	Eigen::Vector2d	point2d_(point3d[0]/point3d[2], point3d[1]/point3d[2]);
	dataVec_[pointId] = PointData(point2d_, point3d[2]);
}


void ImageConstr::setLimits(const Eigen::Vector2d& min, const Eigen::Vector2d& max,
	const double iPercent, const double sPercent, const double damping, const double dampingOffsetPercent)
{
	// precompute avoidance boundaries
	Eigen::Vector2d dist = max - min;
	*iDistMin_ = min + constrDirection_*iPercent*dist;
	*iDistMax_ = max - constrDirection_*iPercent*dist;
	*sDistMin_ = min + constrDirection_*sPercent*dist;
	*sDistMax_ = max - constrDirection_*sPercent*dist;

	damping_ = damping;
	dampingOffset_ = constrDirection_*dampingOffsetPercent*damping;
}


void ImageConstr::computeComponents(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc, const SolverData& data,
	const Eigen::Vector2d& point2d, const double depth, rbd::Jacobian& jac, const int bodyIndex,
	const sva::PTransformd& X_b_p, Eigen::MatrixXd& fullJacobian, Eigen::Vector2d& bCommonTerm)
{
	// compute speed term
	rbd::imagePointJacobian(point2d, depth, *L_img_);
	*surfaceVelocity_ = (jac.velocity(mb, mbc, X_b_p)).vector();
	*speed_ = (*L_img_)*(*surfaceVelocity_);

	// compute norm accel term
	rbd::depthDotJacobian(*speed_, depth, *L_Z_dot_);
	rbd::imagePointJacobianDot(point2d, *speed_, depth, (*L_Z_dot_)*(*surfaceVelocity_), *L_img_dot_);
	*normalAcc_ = (*L_img_)*(jac.normalAcceleration(mb, mbc, data.normalAccB(robotIndex_), X_b_p,
		sva::MotionVecd(Eigen::Vector6d::Zero()))).vector() + (*L_img_dot_)*(*surfaceVelocity_);

	// compute the shortened jacobian
	const auto& shortJacMat = accelFactor_*(*L_img_)*jac.jacobian(mb, mbc, X_b_p*mbc.bodyPosW[bodyIndex]).block(0, 0, 6, jac.dof());

	// fill objects to return for the QP
	jac.fullJacobian(mb, shortJacMat, fullJacobian);
	bCommonTerm = -step_*(*speed_)-accelFactor_*(*normalAcc_);
}


void ImageConstr::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
	nrVars_ = data.nrVars();
	int nrRows = int(2*(dataVec_.size()+dataVecRob_.size()));
	AInEq_.setZero(nrRows, nrVars_);
	bInEq_.setZero(nrRows);
}


void ImageConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs, const SolverData& data)
{
	nrActivated_ = 0;

	const rbd::MultiBodyConfig& mbc = mbcs[robotIndex_];
	const rbd::MultiBody& mb = mbs[robotIndex_];

	// For each independent point
	for(const PointData& ptdata: dataVec_)
	{
		Eigen::Vector2d point2d_ = ptdata.point2d;

		Eigen::Vector2d bCommonTerm;
		computeComponents(mb, mbc, data, ptdata.point2d, ptdata.depthEstimate, jac_, bodyIndex_, X_b_gaze_, jacMat_, bCommonTerm);

		// For x and y
		for(std::size_t i=0; i<2; ++i)
		{
			bool isConstrActive = false;
			std::size_t iOther = (i==0) ? 1 : 0;

			// check occlusion constraint if it is within the 2D image bounds
			if((constrDirection_==1.) || ((point2d_[iOther] > (*iDistMin_)[iOther]) && (point2d_[iOther] < (*iDistMax_)[iOther])))
			{
				if((constrDirection_*point2d_[i] < constrDirection_*(*iDistMin_)[i]) //check min
						&& ((constrDirection_==1.) || (point2d_[i] < 0.))) //handle occlusion constraint ambiquity
				{
					if(constrDirection_*point2d_  [i] > constrDirection_*(*sDistMin_)[i])
					{
						isConstrActive = true;
						ineqInversion_ = -1.*constrDirection_;
						bInEq_(nrActivated_) = constrDirection_*(damping_*((point2d_[i] - (*sDistMin_)[i])/((*iDistMin_)[i] - (*sDistMin_)[i])) - dampingOffset_) + ineqInversion_*bCommonTerm[i];
					}
					else
					{
						//give maxmimum avoidance when in the safety limit
						bInEq_(nrActivated_) = - constrDirection_*dampingOffset_ + ineqInversion_*bCommonTerm[i];
					}
				}
				else if((constrDirection_*point2d_[i] > constrDirection_*(*iDistMax_)[i]) // check max
					&& ((constrDirection_==1.) || (point2d_[i] > 0.))) //handle occlusion constraint ambiquity
				{
					if(constrDirection_*point2d_[i] < constrDirection_*(*sDistMax_)[i])
					{
						isConstrActive = true;
						ineqInversion_ = 1.*constrDirection_;
						bInEq_(nrActivated_) = constrDirection_*(damping_*((point2d_[i] - (*sDistMax_)[i])/((*iDistMax_)[i] - (*sDistMax_)[i])) - dampingOffset_) + ineqInversion_*bCommonTerm[i];
					}
					else
					{
						//give maxmimum avoidance when in the safety limit
						bInEq_(nrActivated_) = - constrDirection_*dampingOffset_ + ineqInversion_*bCommonTerm[i];
					}
				}
			}

			// fill the needed QP constraint data
			if(isConstrActive)
			{
				// update the QP inequality constraint variables
				AInEq_.block(nrActivated_, alphaDBegin_, 1, mb.nrDof()) = ineqInversion_*jacMat_.block(i, 0, 1, mb.nrDof());
				++nrActivated_;
			}
		}
	}
}


std::string ImageConstr::nameInEq() const
{
	return "ImageConstr";
}


std::string ImageConstr::descInEq(const std::vector<rbd::MultiBody>& /* mbs */,
	int line)
{
	int curLine = 0;
	for(std::size_t i=0; i<dataVec_.size(); ++i)
	{
		// For x and y
		for(std::size_t j=0; j<2; ++j)
		{
			if(curLine == line)
			{
				std::stringstream ss;
				ss << "pointId: " << i << std::endl;
				if(j==0)
				{
					ss << "x" << std::endl;
				}
				else
				{
					ss << "y" << std::endl;
				}
				ss << "normalized 2d location: " << std::endl << dataVec_[i].point2d << std::endl;
				return ss.str();
			}
			++curLine;
		}
	}
	return "";
}


int ImageConstr::nrInEq() const
{
	return nrActivated_;
}


int ImageConstr::maxInEq() const
{
	return int(2*(dataVec_.size()+dataVecRob_.size()));
}


const Eigen::MatrixXd& ImageConstr::AInEq() const
{
	return AInEq_;
}


const Eigen::VectorXd& ImageConstr::bInEq() const
{
	return bInEq_;
}



} // namespace qp

} // namespace tasks
