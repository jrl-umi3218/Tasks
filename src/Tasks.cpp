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
#include "Tasks.h"

// includes
// std
#include <set>

// rbd
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

namespace tasks
{


/**
	*													PositionTask
	*/


PositionTask::PositionTask(const rbd::MultiBody& mb, int bodyId,
	const Eigen::Vector3d& pos, const Eigen::Vector3d& bodyPoint):
	pos_(pos),
	point_(bodyPoint),
	bodyIndex_(mb.bodyIndexById(bodyId)),
	jac_(mb, bodyId, bodyPoint),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


void PositionTask::position(const Eigen::Vector3d& pos)
{
	pos_ = pos;
}


const Eigen::Vector3d& PositionTask::position() const
{
	return pos_;
}


void PositionTask::bodyPoint(const Eigen::Vector3d& point)
{
	jac_.point(point);
}


const Eigen::Vector3d& PositionTask::bodyPoint() const
{
	return jac_.point();
}


void PositionTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = pos_ - (point_*mbc.bodyPosW[bodyIndex_]).translation();
	speed_ = jac_.velocity(mb, mbc).linear();
	normalAcc_ = jac_.normalAcceleration(mb, mbc).linear();

	const auto& shortJacMat =
		jac_.jacobian(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void PositionTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = pos_ - (point_*mbc.bodyPosW[bodyIndex_]).translation();
	speed_ = jac_.velocity(mb, mbc).linear();
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB).linear();

	const auto& shortJacMat =
		jac_.jacobian(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void PositionTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	const auto& shortJacMat =
		jac_.jacobianDot(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacDotMat_);
}


const Eigen::VectorXd& PositionTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& PositionTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& PositionTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& PositionTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& PositionTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													OrientationTask
	*/


OrientationTask::OrientationTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Quaterniond& ori):
	ori_(ori.matrix()),
	bodyIndex_(mb.bodyIndexById(bodyId)),
	jac_(mb, bodyId),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


OrientationTask::OrientationTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Matrix3d& ori):
	ori_(ori),
	bodyIndex_(mb.bodyIndexById(bodyId)),
	jac_(mb, bodyId),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


void OrientationTask::orientation(const Eigen::Quaterniond& ori)
{
	ori_ = ori.matrix();
}


void OrientationTask::orientation(const Eigen::Matrix3d& ori)
{
	ori_ = ori;
}


const Eigen::Matrix3d& OrientationTask::orientation() const
{
	return ori_;
}


void OrientationTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = sva::rotationError(mbc.bodyPosW[bodyIndex_].rotation(), ori_, 1e-7);
	speed_ = jac_.velocity(mb, mbc).angular();
	normalAcc_ = jac_.normalAcceleration(mb, mbc).angular();

	const auto& shortJacMat = jac_.jacobian(mb, mbc).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void OrientationTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = sva::rotationError(mbc.bodyPosW[bodyIndex_].rotation(), ori_, 1e-7);
	speed_ = jac_.velocity(mb, mbc).angular();
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB).angular();

	const auto& shortJacMat = jac_.jacobian(mb, mbc).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void OrientationTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	const auto& shortJacMat = jac_.jacobianDot(mb, mbc).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacDotMat_);
}


const Eigen::VectorXd& OrientationTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& OrientationTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& OrientationTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& OrientationTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& OrientationTask::jacDot() const
{
	return jacDotMat_;
}


/**
  *													SurfaceOrientationTask
  */


SurfaceOrientationTask::SurfaceOrientationTask(const rbd::MultiBody& mb,
	int bodyId, const Eigen::Quaterniond& ori, const sva::PTransformd& X_b_s):
	ori_(ori.matrix()),
	bodyIndex_(mb.bodyIndexById(bodyId)),
	jac_(mb, bodyId),
	X_b_s_(X_b_s),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


SurfaceOrientationTask::SurfaceOrientationTask(const rbd::MultiBody& mb,
	int bodyId, const Eigen::Matrix3d& ori, const sva::PTransformd& X_b_s):
	ori_(ori),
	bodyIndex_(mb.bodyIndexById(bodyId)),
	jac_(mb, bodyId),
	X_b_s_(X_b_s),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


void SurfaceOrientationTask::orientation(const Eigen::Quaterniond& ori)
{
	ori_ = ori.matrix();
}


void SurfaceOrientationTask::orientation(const Eigen::Matrix3d& ori)
{
	ori_ = ori;
}


const Eigen::Matrix3d& SurfaceOrientationTask::orientation() const
{
	return ori_;
}


void SurfaceOrientationTask::update(const rbd::MultiBody& mb,
																	const rbd::MultiBodyConfig& mbc)
{
	eval_ = sva::rotationVelocity<double>
		(ori_*mbc.bodyPosW[bodyIndex_].rotation().transpose()*X_b_s_.rotation().transpose(), 1e-7);
	speed_ = jac_.velocity(mb, mbc, X_b_s_).angular();
	normalAcc_ = jac_.normalAcceleration(mb, mbc).angular();

	const auto& shortJacMat =
		jac_.jacobian(mb, mbc, X_b_s_*mbc.bodyPosW[bodyIndex_]).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void SurfaceOrientationTask::update(const rbd::MultiBody& mb,
																	const rbd::MultiBodyConfig& mbc,
																	const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = sva::rotationVelocity<double>
		(ori_*mbc.bodyPosW[bodyIndex_].rotation().transpose()*X_b_s_.rotation().transpose(), 1e-7);
	speed_ = jac_.velocity(mb, mbc, X_b_s_).angular();
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB).angular();

	const auto& shortJacMat =
		jac_.jacobian(mb, mbc, X_b_s_*mbc.bodyPosW[bodyIndex_]).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void SurfaceOrientationTask::updateDot(const rbd::MultiBody& mb,
																		 const rbd::MultiBodyConfig& mbc)
{
	const auto& shortJacMat = jac_.bodyJacobianDot(mb, mbc).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacDotMat_);
}


const Eigen::VectorXd& SurfaceOrientationTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& SurfaceOrientationTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& SurfaceOrientationTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& SurfaceOrientationTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& SurfaceOrientationTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													PostureTask
	*/


PostureTask::PostureTask(const rbd::MultiBody& mb, std::vector<std::vector<double> > q):
	q_(q),
	eval_(mb.nrDof()),
	jacMat_(mb.nrDof(), mb.nrDof()),
	jacDotMat_(mb.nrDof(), mb.nrDof())
{
	eval_.setZero();
	jacMat_.setIdentity();
	jacDotMat_.setZero();

	if(mb.nrDof() > 0 && mb.joint(0).type() == rbd::Joint::Free)
	{
		for(int i = 0; i < 6; ++i)
		{
			jacMat_(i, i) = 0;
		}
	}
}


void PostureTask::posture(std::vector<std::vector<double> > q)
{
	q_ = q;
}


const std::vector<std::vector<double> > PostureTask::posture() const
{
	return q_;
}


void PostureTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	int pos = mb.jointPosInDof(1);

	// we drop the first joint (fixed or free flyier).
	for(int i = 1; i < mb.nrJoints(); ++i)
	{
		// if dof == 1 is a prismatic/revolute joint
		// else if dof == 4 is a spherical one
		// else is a fixed one
		if(mb.joint(i).dof() == 1)
		{
			eval_(pos) = q_[i][0] - mbc.q[i][0];
			++pos;
		}
		else if(mb.joint(i).dof() == 4)
		{
			Matrix3d orid(
				Quaterniond(q_[i][0], q_[i][1], q_[i][2], q_[i][3]).matrix());

			Vector3d err = sva::rotationError(mbc.jointConfig[i].rotation(), orid);

			eval_.segment(pos, 3) = err;
			pos += 3;
		}
	}
}


void PostureTask::updateDot(const rbd::MultiBody& /* mb */, const rbd::MultiBodyConfig& /* mbc */)
{}


const Eigen::VectorXd& PostureTask::eval() const
{
	return eval_;
}


const Eigen::MatrixXd& PostureTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& PostureTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													CoMTask
	*/


CoMTask::CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com):
	com_(com),
	jac_(mb),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{}


CoMTask::CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com,
								std::vector<double> weight):
	com_(com),
	jac_(mb, std::move(weight)),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{}


void CoMTask::com(const Eigen::Vector3d& com)
{
	com_ = com;
}


const Eigen::Vector3d CoMTask::com() const
{
	return com_;
}


void CoMTask::updateInertialParameters(const rbd::MultiBody& mb)
{
	jac_.updateInertialParameters(mb);
}


void CoMTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = com_ - rbd::computeCoM(mb, mbc);

	speed_ = jac_.velocity(mb, mbc);
	normalAcc_ = jac_.normalAcceleration(mb, mbc);
	jacMat_ = jac_.jacobian(mb, mbc);
}


void CoMTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const Eigen::Vector3d& com,
	const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = com_ - com;

	speed_ = jac_.velocity(mb, mbc);
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB);
	jacMat_ = jac_.jacobian(mb, mbc);
}


void CoMTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	jacDotMat_ = jac_.jacobianDot(mb, mbc);
}


const Eigen::VectorXd& CoMTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& CoMTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& CoMTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& CoMTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& CoMTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													MultiCoMTask
	*/


MultiCoMTask::MultiCoMTask(const std::vector<rbd::MultiBody>& mbs,
	std::vector<int> robotIndexes, const Eigen::Vector3d& com):
	com_(com),
	robotIndexes_(std::move(robotIndexes)),
	robotsWeight_(),
	jac_(robotIndexes_.size()),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(robotIndexes_.size())
{
	computeRobotsWeight(mbs);

	// create CoMJacobian and jacobian matrix
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];

		jac_[i] = rbd::CoMJacobian(mb,
			std::vector<double>(mb.nrBodies(), robotsWeight_[i]));
		jacMat_[i].resize(3, mb.nrDof());
	}
}


void MultiCoMTask::com(const Eigen::Vector3d& com)
{
	com_ = com;
}


const Eigen::Vector3d MultiCoMTask::com() const
{
	return com_;
}


const std::vector<int>& MultiCoMTask::robotIndexes() const
{
	return robotIndexes_;
}


void MultiCoMTask::updateInertialParameters(const std::vector<rbd::MultiBody>& mbs)
{
	computeRobotsWeight(mbs);

	// upadte CoMJacobian per body weight
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];

		jac_[i].weight(mb, std::vector<double>(mb.nrBodies(), robotsWeight_[i]));
	}
}


void MultiCoMTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs)
{
	eval_ = com_;
	speed_.setZero();
	normalAcc_.setZero();
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];
		const rbd::MultiBodyConfig& mbc = mbcs[r];

		eval_ -=  rbd::computeCoM(mbs[r], mbcs[r])*robotsWeight_[i];
		speed_ += jac_[i].velocity(mb, mbc);
		normalAcc_ += jac_[i].normalAcceleration(mb, mbc);
		jacMat_[i] = jac_[i].jacobian(mb, mbc);
	}
}


void MultiCoMTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const std::vector<std::vector<sva::MotionVecd>>& normalAccB)
{
	eval_ = com_;
	speed_.setZero();
	normalAcc_.setZero();
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];
		const rbd::MultiBodyConfig& mbc = mbcs[r];

		eval_ -=  rbd::computeCoM(mbs[r], mbcs[r])*robotsWeight_[i];
		speed_ += jac_[i].velocity(mb, mbc);
		normalAcc_ += jac_[i].normalAcceleration(mb, mbc, normalAccB[r]);
		jacMat_[i] = jac_[i].jacobian(mb, mbc);
	}
}


void MultiCoMTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const std::vector<Eigen::Vector3d>& coms,
	const std::vector<std::vector<sva::MotionVecd>>& normalAccB)
{
	eval_ = com_;
	speed_.setZero();
	normalAcc_.setZero();
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];
		const rbd::MultiBodyConfig& mbc = mbcs[r];

		eval_ -= coms[r];
		speed_ += jac_[i].velocity(mb, mbc);
		normalAcc_ += jac_[i].normalAcceleration(mb, mbc, normalAccB[r]);
		jacMat_[i] = jac_[i].jacobian(mb, mbc);
	}
}


void MultiCoMTask::computeRobotsWeight(const std::vector<rbd::MultiBody>& mbs)
{
	double totalMass = 0.;

	robotsWeight_.clear();
	robotsWeight_.reserve(robotIndexes_.size());
	// compute the total mass and the weight of each robot
	for(int r: robotIndexes_)
	{
		double rm = std::accumulate(mbs[r].bodies().begin(),
			mbs[r].bodies().end(), 0., [](double ac, const rbd::Body& b)
				{return ac + b.inertia().mass();}
		);
		robotsWeight_.push_back(rm);
		totalMass += rm;
	}

	// divide all robotsWeight values by the total mass
	std::transform(robotsWeight_.begin(), robotsWeight_.end(),
		robotsWeight_.begin(),
		[totalMass](double w){return w/totalMass;});
}


const Eigen::VectorXd& MultiCoMTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& MultiCoMTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& MultiCoMTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& MultiCoMTask::jac(int index) const
{
	return jacMat_[index];
}


/**
	*													MomentumTask
	*/


MomentumTask::MomentumTask(const rbd::MultiBody& mb, const sva::ForceVecd mom):
	momentum_(mom),
	momentumMatrix_(mb),
	eval_(6),
	speed_(6),
	normalAcc_(6),
	jacMat_(6,mb.nrDof()),
	jacDotMat_(6,mb.nrDof())
{
	speed_.setZero();
}


void MomentumTask::momentum(const sva::ForceVecd& mom)
{
	momentum_ = mom;
}


const sva::ForceVecd MomentumTask::momentum() const
{
	return momentum_;
}


void MomentumTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	Eigen::Vector3d com = rbd::computeCoM(mb, mbc);

	eval_ = momentum_.vector() - rbd::computeCentroidalMomentum(mb, mbc, com).vector();
	normalAcc_ = momentumMatrix_.normalMomentumDot(mb, mbc,
		com,  rbd::computeCoMVelocity(mb, mbc)).vector();

	momentumMatrix_.computeMatrix(mb, mbc, com);
	jacMat_ = momentumMatrix_.matrix();
}


void MomentumTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	const std::vector<sva::MotionVecd>& normalAccB)
{
	Eigen::Vector3d com = rbd::computeCoM(mb, mbc);

	eval_ = momentum_.vector() - rbd::computeCentroidalMomentum(mb, mbc, com).vector();
	normalAcc_ = momentumMatrix_.normalMomentumDot(mb, mbc,
		com,  rbd::computeCoMVelocity(mb, mbc), normalAccB).vector();

	momentumMatrix_.computeMatrix(mb, mbc, com);
	jacMat_ = momentumMatrix_.matrix();
}


void MomentumTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	momentumMatrix_.computeMatrixDot(mb, mbc, rbd::computeCoM(mb, mbc),
				rbd::computeCoMVelocity(mb, mbc));
	jacDotMat_ = momentumMatrix_.matrixDot();
}


const Eigen::VectorXd& MomentumTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& MomentumTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& MomentumTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& MomentumTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& MomentumTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													LinVelocityTask
	*/


LinVelocityTask::LinVelocityTask(const rbd::MultiBody& mb, int bodyId,
	const Eigen::Vector3d& v, const Eigen::Vector3d& bodyPoint):
	vel_(v),
	point_(bodyPoint),
	bodyIndex_(mb.bodyIndexById(bodyId)),
	jac_(mb, bodyId, bodyPoint),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
	// this task don't have any derivative
	speed_.setZero();
}


void LinVelocityTask::velocity(const Eigen::Vector3d& v)
{
	vel_ = v;
}


const Eigen::Vector3d& LinVelocityTask::velocity() const
{
	return vel_;
}


void LinVelocityTask::bodyPoint(const Eigen::Vector3d& point)
{
	jac_.point(point);
}


const Eigen::Vector3d& LinVelocityTask::bodyPoint() const
{
	return jac_.point();
}


void LinVelocityTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = vel_ - jac_.velocity(mb, mbc).linear();
	normalAcc_ = jac_.normalAcceleration(mb, mbc).linear();

	const auto& shortJacMat = jac_.jacobian(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void LinVelocityTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = vel_ - jac_.velocity(mb, mbc).linear();
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB).linear();

	const auto& shortJacMat = jac_.jacobian(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void LinVelocityTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	const auto& shortJacMat = jac_.jacobianDot(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacDotMat_);
}


const Eigen::VectorXd& LinVelocityTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& LinVelocityTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& LinVelocityTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& LinVelocityTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& LinVelocityTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													OrientationTrackingTask
	*/


OrientationTrackingTask::OrientationTrackingTask(const rbd::MultiBody& mb,
	int bodyId, const Eigen::Vector3d& bodyPoint, const Eigen::Vector3d& bodyAxis,
	const std::vector<int>& trackingJointsId,
	const Eigen::Vector3d& trackedPoint):
	bodyIndex_(mb.bodyIndexById(bodyId)),
	bodyPoint_(bodyPoint),
	bodyAxis_(bodyAxis),
	zeroJacIndex_(),
	trackedPoint_(trackedPoint),
	jac_(mb, bodyId),
	eval_(3),
	shortJacMat_(3, jac_.dof()),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
	std::set<int> trackingJointsIndex;
	for(int id: trackingJointsId)
	{
		trackingJointsIndex.insert(mb.jointIndexById(id));
	}

	int jacPos = 0;
	for(int i: jac_.jointsPath())
	{
		const rbd::Joint& curJoint = mb.joint(i);
		if(trackingJointsIndex.find(i) == std::end(trackingJointsIndex))
		{
			for(int j = 0; j < curJoint.dof(); ++j)
			{
				zeroJacIndex_.push_back(jacPos + j);
			}
		}

		jacPos += curJoint.dof();
	}
}


void OrientationTrackingTask::trackedPoint(const Eigen::Vector3d& tp)
{
	trackedPoint_ = tp;
}


const Eigen::Vector3d& OrientationTrackingTask::trackedPoint() const
{
	return trackedPoint_;
}


void OrientationTrackingTask::bodyPoint(const Eigen::Vector3d& bp)
{
	bodyPoint_ = sva::PTransformd(bp);
}


const Eigen::Vector3d& OrientationTrackingTask::bodyPoint() const
{
	return bodyPoint_.translation();
}


void OrientationTrackingTask::bodyAxis(const Eigen::Vector3d& ba)
{
	bodyAxis_ = ba;
}


const Eigen::Vector3d& OrientationTrackingTask::bodyAxis() const
{
	return bodyAxis_;
}


void OrientationTrackingTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;
	const sva::PTransformd& bodyTf = mbc.bodyPosW[bodyIndex_];
	Vector3d desDir(trackedPoint_ - (bodyPoint_*bodyTf).translation());
	Vector3d curDir(bodyTf.rotation().transpose()*bodyAxis_);
	desDir.normalize();
	curDir.normalize();

	Matrix3d targetOri(
		Quaterniond::FromTwoVectors(curDir, desDir).inverse().matrix());

	eval_ = sva::rotationError<double>(mbc.bodyPosW[bodyIndex_].rotation(),
														targetOri*mbc.bodyPosW[bodyIndex_].rotation(), 1e-7);

	shortJacMat_ = jac_.jacobian(mb, mbc).block(0, 0, 3, shortJacMat_.cols());
	zeroJacobian(shortJacMat_);
	jac_.fullJacobian(mb, shortJacMat_, jacMat_);
}


void OrientationTrackingTask::updateDot(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc)
{
	shortJacMat_ = jac_.jacobianDot(mb, mbc).block(0, 0, 3, shortJacMat_.cols());
	zeroJacobian(shortJacMat_);
	jac_.fullJacobian(mb, shortJacMat_, jacDotMat_);
}


const Eigen::MatrixXd& OrientationTrackingTask::jac()
{
	return jacMat_;
}


const Eigen::MatrixXd& OrientationTrackingTask::jacDot()
{
	return jacDotMat_;
}


const Eigen::VectorXd& OrientationTrackingTask::eval()
{
	return eval_;
}


void OrientationTrackingTask::zeroJacobian(Eigen::MatrixXd& jac) const
{
	for(int i: zeroJacIndex_)
	{
		jac.col(i).setZero();
	}
}


} // namespace tasks
