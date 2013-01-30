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
  shortJacMat_(3, jac_.dof()),
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

	shortJacMat_ = jac_.jacobian(mb, mbc).block(3, 0, 3, shortJacMat_.cols());
	jac_.fullJacobian(mb, shortJacMat_, jacMat_);
}


void PositionTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	shortJacMat_ = jac_.jacobianDot(mb, mbc).block(3, 0, 3, shortJacMat_.cols());
	jac_.fullJacobian(mb, shortJacMat_, jacDotMat_);
}


const Eigen::VectorXd& PositionTask::eval() const
{
	return eval_;
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
  shortJacMat_(3, jac_.dof()),
  jacMat_(3, mb.nrDof()),
  jacDotMat_(3, mb.nrDof())
{
}


OrientationTask::OrientationTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Matrix3d& ori):
  ori_(ori),
  bodyIndex_(mb.bodyIndexById(bodyId)),
  jac_(mb, bodyId),
  eval_(3),
  shortJacMat_(3, jac_.dof()),
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

	shortJacMat_ = jac_.jacobian(mb, mbc).block(0, 0, 3, shortJacMat_.cols());
	jac_.fullJacobian(mb, shortJacMat_, jacMat_);
}


void OrientationTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	shortJacMat_ = jac_.jacobianDot(mb, mbc).block(0, 0, 3, shortJacMat_.cols());
	jac_.fullJacobian(mb, shortJacMat_, jacDotMat_);
}


const Eigen::VectorXd& OrientationTask::eval() const
{
	return eval_;
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
	for(std::size_t i = 1; i < mb.nrJoints(); ++i)
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


void CoMTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = com_ - rbd::computeCoM(mb, mbc);

	jacMat_ = jac_.jacobian(mb, mbc).block(3, 0, 3, mb.nrDof());
}


void CoMTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	jacDotMat_ = jac_.jacobianDot(mb, mbc).block(3, 0, 3, mb.nrDof());
}


const Eigen::VectorXd& CoMTask::eval() const
{
	return eval_;
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
	*													LinVelocityTask
	*/


LinVelocityTask::LinVelocityTask(const rbd::MultiBody& mb, int bodyId,
  const Eigen::Vector3d& v, const Eigen::Vector3d& bodyPoint):
  vel_(v),
  point_(bodyPoint),
  bodyIndex_(mb.bodyIndexById(bodyId)),
  jac_(mb, bodyId, bodyPoint),
  eval_(3),
  shortJacMat_(3, jac_.dof()),
  jacMat_(3, mb.nrDof()),
  jacDotMat_(3, mb.nrDof())
{
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
	sva::PTransform E_0_b(mbc.bodyPosW[bodyIndex_].rotation());
	eval_ = vel_ - (E_0_b.invMul(point_*mbc.bodyVelB[bodyIndex_])).linear();

	shortJacMat_ = jac_.jacobian(mb, mbc).block(3, 0, 3, shortJacMat_.cols());
	jac_.fullJacobian(mb, shortJacMat_, jacMat_);
}


void LinVelocityTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	shortJacMat_ = jac_.jacobianDot(mb, mbc).block(3, 0, 3, shortJacMat_.cols());
	jac_.fullJacobian(mb, shortJacMat_, jacDotMat_);
}


const Eigen::VectorXd& LinVelocityTask::eval() const
{
	return eval_;
}


const Eigen::MatrixXd& LinVelocityTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& LinVelocityTask::jacDot() const
{
	return jacDotMat_;
}

} // namespace tasks
