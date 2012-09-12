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
#include "QPTasks.h"

// includes
// std
#include <cmath>
#include <iostream>

// Eigen
#include <Eigen/Geometry>

// RBDyn
#include <MultiBody.h>
#include <MultiBodyConfig.h>

namespace tasks
{

namespace qp
{


/**
	*														SetPointTask
	*/


SetPointTask::SetPointTask(const rbd::MultiBody& mb, HighLevelTask* hlTask,
  double stiffness, double weight):
  Task(weight),
  hlTask_(hlTask),
  stiffness_(stiffness),
  stiffnessSqrt_(2.*std::sqrt(stiffness)),
  Q_(mb.nrDof(), mb.nrDof()),
  C_(mb.nrDof()),
  alphaVec_(mb.nrDof())
{}


void SetPointTask::stiffness(double stiffness)
{
  stiffness_ = stiffness;
  stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}



void SetPointTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	hlTask_->update(mb, mbc);

	const Eigen::MatrixXd& J = hlTask_->jac();
	const Eigen::MatrixXd& JD = hlTask_->jacDot();
	const Eigen::VectorXd& err = hlTask_->eval();
	rbd::paramToVector(mbc.alpha, alphaVec_);

	Q_ = J.transpose()*J;
	C_ = -J.transpose()*(stiffness_*err - stiffnessSqrt_*J*alphaVec_ - JD*alphaVec_);
}


const Eigen::MatrixXd& SetPointTask::Q() const
{
	return Q_;
}


const Eigen::VectorXd& SetPointTask::C() const
{
	return C_;
}


/**
	*												PostureTask
	*/


PostureTask::PostureTask(const rbd::MultiBody& mb,
	std::vector<std::vector<double>> q,
	double stiffness, double weight):
	Task(weight),
	pt_(mb, q),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	Q_(mb.nrDof(), mb.nrDof()),
	C_(mb.nrDof()),
	alphaVec_(mb.nrDof())
{}


void PostureTask::stiffness(double stiffness)
{
	stiffness_ = stiffness;
	stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}


void PostureTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	pt_.update(mb, mbc);
	rbd::paramToVector(mbc.alpha, alphaVec_);

	Q_ = pt_.jac();
	C_.setZero();

	// joint
	C_.segment(mb.jointPosInDof(1), mb.nrDof()) = -stiffness_*pt_.eval() +
		stiffnessSqrt_*alphaVec_.segment(mb.jointPosInDof(1), mb.nrDof());
}

const Eigen::MatrixXd& PostureTask::Q() const
{
	return Q_;
}

const Eigen::VectorXd& PostureTask::C() const
{
	return C_;
}


/**
	*											PositionTask
	*/


PositionTask::PositionTask(const rbd::MultiBody& mb, int bodyId,
  const Eigen::Vector3d& pos, const Eigen::Vector3d& bodyPoint):
  pt_(mb, bodyId, pos, bodyPoint)
{
}


int PositionTask::dim()
{
	return 3;
}


void PositionTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	pt_.update(mb, mbc);
	pt_.updateDot(mb, mbc);
}

const Eigen::MatrixXd& PositionTask::jac()
{
	return pt_.jac();
}

const Eigen::MatrixXd& PositionTask::jacDot()
{
	return pt_.jacDot();
}

const Eigen::VectorXd& PositionTask::eval()
{
	return pt_.eval();
}


/**
	*																OrientationTask
	*/


OrientationTask::OrientationTask(const rbd::MultiBody& mb, int bodyId,
  const Eigen::Quaterniond& ori):
  ot_(mb, bodyId, ori)
{}


OrientationTask::OrientationTask(const rbd::MultiBody& mb, int bodyId,
  const Eigen::Matrix3d& ori):
  ot_(mb, bodyId, ori)
{}


int OrientationTask::dim()
{
	return 3;
}


void OrientationTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	ot_.update(mb, mbc);
	ot_.updateDot(mb, mbc);
}


const Eigen::MatrixXd& OrientationTask::jac()
{
	return ot_.jac();
}


const Eigen::MatrixXd& OrientationTask::jacDot()
{
	return ot_.jacDot();
}


const Eigen::VectorXd& OrientationTask::eval()
{
	return ot_.eval();
}


/**
	*													CoMTask
	*/


CoMTask::CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com):
	ct_(mb, com)
{}


int CoMTask::dim()
{
  return 3;
}


void CoMTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	ct_.update(mb, mbc);
	ct_.updateDot(mb, mbc);
}


const Eigen::MatrixXd& CoMTask::jac()
{
	return ct_.jac();
}


const Eigen::MatrixXd& CoMTask::jacDot()
{
	return ct_.jacDot();
}


const Eigen::VectorXd& CoMTask::eval()
{
	return ct_.eval();
}

} // namespace qp

} // namespace tasks

