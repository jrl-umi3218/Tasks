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
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

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
	dimWeight_(Eigen::VectorXd::Ones(hlTask->dim())),
	Q_(mb.nrDof(), mb.nrDof()),
	C_(mb.nrDof())
{}


SetPointTask::SetPointTask(const rbd::MultiBody& mb, HighLevelTask* hlTask,
	double stiffness, Eigen::VectorXd dimWeight, double weight):
	Task(weight),
	hlTask_(hlTask),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	dimWeight_(dimWeight),
	Q_(mb.nrDof(), mb.nrDof()),
	C_(mb.nrDof())
{}


void SetPointTask::stiffness(double stiffness)
{
	stiffness_ = stiffness;
	stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}


void SetPointTask::dimWeight(Eigen::VectorXd& dim)
{
	dimWeight_ = dim;
}


void SetPointTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	hlTask_->update(mb, mbc);

	const Eigen::MatrixXd& J = hlTask_->jac();
	const Eigen::VectorXd& err = hlTask_->eval();
	const Eigen::VectorXd& speed = hlTask_->speed();
	const Eigen::VectorXd& normalAcc = hlTask_->normalAcc();

	Q_.noalias() = J.transpose()*dimWeight_.asDiagonal()*J;
	C_.noalias() = -J.transpose()*dimWeight_.asDiagonal()*(stiffness_*err -
			stiffnessSqrt_*speed - normalAcc);
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
	*														PIDTask
	*/


PIDTask::PIDTask(const rbd::MultiBody& mb,
	HighLevelTask* hlTask,
	double P, double I, double D, double weight):
	Task(weight),
	hlTask_(hlTask),
	P_(P),
	I_(I),
	D_(D),
	dimWeight_(Eigen::VectorXd::Ones(hlTask->dim())),
	Q_(mb.nrDof(), mb.nrDof()),
	C_(mb.nrDof()),
	error_(hlTask->dim()),
	errorD_(hlTask->dim()),
	errorI_(hlTask->dim())
{}


PIDTask::PIDTask(const rbd::MultiBody& mb,
	HighLevelTask* hlTask,
	double P, double I, double D,
	Eigen::VectorXd dimWeight, double weight):
	Task(weight),
	hlTask_(hlTask),
	P_(P),
	I_(I),
	D_(D),
	dimWeight_(dimWeight),
	Q_(mb.nrDof(), mb.nrDof()),
	C_(mb.nrDof()),
	error_(hlTask->dim()),
	errorD_(hlTask->dim()),
	errorI_(hlTask->dim())
{}


double PIDTask::P() const
{
	return P_;
}


void PIDTask::P(double p)
{
	P_ = p;
}


double PIDTask::I() const
{
	return I_;
}


void PIDTask::I(double i)
{
	I_ = i;
}


double PIDTask::D() const
{
	return D_;
}


void PIDTask::D(double d)
{
	D_ = d;
}


void PIDTask::dimWeight(Eigen::VectorXd& dim)
{
	dimWeight_ = dim;
}


void PIDTask::error(const Eigen::VectorXd& err)
{
	error_ = err;
}


void PIDTask::errorD(const Eigen::VectorXd& errD)
{
	errorD_ = errD;
}


void PIDTask::errorI(const Eigen::VectorXd& errI)
{
	errorI_ = errI;
}


void PIDTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	hlTask_->update(mb, mbc);

	const Eigen::MatrixXd& J = hlTask_->jac();
	const Eigen::VectorXd& normalAcc = hlTask_->normalAcc();

	Q_.noalias() = J.transpose()*dimWeight_.asDiagonal()*J;
	C_.noalias() = -J.transpose()*dimWeight_.asDiagonal()*(P_*error_ -
			D_*errorD_ - I_*errorI_ - normalAcc);
}


const Eigen::MatrixXd& PIDTask::Q() const
{
	return Q_;
}


const Eigen::VectorXd& PIDTask::C() const
{
	return C_;
}


/**
	*														TargetObjectiveTask
	*/


TargetObjectiveTask::TargetObjectiveTask(const rbd::MultiBody& mb,
	HighLevelTask* hlTask, double timeStep, double dur,
	const Eigen::VectorXd& objDot, double weight):
	Task(weight),
	hlTask_(hlTask),
	dt_(timeStep),
	objDot_(objDot),
	dimWeight_(hlTask->dim()),
	phi_(hlTask->dim()),
	psi_(hlTask->dim()),
	Q_(mb.nrDof(), mb.nrDof()),
	C_(mb.nrDof())
{
	duration(dur);
	dimWeight_.setOnes();
}


TargetObjectiveTask::TargetObjectiveTask(const rbd::MultiBody& mb,
	HighLevelTask* hlTask,
	double timeStep, double dur, const Eigen::VectorXd& objDot,
	const Eigen::VectorXd& dimWeight, double weight):
	Task(weight),
	hlTask_(hlTask),
	dt_(timeStep),
	objDot_(objDot),
	dimWeight_(dimWeight),
	phi_(hlTask->dim()),
	psi_(hlTask->dim()),
	Q_(mb.nrDof(), mb.nrDof()),
	C_(mb.nrDof())
{
	duration(dur);
}


double TargetObjectiveTask::duration() const
{
	return (nrIter_ - iter_)*dt_;
}


void TargetObjectiveTask::duration(double d)
{
	nrIter_ = static_cast<int>(std::round(d/dt_));
	iter_ = 0;
}


void TargetObjectiveTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	hlTask_->update(mb, mbc);

	const MatrixXd& J = hlTask_->jac();
	const VectorXd& err = hlTask_->eval();
	const VectorXd& speed = hlTask_->speed();
	const VectorXd& normalAcc = hlTask_->normalAcc();

	// M·[phi, psi]^T = Obj

	// M =
	//  ⎡          2            2⎤
	//  ⎢(-t₀ + tf)   (-t₀ + tf) ⎥
	//  ⎢───────────  ───────────⎥
	//  ⎢     3            6     ⎥
	//  ⎢                        ⎥
	//  ⎢   t₀   tf      t₀   tf ⎥
	//  ⎢ - ── + ──    - ── + ── ⎥
	//  ⎣   2    2       2    2  ⎦

	// M^I =
	//  ⎡          6               2    ⎤
	//  ⎢ ───────────────────   ─────── ⎥
	//  ⎢   2               2   t₀ - tf ⎥
	//  ⎢ t₀  - 2⋅t₀⋅tf + tf            ⎥
	//  ⎢                               ⎥
	//  ⎢          6               4    ⎥
	//  ⎢─────────────────────  ────────⎥
	//  ⎢    2               2  -t₀ + tf⎥
	//  ⎣- t₀  + 2⋅t₀⋅tf - tf           ⎦

	// Obj = [ err - (tf - t₀)·J α, objDot - J·α ]

	double d = (nrIter_ - iter_)*dt_;
	double ds = std::pow(d, 2);

	Matrix2d MI;
	Vector2d Obj;

	MI << 6./ds, 2./(-d),
				6./(-ds), 4./d;

	for(int i = 0; i < hlTask_->dim(); ++i)
	{
		Obj << err(i) - d*speed(i),
					 objDot_(i) - speed(i);
		Vector2d pp(MI*Obj);
		phi_(i) = pp(0);
		psi_(i) = pp(1);
	}

	Q_ = (J.array().colwise()*dimWeight_.array()).matrix().transpose()*J;
	C_ = -(J.array().colwise()*dimWeight_.array()).matrix().transpose()*
		(phi_ - normalAcc);

	++iter_;
}


const Eigen::MatrixXd& TargetObjectiveTask::Q() const
{
	return Q_;
}


const Eigen::VectorXd& TargetObjectiveTask::C() const
{
	return C_;
}


/**
	*														QuadraticTask
	*/


QuadraticTask::QuadraticTask(const rbd::MultiBody& mb, HighLevelTask* hlTask,
  double weight):
  Task(weight),
  hlTask_(hlTask),
  Q_(mb.nrDof(), mb.nrDof()),
  C_(mb.nrDof())
{}


void QuadraticTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	hlTask_->update(mb, mbc);

	const Eigen::MatrixXd& J = hlTask_->jac();
	const Eigen::VectorXd& err = hlTask_->eval();
	const Eigen::VectorXd& normalAcc = hlTask_->normalAcc();

	Q_ = J.transpose()*J;
	C_ = -J.transpose()*(err - normalAcc);
}


const Eigen::MatrixXd& QuadraticTask::Q() const
{
	return Q_;
}


const Eigen::VectorXd& QuadraticTask::C() const
{
	return C_;
}


/**
	*												LinWeightTask
	*/


LinWeightTask::LinWeightTask(Task* t, double step, double objWeight):
  Task(0.),
  task_(t),
  step_(step),
  objWeight_(objWeight)
{
}


void LinWeightTask::weight(double w)
{
	objWeight_ = w;
}


std::pair<int, int> LinWeightTask::begin() const
{
	return task_->begin();
}


void LinWeightTask::updateNrVars(const rbd::MultiBody& mb,
	const SolverData& data)
{
	task_->updateNrVars(mb, data);
}


void LinWeightTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	double curW = Task::weight();
	if(objWeight_ > curW)
	{
		curW = std::min(objWeight_, curW + step_);
	}
	else
	{
		curW = std::max(objWeight_, curW - step_);
	}

	Task::weight(curW);

	task_->update(mb, mbc);
}


const Eigen::MatrixXd& LinWeightTask::Q() const
{
	return task_->Q();
}


const Eigen::VectorXd& LinWeightTask::C() const
{
	return task_->C();
}


/**
	*												PostureTask
	*/


PostureTask::PostureTask(const rbd::MultiBody& mb,
	std::vector<std::vector<double> > q,
	double stiffness, double weight):
	Task(weight),
	pt_(mb, q),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	jointDatas_(),
	Q_(mb.nrDof(), mb.nrDof()),
	C_(mb.nrDof()),
	alphaVec_(mb.nrDof())
{}


void PostureTask::stiffness(double stiffness)
{
	stiffness_ = stiffness;
	stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}


void PostureTask::jointsStiffness(const rbd::MultiBody& mb,
																const std::vector<JointStiffness>& jsv)
{
	jointDatas_.clear();
	jointDatas_.reserve(jsv.size());
	for(const JointStiffness& js: jsv)
	{
		int jointIndex = mb.jointIndexById(js.jointId);
		jointDatas_.push_back({js.stiffness, 2.*std::sqrt(js.stiffness),
													mb.jointPosInDof(jointIndex),
													mb.joint(jointIndex).dof()});
	}
}


void PostureTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	pt_.update(mb, mbc);
	rbd::paramToVector(mbc.alpha, alphaVec_);

	Q_ = pt_.jac();
	C_.setZero();

	int deb = mb.jointPosInDof(1);
	int end = mb.nrDof() - deb;
	// joint
	C_.segment(deb, end) = -stiffness_*pt_.eval().segment(deb, end) +
		stiffnessSqrt_*alphaVec_.segment(deb, end);

	for(const JointData& pjd: jointDatas_)
	{
		C_.segment(pjd.start, pjd.size) +=
				-pjd.stiffness*pt_.eval().segment(pjd.start, pjd.size) +
				pjd.stiffnessSqrt*alphaVec_.segment(pjd.start, pjd.size);
	}
}

const Eigen::MatrixXd& PostureTask::Q() const
{
	return Q_;
}

const Eigen::VectorXd& PostureTask::C() const
{
	return C_;
}

const Eigen::VectorXd& PostureTask::eval() const
{
	return pt_.eval();
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
}

const Eigen::MatrixXd& PositionTask::jac()
{
	return pt_.jac();
}

const Eigen::VectorXd& PositionTask::eval()
{
	return pt_.eval();
}

const Eigen::VectorXd& PositionTask::speed()
{
	return pt_.speed();
}

const Eigen::VectorXd& PositionTask::normalAcc()
{
	return pt_.normalAcc();
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
}


const Eigen::MatrixXd& OrientationTask::jac()
{
	return ot_.jac();
}


const Eigen::VectorXd& OrientationTask::eval()
{
	return ot_.eval();
}


const Eigen::VectorXd& OrientationTask::speed()
{
	return ot_.speed();
}


const Eigen::VectorXd& OrientationTask::normalAcc()
{
	return ot_.normalAcc();
}


/**
	*													CoMTask
	*/


CoMTask::CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com):
	ct_(mb, com)
{}


CoMTask::CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com,
							 std::vector<double> weight):
	ct_(mb, com, std::move(weight))
{}


int CoMTask::dim()
{
	return 3;
}


void CoMTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	ct_.update(mb, mbc);
}


const Eigen::MatrixXd& CoMTask::jac()
{
	return ct_.jac();
}


const Eigen::VectorXd& CoMTask::eval()
{
	return ct_.eval();
}


const Eigen::VectorXd& CoMTask::speed()
{
	return ct_.speed();
}


const Eigen::VectorXd& CoMTask::normalAcc()
{
	return ct_.normalAcc();
}


/**
	*														ContactTask
	*/


void ContactTask::error(const Eigen::Vector3d& error)
{
	error_ = error;
}


void ContactTask::errorD(const Eigen::Vector3d& errorD)
{
	errorD_ = errorD;
}


void ContactTask::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& data)
{
	int nrLambda = 0;
	begin_ = data.lambdaBegin();
	std::vector<FrictionCone> cones;

	for(const UnilateralContact& uc: data.unilateralContacts())
	{
		int curLambda = 0;
		for(std::size_t i = 0; i < uc.points.size(); ++i)
		{
			curLambda += uc.nrLambda(static_cast<int>(i));
		}

		if(uc.bodyId == bodyId_)
		{
			nrLambda = curLambda;
			cones = std::vector<FrictionCone>(uc.points.size(), uc.cone);
			break;
		}

		begin_ += curLambda;
	}

	// if body Id is not unilateral we search in bilateral
	if(nrLambda == 0)
	{
		for(const BilateralContact& uc: data.bilateralContacts())
		{
			int curLambda = 0;
			for(std::size_t i = 0; i < uc.points.size(); ++i)
			{
				curLambda += uc.nrLambda(static_cast<int>(i));
			}

			if(uc.bodyId == bodyId_)
			{
				nrLambda = curLambda;
				cones = uc.cones;
				break;
			}

			begin_ += curLambda;
		}
	}

	conesJac_.resize(3, nrLambda);
	int index = 0;
	for(const FrictionCone& fc: cones)
	{
		for(const Eigen::Vector3d& gen: fc.generators)
		{
			conesJac_.col(index) = gen;
			++index;
		}
	}

	Q_.resize(nrLambda, nrLambda);
	Q_.noalias() = conesJac_.transpose()*conesJac_;
	C_.setZero(nrLambda);
}


void ContactTask::update(const rbd::MultiBody& /* mb */,
	const rbd::MultiBodyConfig& /* mbc */)
{
	C_.noalias() = -conesJac_.transpose()*
			(stiffness_*error_ - stiffnessSqrt_*errorD_);
}


const Eigen::MatrixXd& ContactTask::Q() const
{
	return Q_;
}


const Eigen::VectorXd& ContactTask::C() const
{
	return C_;
}


/**
	*														GripperTorqueTask
	*/


void GripperTorqueTask::updateNrVars(const rbd::MultiBody& /* mb */,
	const SolverData& data)
{
	using namespace Eigen;
	bool found = false;

	begin_ = data.bilateralBegin();
	for(const BilateralContact& bc: data.bilateralContacts())
	{
		int curLambda = 0;
		// compute the number of lambda needed by the current bilateral
		for(std::size_t i = 0; i < bc.points.size(); ++i)
		{
			curLambda += bc.nrLambda(static_cast<int>(i));
		}

		if(bc.bodyId == bodyId_)
		{
			found = true;
			Q_.setZero(curLambda, curLambda);
			C_.resize(curLambda);

			int pos = 0;
			// minimize Torque applied on the gripper motor
			// min Sum_i^nrF  T_i·( p_i^T_o x f_i)
			for(std::size_t i = 0; i < bc.cones.size(); ++i)
			{
				Vector3d T_o_p = bc.points[i] - origin_;
				for(std::size_t j = 0; j < bc.cones[i].generators.size(); ++j)
				{
					// we use abs because the contact force cannot apply
					// negative torque on the gripper
					C_(pos) = std::abs(
						axis_.transpose()*(T_o_p.cross(bc.cones[i].generators[j])));
					++pos;
				}
			}
			break;
		}

		begin_ += curLambda;
	}

	// if no contact was found we don't activate the task
	// (safe position and empty matrix)
	if(!found)
	{
		begin_ = 0;
		Q_.resize(0, 0);
		C_.resize(0);
	}
}


void GripperTorqueTask::update(const rbd::MultiBody& /* mb */,
	const rbd::MultiBodyConfig& /* mbc */)
{ }


const Eigen::MatrixXd& GripperTorqueTask::Q() const
{
	return Q_;
}


const Eigen::VectorXd& GripperTorqueTask::C() const
{
	return C_;
}


/**
	*											LinVelocityTask
	*/


LinVelocityTask::LinVelocityTask(const rbd::MultiBody& mb, int bodyId,
  const Eigen::Vector3d& speed, const Eigen::Vector3d& bodyPoint):
  pt_(mb, bodyId, speed, bodyPoint)
{
}


int LinVelocityTask::dim()
{
	return 3;
}


void LinVelocityTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	pt_.update(mb, mbc);
}


const Eigen::MatrixXd& LinVelocityTask::jac()
{
	return pt_.jac();
}


const Eigen::VectorXd& LinVelocityTask::eval()
{
	return pt_.eval();
}


const Eigen::VectorXd& LinVelocityTask::speed()
{
	return pt_.speed();
}


const Eigen::VectorXd& LinVelocityTask::normalAcc()
{
	return pt_.normalAcc();
}


/**
	*											OrientationTrackingTask
	*/


OrientationTrackingTask::OrientationTrackingTask(const rbd::MultiBody& mb, int bodyId,
	const Eigen::Vector3d& bodyPoint, const Eigen::Vector3d& bodyAxis,
	const std::vector<int>& trackingJointsId,
	const Eigen::Vector3d& trackedPoint):
	ott_(mb, bodyId, bodyPoint, bodyAxis, trackingJointsId, trackedPoint),
	alphaVec_(mb.nrDof()),
	speed_(3),
	normalAcc_(3)
{}


int OrientationTrackingTask::dim()
{
	return 3;
}


void OrientationTrackingTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	ott_.update(mb, mbc);
	rbd::paramToVector(mbc.alpha, alphaVec_);

	speed_.noalias() = ott_.jac()*alphaVec_;
	normalAcc_.noalias() = ott_.jacDot()*alphaVec_;
}


const Eigen::MatrixXd& OrientationTrackingTask::jac()
{
	return ott_.jac();
}


const Eigen::VectorXd& OrientationTrackingTask::eval()
{
	return ott_.eval();
}


const Eigen::VectorXd& OrientationTrackingTask::speed()
{
	return speed_;
}


const Eigen::VectorXd& OrientationTrackingTask::normalAcc()
{
	return normalAcc_;
}

} // namespace qp

} // namespace tasks

