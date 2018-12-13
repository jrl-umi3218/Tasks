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
#include "Tasks/QPTasks.h"

// includes
// std
#include <cmath>
#include <iterator>
#include <set>

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
	*														SetPointTaskCommon
	*/


SetPointTaskCommon::SetPointTaskCommon(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double weight):
	Task(weight),
	hlTask_(hlTask),
	error_(hlTask->dim()),
	dimWeight_(Eigen::VectorXd::Ones(hlTask->dim())),
	robotIndex_(rI),
	alphaDBegin_(0),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	preQ_(hlTask->dim(), mbs[rI].nrDof()),
	preC_(hlTask->dim())
{}


SetPointTaskCommon::SetPointTaskCommon(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	const Eigen::VectorXd& dimWeight, double weight):
	Task(weight),
	hlTask_(hlTask),
	error_(hlTask->dim()),
	dimWeight_(dimWeight),
	robotIndex_(rI),
	alphaDBegin_(0),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	preQ_(hlTask->dim(), mbs[rI].nrDof()),
	preC_(hlTask->dim())
{}


void SetPointTaskCommon::dimWeight(const Eigen::VectorXd& dim)
{
	dimWeight_ = dim;
}


void SetPointTaskCommon::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
}


void SetPointTaskCommon::computeQC(Eigen::VectorXd& error)
{
	const Eigen::MatrixXd& J = hlTask_->jac();
	const Eigen::VectorXd& normalAcc = hlTask_->normalAcc();

	error.noalias() -= normalAcc;
	preC_.noalias() = dimWeight_.asDiagonal()*error;
	C_.noalias() = -J.transpose()*preC_;

	preQ_.noalias() = dimWeight_.asDiagonal()*J;
	Q_.noalias() = J.transpose()*preQ_;
}


const Eigen::MatrixXd& SetPointTaskCommon::Q() const
{
	return Q_;
}


const Eigen::VectorXd& SetPointTaskCommon::C() const
{
	return C_;
}


/**
	*														SetPointTask
	*/


SetPointTask::SetPointTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double stiffness, double weight):
	SetPointTaskCommon(mbs, rI, hlTask, weight),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness))
{}


SetPointTask::SetPointTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double stiffness, const Eigen::VectorXd& dimWeight, double weight):
	SetPointTaskCommon(mbs, rI, hlTask, dimWeight, weight),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness))
{}


void SetPointTask::stiffness(double stiffness)
{
	stiffness_ = stiffness;
	stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}


void SetPointTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	hlTask_->update(mbs, mbcs, data);

	const Eigen::VectorXd& err = hlTask_->eval();
	const Eigen::VectorXd& speed = hlTask_->speed();

	error_.noalias() = stiffness_*err;
	error_.noalias() -= stiffnessSqrt_*speed;
	computeQC(error_);
}


/**
	*														TrackingTask
	*/


TrackingTask::TrackingTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double gainPos, double gainVel, double weight):
	SetPointTaskCommon(mbs, rI, hlTask, weight),
	gainPos_(gainPos),
	gainVel_(gainVel),
	errorPos_(Eigen::VectorXd::Zero(hlTask->dim())),
	errorVel_(Eigen::VectorXd::Zero(hlTask->dim())),
	refAccel_(Eigen::VectorXd::Zero(hlTask->dim()))
{}


TrackingTask::TrackingTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double gainPos, double gainVel, const Eigen::VectorXd& dimWeight, double weight):
	SetPointTaskCommon(mbs, rI, hlTask, dimWeight, weight),
	gainPos_(gainPos),
	gainVel_(gainVel),
	errorPos_(Eigen::VectorXd::Zero(hlTask->dim())),
	errorVel_(Eigen::VectorXd::Zero(hlTask->dim())),
	refAccel_(Eigen::VectorXd::Zero(hlTask->dim()))
{}


void TrackingTask::setGains(double gainPos, double gainVel)
{
	gainPos_ = gainPos;
	gainVel_ = gainVel;
}


void TrackingTask::errorPos(const Eigen::VectorXd& errorPos)
{
	errorPos_ = errorPos;
}


void TrackingTask::errorVel(const Eigen::VectorXd& errorVel)
{
	errorVel_ = errorVel;
}


void TrackingTask::refAccel(const Eigen::VectorXd& refAccel)
{
	refAccel_ = refAccel;
}


void TrackingTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	hlTask_->update(mbs, mbcs, data);

	error_.noalias() = gainPos_*errorPos_;
	error_.noalias() += gainVel_*errorVel_;
	error_.noalias() += refAccel_;
	computeQC(error_);
}


/**
	*														TrajectoryTask
	*/


TrajectoryTask::TrajectoryTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double gainPos, double gainVel, double weight):
	SetPointTaskCommon(mbs, rI, hlTask, weight),
	stiffness_(gainPos*Eigen::VectorXd::Ones(hlTask->dim())),
	damping_(gainVel*Eigen::VectorXd::Ones(hlTask->dim())),
	refVel_(Eigen::VectorXd::Zero(hlTask->dim())),
	refAccel_(Eigen::VectorXd::Zero(hlTask->dim()))
{}


TrajectoryTask::TrajectoryTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double gainPos, double gainVel, const Eigen::VectorXd& dimWeight, double weight):
	SetPointTaskCommon(mbs, rI, hlTask, dimWeight, weight),
	stiffness_(gainPos*Eigen::VectorXd::Ones(hlTask->dim())),
	damping_(gainVel*Eigen::VectorXd::Ones(hlTask->dim())),
	refVel_(Eigen::VectorXd::Zero(hlTask->dim())),
	refAccel_(Eigen::VectorXd::Zero(hlTask->dim()))
{}


void TrajectoryTask::setGains(double gainPos, double gainVel)
{
	stiffness_ = gainPos * Eigen::VectorXd::Ones(hlTask_->dim());
	damping_ = gainVel * Eigen::VectorXd::Ones(hlTask_->dim());
}

void TrajectoryTask::setGains(const Eigen::VectorXd & stiffness,
															const Eigen::VectorXd & damping)
{
	assert(stiffness.size() == stiffness_.size());
	assert(damping.size() == damping_.size());
	stiffness_ = stiffness;
	damping_ = damping;
}

void TrajectoryTask::stiffness(double gainPos)
{
	stiffness_ = gainPos * Eigen::VectorXd::Ones(hlTask_->dim());
}

void TrajectoryTask::stiffness(const Eigen::VectorXd & stiffness)
{
	assert(stiffness.size() == stiffness_.size());
	stiffness_ = stiffness;
}

const Eigen::VectorXd & TrajectoryTask::stiffness() const
{
	return stiffness_;
}

void TrajectoryTask::damping(double gainVel)
{
	damping_ = gainVel * Eigen::VectorXd::Ones(hlTask_->dim());
}

void TrajectoryTask::damping(const Eigen::VectorXd & damping)
{
	assert(damping.size() == damping_.size());
	damping_ = damping;
}


const Eigen::VectorXd & TrajectoryTask::damping() const
{
	return damping_;
}


void TrajectoryTask::refVel(const Eigen::VectorXd& refVel)
{
	refVel_ = refVel;
}

const Eigen::VectorXd & TrajectoryTask::refVel() const
{
	return refVel_;
}

void TrajectoryTask::refAccel(const Eigen::VectorXd& refAccel)
{
	refAccel_ = refAccel;
}

const Eigen::VectorXd & TrajectoryTask::refAccel() const
{
	return refAccel_;
}

void TrajectoryTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	hlTask_->update(mbs, mbcs, data);

	const Eigen::VectorXd& err = hlTask_->eval();
	const Eigen::VectorXd& speed = hlTask_->speed();

	error_.noalias() = stiffness_.asDiagonal()*err;
	error_.noalias() += damping_.asDiagonal()*(refVel_ - speed);
	error_.noalias() += refAccel_;
	computeQC(error_);
}


/**
	*														PIDTask
	*/


PIDTask::PIDTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double P, double I, double D, double weight):
	SetPointTaskCommon(mbs, rI, hlTask, weight),
	P_(P),
	I_(I),
	D_(D),
	error_(Eigen::VectorXd::Zero(hlTask->dim())),
	errorD_(Eigen::VectorXd::Zero(hlTask->dim())),
	errorI_(Eigen::VectorXd::Zero(hlTask->dim()))
{}


PIDTask::PIDTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double P, double I, double D,
	const Eigen::VectorXd& dimWeight, double weight):
	SetPointTaskCommon(mbs, rI, hlTask, dimWeight, weight),
	P_(P),
	I_(I),
	D_(D),
	error_(Eigen::VectorXd::Zero(hlTask->dim())),
	errorD_(Eigen::VectorXd::Zero(hlTask->dim())),
	errorI_(Eigen::VectorXd::Zero(hlTask->dim()))
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


void PIDTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	hlTask_->update(mbs, mbcs, data);

	error_.noalias() = P_*error_;
	error_.noalias() -= D_*errorD_;
	error_.noalias() -= I_*errorI_;
	computeQC(error_);
}


/**
	*														TargetObjectiveTask
	*/


TargetObjectiveTask::TargetObjectiveTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, HighLevelTask* hlTask, double timeStep, double dur,
	const Eigen::VectorXd& objDot, double weight):
	Task(weight),
	hlTask_(hlTask),
	dt_(timeStep),
	objDot_(objDot),
	dimWeight_(Eigen::VectorXd::Ones(hlTask->dim())),
	robotIndex_(rI),
	alphaDBegin_(0),
	phi_(hlTask->dim()),
	psi_(hlTask->dim()),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	preQ_(hlTask->dim(), mbs[rI].nrDof()),
	CVecSum_(hlTask->dim()),
	preC_(hlTask->dim())
{
	duration(dur);
}


TargetObjectiveTask::TargetObjectiveTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, HighLevelTask* hlTask,
	double timeStep, double dur, const Eigen::VectorXd& objDot,
	const Eigen::VectorXd& dimWeight, double weight):
	Task(weight),
	hlTask_(hlTask),
	dt_(timeStep),
	objDot_(objDot),
	dimWeight_(dimWeight),
	robotIndex_(rI),
	alphaDBegin_(0),
	phi_(hlTask->dim()),
	psi_(hlTask->dim()),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	preQ_(hlTask->dim(), mbs[rI].nrDof()),
	CVecSum_(hlTask->dim()),
	preC_(hlTask->dim())
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


void TargetObjectiveTask::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
}


void TargetObjectiveTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs, const SolverData& data)
{
	using namespace Eigen;

	hlTask_->update(mbs, mbcs, data);

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

	preQ_.noalias() = dimWeight_.asDiagonal()*J;
	Q_.noalias() = J.transpose()*preQ_;

	CVecSum_.noalias() = phi_ - normalAcc;
	preC_.noalias() = dimWeight_.asDiagonal()*CVecSum_;
	C_.noalias() = -J.transpose()*preC_;

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
	*												JointsSelector
	*/


JointsSelector JointsSelector::ActiveJoints(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, HighLevelTask* hl, const std::vector<std::string>& activeJointsName,
	const std::map<std::string, std::vector<std::array<int, 2>>>& activeDofs)
{
	return JointsSelector(mbs, robotIndex, hl, activeJointsName, activeDofs);
}


JointsSelector JointsSelector::UnactiveJoints(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, HighLevelTask* hl, const std::vector<std::string>& unactiveJointsName,
	const std::map<std::string, std::vector<std::array<int, 2>>>& unactiveDofs)
{
	using namespace std::placeholders;
	const rbd::MultiBody& mb = mbs[robotIndex];

	std::vector<std::string> activeJointsName;
	std::map<std::string, std::vector<std::array<int, 2>>> activeDofs;
	// sort unactiveJointsName by puting them into a set
	std::set<std::string> unactiveJointsNameSet(unactiveJointsName.begin(),
		unactiveJointsName.end());
	// create a set with all joints id
	std::set<std::string> jointsNameSet;
	std::transform(mb.joints().begin(), mb.joints().end(),
		std::inserter(jointsNameSet, jointsNameSet.begin()),
		std::bind(&rbd::Joint::name, _1));

	// remove unactive joints from the set
	std::set_difference(jointsNameSet.begin(), jointsNameSet.end(),
		unactiveJointsNameSet.begin(), unactiveJointsNameSet.end(),
		std::inserter(activeJointsName, activeJointsName.begin()));

	// Transform unactive dofs into active dofs
	for(const auto & uad : unactiveDofs)
	{
		const auto & jN = uad.first;
		const auto & jDofs = uad.second;
		const auto & j = mb.joint(mb.jointIndexByName(jN));
		int ji = 0;
		for(int i = 0; i < jDofs.size(); ++i)
		{
			auto dofStart = jDofs[i][0];
			if(dofStart > ji)
			{
				activeDofs[jN].push_back({{ji, dofStart - ji}});
			}
			ji = dofStart + jDofs[i][1];
		}
		auto finalDof = jDofs.back()[0] + jDofs.back()[1];
		if(finalDof < j.dof())
		{
			activeDofs[jN].push_back({{finalDof, j.dof() - finalDof}});
		}
	}

	return JointsSelector(mbs, robotIndex, hl, activeJointsName);
}


JointsSelector::JointsSelector(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
	HighLevelTask* hl, const std::vector<std::string>& selectedJointsName,
	const std::map<std::string, std::vector<std::array<int, 2>>>& activeDofs):
	jac_(Eigen::MatrixXd::Zero(hl->dim(), mbs[robotIndex].nrDof())),
	selectedJoints_(),
	hl_(hl)
{
	const rbd::MultiBody& mb = mbs[robotIndex];
	selectedJoints_.reserve(selectedJointsName.size());
	for(const std::string& jName: selectedJointsName)
	{
		int index = mb.jointIndexByName(jName);
		if(activeDofs.count(jName))
		{
			auto pInDof = mb.jointPosInDof(index);
			for(const auto & jdof : activeDofs.at(jName))
			{
				assert(jdof[0] + jdof[1] < mb.joint(index).dof());
				selectedJoints_.push_back({pInDof + jdof[0], jdof[1]});
			}
		}
		else
		{
			selectedJoints_.push_back({mb.jointPosInDof(index), mb.joint(index).dof()});
		}
	}
	// sort data in posInDof order
	std::sort(selectedJoints_.begin(), selectedJoints_.end(),
		[](const SelectedData& s1, const SelectedData& s2)
			{return s1.posInDof < s2.posInDof;});
}


int JointsSelector::dim()
{
	return hl_->dim();
}


void JointsSelector::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	hl_->update(mbs, mbcs, data);
	const Eigen::MatrixXd& jac = hl_->jac();
	for(SelectedData sd: selectedJoints_)
	{
		jac_.block(0, sd.posInDof, jac_.rows(), sd.dof) =
			jac.block(0, sd.posInDof, jac_.rows(), sd.dof);
	}
}


const Eigen::MatrixXd& JointsSelector::jac()
{
	return jac_;
}


const Eigen::VectorXd& JointsSelector::eval()
{
	return hl_->eval();
}


const Eigen::VectorXd& JointsSelector::speed()
{
	return hl_->speed();
}


const Eigen::VectorXd& JointsSelector::normalAcc()
{
	return hl_->normalAcc();
}

/** Torque Task **/
TorqueTask::TorqueTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                       const TorqueBound& tb, double weight):
  Task(weight),
  robotIndex_(robotIndex),
  alphaDBegin_(-1),
  lambdaBegin_(-1),
  motionConstr(mbs, robotIndex, tb),
  jointSelector_(mbs[robotIndex].nrDof()),
  Q_(mbs[robotIndex].nrDof(), mbs[robotIndex].nrDof()),
  C_(mbs[robotIndex].nrDof())
{
  jointSelector_.setOnes();
}

TorqueTask::TorqueTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                       const TorqueBound& tb, const Eigen::VectorXd& jointSelect,
                       double weight):
  Task(weight),
  robotIndex_(robotIndex),
  alphaDBegin_(-1),
  lambdaBegin_(-1),
  motionConstr(mbs, robotIndex, tb),
  jointSelector_(jointSelect),
  Q_(mbs[robotIndex].nrDof(), mbs[robotIndex].nrDof()),
  C_(mbs[robotIndex].nrDof())
{
}

TorqueTask::TorqueTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                       const TorqueBound& tb, const std::string& efName,
                       double weight):
  Task(weight),
  robotIndex_(robotIndex),
  alphaDBegin_(-1),
  lambdaBegin_(-1),
  motionConstr(mbs, robotIndex, tb),
  jointSelector_(mbs[robotIndex].nrDof()),
  Q_(mbs[robotIndex].nrDof(), mbs[robotIndex].nrDof()),
  C_(mbs[robotIndex].nrDof())
{
  rbd::Jacobian jac(mbs[robotIndex], efName);
  jointSelector_.setZero();
  for(auto i : jac.jointsPath())
  {
    //Do not add root joint !
    if(i != 0)
    {
    jointSelector_.segment(mbs[robotIndex].jointPosInDof(i),
                           mbs[robotIndex].joint(i).dof()).setOnes();
    }
  }
}

void TorqueTask::updateNrVars(const std::vector<rbd::MultiBody>& mbs,
                              const SolverData& data)
{
  motionConstr.updateNrVars(mbs, data);
  alphaDBegin_ = data.alphaDBegin(robotIndex_);
  lambdaBegin_ = data.lambdaBegin();
  Q_.resize(data.nrVars(), data.nrVars());
  C_.resize(data.nrVars());
}

void TorqueTask::update(const std::vector<rbd::MultiBody>& mbs,
                        const std::vector<rbd::MultiBodyConfig>& mbcs,
                        const SolverData& data)
{
  motionConstr.update(mbs, mbcs, data);
  Q_.noalias() = motionConstr.matrix().transpose()*jointSelector_.asDiagonal()*motionConstr.matrix();
  C_.noalias() = motionConstr.fd().C().transpose()*jointSelector_.asDiagonal()*motionConstr.matrix();
  //C_.setZero();
}

/**
	*												PostureTask
	*/


PostureTask::PostureTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	std::vector<std::vector<double> > q,
	double stiffness, double weight):
	Task(weight),
	pt_(mbs[rI], q),
	stiffness_(stiffness),
	damping_(2.*std::sqrt(stiffness)),
	robotIndex_(rI),
	alphaDBegin_(0),
	jointDatas_(),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	alphaVec_(mbs[rI].nrDof())
{}


void PostureTask::stiffness(double stiffness)
{
	stiffness_ = stiffness;
	damping_ = 2.*std::sqrt(stiffness);
}

void PostureTask::gains(double stiffness)
{
	stiffness_ = stiffness;
	damping_ = 2.*std::sqrt(stiffness);
}

void PostureTask::gains(double stiffness, double damping)
{
	stiffness_ = stiffness;
	damping_ = damping;
}

void PostureTask::jointsStiffness(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<JointStiffness>& jsv)
{
	jointDatas_.clear();
	jointDatas_.reserve(jsv.size());

	const rbd::MultiBody& mb = mbs[robotIndex_];
	for(const JointStiffness& js: jsv)
	{
		int jointIndex = mb.jointIndexByName(js.jointName);
		jointDatas_.push_back({js.stiffness, 2.*std::sqrt(js.stiffness),
													mb.jointPosInDof(jointIndex),
													mb.joint(jointIndex).dof()});
	}
}

void PostureTask::jointsGains(const std::vector<rbd::MultiBody> &mbs,
	const std::vector<JointGains> &jgv)
{
	jointDatas_.clear();
	jointDatas_.reserve(jgv.size());

	const rbd::MultiBody& mb = mbs[robotIndex_];
	for(const JointGains& jg: jgv)
	{
		int jointIndex = mb.jointIndexByName(jg.jointName);
		jointDatas_.push_back({jg.stiffness, jg.damping, mb.jointPosInDof(jointIndex),
			mb.joint(jointIndex).dof()});
	}
}

void PostureTask::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
}


void PostureTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& /* data */)
{
	const rbd::MultiBody& mb = mbs[robotIndex_];
	const rbd::MultiBodyConfig& mbc = mbcs[robotIndex_];

	pt_.update(mb, mbc);
	rbd::paramToVector(mbc.alpha, alphaVec_);

	Q_ = pt_.jac();
	C_.setZero();

	int deb = mb.jointPosInDof(1);
	int end = mb.nrDof() - deb;
	// joint
	C_.segment(deb, end) = -stiffness_*pt_.eval().segment(deb, end) +
		damping_*alphaVec_.segment(deb, end);

	for(const JointData& pjd: jointDatas_)
	{
		C_.segment(pjd.start, pjd.size) =
				-pjd.stiffness*pt_.eval().segment(pjd.start, pjd.size) +
				pjd.damping*alphaVec_.segment(pjd.start, pjd.size);
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


PositionTask::PositionTask(const std::vector<rbd::MultiBody>& mbs, int rI,
	const std::string& bodyName, const Eigen::Vector3d& pos,
	const Eigen::Vector3d& bodyPoint):
	pt_(mbs[rI], bodyName, pos, bodyPoint),
	robotIndex_(rI)
{
}


int PositionTask::dim()
{
	return 3;
}


void PositionTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	pt_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
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


OrientationTask::OrientationTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, const std::string& bodyName,
	const Eigen::Quaterniond& ori):
	ot_(mbs[rI], bodyName, ori),
	robotIndex_(rI)
{}


OrientationTask::OrientationTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, const std::string& bodyName,
	const Eigen::Matrix3d& ori):
	ot_(mbs[rI], bodyName, ori),
	robotIndex_(rI)
{}


int OrientationTask::dim()
{
	return 3;
}


void OrientationTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	ot_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
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
	*											SurfaceTransformTask
	*/


SurfaceTransformTask::SurfaceTransformTask(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex,
	const std::string& bodyName, const sva::PTransformd& X_0_t,
	const sva::PTransformd& X_b_p):
	TransformTaskCommon(mbs, robotIndex, bodyName, X_0_t, X_b_p)
{
}


void SurfaceTransformTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	tt_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
}


/**
	*											TransformTask
	*/


TransformTask::TransformTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
	const std::string& bodyName, const sva::PTransformd& X_0_t,
	const sva::PTransformd& X_b_p, const Eigen::Matrix3d& E_0_c):
	TransformTaskCommon(mbs, robotIndex, bodyName, X_0_t, X_b_p)
{
	tt_.E_0_c(E_0_c);
}


void TransformTask::E_0_c(const Eigen::Matrix3d& E_0_c)
{
	tt_.E_0_c(E_0_c);
}


const Eigen::Matrix3d& TransformTask::E_0_c() const
{
	return tt_.E_0_c();
}


void TransformTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	tt_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
}


/**
	*																SurfaceOrientationTask
	*/


SurfaceOrientationTask::SurfaceOrientationTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, const std::string& bodyName,
	const Eigen::Quaterniond& ori, const sva::PTransformd& X_b_s):
	ot_(mbs[rI], bodyName, ori, X_b_s),
	robotIndex_(rI)
{}


SurfaceOrientationTask::SurfaceOrientationTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, const std::string& bodyName,
	const Eigen::Matrix3d& ori, const sva::PTransformd& X_b_s):
	ot_(mbs[rI], bodyName, ori, X_b_s),
	robotIndex_(rI)
{}


int SurfaceOrientationTask::dim()
{
	return 3;
}


void SurfaceOrientationTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	ot_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
}


const Eigen::MatrixXd& SurfaceOrientationTask::jac()
{
	return ot_.jac();
}


const Eigen::VectorXd& SurfaceOrientationTask::eval()
{
	return ot_.eval();
}


const Eigen::VectorXd& SurfaceOrientationTask::speed()
{
	return ot_.speed();
}


const Eigen::VectorXd& SurfaceOrientationTask::normalAcc()
{
	return ot_.normalAcc();
}


/**
	*																GazeTask
	*/


GazeTask::GazeTask(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, const std::string& bodyName,
	const Eigen::Vector2d& point2d, double depthEstimate,
	const sva::PTransformd& X_b_gaze,
	const Eigen::Vector2d& point2d_ref):
	gazet_(mbs[robotIndex], bodyName, point2d, depthEstimate, X_b_gaze, point2d_ref),
	robotIndex_(robotIndex)
{}


GazeTask::GazeTask(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, const std::string& bodyName,
	const Eigen::Vector3d& point3d, const sva::PTransformd& X_b_gaze,
	const Eigen::Vector2d& point2d_ref):
	gazet_(mbs[robotIndex], bodyName, point3d, X_b_gaze, point2d_ref),
	robotIndex_(robotIndex)
{}

int GazeTask::dim()
{
	return 2;
}


void GazeTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	gazet_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
}


const Eigen::MatrixXd& GazeTask::jac()
{
	return gazet_.jac();
}


const Eigen::VectorXd& GazeTask::eval()
{
	return gazet_.eval();
}


const Eigen::VectorXd& GazeTask::speed()
{
	return gazet_.speed();
}


const Eigen::VectorXd& GazeTask::normalAcc()
{
	return gazet_.normalAcc();
}


/**
	*																PositionBasedVisServoTask
	*/


PositionBasedVisServoTask::PositionBasedVisServoTask(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, const std::string& bodyName,
	const sva::PTransformd& X_t_s,
	const sva::PTransformd& X_b_s):
	pbvst_(mbs[robotIndex], bodyName, X_t_s, X_b_s),
	robotIndex_(robotIndex)
{}

int PositionBasedVisServoTask::dim()
{
	return 6;
}


void PositionBasedVisServoTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	pbvst_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
}


const Eigen::MatrixXd& PositionBasedVisServoTask::jac()
{
	return pbvst_.jac();
}


const Eigen::VectorXd& PositionBasedVisServoTask::eval()
{
	return pbvst_.eval();
}


const Eigen::VectorXd& PositionBasedVisServoTask::speed()
{
	return pbvst_.speed();
}


const Eigen::VectorXd& PositionBasedVisServoTask::normalAcc()
{
	return pbvst_.normalAcc();
}


/**
	*													CoMTask
	*/


CoMTask::CoMTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, const Eigen::Vector3d& com):
	ct_(mbs[rI], com),
	robotIndex_(rI)
{}


CoMTask::CoMTask(const std::vector<rbd::MultiBody>& mbs, int rI,
	const Eigen::Vector3d& com, std::vector<double> weight):
	ct_(mbs[rI], com, std::move(weight)),
	robotIndex_(rI)
{}


void CoMTask::updateInertialParameters(const std::vector<rbd::MultiBody>& mbs)
{
	ct_.updateInertialParameters(mbs[robotIndex_]);
}


int CoMTask::dim()
{
	return 3;
}


void CoMTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	ct_.update(mbs[robotIndex_], mbcs[robotIndex_],
		rbd::computeCoM(mbs[robotIndex_], mbcs[robotIndex_]),
		data.normalAccB(robotIndex_));
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
	*													MultiCoMTask
	*/


MultiCoMTask::MultiCoMTask(const std::vector<rbd::MultiBody>& mbs,
	std::vector<int> rI, const Eigen::Vector3d& com, double stiffness,
	double weight):
	Task(weight),
	alphaDBegin_(-1),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	dimWeight_(Eigen::Vector3d::Ones()),
	posInQ_(rI.size()),
	mct_(mbs, std::move(rI), com),
	Q_(),
	C_(),
	CSum_(),
	preQ_()
{
	init(mbs);
}


MultiCoMTask::MultiCoMTask(const std::vector<rbd::MultiBody>& mbs,
	std::vector<int> rI, const Eigen::Vector3d& com, double stiffness,
	const Eigen::Vector3d& dimWeight, double weight):
	Task(weight),
	alphaDBegin_(-1),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	dimWeight_(dimWeight),
	posInQ_(rI.size()),
	mct_(mbs, std::move(rI), com),
	Q_(),
	C_(),
	CSum_(),
	preQ_()
{
	init(mbs);
}


void MultiCoMTask::updateInertialParameters(const std::vector<rbd::MultiBody>& mbs)
{
	mct_.updateInertialParameters(mbs);
}


void MultiCoMTask::stiffness(double stiffness)
{
	stiffness_ = stiffness;
	stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}


void MultiCoMTask::dimWeight(const Eigen::Vector3d& dim)
{
	dimWeight_ = dim;
}


void MultiCoMTask::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	auto minMaxIndex =
		std::minmax_element(mct_.robotIndexes().begin(), mct_.robotIndexes().end());
	alphaDBegin_ = data.alphaDBegin(*(minMaxIndex.first));
	int lastBegin = data.alphaDBegin(*(minMaxIndex.second));
	int lastAlphaD = data.alphaD(*(minMaxIndex.second));
	int size = lastBegin + lastAlphaD - alphaDBegin_;

	Q_.setZero(size, size);
	C_.setZero(size);

	posInQ_.clear();
	for(int r: mct_.robotIndexes())
	{
		posInQ_.push_back(data.alphaDBegin(r) - alphaDBegin_);
	}
}


void MultiCoMTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	mct_.update(mbs, mbcs, data.normalAccB());
	CSum_ = stiffness_*mct_.eval();
	CSum_ -= stiffnessSqrt_*mct_.speed();
	CSum_ -= mct_.normalAcc();
	for(int i = 0; i < int(posInQ_.size()); ++i)
	{
		int r = mct_.robotIndexes()[i];
		int begin = posInQ_[i];
		int dof = data.alphaD(r);

		const Eigen::MatrixXd& J = mct_.jac(i);
		preQ_.block(0, 0, 3, dof).noalias() = dimWeight_.asDiagonal()*J;

		Q_.block(begin, begin, dof, dof).noalias() =
			J.transpose()*preQ_.block(0, 0, 3, dof);
		C_.segment(begin, dof).noalias() = -J.transpose()*dimWeight_.asDiagonal()*CSum_;
	}
}


const Eigen::MatrixXd& MultiCoMTask::Q() const
{
	return Q_;
}


const Eigen::VectorXd& MultiCoMTask::C() const
{
	return C_;
}


const Eigen::VectorXd& MultiCoMTask::eval() const
{
	return mct_.eval();
}


const Eigen::VectorXd& MultiCoMTask::speed() const
{
	return mct_.speed();
}


void MultiCoMTask::init(const std::vector<rbd::MultiBody>& mbs)
{
	int maxDof = 0;
	for(int r: mct_.robotIndexes())
	{
		maxDof = std::max(maxDof, mbs[r].nrDof());
	}
	preQ_.resize(3, maxDof);
}


/**
	*													MultiRobotTransformTask
	*/


MultiRobotTransformTask::MultiRobotTransformTask(
	const std::vector<rbd::MultiBody>& mbs,
	int r1Index, int r2Index,
	const std::string& r1BodyName, const std::string& r2BodyName,
	const sva::PTransformd& X_r1b_r1s, const sva::PTransformd& X_r2b_r2s,
	double stiffness, double weight):
	Task(weight),
	alphaDBegin_(-1),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	dimWeight_(Eigen::Vector6d::Ones()),
	posInQ_(2, -1),
	robotIndexes_{{r1Index, r2Index}},
	mrtt_(mbs, r1Index, r2Index, r1BodyName, r2BodyName, X_r1b_r1s, X_r2b_r2s),
	Q_(),
	C_(),
	CSum_(Eigen::Vector6d::Zero()),
	preQ_()
{
	int maxDof = 0;
	for(int r: robotIndexes_)
	{
		maxDof = std::max(maxDof, mbs[r].nrDof());
	}
	preQ_.resize(6, maxDof);
}


void MultiRobotTransformTask::X_r1b_r1s(const sva::PTransformd& X_r1b_r1s)
{
	mrtt_.X_r1b_r1s(X_r1b_r1s);
}


const sva::PTransformd& MultiRobotTransformTask::X_r1b_r1s() const
{
	return mrtt_.X_r1b_r1s();
}


void MultiRobotTransformTask::X_r2b_r2s(const sva::PTransformd& X_r2b_r2s)
{
	mrtt_.X_r2b_r2s(X_r2b_r2s);
}


const sva::PTransformd& MultiRobotTransformTask::X_r2b_r2s() const
{
	return mrtt_.X_r2b_r2s();
}


void MultiRobotTransformTask::stiffness(double stiffness)
{
	stiffness_ = stiffness;
	stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}


void MultiRobotTransformTask::dimWeight(const Eigen::Vector6d& dim)
{
	dimWeight_ = dim;
}


void MultiRobotTransformTask::updateNrVars(
	const std::vector<rbd::MultiBody>& /* mbs */, const SolverData& data)
{
	auto minMaxIndex =
		std::minmax_element(robotIndexes_.begin(), robotIndexes_.end());
	alphaDBegin_ = data.alphaDBegin(*(minMaxIndex.first));
	int lastBegin = data.alphaDBegin(*(minMaxIndex.second));
	int lastAlphaD = data.alphaD(*(minMaxIndex.second));
	int size = lastBegin + lastAlphaD - alphaDBegin_;

	Q_.setZero(size, size);
	C_.setZero(size);

	posInQ_.clear();
	for(int r: robotIndexes_)
	{
		posInQ_.push_back(data.alphaDBegin(r) - alphaDBegin_);
	}
}


void MultiRobotTransformTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	mrtt_.update(mbs, mbcs, data.normalAccB());
	CSum_.noalias() = stiffness_*mrtt_.eval();
	CSum_.noalias() -= stiffnessSqrt_*mrtt_.speed();
	CSum_.noalias() -= mrtt_.normalAcc();

	// first we set to zero used part of Q and C
	for(int i = 0; i < int(posInQ_.size()); ++i)
	{
		int r = robotIndexes_[i];
		int begin = posInQ_[i];
		int dof = data.alphaD(r);
		Q_.block(begin, begin, dof, dof).setZero();
		C_.segment(begin, dof).setZero();
	}

	for(int i = 0; i < int(posInQ_.size()); ++i)
	{
		int r = robotIndexes_[i];
		int begin = posInQ_[i];
		int dof = data.alphaD(r);

		const Eigen::MatrixXd& J = mrtt_.jac(i);
		preQ_.block(0, 0, 6, dof).noalias() = dimWeight_.asDiagonal()*J;

		// scince the two robot index could be the same
		// we had to increment the Q and C matrix
		Q_.block(begin, begin, dof, dof).noalias() +=
			J.transpose()*preQ_.block(0, 0, 6, dof);
		C_.segment(begin, dof).noalias() -= J.transpose()*dimWeight_.asDiagonal()*CSum_;
	}
}


const Eigen::MatrixXd& MultiRobotTransformTask::Q() const
{
	return Q_;
}


const Eigen::VectorXd& MultiRobotTransformTask::C() const
{
	return C_;
}


const Eigen::VectorXd& MultiRobotTransformTask::eval() const
{
	return mrtt_.eval();
}


const Eigen::VectorXd& MultiRobotTransformTask::speed() const
{
	return mrtt_.speed();
}


/**
	*													MomentumTask
	*/


MomentumTask::MomentumTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, const sva::ForceVecd& mom):
	momt_(mbs[rI], mom),
	robotIndex_(rI)
{}


int MomentumTask::dim()
{
	return 6;
}


void MomentumTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	momt_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
}


const Eigen::MatrixXd& MomentumTask::jac()
{
	return momt_.jac();
}


const Eigen::VectorXd& MomentumTask::eval()
{
	return momt_.eval();
}


const Eigen::VectorXd& MomentumTask::speed()
{
	return momt_.speed();
}


const Eigen::VectorXd& MomentumTask::normalAcc()
{
	return momt_.normalAcc();
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


void ContactTask::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	int nrLambda = 0;
	begin_ = data.lambdaBegin();
	std::vector<FrictionCone> cones;
	int curLambda = 0;

	if(nrLambda == 0)
	{
		for(const BilateralContact& uc: data.allContacts())
		{
			curLambda = uc.nrLambda();
			if(uc.contactId == contactId_)
			{
				nrLambda = curLambda;
				cones = uc.r1Cones;
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


void ContactTask::update(const std::vector<rbd::MultiBody>& /* mbs */,
	const std::vector<rbd::MultiBodyConfig>& /* mbcs */,
	const SolverData& /* data */)
{
	/*C_.noalias() = -conesJac_.transpose()*
			(stiffness_*error_ - stiffnessSqrt_*errorD_);*/
          C_.noalias() = -conesJac_.transpose()*error_;
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


void GripperTorqueTask::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	using namespace Eigen;
	bool found = false;

	begin_ = data.bilateralBegin();
	for(const BilateralContact& bc: data.bilateralContacts())
	{
		int curLambda = 0;
		// compute the number of lambda needed by the current bilateral
		for(std::size_t i = 0; i < bc.r1Points.size(); ++i)
		{
			curLambda += bc.nrLambda(static_cast<int>(i));
		}

		if(bc.contactId == contactId_)
		{
			found = true;
			Q_.setZero(curLambda, curLambda);
			C_.resize(curLambda);

			int pos = 0;
			// minimize Torque applied on the gripper motor
			// min Sum_i^nrF  T_i·( p_i^T_o x f_i)
			for(std::size_t i = 0; i < bc.r1Cones.size(); ++i)
			{
				Vector3d T_o_p = bc.r1Points[i] - origin_;
				for(std::size_t j = 0; j < bc.r1Cones[i].generators.size(); ++j)
				{
					// we use abs because the contact force cannot apply
					// negative torque on the gripper
					C_(pos) = std::abs(
						axis_.transpose()*(T_o_p.cross(bc.r1Cones[i].generators[j])));
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


void GripperTorqueTask::update(const std::vector<rbd::MultiBody>& /* mbs */,
	const std::vector<rbd::MultiBodyConfig>& /* mbcs */,
	const SolverData& /* data */)
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


LinVelocityTask::LinVelocityTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, const std::string& bodyName,
	const Eigen::Vector3d& speed, const Eigen::Vector3d& bodyPoint):
	pt_(mbs[rI], bodyName, speed, bodyPoint),
	robotIndex_(rI)
{
}


int LinVelocityTask::dim()
{
	return 3;
}


void LinVelocityTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	pt_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
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


OrientationTrackingTask::OrientationTrackingTask(
	const std::vector<rbd::MultiBody>& mbs, int rI, const std::string& bodyName,
	const Eigen::Vector3d& bodyPoint, const Eigen::Vector3d& bodyAxis,
	const std::vector<std::string>& trackingJointsName,
	const Eigen::Vector3d& trackedPoint):
	robotIndex_(rI),
	ott_(mbs[rI], bodyName, bodyPoint, bodyAxis, trackingJointsName, trackedPoint),
	alphaVec_(mbs[rI].nrDof()),
	speed_(3),
	normalAcc_(3)
{}


int OrientationTrackingTask::dim()
{
	return 3;
}


void OrientationTrackingTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& /* data */)
{
	ott_.update(mbs[robotIndex_], mbcs[robotIndex_]);
	rbd::paramToVector(mbcs[robotIndex_].alpha, alphaVec_);

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

/**
	*											RelativeDistTask
	*/

RelativeDistTask::RelativeDistTask(const std::vector<rbd::MultiBody>& mbs, const int rIndex, const double timestep,
	tasks::RelativeDistTask::rbInfo& rbi1, tasks::RelativeDistTask::rbInfo& rbi2, const Eigen::Vector3d& u1, const Eigen::Vector3d& u2) :
	rIndex_(rIndex),
	rdt_(mbs[rIndex], timestep, rbi1, rbi2, u1, u2)
{
}

int RelativeDistTask::dim()
{
	return 1;
}


void RelativeDistTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	rdt_.update(mbs[rIndex_], mbcs[rIndex_], data.normalAccB(rIndex_));

}


const Eigen::MatrixXd& RelativeDistTask::jac()
{
	return rdt_.jac();
}


const Eigen::VectorXd& RelativeDistTask::eval()
{
	return rdt_.eval();
}


const Eigen::VectorXd& RelativeDistTask::speed()
{
	return rdt_.speed();
}


const Eigen::VectorXd& RelativeDistTask::normalAcc()
{
	return rdt_.normalAcc();
}


/**
	*											VectorOrientationTask
	*/


VectorOrientationTask::VectorOrientationTask(const std::vector<rbd::MultiBody>& mbs, int rI,
	const std::string& bodyName, const Eigen::Vector3d& bodyVector,
	const Eigen::Vector3d& targetVector):
	vot_(mbs[rI], bodyName, bodyVector, targetVector),
	robotIndex_(rI)
{
}


int VectorOrientationTask::dim()
{
	return 3;
}


void VectorOrientationTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	vot_.update(mbs[robotIndex_], mbcs[robotIndex_], data.normalAccB(robotIndex_));
}


const Eigen::MatrixXd& VectorOrientationTask::jac()
{
	return vot_.jac();
}


const Eigen::VectorXd& VectorOrientationTask::eval()
{
	return vot_.eval();
}


const Eigen::VectorXd& VectorOrientationTask::speed()
{
	return vot_.speed();
}


const Eigen::VectorXd& VectorOrientationTask::normalAcc()
{
	return vot_.normalAcc();
}

} // namespace qp

} // namespace tasks

