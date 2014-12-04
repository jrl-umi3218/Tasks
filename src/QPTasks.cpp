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
	*														SetPointTask
	*/


SetPointTask::SetPointTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double stiffness, double weight):
	Task(weight),
	hlTask_(hlTask),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	dimWeight_(Eigen::VectorXd::Ones(hlTask->dim())),
	robotIndex_(rI),
	alphaDBegin_(0),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	preQ_(hlTask->dim(), mbs[rI].nrDof()),
	CVecSum_(hlTask->dim()),
	preC_(hlTask->dim())
{}


SetPointTask::SetPointTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double stiffness, const Eigen::VectorXd& dimWeight, double weight):
	Task(weight),
	hlTask_(hlTask),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	dimWeight_(dimWeight),
	robotIndex_(rI),
	alphaDBegin_(0),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	preQ_(hlTask->dim(), mbs[rI].nrDof()),
	CVecSum_(hlTask->dim()),
	preC_(hlTask->dim())
{}


void SetPointTask::stiffness(double stiffness)
{
	stiffness_ = stiffness;
	stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}


void SetPointTask::dimWeight(const Eigen::VectorXd& dim)
{
	dimWeight_ = dim;
}


void SetPointTask::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
}


void SetPointTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	hlTask_->update(mbs, mbcs, data);

	const Eigen::MatrixXd& J = hlTask_->jac();
	const Eigen::VectorXd& err = hlTask_->eval();
	const Eigen::VectorXd& speed = hlTask_->speed();
	const Eigen::VectorXd& normalAcc = hlTask_->normalAcc();

	preQ_.noalias() = dimWeight_.asDiagonal()*J;
	Q_.noalias() = J.transpose()*preQ_;

	CVecSum_.noalias() = stiffness_*err;
	CVecSum_.noalias() -= stiffnessSqrt_*speed;
	CVecSum_.noalias() -= normalAcc;
	preC_.noalias() = dimWeight_.asDiagonal()*CVecSum_;
	C_.noalias() = -J.transpose()*preC_;
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


PIDTask::PIDTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double P, double I, double D, double weight):
	Task(weight),
	hlTask_(hlTask),
	P_(P),
	I_(I),
	D_(D),
	dimWeight_(Eigen::VectorXd::Ones(hlTask->dim())),
	robotIndex_(rI),
	alphaDBegin_(0),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	error_(hlTask->dim()),
	errorD_(hlTask->dim()),
	errorI_(hlTask->dim()),
	preQ_(hlTask->dim(), mbs[rI].nrDof()),
	CVecSum_(hlTask->dim()),
	preC_(hlTask->dim())
{}


PIDTask::PIDTask(const std::vector<rbd::MultiBody>& mbs,
	int rI,
	HighLevelTask* hlTask,
	double P, double I, double D,
	const Eigen::VectorXd& dimWeight, double weight):
	Task(weight),
	hlTask_(hlTask),
	P_(P),
	I_(I),
	D_(D),
	dimWeight_(dimWeight),
	robotIndex_(rI),
	alphaDBegin_(0),
	Q_(mbs[rI].nrDof(), mbs[rI].nrDof()),
	C_(mbs[rI].nrDof()),
	error_(hlTask->dim()),
	errorD_(hlTask->dim()),
	errorI_(hlTask->dim()),
	preQ_(hlTask->dim(), mbs[rI].nrDof()),
	CVecSum_(hlTask->dim()),
	preC_(hlTask->dim())
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


void PIDTask::dimWeight(const Eigen::VectorXd& dim)
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


void PIDTask::updateNrVars(const std::vector<rbd::MultiBody>& /* mbs */,
	const SolverData& data)
{
	alphaDBegin_ = data.alphaDBegin(robotIndex_);
}


void PIDTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	hlTask_->update(mbs, mbcs, data);

	const Eigen::MatrixXd& J = hlTask_->jac();
	const Eigen::VectorXd& normalAcc = hlTask_->normalAcc();

	preQ_.noalias() = dimWeight_.asDiagonal()*J;
	Q_.noalias() = J.transpose()*preQ_;

	CVecSum_.noalias() = P_*error_;
	CVecSum_.noalias() -= D_*errorD_;
	CVecSum_.noalias() -= I_*errorI_;
	CVecSum_.noalias() -= normalAcc;
	preC_.noalias() = dimWeight_.asDiagonal()*CVecSum_;
	C_.noalias() = -J.transpose()*preC_;
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
	int robotIndex, HighLevelTask* hl, const std::vector<int>& activeJointsId)
{
	return JointsSelector(mbs, robotIndex, hl, activeJointsId);
}


JointsSelector JointsSelector::UnactiveJoints(const std::vector<rbd::MultiBody>& mbs,
	int robotIndex, HighLevelTask* hl, const std::vector<int>& unactiveJointsId)
{
	using namespace std::placeholders;
	const rbd::MultiBody& mb = mbs[robotIndex];

	std::vector<int> activeJointsId;
	// sort unactiveJointsId by puting them into a set
	std::set<int> unactiveJointsIdSet(unactiveJointsId.begin(),
		unactiveJointsId.end());
	// create a set with all joints id
	std::set<int> jointsIdSet;
	std::transform(mb.joints().begin(), mb.joints().end(),
		std::inserter(jointsIdSet, jointsIdSet.begin()),
		std::bind(&rbd::Joint::id, _1));

	// remove unactive joints from the set
	std::set_difference(jointsIdSet.begin(), jointsIdSet.end(),
		unactiveJointsIdSet.begin(), unactiveJointsIdSet.end(),
		std::inserter(activeJointsId, activeJointsId.begin()));

	return JointsSelector(mbs, robotIndex, hl, activeJointsId);
}


JointsSelector::JointsSelector(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
	HighLevelTask* hl, const std::vector<int>& selectedJointsId):
	jac_(Eigen::MatrixXd::Zero(hl->dim(), mbs[robotIndex].nrDof())),
	selectedJoints_(),
	hl_(hl)
{
	const rbd::MultiBody& mb = mbs[robotIndex];
	selectedJoints_.reserve(selectedJointsId.size());
	for(int jId: selectedJointsId)
	{
		int index = mb.jointIndexById(jId);
		selectedJoints_.push_back({mb.jointPosInDof(index), mb.joint(index).dof()});
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
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
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
	stiffnessSqrt_ = 2.*std::sqrt(stiffness);
}


void PostureTask::jointsStiffness(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<JointStiffness>& jsv)
{
	jointDatas_.clear();
	jointDatas_.reserve(jsv.size());

	const rbd::MultiBody& mb = mbs[robotIndex_];
	for(const JointStiffness& js: jsv)
	{
		int jointIndex = mb.jointIndexById(js.jointId);
		jointDatas_.push_back({js.stiffness, 2.*std::sqrt(js.stiffness),
													mb.jointPosInDof(jointIndex),
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


PositionTask::PositionTask(const std::vector<rbd::MultiBody>& mbs, int rI,
	int bodyId, const Eigen::Vector3d& pos, const Eigen::Vector3d& bodyPoint):
	pt_(mbs[rI], bodyId, pos, bodyPoint),
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
	int rI, int bodyId,
	const Eigen::Quaterniond& ori):
	ot_(mbs[rI], bodyId, ori),
	robotIndex_(rI)
{}


OrientationTask::OrientationTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, int bodyId,
	const Eigen::Matrix3d& ori):
	ot_(mbs[rI], bodyId, ori),
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
	*																SurfaceOrientationTask
	*/


SurfaceOrientationTask::SurfaceOrientationTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, int bodyId,
	const Eigen::Quaterniond& ori, const sva::PTransformd& X_b_s):
	ot_(mbs[rI], bodyId, ori, X_b_s),
	robotIndex_(rI)
{}


SurfaceOrientationTask::SurfaceOrientationTask(const std::vector<rbd::MultiBody>& mbs,
	int rI, int bodyId,
	const Eigen::Matrix3d& ori, const sva::PTransformd& X_b_s):
	ot_(mbs[rI], bodyId, ori, X_b_s),
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
	int r1Index, int r2Index, int r1BodyId, int r2BodyId,
	const sva::PTransformd& X_r1b_r1s, const sva::PTransformd& X_r2b_r2s,
	double stiffness, double weight):
	Task(weight),
	alphaDBegin_(-1),
	stiffness_(stiffness),
	stiffnessSqrt_(2.*std::sqrt(stiffness)),
	dimWeight_(Eigen::Vector6d::Ones()),
	posInQ_(2, -1),
	robotIndexes_{{r1Index, r2Index}},
	mrtt_(mbs, r1Index, r2Index, r1BodyId, r2BodyId, X_r1b_r1s, X_r2b_r2s),
	Q_(),
	C_(),
	CSum_(),
	preQ_()
{}


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
	for(int i = 0; i < int(posInQ_.size()); ++i)
	{
		int r = robotIndexes_[i];
		int begin = posInQ_[i];
		int dof = data.alphaD(r);

		const Eigen::MatrixXd& J = mrtt_.jac(i);
		preQ_.block(0, 0, 3, dof).noalias() = dimWeight_.asDiagonal()*J;

		Q_.block(begin, begin, dof, dof).noalias() =
			J.transpose()*preQ_.block(0, 0, 3, dof);
		C_.segment(begin, dof).noalias() = -J.transpose()*dimWeight_.asDiagonal()*CSum_;
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

	if(nrLambda == 0)
	{
		for(const BilateralContact& uc: data.allContacts())
		{
			int curLambda = 0;
			for(std::size_t i = 0; i < uc.r1Points.size(); ++i)
			{
				curLambda += uc.nrLambda(static_cast<int>(i));
			}

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
	int rI, int bodyId,
	const Eigen::Vector3d& speed, const Eigen::Vector3d& bodyPoint):
	pt_(mbs[rI], bodyId, speed, bodyPoint),
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
	const std::vector<rbd::MultiBody>& mbs, int rI, int bodyId,
	const Eigen::Vector3d& bodyPoint, const Eigen::Vector3d& bodyAxis,
	const std::vector<int>& trackingJointsId,
	const Eigen::Vector3d& trackedPoint):
	robotIndex_(rI),
	ott_(mbs[rI], bodyId, bodyPoint, bodyAxis, trackingJointsId, trackedPoint),
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

} // namespace qp

} // namespace tasks

