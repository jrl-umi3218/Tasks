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

#pragma once

// includes
//
#include <array>

// Eigen
#include <Eigen/Core>

// Tasks
#include "Tasks.h"
#include "QPSolver.h"
#include "QPMotionConstr.h"

// forward declaration
// RBDyn
namespace rbd
{
class MultiBody;
struct MultiBodyConfig;
}

namespace tasks
{

namespace qp
{


class TASKS_DLLAPI SetPointTaskCommon : public Task
{
public:
	SetPointTaskCommon(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double weight);
	SetPointTaskCommon(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, const Eigen::VectorXd& dimWeight,
		double weight);

	virtual std::pair<int, int> begin() const override
	{
		return std::make_pair(alphaDBegin_, alphaDBegin_);
	}

	void dimWeight(const Eigen::VectorXd& dim);

	const Eigen::VectorXd& dimWeight() const
	{
		return dimWeight_;
	}

	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& Q() const override;
	virtual const Eigen::VectorXd& C() const override;

protected:
	void computeQC(Eigen::VectorXd& error);

protected:
	HighLevelTask* hlTask_;
	Eigen::VectorXd error_;

private:
	Eigen::VectorXd dimWeight_;
	int robotIndex_, alphaDBegin_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
	// cache
	Eigen::MatrixXd preQ_;
	Eigen::VectorXd preC_;
};


class TASKS_DLLAPI SetPointTask : public SetPointTaskCommon
{
public:
	SetPointTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double stiffness, double weight);

	SetPointTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double stiffness, const Eigen::VectorXd& dimWeight,
		double weight);

	double stiffness() const
	{
		return stiffness_;
	}

	void stiffness(double stiffness);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

private:
	double stiffness_, stiffnessSqrt_;
};


class TASKS_DLLAPI TrackingTask : public SetPointTaskCommon
{
public:
	TrackingTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double gainPos, double gainVel, double weight);

	TrackingTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double gainPos, double gainVel,
		const Eigen::VectorXd& dimWeight, double weight);

	void setGains(double gainPos, double gainVel);

	void errorPos(const Eigen::VectorXd& errorPos);
	void errorVel(const Eigen::VectorXd& errorVel);
	void refAccel(const Eigen::VectorXd& refAccel);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

private:
	double gainPos_, gainVel_;
	Eigen::VectorXd errorPos_, errorVel_, refAccel_;
};


class TASKS_DLLAPI TrajectoryTask : public SetPointTaskCommon
{
public:
	TrajectoryTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double gainPos, double gainVel, double weight);

	TrajectoryTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double gainPos, double gainVel,
		const Eigen::VectorXd& dimWeight, double weight);

	void setGains(double gainPos, double gainVel);
	void setGains(const Eigen::VectorXd & stiffness, const Eigen::VectorXd & damping);
	void stiffness(double gainPos);
	void stiffness(const Eigen::VectorXd & stiffness);
	const Eigen::VectorXd & stiffness() const;
	void damping(double gainVel);
	void damping(const Eigen::VectorXd & damping);
	const Eigen::VectorXd & damping() const;

	void refVel(const Eigen::VectorXd& refVel);
	const Eigen::VectorXd & refVel() const;
	void refAccel(const Eigen::VectorXd& refAccel);
	const Eigen::VectorXd & refAccel() const;

	void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

private:
	Eigen::VectorXd stiffness_, damping_;
	Eigen::VectorXd refVel_, refAccel_;
};


/// @deprecated Must be replace by TrackingTask
class TASKS_DLLAPI PIDTask : public SetPointTaskCommon
{
public:
	PIDTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double P, double I, double D, double weight);

	PIDTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask, double P, double I, double D,
		const Eigen::VectorXd& dimWeight, double weight);

	double P() const;
	void P(double p);
	double I() const;
	void I(double i);
	double D() const;
	void D(double d);

	void error(const Eigen::VectorXd& err);
	void errorD(const Eigen::VectorXd& errD);
	void errorI(const Eigen::VectorXd& errI);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

private:
	double P_, I_, D_;
	Eigen::VectorXd error_, errorD_, errorI_;
};


class TASKS_DLLAPI TargetObjectiveTask : public Task
{
public:
	TargetObjectiveTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask,
		double timeStep, double duration, const Eigen::VectorXd& objDot,
		double weight);

	TargetObjectiveTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hlTask,
		double timeStep, double duration, const Eigen::VectorXd& objDot,
		const Eigen::VectorXd& dimWeight, double weight);

	double duration() const;
	void duration(double d);

	int iter() const
	{
		return iter_;
	}
	void iter(int i)
	{
		iter_ = i;
	}

	int nrIter() const
	{
		return nrIter_;
	}
	void nrIter(int i)
	{
		nrIter_ = i;
	}

	const Eigen::VectorXd& objDot() const
	{
		return objDot_;
	}
	void objDot(const Eigen::VectorXd& o)
	{
		objDot_ = o;
	}

	const Eigen::VectorXd& dimWeight() const
	{
		return dimWeight_;
	}
	void dimWeight(const Eigen::VectorXd& o)
	{
		dimWeight_ = o;
	}

	const Eigen::VectorXd& phi() const
	{
		return phi_;
	}
	const Eigen::VectorXd& psi() const
	{
		return psi_;
	}


	virtual std::pair<int, int> begin() const override
	{
		return std::make_pair(alphaDBegin_, alphaDBegin_);
	}

	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& Q() const override;
	virtual const Eigen::VectorXd& C() const override;

private:
	HighLevelTask* hlTask_;

	int iter_, nrIter_;
	double dt_;
	Eigen::VectorXd objDot_;
	Eigen::VectorXd dimWeight_;
	int robotIndex_, alphaDBegin_;

	Eigen::VectorXd phi_, psi_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
	// cache
	Eigen::MatrixXd preQ_;
	Eigen::VectorXd CVecSum_, preC_;
};



class TASKS_DLLAPI JointsSelector : public HighLevelTask
{
public:
	static JointsSelector ActiveJoints(const std::vector<rbd::MultiBody>& mbs,
		int robotIndex, HighLevelTask* hl,
		const std::vector<std::string>& activeJointsName,
		const std::map<std::string, std::vector<std::array<int, 2>>>& activeDofs = {});
	static JointsSelector UnactiveJoints(const std::vector<rbd::MultiBody>& mbs,
		int robotIndex, HighLevelTask* hl,
		const std::vector<std::string>& unactiveJointsName,
		const std::map<std::string, std::vector<std::array<int, 2>>>& unactiveDofs = {});

public:
	struct SelectedData
	{
		int posInDof, dof;
	};

public:
	JointsSelector(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		HighLevelTask* hl, const std::vector<std::string>& selectedJointsName,
		const std::map<std::string, std::vector<std::array<int, 2>>>& activeDofs = {});

	const std::vector<SelectedData> selectedJoints() const
	{
		return selectedJoints_;
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	Eigen::MatrixXd jac_;
	std::vector<SelectedData> selectedJoints_;
	HighLevelTask* hl_;
};



struct JointStiffness
{
	JointStiffness():
		jointName(),
		stiffness()
	{}
	JointStiffness(const std::string& jName, double stif):
		jointName(jName),
		stiffness(stif)
	{}

	std::string jointName;
	double stiffness;
};

struct JointGains
{
	JointGains():
		jointName(),
		stiffness(),
		damping()
	{}
	JointGains(const std::string& jName, double stif) :
		jointName(jName),
		stiffness(stif)
	{
		damping = 2.*std::sqrt(stif);
	}

	JointGains(const std::string& jName, double stif, double damp):
		jointName(jName),
		stiffness(stif),
		damping(damp)
	{}

	std::string jointName;
	double stiffness, damping;
};

class TASKS_DLLAPI TorqueTask : public Task
{
public:
        TorqueTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                   const TorqueBound& tb, double weight);

        TorqueTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                   const TorqueBound& tb, const Eigen::VectorXd& jointSelect,
                   double weight);

        TorqueTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                   const TorqueBound& tb, const std::string& efName,
                   double weight);

	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual std::pair<int, int> begin() const override
	{
		return std::make_pair(0, 0);
	}

	virtual const Eigen::MatrixXd& Q() const override
        {
          return Q_;
        }

	virtual const Eigen::VectorXd& C() const override
        {
          return C_;
        }

        virtual const Eigen::VectorXd& jointSelect() const
        {
          return jointSelector_;
        }

private:
        int robotIndex_;
        int alphaDBegin_, lambdaBegin_;
        MotionConstr motionConstr;
        Eigen::VectorXd jointSelector_;
        Eigen::MatrixXd Q_;
        Eigen::VectorXd C_;
};

class TASKS_DLLAPI PostureTask : public Task
{
public:
	PostureTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		std::vector<std::vector<double> > q, double stiffness, double weight);

	tasks::PostureTask& task()
	{
		return pt_;
	}

	void posture(std::vector<std::vector<double> > q)
	{
		pt_.posture(q);
	}

	const std::vector<std::vector<double> > posture() const
	{
		return pt_.posture();
	}

	double stiffness() const
	{
		return stiffness_;
	}

	double damping() const
	{
		return damping_;
	}

	void stiffness(double stiffness);

	void gains(double stiffness);
	void gains(double stiffness, double damping);

	void jointsStiffness(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<JointStiffness>& jsv);

	void jointsGains(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<JointGains>& jgv);

	virtual std::pair<int, int> begin() const override
	{
		return std::make_pair(alphaDBegin_, alphaDBegin_);
	}

	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& Q() const override;
	virtual const Eigen::VectorXd& C() const override;

	const Eigen::VectorXd& eval() const;

private:
	struct JointData
	{
		double stiffness, damping;
		int start, size;
	};

private:
	tasks::PostureTask pt_;

	double stiffness_;
	double damping_;
	int robotIndex_, alphaDBegin_;

	std::vector<JointData> jointDatas_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
	Eigen::VectorXd alphaVec_;
};



class TASKS_DLLAPI PositionTask : public HighLevelTask
{
public:
	PositionTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const Eigen::Vector3d& pos,
		const Eigen::Vector3d& bodyPoint=Eigen::Vector3d::Zero());

	tasks::PositionTask& task()
	{
		return pt_;
	}

	void position(const Eigen::Vector3d& pos)
	{
		pt_.position(pos);
	}

	const Eigen::Vector3d& position() const
	{
		return pt_.position();
	}

	void bodyPoint(const Eigen::Vector3d& point)
	{
		pt_.bodyPoint(point);
	}

	const Eigen::Vector3d& bodyPoint() const
	{
		return pt_.bodyPoint();
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mb,
		const std::vector<rbd::MultiBodyConfig>& mbc,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	tasks::PositionTask pt_;
	int robotIndex_;
};



class TASKS_DLLAPI OrientationTask : public HighLevelTask
{
public:
	OrientationTask(const std::vector<rbd::MultiBody>& mbs, int robodIndex,
		const std::string& bodyName, const Eigen::Quaterniond& ori);
	OrientationTask(const std::vector<rbd::MultiBody>& mbs, int robodIndex,
		const std::string& bodyName, const Eigen::Matrix3d& ori);

	tasks::OrientationTask& task()
	{
		return ot_;
	}

	void orientation(const Eigen::Quaterniond& ori)
	{
		ot_.orientation(ori);
	}

	void orientation(const Eigen::Matrix3d& ori)
	{
		ot_.orientation(ori);
	}

	const Eigen::Matrix3d& orientation() const
	{
		return ot_.orientation();
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	tasks::OrientationTask ot_;
	int robotIndex_;
};



template <typename transform_task_t>
class TransformTaskCommon : public HighLevelTask
{
public:
	TransformTaskCommon(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const sva::PTransformd& X_0_t,
		const sva::PTransformd& X_b_p=sva::PTransformd::Identity())
	:  tt_(mbs[robotIndex], bodyName, X_0_t, X_b_p),
	   robotIndex_(robotIndex)
	{
	}

	transform_task_t& task()
	{
		return tt_;
	}

	void target(const sva::PTransformd& X_0_t)
	{
		tt_.target(X_0_t);
	}

	const sva::PTransformd& target() const
	{
		return tt_.target();
	}

	void X_b_p(const sva::PTransformd& X_b_p)
	{
		tt_.X_b_p(X_b_p);
	}

	const sva::PTransformd& X_b_p() const
	{
		return tt_.X_b_p();
	}

	virtual int dim() override
	{
	  return 6;
	}

	virtual const Eigen::MatrixXd& jac() override
	{
	  return tt_.jac();
	}

	virtual const Eigen::VectorXd& eval() override
	{
	  return tt_.eval();
	}

	virtual const Eigen::VectorXd& speed() override
	{
	  return tt_.speed();
	}

	virtual const Eigen::VectorXd& normalAcc() override
	{
	  return tt_.normalAcc();
	}

protected:
	transform_task_t tt_;
	int robotIndex_;
};



/// TransformTask in surface frame.
class TASKS_DLLAPI SurfaceTransformTask : public TransformTaskCommon<tasks::SurfaceTransformTask>
{
public:
	SurfaceTransformTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const sva::PTransformd& X_0_t,
		const sva::PTransformd& X_b_p=sva::PTransformd::Identity());

	virtual void update(const std::vector<rbd::MultiBody>& mb,
		const std::vector<rbd::MultiBodyConfig>& mbc,
		const SolverData& data) override;
};



/// TransformTask in world or user frame.
class TASKS_DLLAPI TransformTask : public TransformTaskCommon<tasks::TransformTask>
{
public:
	TransformTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const sva::PTransformd& X_0_t,
		const sva::PTransformd& X_b_p=sva::PTransformd::Identity(),
		const Eigen::Matrix3d& E_0_c=Eigen::Matrix3d::Identity());

	void E_0_c(const Eigen::Matrix3d& E_0_c);
	const Eigen::Matrix3d& E_0_c() const;

	virtual void update(const std::vector<rbd::MultiBody>& mb,
		const std::vector<rbd::MultiBodyConfig>& mbc,
		const SolverData& data) override;
};



class TASKS_DLLAPI SurfaceOrientationTask : public HighLevelTask
{
public:
	SurfaceOrientationTask(const std::vector<rbd::MultiBody>& mbs, int robodIndex,
			const std::string& bodyName, const Eigen::Quaterniond& ori,
			const sva::PTransformd& X_b_s);
	SurfaceOrientationTask(const std::vector<rbd::MultiBody>& mbs, int robodIndex,
			const std::string& bodyName, const Eigen::Matrix3d& ori,
			const sva::PTransformd& X_b_s);

	tasks::SurfaceOrientationTask& task()
	{
			return ot_;
	}

	void orientation(const Eigen::Quaterniond& ori)
	{
			ot_.orientation(ori);
	}

	void orientation(const Eigen::Matrix3d& ori)
	{
			ot_.orientation(ori);
	}

	const Eigen::Matrix3d& orientation() const
	{
			return ot_.orientation();
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
			const std::vector<rbd::MultiBodyConfig>& mbcs,
			const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	tasks::SurfaceOrientationTask ot_;
	int robotIndex_;
};



class TASKS_DLLAPI GazeTask : public HighLevelTask
{
public:
	GazeTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const Eigen::Vector2d& point2d,
		double depthEstimate, const sva::PTransformd& X_b_gaze,
		const Eigen::Vector2d& point2d_ref = Eigen::Vector2d::Zero());
	GazeTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const Eigen::Vector3d& point3d,
		const sva::PTransformd& X_b_gaze,
		const Eigen::Vector2d& point2d_ref = Eigen::Vector2d::Zero());

	tasks::GazeTask& task()
	{
			return gazet_;
	}

	void error(const Eigen::Vector2d& point2d,
		const Eigen::Vector2d& point2d_ref = Eigen::Vector2d::Zero())
	{
			gazet_.error(point2d, point2d_ref);
	}

	void error(const Eigen::Vector3d& point3d,
		const Eigen::Vector2d& point2d_ref = Eigen::Vector2d::Zero())
	{
			gazet_.error(point3d, point2d_ref);
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
			const std::vector<rbd::MultiBodyConfig>& mbcs, const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	tasks::GazeTask gazet_;
	int robotIndex_;
};



class TASKS_DLLAPI PositionBasedVisServoTask : public HighLevelTask
{
public:
	PositionBasedVisServoTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const sva::PTransformd& X_t_s,
		const sva::PTransformd& X_b_s=sva::PTransformd::Identity());

	tasks::PositionBasedVisServoTask& task()
	{
			return pbvst_;
	}

	void error(const sva::PTransformd& X_t_s)
	{
			pbvst_.error(X_t_s);
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
			const std::vector<rbd::MultiBodyConfig>& mbcs, const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	tasks::PositionBasedVisServoTask pbvst_;
	int robotIndex_;
};



class TASKS_DLLAPI CoMTask : public HighLevelTask
{
public:
	CoMTask(const std::vector<rbd::MultiBody>& mbs,
		int robotIndex, const Eigen::Vector3d& com);
	CoMTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const Eigen::Vector3d& com, std::vector<double> weight);

	tasks::CoMTask& task()
	{
		return ct_;
	}

	void com(const Eigen::Vector3d& com)
	{
		ct_.com(com);
	}

	const Eigen::Vector3d com() const
	{
		return ct_.com();
	}

	void updateInertialParameters(const std::vector<rbd::MultiBody>& mbs);

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	tasks::CoMTask ct_;
	int robotIndex_;
};


class TASKS_DLLAPI MultiCoMTask : public Task
{
public:
	MultiCoMTask(const std::vector<rbd::MultiBody>& mb,
		std::vector<int> robotIndexes, const Eigen::Vector3d& com,
		double stiffness, double weight);
	MultiCoMTask(const std::vector<rbd::MultiBody>& mb,
		std::vector<int> robotIndexes, const Eigen::Vector3d& com,
		double stiffness, const Eigen::Vector3d& dimWeight, double weight);

	tasks::MultiCoMTask& task()
	{
		return mct_;
	}

	void com(const Eigen::Vector3d& com)
	{
		mct_.com(com);
	}

	const Eigen::Vector3d com() const
	{
		return mct_.com();
	}

	void updateInertialParameters(const std::vector<rbd::MultiBody>& mbs);

	double stiffness() const
	{
		return stiffness_;
	}

	void stiffness(double stiffness);

	void dimWeight(const Eigen::Vector3d& dim);

	const Eigen::Vector3d& dimWeight() const
	{
		return dimWeight_;
	}

	virtual std::pair<int, int> begin() const override
	{
		return {alphaDBegin_, alphaDBegin_};
	}

	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& Q() const override;
	virtual const Eigen::VectorXd& C() const override;

	const Eigen::VectorXd& eval() const;
	const Eigen::VectorXd& speed() const;

private:
	void init(const std::vector<rbd::MultiBody>& mbs);

private:
	int alphaDBegin_;
	double stiffness_, stiffnessSqrt_;
	Eigen::Vector3d dimWeight_;
	std::vector<int> posInQ_;
	tasks::MultiCoMTask mct_;
	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
	Eigen::Vector3d CSum_;
	// cache
	Eigen::MatrixXd preQ_;
};


class TASKS_DLLAPI MultiRobotTransformTask : public Task
{
public:
	MultiRobotTransformTask(const std::vector<rbd::MultiBody>& mbs,
		int r1Index, int r2Index,
		const std::string& r1BodyName, const std::string& r2BodyName,
		const sva::PTransformd& X_r1b_r1s, const sva::PTransformd& X_r2b_r2s,
		double stiffness, double weight);

	tasks::MultiRobotTransformTask& task()
	{
		return mrtt_;
	}

	void X_r1b_r1s(const sva::PTransformd& X_r1b_r1s);
	const sva::PTransformd& X_r1b_r1s() const;

	void X_r2b_r2s(const sva::PTransformd& X_r2b_r2s);
	const sva::PTransformd& X_r2b_r2s() const;

	double stiffness() const
	{
		return stiffness_;
	}

	void stiffness(double stiffness);

	void dimWeight(const Eigen::Vector6d& dim);

	const Eigen::VectorXd& dimWeight() const
	{
		return dimWeight_;
	}

	virtual std::pair<int, int> begin() const override
	{
		return {alphaDBegin_, alphaDBegin_};
	}

	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& Q() const override;
	virtual const Eigen::VectorXd& C() const override;

	const Eigen::VectorXd& eval() const;
	const Eigen::VectorXd& speed() const;

private:
	int alphaDBegin_;
	double stiffness_, stiffnessSqrt_;
	Eigen::VectorXd dimWeight_;
	std::vector<int> posInQ_, robotIndexes_;
	tasks::MultiRobotTransformTask mrtt_;
	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
	Eigen::VectorXd CSum_;
	// cache
	Eigen::MatrixXd preQ_;
};


class TASKS_DLLAPI MomentumTask : public HighLevelTask
{
public:
	MomentumTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const sva::ForceVecd& mom);

	tasks::MomentumTask& task()
	{
		return momt_;
	}

	void momentum(const sva::ForceVecd& mom)
	{
		momt_.momentum(mom);
	}

	const sva::ForceVecd momentum()
	{
		return momt_.momentum();
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mb,
		const std::vector<rbd::MultiBodyConfig>& mbc,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	tasks::MomentumTask momt_;
	int robotIndex_;
};


class TASKS_DLLAPI ContactTask : public Task
{
public:
	ContactTask(ContactId contactId, double stiffness, double weight):
		Task(weight),
		contactId_(contactId),
		begin_(0),
		stiffness_(stiffness),
		stiffnessSqrt_(2*std::sqrt(stiffness)),
		conesJac_(),
		error_(Eigen::Vector3d::Zero()),
		errorD_(Eigen::Vector3d::Zero()),
		Q_(),
		C_()
	{}

	virtual std::pair<int, int> begin() const override
	{
		return std::make_pair(begin_, begin_);
	}

	void error(const Eigen::Vector3d& error);
	void errorD(const Eigen::Vector3d& errorD);

	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& Q() const override;
	virtual const Eigen::VectorXd& C() const override;

private:
	ContactId contactId_;
	int begin_;

	double stiffness_, stiffnessSqrt_;
	Eigen::MatrixXd conesJac_;
	Eigen::Vector3d error_, errorD_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
};



class TASKS_DLLAPI GripperTorqueTask : public Task
{
public:
	GripperTorqueTask(ContactId contactId, const Eigen::Vector3d& origin,
		const Eigen::Vector3d& axis, double weight):
		Task(weight),
		contactId_(contactId),
		origin_(origin),
		axis_(axis),
		begin_(0),
		Q_(),
		C_()
	{}

	virtual std::pair<int, int> begin() const override
	{
		return std::make_pair(begin_, begin_);
	}

	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& Q() const override;
	virtual const Eigen::VectorXd& C() const override;

private:
	ContactId contactId_;
	Eigen::Vector3d origin_;
	Eigen::Vector3d axis_;
	int begin_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
};



class TASKS_DLLAPI LinVelocityTask : public HighLevelTask
{
public:
	LinVelocityTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const Eigen::Vector3d& vel,
		const Eigen::Vector3d& bodyPoint=Eigen::Vector3d::Zero());

	tasks::LinVelocityTask& task()
	{
		return pt_;
	}

	void velocity(const Eigen::Vector3d& s)
	{
		pt_.velocity(s);
	}

	const Eigen::Vector3d& velocity() const
	{
		return pt_.velocity();
	}

	void bodyPoint(const Eigen::Vector3d& point)
	{
		pt_.bodyPoint(point);
	}

	const Eigen::Vector3d& bodyPoint() const
	{
		return pt_.bodyPoint();
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	tasks::LinVelocityTask pt_;
	int robotIndex_;
};



class TASKS_DLLAPI OrientationTrackingTask : public HighLevelTask
{
public:
	OrientationTrackingTask(const std::vector<rbd::MultiBody>& mbs,
		int robotIndex, const std::string& bodyName,
		const Eigen::Vector3d& bodyPoint, const Eigen::Vector3d& bodyAxis,
		const std::vector<std::string>& trackingJointsName,
		const Eigen::Vector3d& trackedPoint);

	tasks::OrientationTrackingTask& task()
	{
		return ott_;
	}

	void trackedPoint(const Eigen::Vector3d& tp)
	{
		ott_.trackedPoint(tp);
	}

	const Eigen::Vector3d& trackedPoint() const
	{
		return ott_.trackedPoint();
	}

	void bodyPoint(const Eigen::Vector3d& bp)
	{
		ott_.bodyPoint(bp);
	}

	const Eigen::Vector3d& bodyPoint() const
	{
		return ott_.bodyPoint();
	}

	void bodyAxis(const Eigen::Vector3d& ba)
	{
		ott_.bodyAxis(ba);
	}

	const Eigen::Vector3d& bodyAxis() const
	{
		return ott_.bodyAxis();
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	int robotIndex_;
	tasks::OrientationTrackingTask ott_;
	Eigen::VectorXd alphaVec_;
	Eigen::VectorXd speed_, normalAcc_;
};


class TASKS_DLLAPI RelativeDistTask : public HighLevelTask
{
public:
	RelativeDistTask(const std::vector<rbd::MultiBody>& mbs, const int rIndex, const double timestep,
		tasks::RelativeDistTask::rbInfo& rbi1, tasks::RelativeDistTask::rbInfo& rbi2,
		const Eigen::Vector3d& u1=Eigen::Vector3d::Zero(),
		const Eigen::Vector3d& u2=Eigen::Vector3d::Zero());

	tasks::RelativeDistTask& task()
	{
		return rdt_;
	}

	void robotPoint(const rbd::MultiBody& mb, const std::string& bName,
			const Eigen::Vector3d& point)
	{
		int bIndex = mb.bodyIndexByName(bName);
		rdt_.robotPoint(bIndex, point);
	}
	void envPoint(const rbd::MultiBody& mb, const std::string&  bName,
			const Eigen::Vector3d& point)
	{
		int bIndex = mb.bodyIndexByName(bName);
		rdt_.envPoint(bIndex, point);
	}
	void vector(const rbd::MultiBody& mb, const std::string& bName,
			const Eigen::Vector3d& u)
	{
		int bIndex = mb.bodyIndexByName(bName);
		rdt_.vector(bIndex, u);
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;

private:
	int rIndex_;
	tasks::RelativeDistTask rdt_;
};


class TASKS_DLLAPI VectorOrientationTask : public HighLevelTask
{
public:
	VectorOrientationTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		const std::string& bodyName, const Eigen::Vector3d& bodyVector,
		const Eigen::Vector3d& targetVector);

	tasks::VectorOrientationTask& task()
	{
		return vot_;
	}
	void bodyVector(const Eigen::Vector3d& vector)
	{
		vot_.bodyVector(vector);
	}
	const Eigen::Vector3d bodyVector() const
	{
		return vot_.bodyVector();
	}
	void target(const Eigen::Vector3d& vector)
	{
		vot_.target(vector);
	}
	const Eigen::Vector3d& target() const
	{
		return vot_.target();
	}
	const Eigen::Vector3d& actual() const
	{
		return vot_.actual();
	}

	virtual int dim() override;
	virtual void update(const std::vector<rbd::MultiBody> &mbs,
		const std::vector<rbd::MultiBodyConfig> &mbcs,
		const SolverData &data) override;

	virtual const Eigen::MatrixXd& jac() override;
	virtual const Eigen::VectorXd& eval() override;
	virtual const Eigen::VectorXd& speed() override;
	virtual const Eigen::VectorXd& normalAcc() override;
private:
	tasks::VectorOrientationTask vot_;
	int robotIndex_;
};

} // namespace qp

} // namespace tasks
