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
// Eigen
#include <Eigen/Core>

// Tasks
#include "Tasks.h"
#include "QPSolver.h"

// forward declaration
// RBDyn
namespace rbd
{
class MultiBody;
class MultiBodyConfig;
}

namespace tasks
{

namespace qp
{


class SetPointTask : public Task
{
public:
	SetPointTask(const rbd::MultiBody& mb, HighLevelTask* hlTask,
		double stiffness, double weight);

	double stiffness() const
	{
		return stiffness_;
	}

	void stiffness(double stiffness);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual const Eigen::MatrixXd& Q() const;
	virtual const Eigen::VectorXd& C() const;

private:
	HighLevelTask* hlTask_;

	double stiffness_, stiffnessSqrt_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
	Eigen::VectorXd alphaVec_;
};



class PostureTask : public Task
{
public:
	PostureTask(const rbd::MultiBody& mb, std::vector<std::vector<double>> q,
		double stiffness, double weight);

	tasks::PostureTask& task()
	{
		return pt_;
	}

	void posture(std::vector<std::vector<double>> q)
	{
		pt_.posture(q);
	}

	const std::vector<std::vector<double>> posture() const
	{
		return pt_.posture();
	}

	double stiffness() const
	{
		return stiffness_;
	}

	void stiffness(double stiffness);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual const Eigen::MatrixXd& Q() const;
	virtual const Eigen::VectorXd& C() const;

private:
	tasks::PostureTask pt_;

	double stiffness_;
	double stiffnessSqrt_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;
	Eigen::VectorXd alphaVec_;
};



class PositionTask : public HighLevelTask
{
public:
	PositionTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Vector3d& pos);

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

	virtual int dim();
	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual const Eigen::MatrixXd& jac();
	virtual const Eigen::MatrixXd& jacDot();
	virtual const Eigen::VectorXd& eval();

private:
	tasks::PositionTask pt_;
};



class OrientationTask : public HighLevelTask
{
public:
	OrientationTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Quaterniond& ori);
	OrientationTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Matrix3d& ori);

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

	virtual int dim();
	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual const Eigen::MatrixXd& jac();
	virtual const Eigen::MatrixXd& jacDot();
	virtual const Eigen::VectorXd& eval();

private:
	tasks::OrientationTask ot_;
};



class CoMTask : public HighLevelTask
{
public:
	CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com);

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

	virtual int dim();
	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual const Eigen::MatrixXd& jac();
	virtual const Eigen::MatrixXd& jacDot();
	virtual const Eigen::VectorXd& eval();

private:
	tasks::CoMTask ct_;
};

} // namespace qp

} // namespace tasks
