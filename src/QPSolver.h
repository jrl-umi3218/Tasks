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
// std
#include <vector>

// Eigen
#include <Eigen/Core>


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

class SequenceController3;

class Constraint
{
public:
	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc) = 0;
};



class EqualityConstraint : public Constraint
{
public:
	virtual int nrEqLine() = 0;

	virtual const Eigen::MatrixXd& AEq() const = 0;
	virtual const Eigen::VectorXd& BEq() const = 0;
};



class InequalityConstraint : public Constraint
{
public:
	virtual int nrInEqLine() = 0;

	virtual const Eigen::MatrixXd& AInEq() const = 0;
	virtual const Eigen::VectorXd& BInEq() const = 0;
};



class BoundConstraint : public Constraint
{
public:
	virtual int beginVar() = 0;

	virtual const Eigen::VectorXd& Lower() const = 0;
	virtual const Eigen::VectorXd& Upper() const = 0;
};



class Task
{
public:
	Task(double weight):
		weight_(weight)
	{}

	double weight() const
	{
		return weight_;
	}

	void weight(double w)
	{
		weight_ = w;
	}

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc) = 0;

	virtual const Eigen::MatrixXd& Q() const = 0;
	virtual const Eigen::VectorXd& C() const = 0;

private:
	double weight_;
};



class HighLevelTask
{
public:
	virtual int dim() = 0;

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc) = 0;

	virtual const Eigen::MatrixXd& jac() = 0;
	virtual const Eigen::MatrixXd& jacDot() = 0;

	virtual const Eigen::VectorXd& eval() = 0;
};



class QPSolver
{
public:
	QPSolver();

	bool update(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc);

	void updateEqConstrSize();
	void updateInEqConstrSize();

	void nrVars(int nrVars);
	int nrVars() const;

	void addEqualityConstraint(EqualityConstraint* co);
	void addInequalityConstraint(InequalityConstraint* co);
	void addBoundConstraint(BoundConstraint* co);

	void addTask(Task* task);
	void removeTask(Task* task);
	void resetTasks();
	int nrTasks() const;

private:
	std::vector<Constraint*> constr_;
	std::vector<EqualityConstraint*> eqConstr_;
	std::vector<InequalityConstraint*> inEqConstr_;
	std::vector<BoundConstraint*> boundConstr_;

	std::vector<Task*> tasks_;

	int nrVars_;

	Eigen::MatrixXd A1_;
	Eigen::VectorXd B1_;

	Eigen::MatrixXd A2_;
	Eigen::VectorXd B2_;

	Eigen::VectorXd XL_;
	Eigen::VectorXd XU_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;

	Eigen::VectorXd res_;
};

} // namespace qp

} // namespace tasks

