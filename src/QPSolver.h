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


struct FrictionCone
{
	FrictionCone(){}
	FrictionCone(Eigen::Matrix3d frame, int nrGen, double angle);

	std::vector<Eigen::Vector3d> generators;
};


struct Contact
{
	Contact(){}
	Contact(int bodyId, const std::vector<Eigen::Vector3d>& points,
		Eigen::Matrix3d frame, int nrGen, double angle);

	int bodyId;
	Eigen::Vector3d point;
	FrictionCone cone;
};


class Constraint
{
public:
	virtual ~Constraint() {}
	virtual void updateNrVars(const rbd::MultiBody& mb,
		int alphaD, int lambda, int torque, const std::vector<Contact>& cont) = 0;

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc) = 0;
};



class EqualityConstraint
{
public:
	virtual ~EqualityConstraint() {}
	virtual int nrEqLine() = 0;

	virtual const Eigen::MatrixXd& AEq() const = 0;
	virtual const Eigen::VectorXd& BEq() const = 0;
};



class InequalityConstraint
{
public:
	virtual ~InequalityConstraint() {}
	virtual int nrInEqLine() = 0;

	virtual const Eigen::MatrixXd& AInEq() const = 0;
	virtual const Eigen::VectorXd& BInEq() const = 0;
};



class BoundConstraint
{
public:
	virtual ~BoundConstraint() {}
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
	virtual ~Task() {}

	virtual double weight() const
	{
		return weight_;
	}

	virtual void weight(double w)
	{
		weight_ = w;
	}

	virtual std::pair<int, int> begin() const = 0;

	virtual void updateNrVars(const rbd::MultiBody& mb,
		int alphaD, int lambda, int torque, const std::vector<Contact>& cont) = 0;
	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc) = 0;

	virtual const Eigen::MatrixXd& Q() const = 0;
	virtual const Eigen::VectorXd& C() const = 0;

private:
	double weight_;
};



class HighLevelTask
{
public:
	virtual ~HighLevelTask() {}

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

	void nrVars(const rbd::MultiBody& mb, std::vector<Contact> cont);
	int nrVars() const;

	void addEqualityConstraint(EqualityConstraint* co);
	void removeEqualityConstraint(EqualityConstraint* co);
	int nrEqualityConstraints() const;

	void addInequalityConstraint(InequalityConstraint* co);
	void removeInequalityConstraint(InequalityConstraint* co);
	int nrInequalityConstraints() const;

	void addBoundConstraint(BoundConstraint* co);
	void removeBoundConstraint(BoundConstraint* co);
	int nrBoundConstraints() const;

	void addConstraint(Constraint* co);
	void removeConstraint(Constraint* co);
	int nrConstraints() const;

	void addTask(Task* task);
	void removeTask(Task* task);
	void resetTasks();
	int nrTasks() const;

	const Eigen::VectorXd& result() const;
	Eigen::VectorXd alphaDVec() const;
	Eigen::VectorXd lambdaVec() const;
	Eigen::VectorXd torqueVec() const;

private:
	std::vector<Constraint*> constr_;
	std::vector<EqualityConstraint*> eqConstr_;
	std::vector<InequalityConstraint*> inEqConstr_;
	std::vector<BoundConstraint*> boundConstr_;

	std::vector<Task*> tasks_;

	int alphaD_, lambda_, torque_;
	int nrVars_;
	std::vector<Contact> cont_;

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

