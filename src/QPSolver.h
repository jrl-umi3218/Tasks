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

// Eigen
#include <EigenQP/QLD.h>
#include <EigenQP/LSSOL.h>

// Tasks
#include "QPSolverData.h"
#include "QPContacts.h"


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
class Constraint;
class Equality;
class Inequality;
class Bound;
class Task;



class QPSolver
{
public:
	QPSolver(bool silent=false);

	bool update(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc);
	bool updateQLD(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc);
	bool updateLSSOL(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc);

	void updateEqConstrSize();
	void updateInEqConstrSize();

	void nrVars(const rbd::MultiBody& mb,
		std::vector<UnilateralContact> uni,
		std::vector<BilateralContact> bi);
	int nrVars() const;

	void addEqualityConstraint(Equality* co);
	void removeEqualityConstraint(Equality* co);
	int nrEqualityConstraints() const;

	void addInequalityConstraint(Inequality* co);
	void removeInequalityConstraint(Inequality* co);
	int nrInequalityConstraints() const;

	void addBoundConstraint(Bound* co);
	void removeBoundConstraint(Bound* co);
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

	int contactLambdaPosition(int bodyId) const;

protected:
	void updateSolverSize(int nrVar, int nrEq, int nrIneq);
	void updateQLDSize(int nrVar, int nrEq, int nrIneq);
	void updateLSSOLSize(int nrVar, int nrEq, int nrIneq);

	void preUpdate(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc);
	void postUpdate(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc,
		bool success, const Eigen::VectorXd& result);

private:
	std::vector<Constraint*> constr_;
	std::vector<Equality*> eqConstr_;
	std::vector<Inequality*> inEqConstr_;
	std::vector<Bound*> boundConstr_;

	std::vector<Task*> tasks_;

	SolverData data_;

	int nrEq_;
	Eigen::MatrixXd A1_;
	Eigen::VectorXd B1_;

	int nrInEq_;
	Eigen::MatrixXd A2_;
	Eigen::VectorXd B2_;

	Eigen::VectorXd XL_;
	Eigen::VectorXd XU_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;

	Eigen::VectorXd res_;
	Eigen::VectorXd torqueRes_;

	Eigen::QLD qld_;
	Eigen::StdLSSOL lssol_;

	bool silent_;
};



class Constraint
{
public:
	virtual ~Constraint() {}
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data) = 0;

	virtual void update(const rbd::MultiBody& mb,
		const rbd::MultiBodyConfig& mbc) = 0;
};



template<typename... Fun>
class ConstraintFunction : public Constraint, public Fun...
{
public:
	virtual ~ConstraintFunction() {}

	void addToSolver(QPSolver& sol)
	{
		sol.addConstraint(this);
		addToSolver_p(sol, (&Fun::addToSolver)...);
	}

	void removeFromSolver(QPSolver& sol)
	{
		sol.removeConstraint(this);
		removeFromSolver_p(sol, (&Fun::removeFromSolver)...);
	}

private:
	void addToSolver_p(QPSolver& /* sol */)
	{}

	template<typename F, typename... NextFun>
	void addToSolver_p(QPSolver& sol, F function, NextFun... nFun)
	{
		(this->*function)(sol);
		addToSolver_p(sol, nFun...);
	}

	void removeFromSolver_p(QPSolver& /* sol */)
	{}

	template<typename F, typename... NextFun>
	void removeFromSolver_p(QPSolver& sol, F function, NextFun... nFun)
	{
		(this->*function)(sol);
		removeFromSolver_p(sol, nFun...);
	}
};



class Equality
{
public:
	virtual ~Equality() {}
	virtual int maxEq() const = 0;
	virtual int nrEq() const { return maxEq(); }

	virtual const Eigen::MatrixXd& AEq() const = 0;
	virtual const Eigen::VectorXd& BEq() const = 0;

	void addToSolver(QPSolver& sol)
	{
		sol.addEqualityConstraint(this);
	}

	void removeFromSolver(QPSolver& sol)
	{
		sol.removeEqualityConstraint(this);
	}
};



class Inequality
{
public:
	virtual ~Inequality() {}
	virtual int maxInEq() const = 0;
	virtual int nrInEq() const { return maxInEq(); }

	virtual const Eigen::MatrixXd& AInEq() const = 0;
	virtual const Eigen::VectorXd& BInEq() const = 0;

	void addToSolver(QPSolver& sol)
	{
		sol.addInequalityConstraint(this);
	}

	void removeFromSolver(QPSolver& sol)
	{
		sol.removeInequalityConstraint(this);
	}
};



class Bound
{
public:
	virtual ~Bound() {}
	virtual int beginVar() const = 0;

	virtual const Eigen::VectorXd& Lower() const = 0;
	virtual const Eigen::VectorXd& Upper() const = 0;

	void addToSolver(QPSolver& sol)
	{
		sol.addBoundConstraint(this);
	}

	void removeFromSolver(QPSolver& sol)
	{
		sol.removeBoundConstraint(this);
	}
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
		const SolverData& data) = 0;
	virtual void update(const rbd::MultiBody& mb,
		const rbd::MultiBodyConfig& mbc) = 0;

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

	virtual void update(const rbd::MultiBody& mb,
		const rbd::MultiBodyConfig& mbc) = 0;

	virtual const Eigen::MatrixXd& jac() = 0;
	virtual const Eigen::MatrixXd& jacDot() = 0;

	virtual const Eigen::VectorXd& eval() = 0;
};



} // namespace qp

} // namespace tasks

