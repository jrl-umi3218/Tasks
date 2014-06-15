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

	/*! \brief solve the problem
	 *  \param mb current multibody
	 *  \param mbc result of the solving problem
	 */
	bool solve(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc);
	bool solveLSSOL(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc);

	void updateConstrSize();

	void nrVars(const rbd::MultiBody& mb,
		std::vector<UnilateralContact> uni,
		std::vector<BilateralContact> bi);
	int nrVars() const;

	/// call updateNrVars on all tasks
	void updateTasksNrVars(const rbd::MultiBody& mb) const;
	/// call updateNrVars on all constraints
	void updateConstrsNrVars(const rbd::MultiBody& mb) const;
	/// call updateNrVars on all tasks and constraints
	void updateNrVars(const rbd::MultiBody& mb) const;

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

	const SolverData& data() const;

	const Eigen::VectorXd& result() const;
	Eigen::VectorXd alphaDVec() const;
	Eigen::VectorXd lambdaVec() const;

	int contactLambdaPosition(int bodyId) const;

protected:
	void updateSolverSize(int nrVar, int nrConstr);
	void updateLSSOLSize(int nrVar, int nrConstr);

	void preUpdate(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc);
	void postUpdate(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc,
		bool success, const Eigen::VectorXd& result);

private:
	std::vector<Constraint*> constr_;
	std::vector<Inequality*> inEqConstr_;
	std::vector<Bound*> boundConstr_;

	std::vector<Task*> tasks_;

	SolverData data_;

	int nrALine_;
	Eigen::MatrixXd A_;
	Eigen::VectorXd AL_, AU_;

	Eigen::VectorXd XL_;
	Eigen::VectorXd XU_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;

	Eigen::VectorXd res_;

	Eigen::LSSOL lssol_;

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



class Inequality
{
public:
	virtual ~Inequality() {}
	virtual int maxInEq() const = 0;
	virtual int nrInEq() const { return maxInEq(); }

	virtual const Eigen::MatrixXd& AInEq() const = 0;
	virtual const Eigen::VectorXd& LowerInEq() const = 0;
	virtual const Eigen::VectorXd& UpperInEq() const = 0;

	virtual std::string nameInEq() const = 0;
	virtual std::string descInEq(const rbd::MultiBody& mb, int i) = 0;

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

	virtual std::string nameBound() const = 0;
	virtual std::string descBound(const rbd::MultiBody& mb, int i) = 0;

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
		const rbd::MultiBodyConfig& mbc,
		const SolverData& data) = 0;

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
		const rbd::MultiBodyConfig& mbc,
		const SolverData& data) = 0;

	virtual const Eigen::MatrixXd& jac() = 0;
	virtual const Eigen::VectorXd& eval() = 0;
	virtual const Eigen::VectorXd& speed() = 0;
	virtual const Eigen::VectorXd& normalAcc() = 0;
};



} // namespace qp

} // namespace tasks

