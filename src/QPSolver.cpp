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
#include "QPSolver.h"

// includes
// std
#include <limits>
#include <cmath>

// RBDyn
#include <MultiBody.h>
#include <MultiBodyConfig.h>

// Tasks
#include "QL.h"


namespace tasks
{

namespace qp
{


/**
	*													QPSolver
	*/


QPSolver::QPSolver():
  constr_(),
  eqConstr_(),
  inEqConstr_(),
  boundConstr_(),
  tasks_(),
  nrVars_(0),
  A1_(),B1_(),A2_(),B2_(),
  XL_(),XU_(),
  Q_(),C_(),
  res_()
{ }


bool QPSolver::update(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
	for(std::size_t i = 0; i < constr_.size(); ++i)
	{
		constr_[i]->update(mb, mbc);
	}

	for(std::size_t i = 0; i < tasks_.size(); ++i)
	{
		tasks_[i]->update(mb, mbc);
	}

	A1_.setZero();
	B1_.setZero();
	A2_.setZero();
	B2_.setZero();
	XL_.fill(-std::numeric_limits<double>::infinity());
	XU_.fill(std::numeric_limits<double>::infinity());
	Q_.setZero();
	C_.setZero();

	int count = 0;
	for(std::size_t i = 0; i < eqConstr_.size(); ++i)
	{
		const Eigen::MatrixXd& A1 = eqConstr_[i]->AEq();
		const Eigen::VectorXd& B1 = eqConstr_[i]->BEq();

		A1_.block(count, 0, A1.rows(), nrVars_) = A1;
		B1_.segment(count, A1.rows()) = B1;

		count += A1.rows();
	}

	count = 0;
	for(std::size_t i = 0; i < inEqConstr_.size(); ++i)
	{
		const Eigen::MatrixXd& A2 = inEqConstr_[i]->AInEq();
		const Eigen::VectorXd& B2 = inEqConstr_[i]->BInEq();

		A2_.block(count, 0, A2.rows(), nrVars_) = A2;
		B2_.segment(count, A2.rows()) = B2;

		count += A2.rows();
	}

	for(std::size_t i = 0; i < boundConstr_.size(); ++i)
	{
		const Eigen::VectorXd& XL = boundConstr_[i]->Lower();
		const Eigen::VectorXd& XU = boundConstr_[i]->Upper();
		int bv = boundConstr_[i]->beginVar();

		XL_.segment(bv, XL.size()) = XL;
		XU_.segment(bv, XU.size()) = XU;
	}

	for(std::size_t i = 0; i < tasks_.size(); ++i)
	{
		const Eigen::MatrixXd& Q = tasks_[i]->Q();
		const Eigen::VectorXd& C = tasks_[i]->C();
		int a = Q.rows();
		Q_.block(0, 0, a, a) += tasks_[i]->weight()*Q;
		C_.segment(0, a) += tasks_[i]->weight()*C;
	}

	res_.setZero();
	bool success = false;
	double iter = 1e-8;
	while(!success && iter < 1e-3)
	{
		success = solveQP(A1_.cols(), A1_.rows(), A2_.rows(),
			Q_, C_, A1_, B1_, A2_, B2_, XL_, XU_, res_, iter);
		iter *= 10.;
	}

	if(success)
	{
		rbd::vectorToParam(res_.head(mb.nrDof()), mbc.alphaD);
	}

	return success;
}


void QPSolver::updateEqConstrSize()
{
	int nbEq = 0;
	for(std::size_t i = 0; i < eqConstr_.size(); ++i)
	{
		nbEq += eqConstr_[i]->nrEqLine();
	}

	A1_.resize(nbEq, nrVars_);
	B1_.resize(nbEq);
}


void QPSolver::updateInEqConstrSize()
{
	int nbInEq = 0;
	for(std::size_t i = 0; i < inEqConstr_.size(); ++i)
	{
		nbInEq += inEqConstr_[i]->nrInEqLine();
	}

	A2_.resize(nbInEq, nrVars_);
	B2_.resize(nbInEq);
}


void QPSolver::nrVars(int nrVars)
{
	nrVars_ = nrVars;

	if(XL_.rows() != nrVars_)
	{
		XL_.resize(nrVars_);
		XU_.resize(nrVars_);

		Q_.resize(nrVars_, nrVars_);
		C_.resize(nrVars_);

		res_.resize(nrVars_);
	}
}


int QPSolver::nrVars() const
{
	return nrVars_;
}


void QPSolver::addEqualityConstraint(EqualityConstraint* co)
{
	eqConstr_.push_back(co);
	constr_.push_back(co);
}


void QPSolver::addInequalityConstraint(InequalityConstraint* co)
{
	inEqConstr_.push_back(co);
	constr_.push_back(co);
}


void QPSolver::addBoundConstraint(BoundConstraint* co)
{
	boundConstr_.push_back(co);
	constr_.push_back(co);
}


void QPSolver::addTask(Task* task)
{
	tasks_.push_back(task);
}


void QPSolver::removeTask(Task* task)
{
	std::vector<Task*>::iterator it;
	for(it = tasks_.begin(); it != tasks_.end(); ++it)
	{
		if(*it == task)
		{
			tasks_.erase(it);
			break;
		}
	}
}


int QPSolver::nrTasks() const
{
	return tasks_.size();
}


void QPSolver::resetTasks()
{
	tasks_.clear();
}

} // namespace qp

} // namespace tasks
