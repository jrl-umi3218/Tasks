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
#include <iostream>
#include <limits>
#include <cmath>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>


namespace tasks
{

namespace qp
{


// Value add to the diagonal to ensure positive matrix
static const double DIAG_CONSTANT = 1e-5;


/**
	*													QPSolver
	*/



QPSolver::QPSolver(bool silent):
	constr_(),
	inEqConstr_(),
	boundConstr_(),
	tasks_(),
	nrALine_(0),
	A_(),AL_(),AU_(),
	XL_(),XU_(),
	Q_(),C_(),
	res_(),
	silent_(silent)
{
	lssol_.warm(false);
	lssol_.feasibilityTol(1e-6);
}


bool QPSolver::solve(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
	return solveLSSOL(mb, mbc);
}


bool QPSolver::solveLSSOL(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
	preUpdate(mb, mbc);

	bool success = lssol_.solve(Q_, C_,
		A_.block(0, 0, nrALine_, data_.nrVars_), int(A_.rows()),
		AL_.segment(0, nrALine_), AU_.segment(0, nrALine_), XL_, XU_);

	if(!success)
	{
		std::cerr << "lssol output: " << lssol_.fail() << std::endl;
		const Eigen::VectorXi& istate = lssol_.istate();
		// check bound constraint
		for(int i = 0; i < data_.nrVars_; ++i)
		{
			if(istate(i) < 0)
			{
				for(Bound* b: boundConstr_)
				{
					int start = b->beginVar();
					int end = start + int(b->Lower().rows());
					if(i >= start && i < end)
					{
						int line = i - start;
						std::cerr << b->nameBound() << " violated at line: " << line << std::endl;
						std::cerr << b->descBound(mb, line) << std::endl;
						std::cerr << XL_(i) << " <= " << lssol_.result()(i) << " <= " << XU_(i) << std::endl;
						break;
					}
				}
			}
		}

		// check inequality constraint
		for(int i = 0; i < nrALine_; ++i)
		{
			int iInIstate = i + data_.nrVars_;
			if(istate(iInIstate) < 0)
			{
				int start = 0;
				int end = 0;
				for(Inequality* ie: inEqConstr_)
				{
					end += ie->nrInEq();
					if(i >= start && i < end)
					{
						int line = i - start;
						std::cerr << ie->nameInEq() << " violated at line: " << line << std::endl;
						std::cerr << ie->descInEq(mb, line) << std::endl;
						std::cerr << ie->LowerInEq()(line) << " <= " <<
												 ie->AInEq().row(line)*lssol_.result() <<" <= " <<
												 ie->UpperInEq()(line) << std::endl;
						break;
					}
					start = end;
				}
			}
		}
	}

	postUpdate(mb, mbc, success, lssol_.result());

	return success;
}


void QPSolver::updateConstrSize()
{
	int maxALine = 0;
	for(std::size_t i = 0; i < inEqConstr_.size(); ++i)
	{
		maxALine += inEqConstr_[i]->maxInEq();
	}

	nrALine_ = 0;
	A_.resize(maxALine, data_.nrVars_);
	AL_.resize(maxALine);
	AU_.resize(maxALine);

	updateSolverSize(data_.nrVars_, maxALine);
}


void QPSolver::nrVars(const rbd::MultiBody& mb,
	std::vector<UnilateralContact> uni,
	std::vector<BilateralContact> bi)
{
	data_.alphaD_ = mb.nrDof();
	data_.lambda_ = 0;
	data_.torque_ = (mb.nrDof() - mb.joint(0).dof());
	data_.uniCont_ = uni;
	data_.biCont_ = bi;

	// counting unilateral contact
	for(const UnilateralContact& c: data_.uniCont_)
	{
		for(std::size_t i = 0; i < c.points.size(); ++i)
		{
			data_.lambda_ += c.nrLambda(int(i));
		}
	}
	data_.lambdaUni_ = data_.lambda_;

	// counting bilateral contact
	for(const BilateralContact& c: data_.biCont_)
	{
		for(std::size_t i = 0; i < c.points.size(); ++i)
		{
			data_.lambda_ += c.nrLambda(int(i));
		}
	}
	data_.lambdaBi_ = data_.lambda_ - data_.lambdaUni_;

	data_.nrVars_ = data_.alphaD_ + data_.lambda_;

	if(XL_.rows() != data_.nrVars_)
	{
		XL_.resize(data_.nrVars_);
		XU_.resize(data_.nrVars_);

		Q_.resize(data_.nrVars_, data_.nrVars_);
		C_.resize(data_.nrVars_);

		res_.resize(data_.nrVars_);
	}

	for(Task* t: tasks_)
	{
		t->updateNrVars(mb, data_);
	}

	for(Constraint* c: constr_)
	{
		c->updateNrVars(mb, data_);
	}

	updateSolverSize(data_.nrVars_, int(A_.rows()));
}


int QPSolver::nrVars() const
{
	return data_.nrVars_;
}


void QPSolver::addInequalityConstraint(Inequality* co)
{
	inEqConstr_.push_back(co);
}


void QPSolver::removeInequalityConstraint(Inequality* co)
{
	inEqConstr_.erase(std::find(inEqConstr_.begin(), inEqConstr_.end(), co));
}


int QPSolver::nrInequalityConstraints() const
{
	return static_cast<int>(inEqConstr_.size());
}


void QPSolver::addBoundConstraint(Bound* co)
{
	boundConstr_.push_back(co);
}


void QPSolver::removeBoundConstraint(Bound* co)
{
	boundConstr_.erase(std::find(boundConstr_.begin(), boundConstr_.end(), co));
}


int QPSolver::nrBoundConstraints() const
{
	return static_cast<int>(boundConstr_.size());
}


void QPSolver::addConstraint(Constraint* co)
{
	if(std::find(constr_.begin(), constr_.end(), co) == constr_.end())
	{
		constr_.push_back(co);
	}
}


void QPSolver::removeConstraint(Constraint* co)
{
	auto it = std::find(constr_.begin(), constr_.end(), co);
	if(it != constr_.end())
	{
		constr_.erase(it);
	}
}

int QPSolver::nrConstraints() const
{
	return static_cast<int>(constr_.size());
}


void QPSolver::addTask(Task* task)
{
	if(std::find(tasks_.begin(), tasks_.end(), task) == tasks_.end())
	{
		tasks_.push_back(task);
	}
}


void QPSolver::removeTask(Task* task)
{
	auto it = std::find(tasks_.begin(), tasks_.end(), task);
	if(it != tasks_.end())
	{
		tasks_.erase(it);
	}
}


int QPSolver::nrTasks() const
{
	return static_cast<int>(tasks_.size());
}


void QPSolver::resetTasks()
{
	tasks_.clear();
}


const Eigen::VectorXd& QPSolver::result() const
{
	return res_;
}


Eigen::VectorXd QPSolver::alphaDVec() const
{
	return res_.head(data_.alphaD_);
}


Eigen::VectorXd QPSolver::lambdaVec() const
{
	return res_.segment(data_.alphaD_, data_.lambda_);
}


int QPSolver::contactLambdaPosition(int bodyId) const
{
	int pos = 0;
	for(const UnilateralContact& uc: data_.unilateralContacts())
	{
		if(uc.bodyId == bodyId)
		{
			return pos;
		}

		for(std::size_t i = 0; i < uc.points.size(); ++i)
		{
			pos += uc.nrLambda(int(i));
		}
	}

	for(const BilateralContact& bc: data_.bilateralContacts())
	{
		if(bc.bodyId == bodyId)
		{
			return pos;
		}

		for(std::size_t i = 0; i < bc.points.size(); ++i)
		{
			pos += bc.nrLambda(int(i));
		}
	}

	return -1;
}


void QPSolver::updateSolverSize(int nrVar, int nrConstr)
{
	updateLSSOLSize(nrVar, nrConstr);
}


void QPSolver::updateLSSOLSize(int nrVar, int nrConstr)
{
	lssol_.problem(nrVar, nrConstr);
}


void QPSolver::preUpdate(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
	for(std::size_t i = 0; i < constr_.size(); ++i)
	{
		constr_[i]->update(mb, mbc);
	}

	for(std::size_t i = 0; i < tasks_.size(); ++i)
	{
		tasks_[i]->update(mb, mbc);
	}

	A_.setZero();
	AL_.setZero();
	AU_.setZero();
	XL_.fill(-std::numeric_limits<double>::infinity());
	XU_.fill(std::numeric_limits<double>::infinity());
	Q_.setZero();
	C_.setZero();

	nrALine_ = 0;
	for(std::size_t i = 0; i < inEqConstr_.size(); ++i)
	{
		// ineq constraint can return a matrix with more line
		// than the number of constraint
		int nrConstr = inEqConstr_[i]->nrInEq();
		const Eigen::MatrixXd& A = inEqConstr_[i]->AInEq();
		const Eigen::VectorXd& AL = inEqConstr_[i]->LowerInEq();
		const Eigen::VectorXd& AU = inEqConstr_[i]->UpperInEq();

		A_.block(nrALine_, 0, nrConstr, data_.nrVars_) =
			A.block(0, 0, nrConstr, data_.nrVars_);
		AL_.segment(nrALine_, nrConstr) = AL.head(nrConstr);
		AU_.segment(nrALine_, nrConstr) = AU.head(nrConstr);

		nrALine_ += nrConstr;
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
		std::pair<int, int> b = tasks_[i]->begin();

		int r = static_cast<int>(Q.rows());
		int c = static_cast<int>(Q.cols());

		Q_.block(b.first, b.second, r, c) += tasks_[i]->weight()*Q;
		C_.segment(b.first, r) += tasks_[i]->weight()*C;
	}

	// try to transform Q_ to a positive matrix
	// we just add a small value to the diagonal since
	// the first necessary condition is to have
	// Q_(i,i) > 0
	// may be we can try to check the second
	// condition in a near futur
	// Q_(i,i) + Q_(j,j) > 2·Q_(i,j) for i≠j
	for(int i = 0; i < data_.nrVars_; ++i)
	{
		if(std::abs(Q_(i, i)) < DIAG_CONSTANT)
		{
			Q_(i, i) += DIAG_CONSTANT;
		}
	}
}


void QPSolver::postUpdate(const rbd::MultiBody& /* mb */,
	rbd::MultiBodyConfig& mbc, bool success, const Eigen::VectorXd& result)
{
	if(success)
	{
		res_ = result;

		rbd::vectorToParam(res_.head(data_.alphaD_), mbc.alphaD);

		// don't write contact force to the structure since contact force are used
		// to compute C vector.
	}
}


} // namespace qp

} // namespace tasks
