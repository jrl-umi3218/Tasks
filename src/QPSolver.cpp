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

// Tasks
#include "GenQPSolver.h"


namespace tasks
{

namespace qp
{



/**
	*													QPSolver
	*/



QPSolver::QPSolver(bool silent):
	constr_(),
	eqConstr_(),
	inEqConstr_(),
	genInEqConstr_(),
	boundConstr_(),
	tasks_(),
	maxEqLines_(0),
	maxInEqLines_(0),
	maxGenInEqLines_(0),
	solver_(createQPSolver("LSSOL")),
	silent_(silent)
{
}


// must declare it in cpp because of GenQPSolver fwd declarition
QPSolver::~QPSolver()
{}


bool QPSolver::solve(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
	preUpdate(mb, mbc);

	bool success = solver_->solve();

	/*
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
	*/

	postUpdate(mb, mbc, success);

	return success;
}


template <typename T>
int accumMaxLines(int acc, T* constr)
{
	return acc + constr_traits<T>::maxLines(constr);
}


void QPSolver::updateConstrSize()
{
	maxEqLines_ = std::accumulate(eqConstr_.begin(), eqConstr_.end(), 0,
		accumMaxLines<Equality>);
	maxInEqLines_ = std::accumulate(inEqConstr_.begin(), inEqConstr_.end(), 0,
		accumMaxLines<Inequality>);
	maxGenInEqLines_ = std::accumulate(genInEqConstr_.begin(), genInEqConstr_.end(),
		0, accumMaxLines<GenInequality>);

	solver_->updateSize(data_.nrVars_, maxEqLines_, maxInEqLines_, maxGenInEqLines_);
}


void QPSolver::nrVars(const rbd::MultiBody& mb,
	std::vector<UnilateralContact> uni,
	std::vector<BilateralContact> bi)
{
	data_.normalAccB_.resize(mb.nrBodies());
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

	for(Task* t: tasks_)
	{
		t->updateNrVars(mb, data_);
	}

	for(Constraint* c: constr_)
	{
		c->updateNrVars(mb, data_);
	}

	solver_->updateSize(data_.nrVars_, maxEqLines_, maxInEqLines_, maxGenInEqLines_);
}


int QPSolver::nrVars() const
{
	return data_.nrVars_;
}


void QPSolver::updateTasksNrVars(const rbd::MultiBody& mb) const
{
	for(Task* t: tasks_)
	{
		t->updateNrVars(mb, data_);
	}
}


void QPSolver::updateConstrsNrVars(const rbd::MultiBody& mb) const
{
	for(Constraint* c: constr_)
	{
		c->updateNrVars(mb, data_);
	}
}


void QPSolver::updateNrVars(const rbd::MultiBody& mb) const
{
	updateTasksNrVars(mb);
	updateConstrsNrVars(mb);
}


void QPSolver::addEqualityConstraint(Equality* co)
{
	eqConstr_.push_back(co);
}


void QPSolver::removeEqualityConstraint(Equality* co)
{
	eqConstr_.erase(std::find(eqConstr_.begin(), eqConstr_.end(), co));
}


int QPSolver::nrEqualityConstraints() const
{
	return static_cast<int>(eqConstr_.size());
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


void QPSolver::addGenInequalityConstraint(GenInequality* co)
{
	genInEqConstr_.push_back(co);
}


void QPSolver::removeGenInequalityConstraint(GenInequality* co)
{
	genInEqConstr_.erase(std::find(genInEqConstr_.begin(), genInEqConstr_.end(), co));
}


int QPSolver::nrGenInequalityConstraints() const
{
	return static_cast<int>(genInEqConstr_.size());
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


const SolverData& QPSolver::data() const
{
	return data_;
}


SolverData& QPSolver::data()
{
	return data_;
}


const Eigen::VectorXd& QPSolver::result() const
{
	return solver_->result();
}


Eigen::VectorXd QPSolver::alphaDVec() const
{
	return solver_->result().head(data_.alphaD_);
}


Eigen::VectorXd QPSolver::lambdaVec() const
{
	return solver_->result().segment(data_.alphaD_, data_.lambda_);
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


void QPSolver::preUpdate(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
	data_.computeNormalAccB(mb, mbc);
	for(std::size_t i = 0; i < constr_.size(); ++i)
	{
		constr_[i]->update(mb, mbc, data_);
	}

	for(std::size_t i = 0; i < tasks_.size(); ++i)
	{
		tasks_[i]->update(mb, mbc, data_);
	}

	solver_->updateMatrix(tasks_, eqConstr_, inEqConstr_, genInEqConstr_,
		boundConstr_);
}


void QPSolver::postUpdate(const rbd::MultiBody& /* mb */,
	rbd::MultiBodyConfig& mbc, bool success)
{
	if(success)
	{
		rbd::vectorToParam(solver_->result().head(data_.alphaD_), mbc.alphaD);

		// don't write contact force to the structure since contact force are used
		// to compute C vector.
	}
}


} // namespace qp

} // namespace tasks
