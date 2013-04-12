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
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// Tasks
#include "QL.h"


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
  boundConstr_(),
  tasks_(),
  A1_(),B1_(),A2_(),B2_(),
  XL_(),XU_(),
  Q_(),C_(),
  res_(),
  silent_(silent)
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

		A1_.block(count, 0, A1.rows(), data_.nrVars_) = A1;
		B1_.segment(count, A1.rows()) = B1;

		count += A1.rows();
	}

	count = 0;
	for(std::size_t i = 0; i < inEqConstr_.size(); ++i)
	{
		const Eigen::MatrixXd& A2 = inEqConstr_[i]->AInEq();
		const Eigen::VectorXd& B2 = inEqConstr_[i]->BInEq();

		A2_.block(count, 0, A2.rows(), data_.nrVars_) = A2;
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
		std::pair<int, int> b = tasks_[i]->begin();

		int r = Q.rows();
		int c = Q.cols();

		Q_.block(b.first, b.second, r, c) += tasks_[i]->weight()*Q;
		C_.segment(b.first, r) += tasks_[i]->weight()*C;
	}

	res_.setZero();
	bool success = false;
	double iter = 1e-8;
	while(!success && iter < 1e-3)
	{
		success = solveQP(A1_.cols(), A1_.rows(), A2_.rows(),
			Q_, C_, A1_, B1_, A2_, B2_, XL_, XU_, res_, iter, silent_);
		iter *= 10.;
	}

	if(success)
	{
		rbd::vectorToParam(res_.head(data_.alphaD_), mbc.alphaD);
		rbd::vectorToParam(res_.tail(data_.torque_), mbc.jointTorque);

		// don't write contact force to the structure since there are
		// to compute C vector.
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

	A1_.resize(nbEq, data_.nrVars_);
	B1_.resize(nbEq);
}


void QPSolver::updateInEqConstrSize()
{
	int nbInEq = 0;
	for(std::size_t i = 0; i < inEqConstr_.size(); ++i)
	{
		nbInEq += inEqConstr_[i]->nrInEqLine();
	}

	A2_.resize(nbInEq, data_.nrVars_);
	B2_.resize(nbInEq);
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
		data_.lambda_ += c.nrLambda()*c.points.size();
	}
	data_.lambdaUni_ = data_.lambda_;

	// counting bilateral contact
	for(const BilateralContact& c: data_.biCont_)
	{
		data_.lambda_ += c.nrLambda()*c.points.size();
	}
	data_.lambdaBi_ = data_.lambda_ - data_.lambdaUni_;

	data_.nrVars_ = data_.alphaD_ + data_.lambda_ + data_.torque_;

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
}


int QPSolver::nrVars() const
{
	return data_.nrVars_;
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
	return eqConstr_.size();
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
	return inEqConstr_.size();
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
	return boundConstr_.size();
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
	return constr_.size();
}


void QPSolver::addTask(Task* task)
{
	tasks_.push_back(task);
}


void QPSolver::removeTask(Task* task)
{
	tasks_.erase(std::find(tasks_.begin(), tasks_.end(), task));
}


int QPSolver::nrTasks() const
{
	return tasks_.size();
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


Eigen::VectorXd QPSolver::torqueVec() const
{
	return res_.tail(data_.torque_);
}


} // namespace qp

} // namespace tasks
