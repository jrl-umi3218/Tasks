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

// associated header
#include "QLDQPSolver.h"

// includes
// Tasks
#include "GenQPUtils.h"
#include "QPSolver.h"


namespace tasks
{

namespace qp
{


QLDQPSolver::QLDQPSolver():
	qld_(),
	Aeq_(),Aineq_(),
	beq_(), bineq_(),
	Aeq_full_(),Aineq_full_(),
	XL_(),XU_(),
	XL_full_(),XU_full_(),
	Q_(),C_(),
	Q_full_(),C_full_(),
	nrAeqLines_(0), nrAineqLines_(0)
{
}


void QLDQPSolver::updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq)
{
	int maxAeqLines = nrEq;
	int maxAineqLines = nrInEq + nrGenInEq*2;

	Aeq_full_.resize(maxAeqLines, nrVars);
	Aineq_full_.resize(maxAineqLines, nrVars);

	beq_.resize(maxAeqLines);
	bineq_.resize(maxAineqLines);

	XL_full_.resize(nrVars);
	XU_full_.resize(nrVars);

	Q_full_.resize(nrVars, nrVars);
	C_full_.resize(nrVars);

	if(dependencies_.size())
	{
		int reducedNrVars = nrVars - static_cast<int>(dependencies_.size());

		Aeq_.resize(maxAeqLines, reducedNrVars);
		Aineq_.resize(maxAineqLines, reducedNrVars);

		XL_.resize(reducedNrVars);
		XU_.resize(reducedNrVars);

		Q_.resize(reducedNrVars, reducedNrVars);
		C_.resize(reducedNrVars);

		X_full_.resize(nrVars);

		qld_.problem(reducedNrVars, maxAeqLines, maxAineqLines);
		std::cout << "(before) nrVars " << nrVars << std::endl;
		std::cout << "nrVars " << reducedNrVars << std::endl;
		std::cout << "maxAeqLines " << maxAeqLines << std::endl;
		std::cout << "maxAineqLines " << maxAineqLines << std::endl;
	}
	else
	{
		qld_.problem(nrVars, maxAeqLines, maxAineqLines);
	}

}


void QLDQPSolver::updateMatrix(
	const std::vector<Task*>& tasks,
	const std::vector<Equality*>& eqConstr,
	const std::vector<Inequality*>& inEqConstr,
	const std::vector<GenInequality*>& genInEqConstr,
	const std::vector<Bound*>& boundConstr)
{
	Aeq_full_.setZero();
	Aineq_full_.setZero();
	beq_.setZero();
	bineq_.setZero();
	XL_full_.fill(-std::numeric_limits<double>::infinity());
	XU_full_.fill(std::numeric_limits<double>::infinity());
	Q_full_.setZero();
	C_full_.setZero();

	const int nrVars = int(Q_full_.rows());

	nrAeqLines_ = 0;
	nrAeqLines_ = fillEq(eqConstr, nrVars, nrAeqLines_, Aeq_full_, beq_);
	nrAineqLines_ = 0;
	nrAineqLines_ = fillInEq(inEqConstr, nrVars, nrAineqLines_, Aineq_full_, bineq_);
	nrAineqLines_ = fillGenInEq(genInEqConstr, nrVars, nrAineqLines_, Aineq_full_, bineq_);

	fillBound(boundConstr, XL_full_, XU_full_);
	fillQC(tasks, nrVars, Q_full_, C_full_);
	if(dependencies_.size())
	{
		Aeq_.setZero();
		Aineq_.setZero();
		XL_.fill(-std::numeric_limits<double>::infinity());
		XU_.fill(std::numeric_limits<double>::infinity());
		Q_.setZero();
		C_.setZero();
		reduceA(Aeq_full_, Aeq_, full_to_reduced_, reduced_to_full_, dependencies_);
		reduceA(Aineq_full_, Aineq_, full_to_reduced_, reduced_to_full_, dependencies_);
		reduceBound(XL_full_, XL_, XU_full_, XU_, full_to_reduced_, reduced_to_full_, dependencies_);
		reduceQC(Q_full_, C_full_, Q_, C_, full_to_reduced_, reduced_to_full_, dependencies_);
	}
}


bool QLDQPSolver::solve()
{
	bool success = false;
	if(dependencies_.size())
	{
		success = qld_.solve(Q_, C_,
			Aeq_.block(0, 0, nrAeqLines_, int(Aeq_.cols())), beq_.segment(0, nrAeqLines_),
			Aineq_.block(0, 0, nrAineqLines_, int(Aineq_.cols())), bineq_.segment(0, nrAineqLines_),
			XL_, XU_, 1e-6);
		expandResult(qld_.result(), X_full_,
								 reduced_to_full_,
								 dependencies_);
	}
	else
	{
		success = qld_.solve(Q_full_, C_full_,
			Aeq_full_.block(0, 0, nrAeqLines_, int(Aeq_full_.cols())), beq_.segment(0, nrAeqLines_),
			Aineq_full_.block(0, 0, nrAineqLines_, int(Aineq_full_.cols())), bineq_.segment(0, nrAineqLines_),
			XL_full_, XU_full_, 1e-6);
	}
	return success;
}


const Eigen::VectorXd& QLDQPSolver::result() const
{
	if(dependencies_.size())
	{
		return X_full_;
	}
	else
	{
		return qld_.result();
	}
}


std::ostream& QLDQPSolver::errorMsg(
	const std::vector<rbd::MultiBody>& /* mbs */,
	const std::vector<Task*>& /* tasks */,
	const std::vector<Equality*>& /* eqConstr */,
	const std::vector<Inequality*>& /* inEqConstr */,
	const std::vector<GenInequality*>& /* genInEqConstr */,
	const std::vector<Bound*>& /* boundConstr */,
	std::ostream& out) const
{
	return out;
}


} // namespace qp

} // namespace tasks
