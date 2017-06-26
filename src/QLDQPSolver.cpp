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
#include "Tasks/QPSolver.h"


namespace tasks
{

namespace qp
{


QLDQPSolver::QLDQPSolver():
	qld_(),
	Aeq_(),Aineq_(),
	beq_(), bineq_(),
	AeqFull_(),AineqFull_(),
	XL_(),XU_(),
	XLFull_(),XUFull_(),
	Q_(),C_(),
	QFull_(),CFull_(),
	nrAeqLines_(0), nrAineqLines_(0)
{
}


void QLDQPSolver::updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq)
{
	int maxAeqLines = nrEq;
	int maxAineqLines = nrInEq + nrGenInEq*2;

	AeqFull_.resize(maxAeqLines, nrVars);
	AineqFull_.resize(maxAineqLines, nrVars);

	beq_.resize(maxAeqLines);
	bineq_.resize(maxAineqLines);

	XLFull_.resize(nrVars);
	XUFull_.resize(nrVars);

	QFull_.resize(nrVars, nrVars);
	CFull_.resize(nrVars);

	if(dependencies_.size())
	{
		int reducedNrVars = nrVars - static_cast<int>(dependencies_.size());

		Aeq_.resize(maxAeqLines, reducedNrVars);
		Aineq_.resize(maxAineqLines, reducedNrVars);

		XL_.resize(reducedNrVars);
		XU_.resize(reducedNrVars);

		Q_.resize(reducedNrVars, reducedNrVars);
		C_.resize(reducedNrVars);

		XFull_.resize(nrVars);

		qld_.problem(reducedNrVars, maxAeqLines, maxAineqLines);
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
	AeqFull_.setZero();
	AineqFull_.setZero();
	beq_.setZero();
	bineq_.setZero();
	XLFull_.fill(-std::numeric_limits<double>::infinity());
	XUFull_.fill(std::numeric_limits<double>::infinity());
	QFull_.setZero();
	CFull_.setZero();

	const int nrVars = int(QFull_.rows());

	nrAeqLines_ = 0;
	nrAeqLines_ = fillEq(eqConstr, nrVars, nrAeqLines_, AeqFull_, beq_);
	nrAineqLines_ = 0;
	nrAineqLines_ = fillInEq(inEqConstr, nrVars, nrAineqLines_, AineqFull_, bineq_);
	nrAineqLines_ = fillGenInEq(genInEqConstr, nrVars, nrAineqLines_, AineqFull_, bineq_);

	fillBound(boundConstr, XLFull_, XUFull_);
	fillQC(tasks, nrVars, QFull_, CFull_);
	if(dependencies_.size())
	{
		Aeq_.setZero();
		Aineq_.setZero();
		XL_.fill(-std::numeric_limits<double>::infinity());
		XU_.fill(std::numeric_limits<double>::infinity());
		Q_.setZero();
		C_.setZero();
		reduceA(AeqFull_, Aeq_, fullToReduced_, reducedToFull_, dependencies_);
		reduceA(AineqFull_, Aineq_, fullToReduced_, reducedToFull_, dependencies_);
		reduceBound(XLFull_, XL_, XUFull_, XU_, fullToReduced_, reducedToFull_, dependencies_);
		reduceQC(QFull_, CFull_, Q_, C_, fullToReduced_, reducedToFull_, dependencies_);
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
			XL_, XU_, false, 1e-6);
		expandResult(qld_.result(), XFull_,
								 reducedToFull_,
								 dependencies_);
	}
	else
	{
		success = qld_.solve(QFull_, CFull_,
			AeqFull_.block(0, 0, nrAeqLines_, int(AeqFull_.cols())), beq_.segment(0, nrAeqLines_),
			AineqFull_.block(0, 0, nrAineqLines_, int(AineqFull_.cols())), bineq_.segment(0, nrAineqLines_),
			XLFull_, XUFull_, false, 1e-6);
	}
	return success;
}


const Eigen::VectorXd& QLDQPSolver::result() const
{
	if(dependencies_.size())
	{
		return XFull_;
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
