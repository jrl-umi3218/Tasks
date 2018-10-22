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
#include "LSSOLQPSolver.h"

// includes
// Tasks
#include "GenQPUtils.h"
#include "Tasks/QPSolver.h"


namespace tasks
{

namespace qp
{


LSSOLQPSolver::LSSOLQPSolver():
	lssol_(),
	A_(),AL_(),AU_(),
	AFull_(),
	XL_(),XU_(),
	XLFull_(),XUFull_(),
	Q_(),C_(),
	QFull_(),CFull_(),
	XFull_(),
	nrALines_(0)
{
	lssol_.warm(true);
	lssol_.feasibilityTol(1e-6);
	lssol_.persistance(true);
	lssol_.forceMinSize(false);
	lssol_.feasibilityMaxIter(4*lssol_.feasibilityMaxIter());
	lssol_.optimalityMaxIter(4*lssol_.optimalityMaxIter());
}


void LSSOLQPSolver::updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq)
{
	int maxALines = nrEq + nrInEq + nrGenInEq;
	AFull_.resize(maxALines, nrVars);
	AL_.resize(maxALines);
	AU_.resize(maxALines);

	XLFull_.resize(nrVars);
	XUFull_.resize(nrVars);

	QFull_.resize(nrVars, nrVars);
	CFull_.resize(nrVars);

	if(dependencies_.size())
	{
		int nrReducedVars = nrVars - static_cast<int>(dependencies_.size());

		A_.resize(maxALines, nrReducedVars);
		Q_.resize(nrReducedVars, nrReducedVars);
		C_.resize(nrReducedVars);

		XL_.resize(nrReducedVars);
		XU_.resize(nrReducedVars);
		XFull_.resize(nrVars);

		lssol_.resize(nrReducedVars, maxALines, Eigen::lssol::QP2);
	}
	else
	{
		lssol_.resize(nrVars, maxALines, Eigen::lssol::QP2);
	}
}


void LSSOLQPSolver::updateMatrix(
	const std::vector<Task*>& tasks,
	const std::vector<Equality*>& eqConstr,
	const std::vector<Inequality*>& inEqConstr,
	const std::vector<GenInequality*>& genInEqConstr,
	const std::vector<Bound*>& boundConstr)
{
	AFull_.setZero();
	AL_.setZero();
	AU_.setZero();
	XLFull_.fill(-std::numeric_limits<double>::infinity());
	XUFull_.fill(std::numeric_limits<double>::infinity());
	QFull_.setZero();
	CFull_.setZero();

	const int nrVars = int(QFull_.rows());

	nrALines_ = 0;
	nrALines_ = fillEq(eqConstr, nrVars, nrALines_, AFull_, AL_, AU_);
	nrALines_ = fillInEq(inEqConstr, nrVars, nrALines_, AFull_, AL_, AU_);
	nrALines_ = fillGenInEq(genInEqConstr, nrVars, nrALines_, AFull_, AL_, AU_);

	fillBound(boundConstr, XLFull_, XUFull_);
	fillQC(tasks, nrVars, QFull_, CFull_);

	if(dependencies_.size())
	{
		reduceA(AFull_, A_, fullToReduced_, reducedToFull_, dependencies_);
		reduceBound(XLFull_, XL_, XUFull_, XU_, fullToReduced_, reducedToFull_, dependencies_);
		reduceQC(QFull_, CFull_, Q_, C_, fullToReduced_, reducedToFull_, dependencies_);
	}
}


bool LSSOLQPSolver::solve()
{
	bool success = false;
	if(dependencies_.size())
	{
		success = lssol_.solve(XL_, XU_, Q_, C_,
													 A_.block(0, 0, nrALines_, A_.cols()),
													 AL_.segment(0, nrALines_), AU_.segment(0, nrALines_));
		expandResult(lssol_.result(), XFull_,
									reducedToFull_,
									dependencies_);
	}
	else
	{
		success = lssol_.solve(XLFull_, XUFull_, QFull_, CFull_,
													 AFull_.block(0, 0, nrALines_, AFull_.cols()),
													 AL_.segment(0, nrALines_), AU_.segment(0, nrALines_));
	}
	return success;
}


const Eigen::VectorXd& LSSOLQPSolver::result() const
{
	if(dependencies_.size())
	{
		return XFull_;
	}
	else
	{
		return lssol_.result();
	}
}


std::ostream& LSSOLQPSolver::errorMsg(
	const std::vector<rbd::MultiBody>& mbs,
	const std::vector<Task*>& /* tasks */,
	const std::vector<Equality*>& eqConstr,
	const std::vector<Inequality*>& inEqConstr,
	const std::vector<GenInequality*>& genInEqConstr,
	const std::vector<Bound*>& boundConstr,
	std::ostream& out) const
{
	const int nrVars = int(Q_.rows());

	out << "lssol output (" << lssol_.inform() << "): ";
	out << std::endl;
	lssol_.inform(out);

	const Eigen::VectorXi& istate = lssol_.istate();
	// check bound constraint
	for(int i = 0; i < nrVars; ++i)
	{
		if(istate(i) < 0)
		{
			for(Bound* b: boundConstr)
			{
				int start = b->beginVar();
				int end = start + int(b->Lower().rows());
				if(i >= start && i < end)
				{
					int line = i - start;
					out << b->nameBound() << " violated at line: " << line << std::endl;
					out << b->descBound(mbs, line) << std::endl;
					out << XL_(i) << " <= " << lssol_.result()(i) << " <= " << XU_(i) << std::endl;
					break;
				}
			}
		}
	}

	// check inequality constraint
	for(int i = 0; i < nrALines_; ++i)
	{
		int iInIstate = i + nrVars;
		if(istate(iInIstate) < 0)
		{
			int start = 0;
			int end = 0;

			constrErrorMsg(mbs, lssol_.result(), i, eqConstr, start, end, out);
			constrErrorMsg(mbs, lssol_.result(), i, inEqConstr, start, end, out);
			constrErrorMsg(mbs, lssol_.result(), i, genInEqConstr, start, end, out);
			out << std::endl;
		}
	}
	return out;
}


} // namespace qp

} // namespace tasks
