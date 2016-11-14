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
#include "QPSolver.h"


namespace tasks
{

namespace qp
{


LSSOLQPSolver::LSSOLQPSolver():
	lssol_(),
	A_(),AL_(),AU_(),
	A_full_(),
	XL_(),XU_(),
	XL_full_(),XU_full_(),
	Q_(),C_(),
	Q_full_(),C_full_(),
	X_full_(),
	nrALines_(0)
{
	lssol_.warm(true);
	lssol_.feasibilityTol(1e-6);
}


void LSSOLQPSolver::updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq)
{
	int maxALines = nrEq + nrInEq + nrGenInEq;
	A_full_.resize(maxALines, nrVars);
	AL_.resize(maxALines);
	AU_.resize(maxALines);

	XL_full_.resize(nrVars);
	XU_full_.resize(nrVars);

	Q_full_.resize(nrVars, nrVars);
	C_full_.resize(nrVars);

	if(dependencies_.size())
	{
		int nrReducedVars = nrVars - static_cast<int>(dependencies_.size());

		A_.resize(maxALines, nrReducedVars);
		Q_.resize(nrReducedVars, nrReducedVars);
		C_.resize(nrReducedVars);

		XL_.resize(nrReducedVars);
		XU_.resize(nrReducedVars);
		X_full_.resize(nrVars);

		lssol_.problem(nrReducedVars, maxALines);
	}
	else
	{
		lssol_.problem(nrVars, maxALines);
	}
}


void LSSOLQPSolver::updateMatrix(
	const std::vector<Task*>& tasks,
	const std::vector<Equality*>& eqConstr,
	const std::vector<Inequality*>& inEqConstr,
	const std::vector<GenInequality*>& genInEqConstr,
	const std::vector<Bound*>& boundConstr)
{
	A_full_.setZero();
	AL_.setZero();
	AU_.setZero();
	XL_full_.fill(-std::numeric_limits<double>::infinity());
	XU_full_.fill(std::numeric_limits<double>::infinity());
	Q_full_.setZero();
	C_full_.setZero();

	const int nrVars = int(Q_full_.rows());

	nrALines_ = 0;
	nrALines_ = fillEq(eqConstr, nrVars, nrALines_, A_full_, AL_, AU_);
	nrALines_ = fillInEq(inEqConstr, nrVars, nrALines_, A_full_, AL_, AU_);
	nrALines_ = fillGenInEq(genInEqConstr, nrVars, nrALines_, A_full_, AL_, AU_);

	fillBound(boundConstr, XL_full_, XU_full_);
	fillQC(tasks, nrVars, Q_full_, C_full_);

	if(dependencies_.size())
	{
		reduceA(A_full_, A_, full_to_reduced_, reduced_to_full_, dependencies_);
		reduceBound(XL_full_, XL_, XU_full_, XU_, full_to_reduced_, reduced_to_full_, dependencies_);
		reduceQC(Q_full_, C_full_, Q_, C_, full_to_reduced_, reduced_to_full_, dependencies_);
	}
}


bool LSSOLQPSolver::solve()
{
	bool success = false;
	if(dependencies_.size())
	{
		success = lssol_.solve(Q_, C_,
			A_.block(0, 0, nrALines_, int(A_.cols())), int(A_.rows()),
			AL_.segment(0, nrALines_), AU_.segment(0, nrALines_), XL_, XU_);
		expandResult(lssol_.result(), X_full_,
									reduced_to_full_,
									dependencies_);
	}
	else
	{
		success = lssol_.solve(Q_full_, C_full_,
			A_full_.block(0, 0, nrALines_, int(A_full_.cols())), int(A_full_.rows()),
			AL_.segment(0, nrALines_), AU_.segment(0, nrALines_), XL_full_, XU_full_);
	}
	return success;
}


const Eigen::VectorXd& LSSOLQPSolver::result() const
{
	if(dependencies_.size())
	{
		return X_full_;
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

	out << "lssol output: " << lssol_.fail() << std::endl;
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
