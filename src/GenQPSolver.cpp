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
#include "GenQPSolver.h"

// includes
// tasks
#include "QPSolver.h"


namespace tasks
{

namespace qp
{



/**
	*													Utility functions
	*/



// Value add to the diagonal to ensure positive matrix
static const double DIAG_CONSTANT = 1e-4;


void fillQC(const std::vector<Task*>& tasks, int nrVars,
	Eigen::MatrixXd& Q, Eigen::VectorXd& C)
{
	for(std::size_t i = 0; i < tasks.size(); ++i)
	{
		const Eigen::MatrixXd& Qi = tasks[i]->Q();
		const Eigen::VectorXd& Ci = tasks[i]->C();
		std::pair<int, int> b = tasks[i]->begin();

		int r = static_cast<int>(Qi.rows());
		int c = static_cast<int>(Qi.cols());

		Q.block(b.first, b.second, r, c) += tasks[i]->weight()*Qi;
		C.segment(b.first, r) += tasks[i]->weight()*Ci;
	}

	// try to transform Q_ to a positive matrix
	// we just add a small value to the diagonal since
	// the first necessary condition is to have
	// Q_(i,i) > 0
	// may be we can try to check the second
	// condition in a near futur
	// Q_(i,i) + Q_(j,j) > 2·Q_(i,j) for i≠j
	for(int i = 0; i < nrVars; ++i)
	{
		if(std::abs(Q(i, i)) < DIAG_CONSTANT)
		{
			Q(i, i) += DIAG_CONSTANT;
		}
	}
}


int fillEq(const std::vector<Equality*>& eq, int nrVars,
	int nrALines, Eigen::MatrixXd& A, Eigen::VectorXd& AL, Eigen::VectorXd& AU)
{
	for(std::size_t i = 0; i < eq.size(); ++i)
	{
		// ineq constraint can return a matrix with more line
		// than the number of constraint
		int nrConstr = eq[i]->nrEq();
		const Eigen::MatrixXd& Ai = eq[i]->AEq();
		const Eigen::VectorXd& bi = eq[i]->bEq();

		A.block(nrALines, 0, nrConstr, nrVars) =
			Ai.block(0, 0, nrConstr, nrVars);
		AL.segment(nrALines, nrConstr) = bi.head(nrConstr);
		AU.segment(nrALines, nrConstr) = bi.head(nrConstr);

		nrALines += nrConstr;
	}

	return nrALines;
}


int fillInEq(const std::vector<Inequality*>& inEq, int nrVars,
	int nrALines, Eigen::MatrixXd& A, Eigen::VectorXd& AL, Eigen::VectorXd& AU)
{
	for(std::size_t i = 0; i < inEq.size(); ++i)
	{
		// ineq constraint can return a matrix with more line
		// than the number of constraint
		int nrConstr = inEq[i]->nrInEq();
		const Eigen::MatrixXd& Ai = inEq[i]->AInEq();
		const Eigen::VectorXd& bi = inEq[i]->bInEq();

		A.block(nrALines, 0, nrConstr, nrVars) =
			Ai.block(0, 0, nrConstr, nrVars);
		AL.segment(nrALines, nrConstr).fill(-std::numeric_limits<double>::infinity());
		AU.segment(nrALines, nrConstr) = bi.head(nrConstr);

		nrALines += nrConstr;
	}

	return nrALines;
}


int fillGenInEq(const std::vector<GenInequality*>& genInEq, int nrVars,
	int nrALines, Eigen::MatrixXd& A, Eigen::VectorXd& AL, Eigen::VectorXd& AU)
{
	for(std::size_t i = 0; i < genInEq.size(); ++i)
	{
		// ineq constraint can return a matrix with more line
		// than the number of constraint
		int nrConstr = genInEq[i]->nrGenInEq();
		const Eigen::MatrixXd& Ai = genInEq[i]->AGenInEq();
		const Eigen::VectorXd& ALi = genInEq[i]->LowerGenInEq();
		const Eigen::VectorXd& AUi = genInEq[i]->UpperGenInEq();

		A.block(nrALines, 0, nrConstr, nrVars) =
			Ai.block(0, 0, nrConstr, nrVars);
		AL.segment(nrALines, nrConstr) = ALi.head(nrConstr);
		AU.segment(nrALines, nrConstr) = AUi.head(nrConstr);

		nrALines += nrConstr;
	}

	return nrALines;
}


void fillBound(const std::vector<Bound*>& bounds,
	Eigen::VectorXd& XL, Eigen::VectorXd& XU)
{
	for(std::size_t i = 0; i < bounds.size(); ++i)
	{
		const Eigen::VectorXd& XLi = bounds[i]->Lower();
		const Eigen::VectorXd& XUi = bounds[i]->Upper();
		int bv = bounds[i]->beginVar();

		XL.segment(bv, XLi.size()) = XLi;
		XU.segment(bv, XUi.size()) = XUi;
	}
}



/**
	*																LSSOLQPSolver
	*/



LSSOLQPSolver::LSSOLQPSolver():
	A_(),AL_(),AU_(),
	XL_(),XU_(),
	Q_(),C_(),
	nrALines_(0)
{
	lssol_.warm(false);
	lssol_.feasibilityTol(1e-6);
}


void LSSOLQPSolver::updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq)
{
	int maxALines = nrEq + nrInEq + nrGenInEq;
	A_.resize(maxALines, nrVars);
	AL_.resize(maxALines);
	AU_.resize(maxALines);

	XL_.resize(nrVars);
	XU_.resize(nrVars);

	Q_.resize(nrVars, nrVars);
	C_.resize(nrVars);

	lssol_.problem(nrVars, maxALines);
}


void LSSOLQPSolver::updateMatrix(
	const std::vector<Task*>& tasks,
	const std::vector<Equality*>& eqConstr,
	const std::vector<Inequality*>& inEqConstr,
	const std::vector<GenInequality*>& genInEqConstr,
	const std::vector<Bound*>& boundConstr)
{
	A_.setZero();
	AL_.setZero();
	AU_.setZero();
	XL_.fill(-std::numeric_limits<double>::infinity());
	XU_.fill(std::numeric_limits<double>::infinity());
	Q_.setZero();
	C_.setZero();

	const int nrVars = int(Q_.rows());

	nrALines_ = 0;
	nrALines_ = fillEq(eqConstr, nrVars, nrALines_, A_, AL_, AU_);
	nrALines_ = fillInEq(inEqConstr, nrVars, nrALines_, A_, AL_, AU_);
	nrALines_ = fillGenInEq(genInEqConstr, nrVars, nrALines_, A_, AL_, AU_);

	fillBound(boundConstr, XL_, XU_);
	fillQC(tasks, nrVars, Q_, C_);
}


bool LSSOLQPSolver::solve()
{
	bool success = lssol_.solve(Q_, C_,
		A_.block(0, 0, nrALines_, int(A_.cols())), int(A_.rows()),
		AL_.segment(0, nrALines_), AU_.segment(0, nrALines_), XL_, XU_);
	return success;
}


const Eigen::VectorXd& LSSOLQPSolver::result() const
{
	return lssol_.result();
}



} // namespace qp

} // namespace tasks
