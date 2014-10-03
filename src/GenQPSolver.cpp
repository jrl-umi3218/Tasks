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
// std
#include <map>

// tasks
#include "QPSolver.h"


namespace tasks
{

namespace qp
{



/**
	*													Factory
	*/



template <typename T>
T* allocateQP()
{
	return new T;
}


static const
std::map<std::string, std::function<GenQPSolver*(void)> > qpFactory = {
	{"LSSOL", allocateQP<LSSOLQPSolver>},
	{"QLD", allocateQP<QLDQPSolver>}
};


GenQPSolver* createQPSolver(const std::string& name)
{
	return qpFactory.at(name)();
}



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


// general qp form


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


// standard qp form


int fillEq(const std::vector<Equality*>& eq, int nrVars,
	int nrALines, Eigen::MatrixXd& A, Eigen::VectorXd& b)
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
		b.segment(nrALines, nrConstr) = bi.head(nrConstr);

		nrALines += nrConstr;
	}

	return nrALines;
}


int fillInEq(const std::vector<Inequality*>& inEq, int nrVars,
	int nrALines, Eigen::MatrixXd& A, Eigen::VectorXd& b)
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
		b.segment(nrALines, nrConstr) = bi.head(nrConstr);

		nrALines += nrConstr;
	}

	return nrALines;
}


int fillGenInEq(const std::vector<GenInequality*>& genInEq, int nrVars,
	int nrALines, Eigen::MatrixXd& A, Eigen::VectorXd& b)
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
			-Ai.block(0, 0, nrConstr, nrVars);
		b.segment(nrALines, nrConstr) = -ALi.head(nrConstr);

		nrALines += nrConstr;

		A.block(nrALines, 0, nrConstr, nrVars) =
			Ai.block(0, 0, nrConstr, nrVars);
		b.segment(nrALines, nrConstr) = AUi.head(nrConstr);

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


// print of a constraint at a given line
template<typename T>
std::ostream& printConstr(const Eigen::VectorXd& result, T* constr, int line,
	std::ostream& out);

template<>
std::ostream& printConstr(const Eigen::VectorXd& result, Equality* constr,
	int line, std::ostream& out)
{
	out << constr->AEq().row(line)*result <<" = " <<
				 constr->bEq()(line);
	return out;
}

template<>
std::ostream& printConstr(const Eigen::VectorXd& result, Inequality* constr,
	int line, std::ostream& out)
{
	out << constr->AInEq().row(line)*result <<" <= " <<
				 constr->bInEq()(line);
	return out;
}

template<>
std::ostream& printConstr(const Eigen::VectorXd& result, GenInequality* constr,
	int line, std::ostream& out)
{
	out << constr->LowerGenInEq()(line) << " <= " <<
				 constr->AGenInEq().row(line)*result <<" <= " <<
				 constr->UpperGenInEq()(line);
	return out;
}


template <typename T>
std::ostream& constrErrorMsg(const rbd::MultiBody& mb, const Eigen::VectorXd& result,
	int ALine, const std::vector<T*>& constr, int& start, int& end,
	std::ostream& out)
{
	for(T* e: constr)
	{
		end += constr_traits<T>::nrLines(e);
		if(ALine >= start && ALine < end)
		{
			int line = ALine - start;
			out << constr_traits<T>::name(e) << " violated at line: " << line << std::endl;
			out << constr_traits<T>::desc(e, mb, line) << std::endl;
			printConstr(result, e, line, out) << std::endl;
			start = end;
			break;
		}
		start = end;
	}
	return out;
}



/**
	*																LSSOLQPSolver
	*/



LSSOLQPSolver::LSSOLQPSolver():
	lssol_(),
	A_(),AL_(),AU_(),
	XL_(),XU_(),
	Q_(),C_(),
	nrALines_(0)
{
	lssol_.warm(true);
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
					out << b->descBound(mb, line) << std::endl;
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

			constrErrorMsg(mb, lssol_.result(), i, eqConstr, start, end, out);
			constrErrorMsg(mb, lssol_.result(), i, inEqConstr, start, end, out);
			constrErrorMsg(mb, lssol_.result(), i, genInEqConstr, start, end, out);
			out << std::endl;
		}
	}
	return out;
}



/**
	*																QLDQPSolver
	*/



QLDQPSolver::QLDQPSolver():
	qld_(),
	Aeq_(),Aineq_(),
	beq_(), bineq_(),
	XL_(),XU_(),
	Q_(),C_(),
	nrAeqLines_(0), nrAineqLines_(0)
{
}


void QLDQPSolver::updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq)
{
	int maxAeqLines = nrEq;
	int maxAineqLines = nrInEq + nrGenInEq*2;

	Aeq_.resize(maxAeqLines, nrVars);
	Aineq_.resize(maxAineqLines, nrVars);

	beq_.resize(maxAeqLines);
	bineq_.resize(maxAineqLines);

	XL_.resize(nrVars);
	XU_.resize(nrVars);

	Q_.resize(nrVars, nrVars);
	C_.resize(nrVars);

	qld_.problem(nrVars, maxAeqLines, maxAineqLines);
}


void QLDQPSolver::updateMatrix(
	const std::vector<Task*>& tasks,
	const std::vector<Equality*>& eqConstr,
	const std::vector<Inequality*>& inEqConstr,
	const std::vector<GenInequality*>& genInEqConstr,
	const std::vector<Bound*>& boundConstr)
{
	Aeq_.setZero();
	Aineq_.setZero();
	beq_.setZero();
	bineq_.setZero();
	XL_.fill(-std::numeric_limits<double>::infinity());
	XU_.fill(std::numeric_limits<double>::infinity());
	Q_.setZero();
	C_.setZero();

	const int nrVars = int(Q_.rows());

	nrAeqLines_ = 0;
	nrAeqLines_ = fillEq(eqConstr, nrVars, nrAeqLines_, Aeq_, beq_);
	nrAineqLines_ = 0;
	nrAineqLines_ = fillInEq(inEqConstr, nrVars, nrAineqLines_, Aineq_, bineq_);
	nrAineqLines_ = fillGenInEq(genInEqConstr, nrVars, nrAineqLines_, Aineq_, bineq_);

	fillBound(boundConstr, XL_, XU_);
	fillQC(tasks, nrVars, Q_, C_);
}


bool QLDQPSolver::solve()
{
	bool success = qld_.solve(Q_, C_,
		Aeq_.block(0, 0, nrAeqLines_, int(Aeq_.cols())), beq_.segment(0, nrAeqLines_),
		Aineq_.block(0, 0, nrAineqLines_, int(Aineq_.cols())), bineq_.segment(0, nrAineqLines_),
		XL_, XU_, 1e-6);
	return success;
}


const Eigen::VectorXd& QLDQPSolver::result() const
{
	return qld_.result();
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
