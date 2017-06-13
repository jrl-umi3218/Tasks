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

#pragma once

// includes
// std
#include <vector>

// Eigen
#include <Eigen/Core>

// Tasks
#include "Tasks/QPSolver.h"


namespace tasks
{

namespace qp
{


// Value add to the diagonal to ensure positive matrix
static const double DIAG_CONSTANT = 1e-4;


/**
	* Fill the \f$ Q \f$ matrix and the \f$ c \f$ vector based on the
	* task list.
	*/
inline void fillQC(const std::vector<Task*>& tasks, int nrVars,
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

/**
	* Reduce \f$ Q \f$ matrix and the \f$ c \f$ vector based on the dependencies list
	*/
inline void reduceQC(const Eigen::MatrixXd & QFull, const Eigen::VectorXd & CFull,
										 Eigen::MatrixXd & Q, Eigen::VectorXd & C,
										 const std::vector<int> & fullToReduced,
										 const std::vector<int> & reducedToFull,
										 const std::vector<std::tuple<int, int, double>> & dependencies)
{
	/* Start by moving the non-reduced variables to their new location */
	for(size_t i = 0; i < reducedToFull.size(); ++i)
	{
		for(size_t j = i; j < reducedToFull.size(); ++j)
		{
			Q(i,j) = QFull(reducedToFull[i], reducedToFull[j]);
		}
		C(i) = CFull(reducedToFull[i]);
	}
	/* Apply the reduction */
	for(const auto & d : dependencies)
	{
		const int & primaryFullI = std::get<0>(d);
		const int & primaryReducedI = fullToReduced[primaryFullI];
		const int & replicaFullI = std::get<1>(d);
		const double & alpha = std::get<2>(d);
		/* Apply reduction to C */
		C(primaryReducedI) += alpha*CFull(replicaFullI);
		/* Add cross-terms to Q */
		for(size_t i = 0; i < reducedToFull.size(); ++i)
		{
			Q(i, primaryReducedI) += alpha*QFull(reducedToFull[i], replicaFullI);
		}
		/* Add diagonal element to Q */
		Q(primaryReducedI, primaryReducedI) += alpha*alpha*QFull(replicaFullI, replicaFullI);
	}
	/* Update the lower triangular part of Q */
	Q = Q.selfadjointView<Eigen::Upper>();
}

// general qp form


/**
	* Fill the \f$ A \f$ matrix and the \f$ L \f$ and \f$ U \f$ bounds vectors
	* based on the equality constaint list.
	*/
inline int fillEq(const std::vector<Equality*>& eq, int nrVars,
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


/**
	* Fill the \f$ A \f$ matrix and the \f$ L \f$ and \f$ U \f$ bounds vectors
	* based on the inequality constaint list.
	*/
inline int fillInEq(const std::vector<Inequality*>& inEq, int nrVars,
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


/**
	* Fill the \f$ A \f$ matrix and the \f$ L \f$ and \f$ U \f$ bounds vectors
	* based on the general inequality constaint list.
	*/
inline int fillGenInEq(const std::vector<GenInequality*>& genInEq, int nrVars,
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


/**
	* Fill the \f$ A \f$ matrix and the \f$ b \f$ vectors
	* based on the equality constaint list.
	*/
inline int fillEq(const std::vector<Equality*>& eq, int nrVars,
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


/**
	* Fill the \f$ A \f$ matrix and the \f$ b \f$ vectors
	* based on the inequality constaint list.
	*/
inline int fillInEq(const std::vector<Inequality*>& inEq, int nrVars,
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


/**
	* Fill the \f$ A \f$ matrix and the \f$ b \f$ vectors
	* based on the general inequality constaint list.
	*/
inline int fillGenInEq(const std::vector<GenInequality*>& genInEq, int nrVars,
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


/**
	* Fill the \f$ L \f$  and \f$ U \f$ bounds vectors
	* based on the bound constaint list.
	*/
inline void fillBound(const std::vector<Bound*>& bounds,
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
	* Reduce \f$ A \f$ matrix based on the dependencies list
	*/
inline void reduceA(const Eigen::MatrixXd & AFull,
										Eigen::MatrixXd & A,
										const std::vector<int> & fullToReduced,
										const std::vector<int> & reducedToFull,
										const std::vector<std::tuple<int, int, double>> & dependencies)
{
	for(size_t i = 0; i <reducedToFull.size(); ++i)
	{
		A.col(i) = AFull.col(reducedToFull[i]);
	}
	for(const auto & d : dependencies)
	{
		const int & primaryFullI = std::get<0>(d);
		const int & primaryReducedI = fullToReduced[primaryFullI];
		const int & replicaFullI = std::get<1>(d);
		const double & alpha = std::get<2>(d);
		A.col(primaryReducedI) += alpha*AFull.col(replicaFullI);
	}
}

/**
	* Reduce bounds vector based on the dependencies list
	*/
inline void reduceBound(const Eigen::VectorXd & XLFull,
												Eigen::VectorXd & XL,
												const Eigen::VectorXd & XUFull,
												Eigen::VectorXd & XU,
												const std::vector<int> & fullToReduced,
												const std::vector<int> & reducedToFull,
												const std::vector<std::tuple<int, int, double>> & dependencies)
{
	for(size_t i = 0; i < static_cast<size_t>(XL.rows()); ++i)
	{
		XL(i) = XLFull(reducedToFull[i]);
		XU(i) = XUFull(reducedToFull[i]);
	}
	for(const auto & d : dependencies)
	{
		const int & primaryFullI = std::get<0>(d);
		const int & primaryReducedI = fullToReduced[primaryFullI];
		const int & replicaFullI = std::get<1>(d);
		const double & alpha = std::get<2>(d);
		if(alpha != 0)
		{
			/** If alpha is negative, the upper/lower bounds should be inverted */
			if(alpha < 0)
			{
				XL(primaryReducedI) = std::max(XL(primaryReducedI), XUFull(replicaFullI)/alpha);
				XU(primaryReducedI) = std::min(XU(primaryReducedI), XLFull(replicaFullI)/alpha);
			}
			else
			{
				XL(primaryReducedI) = std::max(XL(primaryReducedI), XLFull(replicaFullI)/alpha);
				XU(primaryReducedI) = std::min(XU(primaryReducedI), XUFull(replicaFullI)/alpha);
			}
		}
	}
}

/**
	* Return the full variable from the reduced variable
	*/
inline void expandResult(const Eigen::VectorXd & result,
												 Eigen::VectorXd & resultFull,
												 const std::vector<int> & reducedToFull,
												 const std::vector<std::tuple<int, int, double>> & dependencies)
{
	for(size_t i = 0; i < reducedToFull.size(); ++i)
	{
		resultFull(reducedToFull[i]) = result(i);
	}
	for(const auto & d : dependencies)
	{
		const int & primaryFullI = std::get<0>(d);
		const int & replicaFullI = std::get<1>(d);
		const double & alpha = std::get<2>(d);
		resultFull(replicaFullI) = alpha*resultFull(primaryFullI);
	}
}

// print of a constraint at a given line
template<typename T>
std::ostream& printConstr(const Eigen::VectorXd& result, T* constr, int line,
	std::ostream& out);

template<>
inline std::ostream& printConstr(const Eigen::VectorXd& result, Equality* constr,
	int line, std::ostream& out)
{
	out << constr->AEq().row(line)*result <<" = " <<
				 constr->bEq()(line);
	return out;
}

template<>
inline std::ostream& printConstr(const Eigen::VectorXd& result, Inequality* constr,
	int line, std::ostream& out)
{
	out << constr->AInEq().row(line)*result <<" <= " <<
				 constr->bInEq()(line);
	return out;
}

template<>
inline std::ostream& printConstr(const Eigen::VectorXd& result, GenInequality* constr,
	int line, std::ostream& out)
{
	out << constr->LowerGenInEq()(line) << " <= " <<
				 constr->AGenInEq().row(line)*result <<" <= " <<
				 constr->UpperGenInEq()(line);
	return out;
}


template <typename T>
inline std::ostream& constrErrorMsg(const std::vector<rbd::MultiBody>& mbs,
	const Eigen::VectorXd& result,
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
			out << constr_traits<T>::desc(e, mbs, line) << std::endl;
			printConstr(result, e, line, out) << std::endl;
			start = end;
			break;
		}
		start = end;
	}
	return out;
}


} // namespace qp

} // namespace tasks
