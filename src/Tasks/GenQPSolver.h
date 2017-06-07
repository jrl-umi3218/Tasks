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

#include <tasks/config.hh>

// forward declaration
// RBDyn
namespace rbd
{
class MultiBody;
}


namespace tasks
{

namespace qp
{
// forward declarition
class Task;
class Equality;
class Inequality;
class GenInequality;
class Bound;
class GenQPSolver;


/**
	* Factory to create GenQPSolver implementation.
	* Two argument are supported QLD and LSSOL.
	*/
TASKS_DLLAPI GenQPSolver* createQPSolver(const std::string& name);


/**
	* Generic QP solver abstract interface.
	* Solve the following problem:
	* \f{align}
	* \underset{x}{\text{minimize }} & \frac{1}{2} x^T Q x + x^T c\\
	* \text{s.t. } & L \leq \left\{ \begin{array}{c} x \\ A x \end{array} \right\} \leq U
	* \f}
	*/
class TASKS_DLLAPI GenQPSolver
{
public:
	/// Default QP solver.
	static const std::string default_qp_solver;

public:
	virtual ~GenQPSolver() {}

	/**
		* Update the problem size.
		* @param nrVars Variable number.
		* @param nrEq maximum number of equality.
		* @param nrInEq maximum number of inequality.
		* @param nrGenInEq maximum number of general inequality.
		*/
	virtual void updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq) = 0;

	/**
	* Setup dependent variables, only linear dependencies are supported
	* @param nrVars Variable number.
	* @param dependencies List of tuple {primary, replica, factor}
	*/
	virtual void setDependencies(int nrVars, std::vector<std::tuple<int, int, double>> dependencies);

	/**
		* Construct the QP matrices.
		* @param tasks Build \f$ Q \f$ and \f$ c \f$.
		* @param eqConstr Build \f$ A x = b \f$ constraints.
		* @param inEqConstr Build \f$ A x \geq b \f$ constraints.
		* @param genInEqConstr Build \f$ L \leq A x \leq U \f$ constraints.
		* @param boundConstr Build \f$ L \leq x \leq U \f$ constraints.
		*/
	virtual void updateMatrix(const std::vector<Task*>& tasks,
		const std::vector<Equality*>& eqConstr,
		const std::vector<Inequality*>& inEqConstr,
		const std::vector<GenInequality*>& genInEqConstr,
		const std::vector<Bound*>& boundConstr) = 0;

	/**
		* Solve the quadratic program.
		* @return true of success false on failure.
		*/
	virtual bool solve() = 0;

	/// @return Optimal \f$ x \f$ vector.
	virtual const Eigen::VectorXd& result() const = 0;

	/// @return Error message if GenQPSolver::solve has returned false.
	virtual std::ostream& errorMsg(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<Task*>& tasks,
		const std::vector<Equality*>& eqConstr,
		const std::vector<Inequality*>& inEqConstr,
		const std::vector<GenInequality*>& genInEqConstr,
		const std::vector<Bound*>& boundConstr,
		std::ostream& out) const = 0;
protected:
	/** Correspondence between full variable indices and reduced variables */
	std::vector<int> fullToReduced_;
	/** Correspondence between reduced variable indices and full variable indices */
	std::vector<int> reducedToFull_;

	/** Variable dependencies, each tuple gives the primary variable index in the
	 * full variable, replica variable index in the full variable and the factor
	 * in the dependency equation: replica = factor * primary */
	std::vector<std::tuple<int, int, double>> dependencies_;
};


} // namespace qp

} // namespace tasks
