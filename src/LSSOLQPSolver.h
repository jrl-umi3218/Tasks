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
// eigen-lssol
#include <eigen-lssol/LSSOL.h>

// Tasks
#include "GenQPSolver.h"


namespace tasks
{

namespace qp
{


/**
	* GenQPSolver interface implementation with the LSSOL QP solver.
	*/
class TASKS_DLLAPI LSSOLQPSolver : public GenQPSolver
{
public:
	LSSOLQPSolver();

	virtual void updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq);
	virtual void updateMatrix(const std::vector<Task*>& tasks,
		const std::vector<Equality*>& eqConstr,
		const std::vector<Inequality*>& inEqConstr,
		const std::vector<GenInequality*>& genInEqConstr,
		const std::vector<Bound*>& boundConstr);
	virtual bool solve();
	virtual const Eigen::VectorXd& result() const;
	virtual std::ostream& errorMsg(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<Task*>& tasks,
		const std::vector<Equality*>& eqConstr,
		const std::vector<Inequality*>& inEqConstr,
		const std::vector<GenInequality*>& genInEqConstr,
		const std::vector<Bound*>& boundConstr,
		std::ostream& out) const override;

private:
	Eigen::LSSOL lssol_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd AL_, AU_;

	Eigen::VectorXd XL_;
	Eigen::VectorXd XU_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;

	int nrALines_;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace qp

} // namespace tasks
