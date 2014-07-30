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

// Eigen
#include <EigenQP/LSSOL.h>


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



class GenQPSolver
{
public:
	virtual ~GenQPSolver() {}

	virtual void updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq) = 0;
	virtual void updateMatrix(const std::vector<Task*>& tasks,
		const std::vector<Equality*>& eqConstr,
		const std::vector<Inequality*>& inEqConstr,
		const std::vector<GenInequality*>& genInEqConstr,
		const std::vector<Bound*>& boundConstr) = 0;
	virtual bool solve() = 0;
	virtual const Eigen::VectorXd& result() const = 0;
};



class LSSOLQPSolver : public GenQPSolver
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

private:
	Eigen::LSSOL lssol_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd AL_, AU_;

	Eigen::VectorXd XL_;
	Eigen::VectorXd XU_;

	Eigen::MatrixXd Q_;
	Eigen::VectorXd C_;

	int nrALines_;
};


} // namespace qp

} // namespace tasks
