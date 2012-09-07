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
// Eigen
#include <Eigen/Core>

// RBDyn
#include <FD.h>
#include <Jacobian.h>

// Tasks
#include "QPSolver.h"


namespace tasks
{

namespace qp
{

class MotionConstr : public EqualityConstraint, public BoundConstraint, public Constraint
{
public:
	MotionConstr(const rbd::MultiBody& mb);

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		int alphaD, int lambda, int torque, const std::vector<Contact>& cont);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// Equality Constraint
	virtual int nrEqLine();

	virtual const Eigen::MatrixXd& AEq() const;
	virtual const Eigen::VectorXd& BEq() const;

	// Bound Constraint
	virtual int beginVar();

	virtual const Eigen::VectorXd& Lower() const;
	virtual const Eigen::VectorXd& Upper() const;

private:
	struct ContactData
	{
		rbd::Jacobian jac;
		Eigen::MatrixXd transJac;
		std::vector<Eigen::Vector3d> points;
		std::vector<Eigen::Vector3d> normals;
	};

private:
	rbd::ForwardDynamics fd_;

	std::vector<ContactData> cont_;
	Eigen::MatrixXd fullJac_;

	Eigen::MatrixXd AEq_;
	Eigen::VectorXd BEq_;

	Eigen::VectorXd XL_, XU_;

	int nrDof_, nrFor_, nrTor_;
};



class ContactAccConstr : public EqualityConstraint, public Constraint
{
public:
	ContactAccConstr(const rbd::MultiBody& mb);

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		int alphaD, int lambda, int torque, const std::vector<Contact>& cont);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// Equality Constraint
	virtual int nrEqLine();

	virtual const Eigen::MatrixXd& AEq() const;
	virtual const Eigen::VectorXd& BEq() const;

private:
	struct ContactData
	{
		rbd::Jacobian jac;
	};

private:
	std::vector<ContactData> cont_;

	Eigen::MatrixXd fullJac_;
	Eigen::VectorXd alphaVec_;

	Eigen::MatrixXd AEq_;
	Eigen::VectorXd BEq_;

	int nrDof_, nrFor_, nrTor_;
};

} // namespace qp

} // namespace tasks

