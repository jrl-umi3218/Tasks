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
#include <RBDyn/FD.h>
#include <RBDyn/Jacobian.h>

// Tasks
#include "QPSolver.h"


namespace tasks
{

namespace qp
{

class MotionConstrCommon : public ConstraintFunction<GenInequality, Bound>
{
public:
	MotionConstrCommon(const rbd::MultiBody& mb);

	void computeTorque(const Eigen::VectorXd& alphaD,
										const Eigen::VectorXd& lambda);
	const Eigen::VectorXd& torque() const;
	void torque(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc) const;

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	void computeMatrix(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// Description
	virtual std::string nameGenInEq() const;
	virtual std::string descGenInEq(const rbd::MultiBody& mb, int line);
	virtual std::string nameBound() const;
	virtual std::string descBound(const rbd::MultiBody& mb, int line);

	// Inequality Constraint
	virtual int maxGenInEq() const;

	virtual const Eigen::MatrixXd& AGenInEq() const;
	virtual const Eigen::VectorXd& LowerGenInEq() const;
	virtual const Eigen::VectorXd& UpperGenInEq() const;

	// Bound Constraint
	virtual int beginVar() const;

	virtual const Eigen::VectorXd& Lower() const;
	virtual const Eigen::VectorXd& Upper() const;

protected:
	struct ContactData
	{
		ContactData() {}
		ContactData(const rbd::MultiBody& mb, int body,
			std::vector<Eigen::Vector3d> points,
			const std::vector<FrictionCone>& cones);


		rbd::Jacobian jac;
		int body;
		std::vector<Eigen::Vector3d> points;
		std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic> > generators;
		// Hold the translated jacobian
		Eigen::MatrixXd jacTrans;
		// Hold the generator in world frame
		std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic> > generatorsComp;
	};

protected:
	rbd::ForwardDynamics fd_;

	std::vector<ContactData> cont_;
	Eigen::MatrixXd fullJac_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd AL_, AU_;

	Eigen::VectorXd XL_, XU_;

	Eigen::VectorXd curTorque_;

	int nrDof_, nrFor_, nrTor_;
};



class MotionConstr : public MotionConstrCommon
{
public:
	MotionConstr(const rbd::MultiBody& mb,
							 std::vector<std::vector<double>> lTorqueBounds,
							 std::vector<std::vector<double>> uTorqueBounds);

	// Constraint
	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

private:
	Eigen::VectorXd torqueL_, torqueU_;
};


struct SpringJoint
{
	SpringJoint(){}
	SpringJoint(int jId, double K, double C, double O):
		jointId(jId),K(K),C(C),O(O)
	{}

	int jointId;
	double K, C, O;
};


class MotionSpringConstr : public MotionConstrCommon
{
public:
	MotionSpringConstr(const rbd::MultiBody& mb,
										std::vector<std::vector<double>> lTorqueBounds,
										std::vector<std::vector<double>> uTorqueBounds,
										const std::vector<SpringJoint>& springs);

	// Constraint
	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);
private:
	struct SpringJointData
	{
		int index;
		int posInDof;
		double K;
		double C;
		double O;
	};

private:
	Eigen::VectorXd torqueL_, torqueU_;
	std::vector<SpringJointData> springs_;
};



/**
 * @brief Use polynome in function of q to compute torque limits.
 * BEWARE: Only work with 1 dof/param joint
 */
class MotionPolyConstr : public MotionConstrCommon
{
public:
	MotionPolyConstr(const rbd::MultiBody& mb,
									const std::vector<std::vector<Eigen::VectorXd>>& lTorqueBounds,
									const std::vector<std::vector<Eigen::VectorXd>>& uTorqueBounds);

	// Constraint
	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

private:
	std::vector<Eigen::VectorXd> torqueL_, torqueU_;
	std::vector<int> jointIndex_;
};


} // namespace qp

} // namespace tasks
