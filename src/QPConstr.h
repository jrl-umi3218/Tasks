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

// SCD
#include <SCD/Matrix/SCD_Types.h>

// Tasks
#include "QPSolver.h"

// forward declaration
// SCD
namespace SCD
{
class S_Object;
class CD_Pair;
}


namespace tasks
{

namespace qp
{

SCD::Matrix4x4 toSCD(const sva::PTransform& t);

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
		int body;
		rbd::Jacobian jac;
		Eigen::Vector3d point;
		FrictionCone cone;
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



class JointLimitsConstr : public BoundConstraint, public Constraint
{
public:
	JointLimitsConstr(const rbd::MultiBody& mb,
		std::vector<std::vector<double>> lBound,
		std::vector<std::vector<double>> uBound,
		double step);

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		int alphaD, int lambda, int torque, const std::vector<Contact>& cont);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// Bound Constraint
	virtual int beginVar();

	virtual const Eigen::VectorXd& Lower() const;
	virtual const Eigen::VectorXd& Upper() const;

private:
	Eigen::VectorXd lower_, upper_;
	Eigen::VectorXd qMin_, qMax_;
	Eigen::VectorXd qVec_, alphaVec_;
	int begin_;
	double step_;
};



class TorqueLimitsConstr : public BoundConstraint, public Constraint
{
public:
	TorqueLimitsConstr(const rbd::MultiBody& mb,
		std::vector<std::vector<double>> lBound,
		std::vector<std::vector<double>> uBound);

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		int alphaD, int lambda, int torque, const std::vector<Contact>& cont);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// Bound Constraint
	virtual int beginVar();

	virtual const Eigen::VectorXd& Lower() const;
	virtual const Eigen::VectorXd& Upper() const;

private:
	Eigen::VectorXd lower_, upper_;
	int begin_;
	double step_;
};



class SelfCollisionConstr : public InequalityConstraint, public Constraint
{
public:
	SelfCollisionConstr(const rbd::MultiBody& mb, double step);

	void addCollision(const rbd::MultiBody& mb,
		int body1Id, SCD::S_Object* body1, const sva::PTransform& body1T,
		int body2Id, SCD::S_Object* body2, const sva::PTransform& body2T,
		double di, double ds, double damping);
	void rmCollision(int body1Id, int body2Id);
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		int alphaD, int lambda, int torque, const std::vector<Contact>& cont);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// InEquality Constraint
	virtual int nrInEqLine();

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& BInEq() const;

private:
	struct CollData
	{
		CollData(const rbd::MultiBody& mb,
			int body1Id, SCD::S_Object* body1, const sva::PTransform& body1T,
			int body2Id, SCD::S_Object* body2, const sva::PTransform& body2T,
			double di, double ds, double damping);

		SCD::CD_Pair* pair;
		sva::PTransform body1T, body2T;
		Eigen::Vector3d normVecDist;
		rbd::Jacobian jacB1, jacB2;
		double di, ds;
		double damping;
		int body1Id, body2Id;
		int body1, body2;
	};

private:
	std::vector<CollData> dataVec_;
	double step_;
	int nrVars_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd BInEq_;

	Eigen::MatrixXd fullJac_;
	Eigen::MatrixXd fullJacDot_;
	Eigen::VectorXd alphaVec_;
	Eigen::VectorXd calcVec_;
};



class StaticEnvCollisionConstr : public InequalityConstraint, public Constraint
{
public:
	StaticEnvCollisionConstr(const rbd::MultiBody& mb, double step);

	void addCollision(const rbd::MultiBody& mb,
		int bodyId, SCD::S_Object* body, const sva::PTransform& bodyT,
		int envId, SCD::S_Object* env,
		double di, double ds, double damping);
	void rmCollision(int bodyId, int envId);
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		int alphaD, int lambda, int torque, const std::vector<Contact>& cont);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	// InEquality Constraint
	virtual int nrInEqLine();

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& BInEq() const;

private:
	struct CollData
	{
		CollData(const rbd::MultiBody& mb,
			int bodyId, SCD::S_Object* body, const sva::PTransform& bodyT,
			int envId, SCD::S_Object* env,
			double di, double ds, double damping);

		SCD::CD_Pair* pair;
		sva::PTransform bodyT;
		Eigen::Vector3d normVecDist;
		rbd::Jacobian jacB1;
		double di, ds;
		double damping;
		int bodyId, envId;
		int body;
	};

private:
	std::vector<CollData> dataVec_;
	double step_;
	int nrVars_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd BInEq_;

	Eigen::MatrixXd fullJac_;
	Eigen::MatrixXd fullJacDot_;
	Eigen::VectorXd alphaVec_;
	Eigen::VectorXd calcVec_;
};

} // namespace qp

} // namespace tasks

