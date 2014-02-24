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
#include <RBDyn/CoM.h>

// SCD
#include <SCD/Matrix/SCD_Types.h>
#include <SCD/S_Object/S_Sphere.h>

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

SCD::Matrix4x4 toSCD(const sva::PTransformd& t);


class ContactConstrCommon
{
public:
	bool addVirtualContact(int bodyId);
	bool removeVirtualContact(int bodyId);
	void resetVirtualContacts();

protected:
	std::set<int> bodyIdInContact(const rbd::MultiBody& mb,
		const SolverData& data);

protected:
	std::set<int> virtualContacts_;
};


class ContactAccConstr : public ConstraintFunction<Inequality>,
	public ContactConstrCommon
{
public:
	ContactAccConstr(const rbd::MultiBody& mb);

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// Inequality Constraint
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& LowerInEq() const;
	virtual const Eigen::VectorXd& UpperInEq() const;

private:
	struct ContactData
	{
		ContactData(rbd::Jacobian j):
			jac(j)
		{}


		rbd::Jacobian jac;
	};

private:
	std::vector<ContactData> cont_;

	Eigen::MatrixXd fullJac_;
	Eigen::VectorXd alphaVec_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd ALU_;

	int nrDof_, nrFor_, nrTor_;
};



class ContactSpeedConstr : public ConstraintFunction<Inequality>,
	public ContactConstrCommon
{
public:
	ContactSpeedConstr(const rbd::MultiBody& mb, double timeStep);

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// Inequality Constraint
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& LowerInEq() const;
	virtual const Eigen::VectorXd& UpperInEq() const;

private:
	struct ContactData
	{
		ContactData(rbd::Jacobian j):
			jac(j),
			body(j.jointsPath().back())
		{}

		rbd::Jacobian jac;
		int body;
	};

private:
	std::vector<ContactData> cont_;

	Eigen::MatrixXd fullJac_;
	Eigen::VectorXd alphaVec_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd ALU_;

	int nrDof_, nrFor_, nrTor_;
	double timeStep_;
};



class JointLimitsConstr : public ConstraintFunction<Bound>
{
public:
	JointLimitsConstr(const rbd::MultiBody& mb,
		std::vector<std::vector<double> > lBound,
		std::vector<std::vector<double> > uBound,
		double step);

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual std::string nameBound() const;
	virtual std::string descBound(const rbd::MultiBody& mb, int line);

	// Bound Constraint
	virtual int beginVar() const;

	virtual const Eigen::VectorXd& Lower() const;
	virtual const Eigen::VectorXd& Upper() const;

private:
	Eigen::VectorXd lower_, upper_;
	Eigen::VectorXd qMin_, qMax_;
	Eigen::VectorXd qVec_, alphaVec_;
	int begin_;
	double step_;
};



class DamperJointLimitsConstr : public ConstraintFunction<Bound>
{
public:
	DamperJointLimitsConstr(const rbd::MultiBody& mb,
		const std::vector<std::vector<double> >& lBound,
		const std::vector<std::vector<double> >& uBound,
		std::vector<std::vector<double> > lVel,
		std::vector<std::vector<double> > uVel,
		double interPercent, double securityPercent, double damperOffset, double step);

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual std::string nameBound() const;
	virtual std::string descBound(const rbd::MultiBody& mb, int line);

	// Bound Constraint
	virtual int beginVar() const;

	virtual const Eigen::VectorXd& Lower() const;
	virtual const Eigen::VectorXd& Upper() const;

	/// compute damping that avoid speed jump
	double computeDamping(double alpha, double dist, double iDist, double sDist);
	double computeDamper(double dist, double iDist, double sDist, double damping);

private:
	struct DampData
	{
		enum State {Low, Upp, Free};

		DampData(double mi, double ma, double miV, double maV,
						 double idi, double sdi, int vp, int i):
			min(mi), max(ma), minVel(miV), maxVel(maV), iDist(idi), sDist(sdi),
			index(i), vecPos(vp), damping(0.), state(Free)
		{}

		double min, max;
		double minVel, maxVel;
		double iDist, sDist;
		int index;
		int vecPos;
		double damping;
		State state;
	};

private:
	std::vector<DampData> data_;

	Eigen::VectorXd lower_, upper_;
	int begin_;
	double step_;
	double damperOff_;
};



class SelfCollisionConstr : public ConstraintFunction<Inequality>
{
public:
	SelfCollisionConstr(const rbd::MultiBody& mb, double step);

	void addCollision(const rbd::MultiBody& mb, int collId,
		int body1Id, SCD::S_Object* body1, const sva::PTransformd& body1T,
		int body2Id, SCD::S_Object* body2, const sva::PTransformd& body2T,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmCollision(int collId);
	std::size_t nrCollisions() const;
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int nrInEq() const;
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& LowerInEq() const;
	virtual const Eigen::VectorXd& UpperInEq() const;

private:
	struct CollData
	{
		enum class DampingType {Hard, Soft, Free};
		CollData(const rbd::MultiBody& mb, int collId,
			int body1Id, SCD::S_Object* body1, const sva::PTransformd& body1T,
			int body2Id, SCD::S_Object* body2, const sva::PTransformd& body2T,
			double di, double ds, double damping, double dampingOff);

		SCD::CD_Pair* pair;
		sva::PTransformd body1T, body2T;
		Eigen::Vector3d normVecDist;
		rbd::Jacobian jacB1, jacB2;
		double di, ds;
		double damping;
		int collId;
		int body1Id, body2Id;
		int body1, body2;
		DampingType dampingType;
		double dampingOff;
	};

private:
	std::vector<CollData> dataVec_;
	double step_;
	int nrVars_;
	int nrActivated_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd AL_, AU_;

	Eigen::MatrixXd fullJac_;
	Eigen::MatrixXd fullJacDot_;
	Eigen::VectorXd alphaVec_;
	Eigen::VectorXd calcVec_;
};



class StaticEnvCollisionConstr : public ConstraintFunction<Inequality>
{
public:
	StaticEnvCollisionConstr(const rbd::MultiBody& mb, double step);

	void addCollision(const rbd::MultiBody& mb, int collId,
		int bodyId, SCD::S_Object* body, const sva::PTransformd& bodyT,
		int envId, SCD::S_Object* env,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmCollision(int collId);
	std::size_t nrCollisions() const;
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int nrInEq() const;
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& LowerInEq() const;
	virtual const Eigen::VectorXd& UpperInEq() const;

private:
	struct CollData
	{
		enum class DampingType {Hard, Soft, Free};
		CollData(const rbd::MultiBody& mb, int collId,
			int bodyId, SCD::S_Object* body, const sva::PTransformd& bodyT,
			int envId, SCD::S_Object* env,
			double di, double ds, double damping, double dampingOff);

		SCD::CD_Pair* pair;
		sva::PTransformd bodyT;
		Eigen::Vector3d normVecDist;
		rbd::Jacobian jacB1;
		double di, ds;
		double damping;
		int collId;
		int bodyId, envId;
		int body;
		DampingType dampingType;
		double dampingOff;
	};

private:
	std::vector<CollData> dataVec_;
	double step_;
	int nrVars_;
	int nrActivated_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd AL_, AU_;

	Eigen::MatrixXd fullJac_;
	Eigen::MatrixXd fullJacDot_;
	Eigen::VectorXd alphaVec_;
	Eigen::VectorXd calcVec_;
};

class CoMCollisionConstr : public ConstraintFunction<Inequality>
{
public:
	CoMCollisionConstr(const rbd::MultiBody& mb, double step);

	void addCollision(const rbd::MultiBody& mb,
		int collId, SCD::S_Object* env,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmCollision(int collId);
	std::size_t nrCollisions() const;
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int nrInEq() const;
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& LowerInEq() const;
	virtual const Eigen::VectorXd& UpperInEq() const;

private:
	struct CollData
	{
		enum class DampingType {Hard, Soft, Free};
		CollData(const rbd::MultiBody& mb, int collId,
			SCD::S_Object* env,
			double di, double ds, double damping, double dampingOff);
		SCD::S_Sphere comSphere_;
		SCD::CD_Pair* pair;
		Eigen::Vector3d normVecDist;
		rbd::CoMJacobian jacCoM;
		double di, ds;
		double damping;
		int collId;
		DampingType dampingType;
		double dampingOff;
	};

private:
	std::vector<CollData> dataVec_;
	double step_;
	int nrVars_;
	int nrActivated_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd AL_, AU_;

	Eigen::MatrixXd fullJac_;
	Eigen::MatrixXd fullJacDot_;
	Eigen::VectorXd alphaVec_;
	Eigen::VectorXd calcVec_;
};



class GripperTorqueConstr : public ConstraintFunction<Inequality>
{
public:
	GripperTorqueConstr();

	void addGripper(int bodyId, double torqueLimit,
		const Eigen::Vector3d& origin, const Eigen::Vector3d& axis);
	bool rmGripper(int bodyId);
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& LowerInEq() const;
	virtual const Eigen::VectorXd& UpperInEq() const;

private:
	struct GripperData
	{
		GripperData(int bId, double tl,
			const Eigen::Vector3d& o, const Eigen::Vector3d& a);

		int bodyId;
		double torqueLimit;
		Eigen::Vector3d origin;
		Eigen::Vector3d axis;
	};

private:
	std::vector<GripperData> dataVec_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd AL_, AU_;
};

} // namespace qp

} // namespace tasks

