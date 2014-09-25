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
#include <map>
#include <set>

// Eigen
#include <Eigen/Core>

// RBDyn
#include <RBDyn/FD.h>
#include <RBDyn/Jacobian.h>
#include <RBDyn/CoM.h>

// sch
#include <sch/Matrix/SCH_Types.h>
#include <sch/S_Object/S_Sphere.h>

// Tasks
#include "QPSolver.h"

// forward declaration
// sch
namespace sch
{
class S_Object;
class CD_Pair;
}


namespace tasks
{

namespace qp
{

sch::Matrix4x4 tosch(const sva::PTransformd& t);


class ContactConstrCommon
{
public:
	bool addVirtualContact(int bodyId);
	bool removeVirtualContact(int bodyId);
	void resetVirtualContacts();

	bool addDofContact(int bodyId, const Eigen::MatrixXd& dof);
	bool removeDofContact(int bodyId);
	void resetDofContacts();

protected:
	std::set<int> bodyIdInContact(const rbd::MultiBody& mb,
		const SolverData& data);

protected:
	std::set<int> virtualContacts_;
	std::map<int, Eigen::MatrixXd> dofContacts_;
};


class ContactAccConstr : public ConstraintFunction<Equality>,
	public ContactConstrCommon
{
public:
	ContactAccConstr(const rbd::MultiBody& mb);

	void updateDofContacts();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

	virtual std::string nameEq() const;
	virtual std::string descEq(const rbd::MultiBody& mb, int line);

	// Inequality Constraint
	virtual int nrEq() const;
	virtual int maxEq() const;

	virtual const Eigen::MatrixXd& AEq() const;
	virtual const Eigen::VectorXd& bEq() const;

private:
	struct ContactData
	{
		ContactData(rbd::Jacobian j, const Eigen::MatrixXd& d, int bId):
			jac(j),
			dof(d),
			bodyId(bId)
		{}

		rbd::Jacobian jac;
		Eigen::MatrixXd dof;
		int bodyId;
	};

private:
	void updateNrEq();

private:
	std::vector<ContactData> cont_;

	Eigen::MatrixXd fullJac_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd b_;

	int nrEq_;
	int nrDof_, nrFor_, nrTor_;
};



class ContactSpeedConstr : public ConstraintFunction<Equality>,
	public ContactConstrCommon
{
public:
	ContactSpeedConstr(const rbd::MultiBody& mb, double timeStep);

	void updateDofContacts();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

	virtual std::string nameEq() const;
	virtual std::string descEq(const rbd::MultiBody& mb, int line);

	// Inequality Constraint
	virtual int nrEq() const;
	virtual int maxEq() const;

	virtual const Eigen::MatrixXd& AEq() const;
	virtual const Eigen::VectorXd& bEq() const;

private:
	struct ContactData
	{
		ContactData(rbd::Jacobian j, const Eigen::MatrixXd& d, int bId):
			jac(j),
			dof(d),
			body(j.jointsPath().back()),
			bodyId(bId)
		{}

		rbd::Jacobian jac;
		Eigen::MatrixXd dof;
		int body;
		int bodyId;
	};

private:
	void updateNrEq();

private:
	std::vector<ContactData> cont_;

	Eigen::MatrixXd fullJac_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd b_;

	int nrEq_;
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

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

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

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

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
		int body1Id, sch::S_Object* body1, const sva::PTransformd& body1T,
		int body2Id, sch::S_Object* body2, const sva::PTransformd& body2T,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmCollision(int collId);
	std::size_t nrCollisions() const;
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int nrInEq() const;
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& bInEq() const;

private:
	struct CollData
	{
		enum class DampingType {Hard, Soft, Free};
		CollData(const rbd::MultiBody& mb, int collId,
			int body1Id, sch::S_Object* body1, const sva::PTransformd& body1T,
			int body2Id, sch::S_Object* body2, const sva::PTransformd& body2T,
			double di, double ds, double damping, double dampingOff);

		sch::CD_Pair* pair;
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
	Eigen::VectorXd bInEq_;

	Eigen::MatrixXd fullJac_;
	Eigen::VectorXd calcVec_;
};



class StaticEnvCollisionConstr : public ConstraintFunction<Inequality>
{
public:
	StaticEnvCollisionConstr(const rbd::MultiBody& mb, double step);

	void addCollision(const rbd::MultiBody& mb, int collId,
		int bodyId, sch::S_Object* body, const sva::PTransformd& bodyT,
		int envId, sch::S_Object* env,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmCollision(int collId);
	std::size_t nrCollisions() const;
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int nrInEq() const;
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& bInEq() const;

private:
	struct CollData
	{
		enum class DampingType {Hard, Soft, Free};
		CollData(const rbd::MultiBody& mb, int collId,
			int bodyId, sch::S_Object* body, const sva::PTransformd& bodyT,
			int envId, sch::S_Object* env,
			double di, double ds, double damping, double dampingOff);

		sch::CD_Pair* pair;
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
	Eigen::VectorXd bInEq_;

	Eigen::MatrixXd fullJac_;
	Eigen::VectorXd calcVec_;
};



class CoMCollisionConstr : public ConstraintFunction<Inequality>
{
public:
	CoMCollisionConstr(const rbd::MultiBody& mb, double step);

	void addCollision(const rbd::MultiBody& mb,
		int collId, sch::S_Object* env,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmCollision(int collId);
	std::size_t nrCollisions() const;
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int nrInEq() const;
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& bInEq() const;

private:
	struct CollData
	{
		enum class DampingType {Hard, Soft, Free};
		CollData(const rbd::MultiBody& mb, int collId,
			sch::S_Object* env,
			double di, double ds, double damping, double dampingOff);
		sch::S_Sphere comSphere_;
		sch::CD_Pair* pair;
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
	Eigen::VectorXd bInEq_;

	Eigen::VectorXd calcVec_;
};



class CoMIncPlaneConstr : public ConstraintFunction<Inequality>
{
public:
	CoMIncPlaneConstr(const rbd::MultiBody& mb, double step);

	void addPlane(
		int planeId, const Eigen::Vector3d& normal, double offset,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmPlane(int planeId);
	std::size_t nrPlanes() const;
	void reset();

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int nrInEq() const;
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& bInEq() const;

private:
	struct PlaneData
	{
		enum class DampingType {Hard, Soft, Free};
		PlaneData(int planeId,
			const Eigen::Vector3d& normal, double offset,
			double di, double ds, double damping, double dampingOff);
		Eigen::Vector3d normal;
		double offset;
		double dist;
		double di, ds;
		double damping;
		int planeId;
		DampingType dampingType;
		double dampingOff;
	};

private:
	std::vector<PlaneData> dataVec_;
	double step_;
	int nrVars_;
	int nrActivated_;
	std::vector<std::size_t> activated_;

	rbd::CoMJacobian jacCoM_;
	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd bInEq_;

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

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const rbd::MultiBody& mb, int line);

	// InInequality Constraint
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& bInEq() const;

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
	Eigen::VectorXd bInEq_;
};



class ConstantSpeedConstr : public ConstraintFunction<Equality>
{
public:
	ConstantSpeedConstr(const rbd::MultiBody& mb, double timeStep);

	void addConstantSpeed(const rbd::MultiBody& mb, int bodyId,
											const Eigen::Vector3d& bodyPoint,
											const Eigen::MatrixXd& dof,
											const Eigen::VectorXd& speed);
	bool removeConstantSpeed(int bodyId);
	void resetConstantSpeed();
	std::size_t nrConstantSpeed() const;

	// Constraint
	virtual void updateNrVars(const rbd::MultiBody& mb,
		const SolverData& data);

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const SolverData& data);

	virtual std::string nameEq() const;
	virtual std::string descEq(const rbd::MultiBody& mb, int line);

	// Inequality Constraint
	virtual int maxEq() const;

	virtual const Eigen::MatrixXd& AEq() const;
	virtual const Eigen::VectorXd& bEq() const;

private:
	struct ConstantSpeedData
	{
		ConstantSpeedData(rbd::Jacobian j, const Eigen::MatrixXd& d,
										 const Eigen::VectorXd& s, int bId):
			jac(j),
			bodyPoint(j.point()),
			dof(d),
			speed(s),
			body(j.jointsPath().back()),
			bodyId(bId)
		{}

		rbd::Jacobian jac;
		sva::PTransformd bodyPoint;
		Eigen::MatrixXd dof;
		Eigen::VectorXd speed;
		int body;
		int bodyId;
	};

private:
	void updateNrEq();

private:
	std::vector<ConstantSpeedData> cont_;

	Eigen::MatrixXd fullJac_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd b_;

	int nrVars_;
	double timeStep_;
};

} // namespace qp

} // namespace tasks

