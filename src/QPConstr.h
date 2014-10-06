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
class QBound;
class AlphaBound;

namespace qp
{

sch::Matrix4x4 tosch(const sva::PTransformd& t);


class ContactConstrCommon
{
public:
	bool addVirtualContact(const ContactId& contactId);
	bool removeVirtualContact(const ContactId& contactId);
	void resetVirtualContacts();

	bool addDofContact(const ContactId& contactId, const Eigen::MatrixXd& dof);
	bool removeDofContact(const ContactId& contactId);
	void resetDofContacts();

protected:
	struct ContactCommon
	{
		ContactId cId;
		sva::PTransformd X_b1_b2;

		bool operator==(const ContactCommon& cc) const;
		bool operator<(const ContactCommon& cc) const;
	};

protected:
	std::set<ContactCommon> contactCommonInContact(
		const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

protected:
	std::set<ContactId> virtualContacts_;
	std::map<ContactId, Eigen::MatrixXd> dofContacts_;
};


class ContactAccConstr : public ConstraintFunction<Equality>,
	public ContactConstrCommon
{
public:
	ContactAccConstr();

	void updateDofContacts();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameEq() const;
	virtual std::string descEq(const std::vector<rbd::MultiBody>& mb, int line);

	// Inequality Constraint
	virtual int nrEq() const;
	virtual int maxEq() const;

	virtual const Eigen::MatrixXd& AEq() const;
	virtual const Eigen::VectorXd& bEq() const;

private:
	struct ContactSideData
	{
		ContactSideData(int rI, int aDB, double s, const rbd::Jacobian& j):
			robotIndex(rI), alphaDBegin(aDB), sign(s), jac(j)
		{}

		int robotIndex, alphaDBegin;
		double sign;
		rbd::Jacobian jac;
	};

	struct ContactData
	{
		ContactData(std::vector<ContactSideData> csds,
			const Eigen::MatrixXd& d, const ContactId& cId):
			contacts(std::move(csds)),
			dof(d),
			contactId(cId)
		{}

		std::vector<ContactSideData> contacts;
		Eigen::MatrixXd dof;
		ContactId contactId;
	};

private:
	void updateNrEq();

private:
	std::vector<ContactData> cont_;

	std::vector<Eigen::MatrixXd> fullJac_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd b_;

	int nrEq_;
};



class ContactSpeedConstr : public ConstraintFunction<Equality>,
	public ContactConstrCommon
{
public:
	ContactSpeedConstr(double timeStep);

	void updateDofContacts();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameEq() const;
	virtual std::string descEq(const std::vector<rbd::MultiBody>& mbs, int line);

	// Inequality Constraint
	virtual int nrEq() const;
	virtual int maxEq() const;

	virtual const Eigen::MatrixXd& AEq() const;
	virtual const Eigen::VectorXd& bEq() const;

private:
	struct ContactSideData
	{
		ContactSideData(int rI, int aDB, double s, const rbd::Jacobian& j):
			robotIndex(rI), alphaDBegin(aDB), bodyIndex(j.jointsPath().back()),
			sign(s), jac(j)
		{}

		int robotIndex, alphaDBegin, bodyIndex;
		double sign;
		rbd::Jacobian jac;
	};

	struct ContactData
	{
		ContactData(std::vector<ContactSideData> csds,
			const Eigen::MatrixXd& d, const ContactId& cId):
			contacts(std::move(csds)),
			dof(d),
			contactId(cId)
		{}

		std::vector<ContactSideData> contacts;
		Eigen::MatrixXd dof;
		ContactId contactId;
	};

private:
	void updateNrEq();

private:
	std::vector<ContactData> cont_;

	std::vector<Eigen::MatrixXd> fullJac_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd b_;

	int nrEq_;
	double timeStep_;
};



class JointLimitsConstr : public ConstraintFunction<Bound>
{
public:
	JointLimitsConstr(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		QBound bound, double step);

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameBound() const;
	virtual std::string descBound(const std::vector<rbd::MultiBody>& mbs, int line);

	// Bound Constraint
	virtual int beginVar() const;

	virtual const Eigen::VectorXd& Lower() const;
	virtual const Eigen::VectorXd& Upper() const;

private:
	int robotIndex_, alphaDBegin_, alphaDOffset_;
	double step_;
	Eigen::VectorXd qMin_, qMax_;
	Eigen::VectorXd qVec_, alphaVec_;
	Eigen::VectorXd lower_, upper_;
};



class DamperJointLimitsConstr : public ConstraintFunction<Bound>
{
public:
	DamperJointLimitsConstr(const std::vector<rbd::MultiBody>& mbs,
		int robotIndex, const QBound& qBound, const AlphaBound& aBound,
		double interPercent, double securityPercent, double damperOffset, double step);

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameBound() const;
	virtual std::string descBound(const std::vector<rbd::MultiBody>& mbs,
		int line);

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
						 double idi, double sdi, int aDB, int i):
			min(mi), max(ma), minVel(miV), maxVel(maV), iDist(idi), sDist(sdi),
			jointIndex(i), alphaDBegin(aDB), damping(0.), state(Free)
		{}

		double min, max;
		double minVel, maxVel;
		double iDist, sDist;
		int jointIndex;
		int alphaDBegin;
		double damping;
		State state;
	};

private:
	int robotIndex_, alphaDBegin_;
	std::vector<DampData> data_;

	Eigen::VectorXd lower_, upper_;
	double step_;
	double damperOff_;
};



class CollisionConstr : public ConstraintFunction<Inequality>
{
public:
	CollisionConstr(const std::vector<rbd::MultiBody>& mbs, double step);

	void addCollision(const std::vector<rbd::MultiBody>& mbs, int collId,
		int r1Index, int r1BodyId,
		sch::S_Object* body1, const sva::PTransformd& X_op1_o1,
		int r2Index, int r2BodyId,
		sch::S_Object* body2, const sva::PTransformd& X_op2_o2,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmCollision(int collId);
	std::size_t nrCollisions() const;
	void reset();

	void updateNrCollisions();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const std::vector<rbd::MultiBody>& mbs, int line);

	// InInequality Constraint
	virtual int nrInEq() const;
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& bInEq() const;

private:
	struct BodyCollData
	{
		BodyCollData(const rbd::MultiBody& mb,
			int rIndex, int bodyId, sch::S_Object* hull,
			const sva::PTransformd& X_op_o);

		sch::S_Object* hull;
		rbd::Jacobian jac;
		sva::PTransformd X_op_o;
		int rIndex, bIndex, bodyId;
	};

	struct CollData
	{
		enum class DampingType {Hard, Soft, Free};
		CollData(std::vector<BodyCollData> bcds, int collId,
			sch::S_Object* body1, sch::S_Object* body2,
			double di, double ds, double damping, double dampingOff);

		sch::CD_Pair* pair;
		Eigen::Vector3d normVecDist;
		double di, ds;
		double damping;
		std::vector<BodyCollData> bodies;

		DampingType dampingType;
		double dampingOff;
		int collId;
	};

private:
	double computeDamping(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs, const CollData& cd,
		const Eigen::Vector3d& normalVecDist, double dist) const;

private:
	std::vector<CollData> dataVec_;
	double step_;
	int nrActivated_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd bInEq_;

	std::vector<Eigen::MatrixXd> fullJac_;

	int nrVars_;
};



class CoMIncPlaneConstr : public ConstraintFunction<Inequality>
{
public:
	CoMIncPlaneConstr(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		double step);

	void addPlane(
		int planeId, const Eigen::Vector3d& normal, double offset,
		double di, double ds, double damping, double dampingOff=0.);
	bool rmPlane(int planeId);
	std::size_t nrPlanes() const;
	void reset();

	void updateNrPlanes();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const std::vector<rbd::MultiBody>& mbs, int line);

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
	int robotIndex_, alphaDBegin_;
	std::vector<PlaneData> dataVec_;
	double step_;
	int nrVars_;
	int nrActivated_;
	std::vector<std::size_t> activated_;

	rbd::CoMJacobian jacCoM_;
	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd bInEq_;
};



class GripperTorqueConstr : public ConstraintFunction<Inequality>
{
public:
	GripperTorqueConstr();

	void addGripper(const ContactId& cId, double torqueLimit,
		const Eigen::Vector3d& origin, const Eigen::Vector3d& axis);
	bool rmGripper(const ContactId& cId);
	void reset();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mb,
		const SolverData& data);

	virtual void update(const std::vector<rbd::MultiBody>& mb,
		const std::vector<rbd::MultiBodyConfig>& mbc,
		const SolverData& data);

	virtual std::string nameInEq() const;
	virtual std::string descInEq(const std::vector<rbd::MultiBody>& mb, int line);

	// InInequality Constraint
	virtual int maxInEq() const;

	virtual const Eigen::MatrixXd& AInEq() const;
	virtual const Eigen::VectorXd& bInEq() const;

private:
	struct GripperData
	{
		GripperData(const ContactId& cId, double tl,
			const Eigen::Vector3d& o, const Eigen::Vector3d& a);

		ContactId contactId;
		double torqueLimit;
		Eigen::Vector3d origin;
		Eigen::Vector3d axis;
	};

private:
	std::vector<GripperData> dataVec_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd bInEq_;
};



//class ConstantSpeedConstr : public ConstraintFunction<Equality>
//{
//public:
//	ConstantSpeedConstr(const rbd::MultiBody& mb, double timeStep);

//	void addConstantSpeed(const rbd::MultiBody& mb, int bodyId,
//											const Eigen::Vector3d& bodyPoint,
//											const Eigen::MatrixXd& dof,
//											const Eigen::VectorXd& speed);
//	bool removeConstantSpeed(int bodyId);
//	void resetConstantSpeed();
//	std::size_t nrConstantSpeed() const;

//	// Constraint
//	virtual void updateNrVars(const rbd::MultiBody& mb,
//		const SolverData& data);

//	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
//		const SolverData& data);

//	virtual std::string nameEq() const;
//	virtual std::string descEq(const rbd::MultiBody& mb, int line);

//	// Inequality Constraint
//	virtual int maxEq() const;

//	virtual const Eigen::MatrixXd& AEq() const;
//	virtual const Eigen::VectorXd& bEq() const;

//private:
//	struct ConstantSpeedData
//	{
//		ConstantSpeedData(rbd::Jacobian j, const Eigen::MatrixXd& d,
//										 const Eigen::VectorXd& s, int bId):
//			jac(j),
//			bodyPoint(j.point()),
//			dof(d),
//			speed(s),
//			body(j.jointsPath().back()),
//			bodyId(bId)
//		{}

//		rbd::Jacobian jac;
//		sva::PTransformd bodyPoint;
//		Eigen::MatrixXd dof;
//		Eigen::VectorXd speed;
//		int body;
//		int bodyId;
//	};

//private:
//	void updateNrEq();

//private:
//	std::vector<ConstantSpeedData> cont_;

//	Eigen::MatrixXd fullJac_;

//	Eigen::MatrixXd A_;
//	Eigen::VectorXd b_;

//	int nrVars_;
//	double timeStep_;
//};

} // namespace qp

} // namespace tasks

