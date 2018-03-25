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
// Eigen
#include <Eigen/Core>
#include <Eigen/StdVector>

// RBDyn
#include <RBDyn/Jacobian.h>
#include <RBDyn/CoM.h>

// sch
#include <sch/Matrix/SCH_Types.h>

// Tasks
#include "QPSolver.h"

// unique_ptr
#include <memory>

// forward declaration
// sch
namespace sch
{
class S_Object;
class CD_Pair;
}


namespace tasks
{
struct QBound;
struct AlphaBound;

namespace qp
{


/// Convert a sch-core transformation matrix to a sva::PTransformd matrix.
TASKS_DLLAPI sch::Matrix4x4 tosch(const sva::PTransformd& t);


/**
	* Avoid to reach articular position limits based on direct integration.
	* \f[
	* \underline{q} \leq
	* q + \alpha \Delta_{step} + \frac{1}{2} \dot{\alpha} \Delta_{step}^2
	* \leq \overline{q}
	* \f]
	* This constraint can be impossible to fulfill when articulation velocity
	* is really high. Always prefer to use DamperJointLimitsConstr.
	*/
class TASKS_DLLAPI JointLimitsConstr : public ConstraintFunction<Bound>
{
public:
	/**
		* @param mbs Multi-robot system.
		* @param robotIndex Constrained robot Index in mbs.
		* @param bound Articular position bounds.
		* @param step Time step in second.
		*/
	JointLimitsConstr(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		QBound bound, double step);

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual std::string nameBound() const override;
	virtual std::string descBound(const std::vector<rbd::MultiBody>& mbs, int line) override;

	// Bound Constraint
	virtual int beginVar() const override;

	virtual const Eigen::VectorXd& Lower() const override;
	virtual const Eigen::VectorXd& Upper() const override;

private:
	int robotIndex_, alphaDBegin_, alphaDOffset_;
	double step_;
	Eigen::VectorXd qMin_, qMax_;
	Eigen::VectorXd qVec_, alphaVec_;
	Eigen::VectorXd lower_, upper_;
};



/**
	* Avoid to reach articular position and velocity limits based on
	* a velocity damper.
	* For each articulation \f$ j \f$:
	* \f[
	* \frac{\max(-\xi \frac{\underline{d} - d_s}{d_i - d_s},
	*            \underline{\alpha}_{j})
	* - \alpha_{j}}{\Delta_{dt}}
	* \leq \dot{\alpha}_{j} \leq
	* \frac{\min(\xi \frac{\overline{d} - d_s}{d_i - d_s},
	*            \overline{\alpha}_{j})
	* - \alpha_{j}}{\Delta_{dt}}
	* \f]
	* with \f$ \underline{d} \f$ and \f$ \overline{d} \f$ the distance to the
	* lower and upper articular position bound, \f$ d_i \f$ the interactive
	* distance, \f$ d_s \f$ the security distance and \f$ \xi \f$ the damper.
	*
	* The damper \f$ \xi \f$ is calculated automatically each time
	* the distance \f$ \underline{d} \f$ or \f$ \overline{d} \f$ go below
	* the interactive distance \f$ d_i \f$ with
	* the following formula:
	* \f[ \xi = -\frac{d_i - d_s}{d - d_s}\alpha + \xi_{\text{off}} \f]
	*
	* This behavior give the better stability and we usually use
	* the following values:
	* - \f$ d_i = 0.1 (\overline{q} - \underline{q}) \f$
	* - \f$ d_s = 0.01 (\overline{q} - \underline{q}) \f$
	* - \f$ \xi_{\text{off}} = 0.5 \f$
	*/
class TASKS_DLLAPI DamperJointLimitsConstr : public ConstraintFunction<Bound>
{
public:
	/**
		* @param mbs Multi-robot system.
		* @param robotIndex Constrained robot Index in mbs.
		* @param qBound Articular position bounds.
		* @param aBound Articular velocity bounds.
		* @param interPercent \f$ interPercent (\overline{q} - \underline{q}) \f$
		* @param securityPercent \f$ securityPercent (\overline{q} - \underline{q}) \f$
		* @param damperOffset \f$ \xi_{\text{off}} \f$
		* @param step Time step in second.
		*/
	DamperJointLimitsConstr(const std::vector<rbd::MultiBody>& mbs,
		int robotIndex, const QBound& qBound, const AlphaBound& aBound,
		double interPercent, double securityPercent, double damperOffset, double step);

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual std::string nameBound() const override;
	virtual std::string descBound(const std::vector<rbd::MultiBody>& mbs,
		int line) override;

	// Bound Constraint
	virtual int beginVar() const override;

	virtual const Eigen::VectorXd& Lower() const override;
	virtual const Eigen::VectorXd& Upper() const override;

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

/**
	* Avoid that too robot links enter into collision based on a velocity damper.
	* For each collision pair:
	* \f[
	* \dot{d} + \ddot{d}\Delta_{dt} \geq -\xi \frac{d - d_s}{d_i - d_s}
	* \f]
	* with \f$ d \f$ the minimal distance between the two links,
	* \f$ d_i \f$ the interactive distance, \f$ d_s \f$ the security distance
	* and \f$ \xi \f$ the damper.
	*
	* The damper \f$ \xi \f$ can be calculated automatically each time
	* the distance \f$ d \f$ go below the interactive distance \f$ d_i \f$ with
	* the following formula:
	* \f[ \xi = -\frac{d_i - d_s}{d - d_s}\alpha + \xi_{\text{off}} \f]
	*/
class TASKS_DLLAPI CollisionConstr : public ConstraintFunction<Inequality>
{
public:
	/**
		* @param mbs Multi-robot system.
		* @param step Time step in second.
		*/
	CollisionConstr(const std::vector<rbd::MultiBody>& mbs, double step);

	/**
		* Add a collision avoidance constraint.
		* Don't forget to call updateNrCollisions and QPSolver::updateConstrSize.
		* You can also only call QPSolver::nrVars or QPSolver::updateConstrsNrVars
		* or QPSolver::updateNrVars.
		*
		* @param mbs Multi-robot system (must be the same given in the constructor.
		* @param collId Id of this collision, must be unique.
		* @param r1Index First constrained robot Index in mbs.
		* @param r1BodyId Constrained body id in mbs[r1Index].
		* @param body1 sch-core hull associated to the r1BodyId link.
		* @param X_op1_o1 body1 position will be set at each iteration to
		* \f$ {}^{o1}X_{op1} {}^{r1BodyId}X_O \f$.
		* @param r2Index Second constrained robot Index in mbs
		* (can be equal to r1Index).
		* @param r2BodyId Constrained body id in mbs[r2Index].
		* @param body2 sch-core hull associated to the r2BodyId link.
		* @param X_op2_o2 body2 position will be set at each iteration to
		* \f$ {}^{o2}X_{op2} {}^{r2BodyId}X_O \f$.
		* @param di \f$ d_i \f$.
		* @param ds \f$ d_s \f$.
		* @param damping \f$ \xi \f$, if set to 0 the damping is computed automatically.
		* @param dampingOff \f$ \xi_{\text{off}} \f$.
		*/
	void addCollision(const std::vector<rbd::MultiBody>& mbs, int collId,
		int r1Index, const std::string& r1BodyName,
		sch::S_Object* body1, const sva::PTransformd& X_op1_o1,
		int r2Index, const std::string& r2BodyName,
		sch::S_Object* body2, const sva::PTransformd& X_op2_o2,
		double di, double ds, double damping, double dampingOff=0.);

	/**
		* Remove a collision avoidance constraint.
		* @param collId Collision id to remove.
		* @return true if the collision as been removed false if the collision id
		* was associated with no collision.
		*/
	bool rmCollision(int collId);

	/// @return Number of collision constraint.
	std::size_t nrCollisions() const;

	/// Remove all collision constraints.
	void reset();

	/// Reallocate A and b matrix.
	void updateNrCollisions();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual std::string nameInEq() const override;
	virtual std::string descInEq(const std::vector<rbd::MultiBody>& mbs, int line) override;

	// In Inequality Constraint
	virtual int nrInEq() const override;
	virtual int maxInEq() const override;

	virtual const Eigen::MatrixXd& AInEq() const override;
	virtual const Eigen::VectorXd& bInEq() const override;

private:
	struct BodyCollData
	{
		BodyCollData(const rbd::MultiBody& mb,
			int rIndex, const std::string& bodyName,
			sch::S_Object* hull, const sva::PTransformd& X_op_o);

		sch::S_Object* hull;
		rbd::Jacobian jac;
		sva::PTransformd X_op_o;
		int rIndex, bIndex;
		std::string bodyName;
	};

	struct CollData
	{
		enum class DampingType {Hard, Soft, Free};
		CollData(std::vector<BodyCollData> bcds, int collId,
			sch::S_Object* body1, sch::S_Object* body2,
			double di, double ds, double damping, double dampingOff);
		std::unique_ptr<sch::CD_Pair> pair;
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
	int nrActivated_, totalAlphaD_;

	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd bInEq_;

	Eigen::MatrixXd fullJac_, distJac_;

	int nrVars_;

  CollisionConstr(const CollisionConstr&) = delete;
  CollisionConstr& operator=(const CollisionConstr&) = delete;
};



/**
	* Prevent robot CoM to go out of a convex hull.
	* For each plane that compose the convex hull:
	* \f[
	* \dot{d} + \ddot{d}\Delta_{dt} \geq -\xi \frac{d - d_s}{d_i - d_s}
	* \f]
	* with \f$ d \f$ the CoM and the plane,
	* \f$ d_i \f$ the interactive distance, \f$ d_s \f$ the security distance
	* and \f$ \xi \f$ the damper.
	*
	* The damper \f$ \xi \f$ can be calculated automatically each time
	* the distance \f$ d \f$ go below the interactive distance \f$ d_i \f$ with
	* the following formula:
	* \f[ \xi = -\frac{d_i - d_s}{d - d_s}\alpha + \xi_{\text{off}} \f]
	*/
class TASKS_DLLAPI CoMIncPlaneConstr : public ConstraintFunction<Inequality>
{
public:
	/**
		* @param mbs Multi-robot system.
		* @param robotIndex Constrained robot Index in mbs.
		* @param step Time step in second.
		*/
	CoMIncPlaneConstr(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
		double step);

	/**
		* Add a plane to the convex hull.
		* Don't forget to call updateNrPlanes and QPSolver::updateConstrSize.
		* You can also only call QPSolver::nrVars or QPSolver::updateConstrsNrVars
		* or QPSolver::updateNrVars.
		*
		* @param planeId Id of this plane, must be unique.
		* @param normal Normal of the plane, must be unitary and in the
		* allowed zone direction.
		* @param offset Origin of the normal.
		* @param di \f$ d_i \f$.
		* @param ds \f$ d_s \f$.
		* @param damping \f$ \xi \f$, if set to 0 the damping is computed automatically.
		* @param dampingOff \f$ \xi_{\text{off}} \f$.
		*/
	void addPlane(
		int planeId, const Eigen::Vector3d& normal, double offset,
		double di, double ds, double damping, double dampingOff=0.);

	void addPlane(
		int planeId, const Eigen::Vector3d& normal, double offset,
		double di, double ds, double damping,
                const Eigen::Vector3d& speed,
                const Eigen::Vector3d& normalDot, double dampingOff=0.);

	/**
		* Remove a plane.
		* @param planeId Plane id to remove.
		* @return true if the plane as been removed false if plane id
		* was associated with no collision.
		*/
	bool rmPlane(int planeId);

	/// @return Number of plane.
	std::size_t nrPlanes() const;

	/// Remove all plane.
	void reset();

	/// Reallocate A and b matrix.
	void updateNrPlanes();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual std::string nameInEq() const override;
	virtual std::string descInEq(const std::vector<rbd::MultiBody>& mbs, int line) override;

	// In Inequality Constraint
	virtual int nrInEq() const override;
	virtual int maxInEq() const override;

	virtual const Eigen::MatrixXd& AInEq() const override;
	virtual const Eigen::VectorXd& bInEq() const override;

private:
	struct PlaneData
	{
		enum class DampingType {Hard, Soft, Free};
		PlaneData(int planeId,
			const Eigen::Vector3d& normal, double offset,
			double di, double ds, double damping, double dampingOff,
                        const Eigen::Vector3d& speed,
                        const Eigen::Vector3d& normalDot);
		Eigen::Vector3d normal;
		Eigen::Vector3d normalDot;
		double offset;
		double dist;
		double di, ds;
		double damping;
		int planeId;
		DampingType dampingType;
		double dampingOff;
                Eigen::Vector3d speed;
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



class TASKS_DLLAPI GripperTorqueConstr : public ConstraintFunction<Inequality>
{
public:
	GripperTorqueConstr();

	void addGripper(const ContactId& cId, double torqueLimit,
		const Eigen::Vector3d& origin, const Eigen::Vector3d& axis);
	bool rmGripper(const ContactId& cId);
	void reset();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mb,
		const SolverData& data) override;

	virtual void update(const std::vector<rbd::MultiBody>& mb,
		const std::vector<rbd::MultiBodyConfig>& mbc,
		const SolverData& data) override;

	virtual std::string nameInEq() const override;
	virtual std::string descInEq(const std::vector<rbd::MultiBody>& mb, int line) override;

	// In Inequality Constraint
	virtual int maxInEq() const override;

	virtual const Eigen::MatrixXd& AInEq() const override;
	virtual const Eigen::VectorXd& bInEq() const override;

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


/**
	* Constraint some link velocity in link frame by direct integration.
	* \f[
	* S (\underline{v} - v)
	* \leq S a \Delta_{dt} \leq
	* S (\overline{v} - v)
	* \f]
	* with \f$ \underline{v} \f$ and \f$ \overline{v} \f$ the velocity lower and
	* upper bound, \f$ v \f$ the link velocity in link frame,
	* \f$ a \f$ the link acceleration in link frame
	* and \f$ S \f$ the axis selection matrix.
	*/
class TASKS_DLLAPI BoundedSpeedConstr : public ConstraintFunction<GenInequality>
{
public:
	/**
		* @param mbs Multi-robot system.
		* @param robotIndex Constrained robot Index in mbs.
		* @param timeStep Time step in second.
		*/
	BoundedSpeedConstr(const std::vector<rbd::MultiBody>& mbs,
		int robotIndex, double timeStep);

	/**
		* Add a targeted speed constraint.
		* Don't forget to call updateBoundedSpeeds and QPSolver::updateConstrSize.
		* You can also only call QPSolver::nrVars or QPSolver::updateConstrsNrVars
		* or QPSolver::updateNrVars.
		* @param mbs Multi-robot system (must be the same given in the constructor.
		* @param bodyId Constrained body id in mbs[robotIndex].
		* @param bodyPoint static translation applied on the link
		* \f$ v = xlt(bodyPoint) v_{bodyId} \f$.
		* @param dof \f$ S \in \mathbb{R}^{n \times 6} \f$.
		* @param speed Targeted velocity
		* \f$ \underline{v} \in \mathbb{R}^{n} \f$ and
		* \f$ \overline{v} \in \mathbb{R}^{n} \f$.
		*/
	void addBoundedSpeed(const std::vector<rbd::MultiBody>& mbs, 
		const std::string& bodyName,
		const Eigen::Vector3d& bodyPoint, const Eigen::MatrixXd& dof,
		const Eigen::VectorXd& speed);

	/**
		* Add a bounded speed constraint.
		* Don't forget to call updateBoundedSpeeds and QPSolver::updateConstrSize.
		* You can also only call QPSolver::nrVars or QPSolver::updateConstrsNrVars
		* or QPSolver::updateNrVars.
		* @param mbs Multi-robot system (must be the same given in the constructor.
		* @param bodyId Constrained body id in mbs[robotIndex].
		* @param bodyPoint static translation applied on the link
		* \f$ v = xlt(bodyPoint) v_{bodyId} \f$.
		* @param dof \f$ S \in \mathbb{R}^{n \times 6} \f$.
		* @param lowerSpeed Lower velocity \f$ \underline{v} \in \mathbb{R}^{n} \f$.
		* @param upperSpeed Upper velocity \f$ \overline{v} \in \mathbb{R}^{n} \f$.
		*/
	void addBoundedSpeed(const std::vector<rbd::MultiBody>& mbs,
		const std::string& bodyName,
		const Eigen::Vector3d& bodyPoint, const Eigen::MatrixXd& dof,
		const Eigen::VectorXd& lowerSpeed, const Eigen::VectorXd& upperSpeed);

	/**
		* Remove a bounded speed constraint.
		* @param bodyName Remove the constraint applied on bodyName.
		* @return true if the constraint has been removed false if bodyId
		* was associated with no constraint.
		*/
	bool removeBoundedSpeed(const std::string& bodyName);

	/// Remove all bounded speed constraint.
	void resetBoundedSpeeds();

	/// @return Number of bounded speed constraint.
	std::size_t nrBoundedSpeeds() const;

	/// Reallocate A and b matrix.
	void updateBoundedSpeeds();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbc,
		const SolverData& data) override;

	virtual std::string nameGenInEq() const override;
	virtual std::string descGenInEq(const std::vector<rbd::MultiBody>& mb, int line) override;

	// Inequality Constraint
	virtual int maxGenInEq() const override;

	virtual const Eigen::MatrixXd& AGenInEq() const override;
	virtual const Eigen::VectorXd& LowerGenInEq() const override;
	virtual const Eigen::VectorXd& UpperGenInEq() const override;

private:
	struct BoundedSpeedData
	{
		BoundedSpeedData(rbd::Jacobian j, const Eigen::MatrixXd& d,
			const Eigen::VectorXd& ls, const Eigen::VectorXd& us,
			const std::string& bName):
			jac(j),
			bodyPoint(j.point()),
			dof(d),
			lSpeed(ls),
			uSpeed(us),
			body(j.jointsPath().back()),
			bodyName(bName)
		{}

		rbd::Jacobian jac;
		sva::PTransformd bodyPoint;
		Eigen::MatrixXd dof;
		Eigen::VectorXd lSpeed, uSpeed;
		int body;
		std::string bodyName;
	};

private:
	void updateNrEq();

private:
	int robotIndex_, alphaDBegin_;
	std::vector<BoundedSpeedData> cont_;

	Eigen::MatrixXd fullJac_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd lower_, upper_;

	int nrVars_;
	double timeStep_;
};



/**
	* Generalized Image Point Constraint
	* This can be used as either a Field of View constraint
	* or an Occlusion avoidance constraint. Both with damping
	* Note that the point must start within the feasible boundary.
	*/
class TASKS_DLLAPI ImageConstr : public ConstraintFunction<Inequality>
{
public:
	/**
		* @param mbs Multi-robot system.
		* @param robotIndex the index of the robot the constraint is used on
		* @param bName name of the body containing the camera
		* @param X_b_gaze pose of the camera frame relative to the body containing it
		* @param step Time step in second.
		* @param constrDirection used to switch between a field of view constraint (point
		* 	inside the image limits) and occlusion avoidance (point outside the limits) 1 by
		*		default corresponds to FoV, set it to -1 for occlusion avoidance
		*/
	ImageConstr(const std::vector<rbd::MultiBody>& mbs, int robotIndex, const std::string &bName,
		const sva::PTransformd& X_b_gaze, double step, double constrDirection=1.);

	/** Copy constructor */
	ImageConstr(const ImageConstr& rhs);

	/** Move constructor */
	ImageConstr(ImageConstr&&) = default;

	/** Copy assignment operator */
	ImageConstr& operator=(const ImageConstr& rhs);

	/** Move assignment operator */
	ImageConstr& operator=(ImageConstr&&) = default;

	/**
	 * @brief setLimits
	 * @param min
	 * @param max
	 * @param iPercent
	 * @param sPercent
	 * @param damping
	 */
	void setLimits(const Eigen::Vector2d& min, const Eigen::Vector2d& max, const double iPercent,
		const double sPercent, const double damping, const double dampingOffsetPercent);

	/**
	 * @brief addPoint
	 * @param point2d
	 * @param depthEstimate
	 * @return
	 */
	int addPoint(const Eigen::Vector2d& point2d, const double depthEstimate);
	int addPoint(const Eigen::Vector3d& point3d);

	/**
	 * @brief addPoint - overload for adding a self point
	 * @param mbs
	 * @param bodyId
	 * @param X_b_p
	 */
	void addPoint(const std::vector<rbd::MultiBody>& mbs, const std::string &bName,
		const sva::PTransformd& X_b_p=sva::PTransformd::Identity());
	void reset();
	void updatePoint(const int pointId, const Eigen::Vector2d& point2d);
	void updatePoint(const int pointId, const Eigen::Vector2d& point2d, const double depthEstimate);
	void updatePoint(const int pointId, const Eigen::Vector3d& point3d);

	void computeComponents(const rbd::MultiBody &mb, const rbd::MultiBodyConfig &mbc,
		const SolverData &data, const Eigen::Vector2d& point2d, const double depth,
		rbd::Jacobian &jac, const int bodyIndex, const sva::PTransformd& X_b_p,
		Eigen::MatrixXd& fullJacobian, Eigen::Vector2d& bCommonTerm);

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	// In Inequality Constraint
	virtual std::string nameInEq() const override;
	virtual std::string descInEq(const std::vector<rbd::MultiBody>& mbs, int line) override;
	virtual int nrInEq() const override;
	virtual int maxInEq() const override;

	virtual const Eigen::MatrixXd& AInEq() const override;
	virtual const Eigen::VectorXd& bInEq() const override;

private:
	struct PointData
	{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		PointData(const Eigen::Vector2d& pt, const double d);
		Eigen::Vector2d point2d;
		double depthEstimate;
	};
	struct RobotPointData
	{
		RobotPointData(const std::string &bName, const sva::PTransformd& X, const rbd::Jacobian& j);
		std::string bName;
		sva::PTransformd X_b_p;
		rbd::Jacobian jac;
	};

private:
	std::vector<PointData, Eigen::aligned_allocator<PointData>> dataVec_;
	std::vector<RobotPointData> dataVecRob_;
	int robotIndex_, bodyIndex_, alphaDBegin_;
	int nrVars_;
	double step_, accelFactor_;
	int nrActivated_;

	rbd::Jacobian jac_;
	sva::PTransformd X_b_gaze_;
	std::unique_ptr<Eigen::Matrix<double, 2, 6>> L_img_;
	std::unique_ptr<Eigen::Matrix<double, 6, 1>> surfaceVelocity_;
	std::unique_ptr<Eigen::Matrix<double, 1, 6>> L_Z_dot_;
	std::unique_ptr<Eigen::Matrix<double, 2, 6>> L_img_dot_;
	std::unique_ptr<Eigen::Vector2d> speed_;
	std::unique_ptr<Eigen::Vector2d> normalAcc_;
	Eigen::MatrixXd jacMat_;
	std::unique_ptr<Eigen::Vector2d> iDistMin_, iDistMax_, sDistMin_, sDistMax_;
	double damping_, dampingOffset_, ineqInversion_, constrDirection_;
	Eigen::MatrixXd AInEq_;
	Eigen::VectorXd bInEq_;
};


} // namespace qp

} // namespace tasks

