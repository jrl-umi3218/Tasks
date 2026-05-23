/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
// std
#include <vector>

// Eigen
#include <Eigen/Core>

// Tasks
#include "QPSolver.h"

namespace tasks
{

namespace qp
{

/**
 * Base class for coincidence equality constraints.
 *
 * This constraint enforces acceleration-level coincidence between two bodies
 * of the same robot using their Jacobians and normal accelerations.
 *
 * A proportional-derivative (PD) stabilization term is added in the update step
 * to reduce position and velocity errors between the bodies.
 *
 * Derived classes define whether the coincidence applies to translation,
 * rotation, or both.
 */
class TASKS_DLLAPI CoincidenceConstr : public ConstraintFunction<Equality>
{
public:
  /**
   * @param robotIndex Constrained robot Index in mbs.
   * @param body1Name Name of the first body.
   * @param body2Name Name of the second body.
   * @param jointSelector Joint selection vector.
   */
  CoincidenceConstr(int robotIndex,
                    const std::string & body1Name,
                    const std::string & body2Name,
                    const Eigen::VectorXd & jointSelector)
  : robotIndex_(robotIndex), body1Name_(body1Name), body2Name_(body2Name), body1Index_(-1), body2Index_(-1),
    jointSelector_(jointSelector)
  {
  }

  // Equality constraint
  virtual int maxEq() const override;

  virtual const Eigen::MatrixXd & AEq() const override;

  virtual const Eigen::VectorXd & bEq() const override;

  /**
   * Set the joint selector.
   * @param selector Joint selection vector.
   */
  void setJointSelector(const Eigen::VectorXd & selector);

protected:
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;
  int robotIndex_;
  std::string body1Name_, body2Name_;
  int body1Index_, body2Index_;
  Eigen::VectorXd jointSelector_;
};

/**
 * Enforce translational coincidence between two bodies.
 *
 * This constraint enforces equality of linear accelerations between two bodies
 * (6D formulation), using the difference of their full Jacobians.
 *
 * A PD stabilization term is added based on the difference of body world
 * positions and linear velocities.
 */
class TASKS_DLLAPI FixedCoincidenceConstr : public CoincidenceConstr
{
public:
  /**
   * @param robotIndex Constrained robot Index in mbs.
   * @param body1Name Name of the first body.
   * @param body2Name Name of the second body.
   * @param point1 Unused (kept for API compatibility).
   * @param point2 Unused (kept for API compatibility).
   * @param jointSelector Joint selection vector.
   */
  FixedCoincidenceConstr(int robotIndex = 0,
                         const std::string & body1Name = "",
                         const std::string & body2Name = "",
                         const Eigen::Vector3d & point1 = Eigen::Vector3d::Zero(),
                         const Eigen::Vector3d & point2 = Eigen::Vector3d::Zero(),
                         const Eigen::VectorXd & jointSelector = Eigen::VectorXd());

  // Equality constraint
  virtual std::string nameEq() const override;

  virtual std::string descEq(const std::vector<rbd::MultiBody> & mbs, int line) override;

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

private:
  Eigen::Vector3d point1_, point2_;
};

/**
 * Enforce rotational coincidence between two bodies.
 *
 * This constraint enforces equality of angular accelerations between two bodies
 * (3D formulation), using the angular part of their Jacobians.
 *
 * A PD stabilization term is added based on the difference of body world
 * orientations and angular velocities.
 */
class TASKS_DLLAPI RotationalCoincidenceConstr : public CoincidenceConstr
{
public:
  /**
   * @param robotIndex Constrained robot Index in mbs.
   * @param body1Name Name of the first body.
   * @param body2Name Name of the second body.
   * @param jointSelector Joint selection vector.
   */
  RotationalCoincidenceConstr(int robotIndex = 0,
                              const std::string & body1Name = "",
                              const std::string & body2Name = "",
                              const Eigen::VectorXd & jointSelector = Eigen::VectorXd());

  // Equality constraint
  virtual std::string nameEq() const override;

  virtual std::string descEq(const std::vector<rbd::MultiBody> & mbs, int line) override;

  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;
};

} // namespace qp

} // namespace tasks
