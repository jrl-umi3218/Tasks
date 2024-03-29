/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
// std
#include <map>
#include <set>

// Eigen
#include <Eigen/Core>

// RBDyn
#include <RBDyn/Jacobian.h>

// Tasks
#include "QPSolver.h"

namespace tasks
{

namespace qp
{

/**
 * Manage modifier on contact constraint.
 */
class TASKS_DLLAPI ContactConstrCommon
{
public:
  /**
   * Add a virtual contact.
   * A virtual contact don't have any motion constraint apply on it (move freely).
   * It can be seen like a Dof contact with \f$ S \in \mathbb{R}^{0 \times 6} \f$
   */
  bool addVirtualContact(const ContactId & contactId);

  /// Remove a virtual contact.
  bool removeVirtualContact(const ContactId & contactId);

  /// Remove all virtual contact.
  void resetVirtualContacts();

  /**
   * Free some degree of freedom of a contact like in the BoundedSpeedConstr.
   * \f[ \bar{v} = S v \f]
   * with \f$ v \f$ the velocity in the UnilateralContact and BilateralContact
   * \f$ cf \f$ frame.
   * @param dof \f$ S \in \mathbb{R}^{n \times 6} \f$
   * @see BoundedSpeedConstr
   * @see ContactConstr::updateDofContacts
   */
  bool addDofContact(const ContactId & contactId, const Eigen::MatrixXd & dof);

  /**
   * Remove a Dof contact.
   * @see ContactConstr::updateDofContacts
   */
  bool removeDofContact(const ContactId & contactId);

  /**
   * Check if a DoF has specific DoF
   */
  bool hasDoFContact(const ContactId & id) const;

  /**
   * Get the DoF of a contact
   * @see addDofContact
   */
  const Eigen::MatrixXd & dofContact(const ContactId & contactId);

  /**
   * Remove all Dof contact.
   * @see ContactConstr::updateDofContacts
   */
  void resetDofContacts();

protected:
  struct ContactCommon
  {
    ContactId cId;
    sva::PTransformd X_b1_cf;
    sva::PTransformd X_b1_b2;

    bool operator==(const ContactCommon & cc) const;
    bool operator<(const ContactCommon & cc) const;
  };

protected:
  std::set<ContactCommon> contactCommonInContact(const std::vector<rbd::MultiBody> & mbs, const SolverData & data);

protected:
  std::set<ContactId> virtualContacts_;
  std::map<ContactId, Eigen::MatrixXd> dofContacts_;
};

/**
 * Common contact constraint computation.
 */
class TASKS_DLLAPI ContactConstr : public ConstraintFunction<Equality>, public ContactConstrCommon
{
public:
  ContactConstr();

  /**
   * Update \f$ S \f$ matrix based on Dof contact.
   * You must call this method after calling ContactConstrCommon::addDofContact,
   * ContactConstrCommon::removeDofContact
   * and ContactConstrCommon::resetDofContacts.
   */
  void updateDofContacts();

  // Constraint
  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;

  virtual std::string descEq(const std::vector<rbd::MultiBody> & mbs, int line) override;

  // Inequality Constraint
  virtual int nrEq() const override;
  virtual int maxEq() const override;

  virtual const Eigen::MatrixXd & AEq() const override;
  virtual const Eigen::VectorXd & bEq() const override;

protected:
  struct ContactSideData
  {
    ContactSideData(int rI, int aDB, double s, const rbd::Jacobian & j, const sva::PTransformd & Xbp)
    : robotIndex(rI), alphaDBegin(aDB), bodyIndex(j.jointsPath().back()), sign(s), jac(j), X_b_p(Xbp)
    {
    }

    int robotIndex, alphaDBegin, bodyIndex;
    double sign;
    rbd::Jacobian jac;
    sva::PTransformd X_b_p;
  };

  struct ContactData
  {
    ContactData(std::vector<ContactSideData> csds,
                const Eigen::MatrixXd & d,
                int r1,
                int r2,
                int b1,
                int b2,
                const sva::PTransformd & X_bb,
                const sva::PTransformd & X_bcf,
                const ContactId & cId)
    : contacts(std::move(csds)), dof(d), r1Index(r1), r2Index(r2), b1Index(b1), b2Index(b2), X_b1_b2(X_bb),
      X_b1_cf(X_bcf), contactId(cId)
    {
      revDof = Eigen::Matrix6d::Zero();
      // Find which dofs are selected and reverse that selection
      for(int j = 0; j < 6; ++j)
      {
        bool dof_j_active = false;
        for(int r = 0; r < dof.rows(); ++r) { dof_j_active = dof_j_active || dof(r, j) == 1; }
        if(!dof_j_active) { revDof(j, j) = 1; }
      }
    }

    void update(const std::vector<rbd::MultiBodyConfig> & mbcs);

    std::vector<ContactSideData> contacts;
    Eigen::MatrixXd dof;
    Eigen::MatrixXd revDof;
    int r1Index, r2Index;
    int b1Index, b2Index;
    sva::PTransformd X_b1_b2;
    sva::PTransformd X_b1_cf;
    ContactId contactId;
  };

protected:
  void updateNrEq();

protected:
  std::vector<ContactData> cont_;

  Eigen::MatrixXd fullJac_, dofJac_;

  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;

  int nrEq_, totalAlphaD_;
  double timeStep_;
};

/**
 * Contact constraint by targeting a null acceleration.
 * \f[ a = 0 \f]
 * The contact velocity must be null to make this constraint work.
 *
 * This constraint formulation is usually unstable, prefer
 * ContactSpeedConstr and ContactPosConstr.
 */
class TASKS_DLLAPI ContactAccConstr : public ContactConstr
{
public:
  ContactAccConstr();

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual std::string nameEq() const override;
};

/**
 * Contact constraint by targeting a null velocity.
 * \f[ v + a \Delta_{dt} = 0 \f]
 * This constraint formulation is stable for robot/environment contact.
 * For robot/robot contact prefer ContactPosConstr.
 */
class TASKS_DLLAPI ContactSpeedConstr : public ContactConstr
{
public:
  /// @param timeStep Time step in second.
  ContactSpeedConstr(double timeStep);

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual std::string nameEq() const override;

private:
  double timeStep_;
};

/**
 * Contact constraint by targeting a constant frame.
 * \f[ v + a \Delta_{dt} = \frac{\epsilon(q)}{\Delta_{dt}} \f]
 * Where \f$ \epsilon \f$ is the error between the initial
 * and current contact frame.
 *
 * This constraint formulation is stable in all case.
 */
class TASKS_DLLAPI ContactPosConstr : public ContactConstr
{
public:
  /// @param timeStep Time step in second.
  ContactPosConstr(double timeStep);

  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

  virtual std::string nameEq() const override;

private:
  double timeStep_;
};

} // namespace qp

} // namespace tasks
