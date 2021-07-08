/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// associated header
#include "Tasks/QPContactConstr.h"

// includes
// std
#include <set>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// Tasks
#include "utils.h"

namespace tasks
{

namespace qp
{

/**
 *															ContactCommon
 */

bool ContactConstrCommon::ContactCommon::operator==(const ContactCommon & cc) const
{
  return cId == cc.cId;
}

bool ContactConstrCommon::ContactCommon::operator<(const ContactCommon & cc) const
{
  return cId < cc.cId;
}

bool ContactConstrCommon::addVirtualContact(const ContactId & cId)
{
  return virtualContacts_.insert(cId).second;
}

bool ContactConstrCommon::removeVirtualContact(const ContactId & cId)
{
  return virtualContacts_.erase(cId) == 1;
}

void ContactConstrCommon::resetVirtualContacts()
{
  virtualContacts_.clear();
}

bool ContactConstrCommon::addDofContact(const ContactId & cId, const Eigen::MatrixXd & dof)
{
  return dofContacts_.insert({cId, dof}).second;
}

bool ContactConstrCommon::removeDofContact(const ContactId & cId)
{
  return dofContacts_.erase(cId) == 1;
}

void ContactConstrCommon::resetDofContacts()
{
  dofContacts_.clear();
}

std::set<ContactConstrCommon::ContactCommon> ContactConstrCommon::contactCommonInContact(
    const std::vector<rbd::MultiBody> & mbs,
    const SolverData & data)
{
  std::set<ContactCommon> ret;
  auto isValid = [&mbs, this](const ContactId & contactId) {
    // if is virtualContacts we don't add it
    return (virtualContacts_.find(contactId) == virtualContacts_.end());
  };

  for(const BilateralContact & c : data.allContacts())
  {
    if(isValid(c.contactId))
    {
      ret.insert({c.contactId, c.X_b1_cf, c.X_b1_b2});
    }
  }

  return std::move(ret);
}

/**
 *															ContactConstr
 */

void ContactConstr::ContactData::update(const std::vector<rbd::MultiBodyConfig> & mbcs)
{
  if(contacts.size() != 2 && contacts[0].sign > 0)
  {
    // X_b1_b2 is only encoded in the constraint if there is two active robots or if the active robot is second
    return;
  }
  auto & X_b2_cf = contacts.size() == 2 ? contacts[1].X_b_p : contacts[0].X_b_p;
  const auto & mbc1 = mbcs[r1Index];
  const auto & X_0_b1 = mbc1.bodyPosW[b1Index];
  const auto & mbc2 = mbcs[r2Index];
  const auto & X_0_b2 = mbc2.bodyPosW[b2Index];
  auto X_b1_b2_current = X_0_b2 * X_0_b1.inv();
  auto X_b2_cf_current = X_b1_cf * X_b1_b2_current.inv();
  // Only apply the motion allowed by the DoF selection
  Eigen::Vector6d error = revDof * sva::transformError(X_b2_cf_current.inv(), X_b2_cf.inv()).vector();
  auto offset = sva::PTransformd(sva::RotX(error(0)) * sva::RotY(error(1)) * sva::RotZ(error(2)),
                                 Eigen::Vector3d(error(3), error(4), error(5)));
  X_b2_cf = offset * X_b2_cf;
  X_b1_b2 = X_b2_cf.inv() * X_b1_cf;
}

ContactConstr::ContactConstr() : cont_(), fullJac_(), dofJac_(), A_(), b_(), nrEq_(0), totalAlphaD_(0) {}

void ContactConstr::updateDofContacts()
{
  for(ContactData & c : cont_)
  {
    auto it = dofContacts_.find(c.contactId);
    if(it != dofContacts_.end())
    {
      c.dof = it->second;
    }
    else
    {
      c.dof.setIdentity(6, 6);
    }
  }
  updateNrEq();
}

void ContactConstr::updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data)
{
  cont_.clear();
  totalAlphaD_ = data.totalAlphaD();

  int maxDof = std::max_element(mbs.begin(), mbs.end(), compareDof)->nrDof();
  fullJac_.resize(6, maxDof);
  dofJac_.resize(6, maxDof);

  std::set<ContactCommon> contactCSet = contactCommonInContact(mbs, data);
  for(const ContactCommon & cC : contactCSet)
  {
    Eigen::MatrixXd dof(Eigen::MatrixXd::Identity(6, 6));
    auto it = dofContacts_.find(cC.cId);
    if(it != dofContacts_.end())
    {
      dof = it->second;
    }
    std::vector<ContactSideData> contacts;
    auto addContact = [&mbs, &data, &contacts](int rIndex, const std::string & bName, double sign,
                                               const sva::PTransformd & point) {
      if(mbs[rIndex].nrDof() > 0)
      {
        contacts.emplace_back(rIndex, data.alphaDBegin(rIndex), sign, rbd::Jacobian(mbs[rIndex], bName), point);
        return contacts.back().jac.jointsPath().back();
      }
      return mbs[rIndex].bodyIndexByName(bName);
    };
    int r1Index = cC.cId.r1Index;
    int b1Index = addContact(r1Index, cC.cId.r1BodyName, 1., cC.X_b1_cf);
    int r2Index = cC.cId.r2Index;
    int b2Index = addContact(r2Index, cC.cId.r2BodyName, -1., cC.X_b1_cf * cC.X_b1_b2.inv());

    cont_.emplace_back(std::move(contacts), dof, r1Index, r2Index, b1Index, b2Index, cC.X_b1_b2, cC.X_b1_cf, cC.cId);
  }
  updateNrEq();

  A_.setZero(cont_.size() * 6, data.nrVars());
  b_.setZero(cont_.size() * 6);
}

int ContactConstr::nrEq() const
{
  return nrEq_;
}

std::string ContactConstr::descEq(const std::vector<rbd::MultiBody> & mbs, int line)
{
  std::ostringstream oss;
  int contact = line / 6;
  for(const ContactSideData & csd : cont_[contact].contacts)
  {
    int body = csd.jac.jointsPath().back();
    oss << "Contact: " << mbs[csd.robotIndex].body(body).name() << std::endl;
  }
  return oss.str();
}

int ContactConstr::maxEq() const
{
  return int(A_.rows());
}

const Eigen::MatrixXd & ContactConstr::AEq() const
{
  return A_;
}

const Eigen::VectorXd & ContactConstr::bEq() const
{
  return b_;
}

void ContactConstr::updateNrEq()
{
  nrEq_ = 0;
  for(const ContactData & c : cont_)
  {
    nrEq_ += int(c.dof.rows());
  }
}

/**
 *															ContactAccConstr
 */

ContactAccConstr::ContactAccConstr() : ContactConstr() {}

void ContactAccConstr::update(const std::vector<rbd::MultiBody> & mbs,
                              const std::vector<rbd::MultiBodyConfig> & mbcs,
                              const SolverData & data)
{
  using namespace Eigen;

  A_.setZero();
  b_.setZero();
  // J_i*alphaD + JD_i*alpha = 0

  int index = 0;
  for(std::size_t i = 0; i < cont_.size(); ++i)
  {
    ContactData & cd = cont_[i];
    int rows = int(cd.dof.rows());

    cd.update(mbcs);

    for(std::size_t j = 0; j < cd.contacts.size(); ++j)
    {
      ContactSideData & csd = cd.contacts[j];
      const rbd::MultiBody & mb = mbs[csd.robotIndex];
      const rbd::MultiBodyConfig & mbc = mbcs[csd.robotIndex];

      // AEq = J_i
      sva::PTransformd X_0_p = csd.X_b_p * mbc.bodyPosW[csd.bodyIndex];
      const MatrixXd & jacMat = csd.jac.jacobian(mb, mbc, X_0_p);
      dofJac_.block(0, 0, rows, csd.jac.dof()).noalias() = csd.sign * cd.dof * jacMat;
      csd.jac.fullJacobian(mb, dofJac_.block(0, 0, rows, csd.jac.dof()), fullJac_);
      A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() += fullJac_.block(0, 0, rows, mb.nrDof());

      // BEq = -JD_i*alpha
      Vector6d normalAcc = csd.jac
                               .normalAcceleration(mb, mbc, data.normalAccB(csd.robotIndex), csd.X_b_p,
                                                   sva::MotionVecd(Vector6d::Zero()))
                               .vector();
      b_.segment(index, rows).noalias() -= csd.sign * cd.dof * normalAcc;
    }
    index += rows;
  }
}

std::string ContactAccConstr::nameEq() const
{
  return "ContactAccConstr";
}

/**
 *															ContactSpeedConstr
 */

ContactSpeedConstr::ContactSpeedConstr(double timeStep) : ContactConstr(), timeStep_(timeStep) {}

void ContactSpeedConstr::update(const std::vector<rbd::MultiBody> & mbs,
                                const std::vector<rbd::MultiBodyConfig> & mbcs,
                                const SolverData & data)
{
  using namespace Eigen;

  A_.block(0, 0, nrEq_, totalAlphaD_).setZero();
  b_.head(nrEq_).setZero();
  // J_i*alphaD + JD_i*alpha = 0

  int index = 0;
  for(std::size_t i = 0; i < cont_.size(); ++i)
  {
    ContactData & cd = cont_[i];
    int rows = int(cd.dof.rows());

    cd.update(mbcs);

    for(std::size_t j = 0; j < cd.contacts.size(); ++j)
    {
      ContactSideData & csd = cd.contacts[j];
      const rbd::MultiBody & mb = mbs[csd.robotIndex];
      const rbd::MultiBodyConfig & mbc = mbcs[csd.robotIndex];

      // AEq = J_i
      sva::PTransformd X_0_p = csd.X_b_p * mbc.bodyPosW[csd.bodyIndex];
      const MatrixXd & jacMat = csd.jac.jacobian(mb, mbc, X_0_p);
      dofJac_.block(0, 0, rows, csd.jac.dof()).noalias() = csd.sign * cd.dof * jacMat;
      csd.jac.fullJacobian(mb, dofJac_.block(0, 0, rows, csd.jac.dof()), fullJac_);
      A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() += fullJac_.block(0, 0, rows, mb.nrDof());

      // BEq = -JD_i*alpha
      Vector6d normalAcc = csd.jac
                               .normalAcceleration(mb, mbc, data.normalAccB(csd.robotIndex), csd.X_b_p,
                                                   sva::MotionVecd(Vector6d::Zero()))
                               .vector();
      Vector6d velocity = csd.jac.velocity(mb, mbc, csd.X_b_p).vector();
      b_.segment(index, rows).noalias() -= csd.sign * cd.dof * (normalAcc + velocity / timeStep_);
    }

    index += rows;
  }
}

std::string ContactSpeedConstr::nameEq() const
{
  return "ContactSpeedConstr";
}

/**
 *															ContactPosConstr
 */

ContactPosConstr::ContactPosConstr(double timeStep) : ContactConstr(), timeStep_(timeStep) {}

void ContactPosConstr::update(const std::vector<rbd::MultiBody> & mbs,
                              const std::vector<rbd::MultiBodyConfig> & mbcs,
                              const SolverData & data)
{
  using namespace Eigen;

  A_.block(0, 0, nrEq_, totalAlphaD_).setZero();
  b_.head(nrEq_).setZero();
  // J_i*alphaD + JD_i*alpha = 0

  int index = 0;
  for(std::size_t i = 0; i < cont_.size(); ++i)
  {
    ContactData & cd = cont_[i];
    int rows = int(cd.dof.rows());

    cd.update(mbcs);

    for(std::size_t j = 0; j < cd.contacts.size(); ++j)
    {
      ContactSideData & csd = cd.contacts[j];
      const rbd::MultiBody & mb = mbs[csd.robotIndex];
      const rbd::MultiBodyConfig & mbc = mbcs[csd.robotIndex];

      // AEq = J_i
      sva::PTransformd X_0_p = csd.X_b_p * mbc.bodyPosW[csd.bodyIndex];
      const MatrixXd & jacMat = csd.jac.jacobian(mb, mbc, X_0_p);
      dofJac_.block(0, 0, rows, csd.jac.dof()).noalias() = csd.sign * cd.dof * jacMat;
      csd.jac.fullJacobian(mb, dofJac_.block(0, 0, rows, csd.jac.dof()), fullJac_);
      A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() += fullJac_.block(0, 0, rows, mb.nrDof());

      // BEq = -JD_i*alpha
      Vector6d normalAcc = csd.jac
                               .normalAcceleration(mb, mbc, data.normalAccB(csd.robotIndex), csd.X_b_p,
                                                   sva::MotionVecd(Vector6d::Zero()))
                               .vector();
      Vector6d velocity = csd.jac.velocity(mb, mbc, csd.X_b_p).vector();
      b_.segment(index, rows).noalias() -= csd.sign * cd.dof * (normalAcc + velocity / timeStep_);
    }

    // target the derivative of the position error
    sva::PTransformd X_0_b1cf = cd.X_b1_cf * mbcs[cd.contactId.r1Index].bodyPosW[cd.b1Index];
    sva::PTransformd X_0_b2cf = cd.X_b1_cf * cd.X_b1_b2.inv() * mbcs[cd.contactId.r2Index].bodyPosW[cd.b2Index];

    sva::PTransformd X_b1cf_b2cf = X_0_b2cf * X_0_b1cf.inv();
    Eigen::Vector6d error;
    error.head<3>() = sva::rotationVelocity(X_b1cf_b2cf.rotation());
    error.tail<3>() = X_b1cf_b2cf.translation();
    b_.segment(index, rows) += cd.dof * (error / timeStep_);

    index += rows;
  }
}

std::string ContactPosConstr::nameEq() const
{
  return "ContactPosConstr";
}

/**
 *															TorqueFbTermContactPDConstr
 */


TorqueFbTermContactPDConstr::TorqueFbTermContactPDConstr(Eigen::Vector6d stiffness,
                                                         Eigen::Vector6d damping,
                                                         int mainRobotIndex,
                                                         const std::shared_ptr<torque_control::TorqueFeedbackTerm> fbTerm):
        ContactConstr(),
        stiffness_default_(stiffness),
        damping_default_(damping),
        mainRobotIndex_(mainRobotIndex),
        fbTerm_(fbTerm)
{}


void TorqueFbTermContactPDConstr::setPDgainsForContact(const ContactId& cId, const Eigen::Vector6d& stiff,
                                                       const Eigen::Vector6d& damp)
{
        if (contPD_.find(cId) == contPD_.end())
                contPD_.insert(std::make_pair(cId, PDgains(stiff, damp)));
        else
                contPD_.at(cId) = PDgains(stiff, damp);
}


void TorqueFbTermContactPDConstr::update(const std::vector<rbd::MultiBody>& mbs,
                                         const std::vector<rbd::MultiBodyConfig>& mbcs,
                                         const SolverData& data)
{
        using namespace Eigen;
  
        A_.block(0, 0, nrEq_, totalAlphaD_).setZero();
	b_.head(nrEq_).setZero();

        int index = 0;
        
        for(std::size_t i = 0; i < cont_.size(); ++i)
        {
                ContactData& cd = cont_[i];
                int rows = int(cd.dof.rows());

                Vector6d stiffness, damping;

                if (contPD_.find(cd.contactId) != contPD_.end())
                {
                        stiffness = contPD_.at(cd.contactId).stiffness;
                        damping = contPD_.at(cd.contactId).damping;
                }
                else
                {
                        stiffness = stiffness_default_;
                        damping = damping_default_;
                }

                for(std::size_t j = 0; j < cd.contacts.size(); ++j)
                {
                        ContactSideData& csd = cd.contacts[j];
                        const rbd::MultiBody& mb = mbs[csd.robotIndex];
                        const rbd::MultiBodyConfig& mbc = mbcs[csd.robotIndex];

                        // AEq = J_i
                        sva::PTransformd X_0_p = csd.X_b_p*mbc.bodyPosW[csd.bodyIndex];
                        const MatrixXd& jacMat = csd.jac.jacobian(mb, mbc, X_0_p);
                        dofJac_.block(0, 0, rows, csd.jac.dof()).noalias() =
                                csd.sign*cd.dof*jacMat;
                        csd.jac.fullJacobian(mb, dofJac_.block(0, 0, rows, csd.jac.dof()),
                                fullJac_);
                        A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() +=
                                fullJac_.block(0, 0, rows, mb.nrDof());

                        // BEq = dVSurf_obj - JD_i*alpha - J_i*gammaD
                        Vector6d normalAcc = csd.jac.normalAcceleration(
                                mb, mbc, data.normalAccB(csd.robotIndex), csd.X_b_p,
                                sva::MotionVecd(Vector6d::Zero())).vector();
                        Vector6d velocity = csd.jac.velocity(mb, mbc, csd.X_b_p).vector();
                        b_.segment(index, rows).noalias() -=
                          csd.sign*cd.dof*(normalAcc + damping.asDiagonal() * velocity);
                        if (csd.robotIndex == mainRobotIndex_)
                        {
                                b_.segment(index, rows).noalias() -=
                                        csd.sign * cd.dof * fullJac_.block(0, 0, rows, mb.nrDof()) * fbTerm_->gammaD();
                        }
                }

                sva::PTransformd X_0_b1cf =
                        cd.X_b1_cf*mbcs[cd.contactId.r1Index].bodyPosW[cd.b1Index];
                sva::PTransformd X_0_b2cf =
                        cd.X_b1_cf*cd.X_b1_b2.inv()*mbcs[cd.contactId.r2Index].bodyPosW[cd.b2Index];

                sva::PTransformd X_b1cf_b2cf = X_0_b2cf*X_0_b1cf.inv();
                Eigen::Vector6d error;
                error.head<3>() = sva::rotationVelocity(X_b1cf_b2cf.rotation());
                error.tail<3>() = X_b1cf_b2cf.translation();
                b_.segment(index, rows) += cd.dof * stiffness.asDiagonal() * error;
                
                index += rows;
        }
}


std::string TorqueFbTermContactPDConstr::nameEq() const
{
	return "TorqueFbTermContactPDConstr";
}


/**
	*															TorqueFbTermContactHybridConstr
	*/


TorqueFbTermContactHybridConstr::TorqueFbTermContactHybridConstr(Eigen::Vector6d proportional,
                                                                 Eigen::Vector6d derivative,
                                                                 Eigen::Vector6d damping,
                                                                 int mainRobotIndex,
                                                                 const std::shared_ptr<torque_control::TorqueFeedbackTerm> fbTerm):
        ContactConstr(),
        proportional_default_(proportional),
        derivative_default_(derivative),
        damping_default_(damping),
        mainRobotIndex_(mainRobotIndex),
        fbTerm_(fbTerm)
{}
  

void TorqueFbTermContactHybridConstr::setGainsForContact(const ContactId& cId, const Eigen::Vector6d& prop,
                                                         const Eigen::Vector6d& deriv, const Eigen::Vector6d& damp)
{
        if (contGains_.find(cId) == contGains_.end())
                contGains_.insert(std::make_pair(cId, Gains(prop, deriv, damp)));
        else
                contGains_.at(cId) = Gains(prop, deriv, damp);
}


void TorqueFbTermContactHybridConstr::update(const std::vector<rbd::MultiBody>& mbs,
                                             const std::vector<rbd::MultiBodyConfig>& mbcs,
                                             const SolverData& data)
{
        using namespace Eigen;
  
        A_.block(0, 0, nrEq_, totalAlphaD_).setZero();
	b_.head(nrEq_).setZero();

        int index = 0;
        
        for(std::size_t i = 0; i < cont_.size(); ++i)
        {
                ContactData& cd = cont_[i];
                int rows = int(cd.dof.rows());

                Vector6d proportional, derivative, damping;

                if (contGains_.find(cd.contactId) != contGains_.end())
                {
                        proportional = contGains_.at(cd.contactId).proportional;
                        derivative = contGains_.at(cd.contactId).derivative;
                        damping = contGains_.at(cd.contactId).damping;
                }
                else
                {
                        proportional = proportional_default_;
                        derivative = derivative_default_;
                        damping = damping_default_;
                }

                for(std::size_t j = 0; j < cd.contacts.size(); ++j)
                {
                        ContactSideData& csd = cd.contacts[j];
                        const rbd::MultiBody& mb = mbs[csd.robotIndex];
                        const rbd::MultiBodyConfig& mbc = mbcs[csd.robotIndex];

                        // AEq = J_i
                        sva::PTransformd X_0_p = csd.X_b_p*mbc.bodyPosW[csd.bodyIndex];
                        const MatrixXd& jacMat = csd.jac.jacobian(mb, mbc, X_0_p);
                        dofJac_.block(0, 0, rows, csd.jac.dof()).noalias() =
                                csd.sign*cd.dof*jacMat;
                        csd.jac.fullJacobian(mb, dofJac_.block(0, 0, rows, csd.jac.dof()),
                                fullJac_);
                        A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() +=
                                fullJac_.block(0, 0, rows, mb.nrDof());

                        // BEq = dVSurf_obj - JD_i*alpha - J_i*gammaD
                        Vector6d normalAcc = csd.jac.normalAcceleration(
                                mb, mbc, data.normalAccB(csd.robotIndex), csd.X_b_p,
                                sva::MotionVecd(Vector6d::Zero())).vector();
                        Vector6d velocity = csd.jac.velocity(mb, mbc, csd.X_b_p).vector();
                        b_.segment(index, rows).noalias() -=
                          csd.sign*cd.dof*(normalAcc + damping.asDiagonal() * velocity);
                        if (csd.robotIndex == mainRobotIndex_)
                        {
                                b_.segment(index, rows).noalias() -=
                                        csd.sign * cd.dof * fullJac_.block(0, 0, rows, mb.nrDof()) * fbTerm_->gammaD();
                        }
                }

                // Pending to add the admittance component

                index += rows;
        }
}


std::string TorqueFbTermContactHybridConstr::nameEq() const
{
	return "TorqueFbTermContactHybridConstr";
}


} // namespace qp

} // namespace tasks
