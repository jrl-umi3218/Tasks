/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// Tasks
#include "QPContacts.h"

// forward declaration
namespace rbd
{
class MultiBody;
struct MultiBodyConfig;
} // namespace rbd

namespace tasks
{

namespace qp
{

class TASKS_DLLAPI SolverData
{
public:
  friend class QPSolver;
  friend class PassivityPIDTerm_QPSolver;

  SolverData();

  int nrVars() const
  {
    return nrVars_;
  }

  int totalAlphaD() const
  {
    return totalAlphaD_;
  }

  int totalLambda() const
  {
    return totalLambda_;
  }

  int alphaD(int robotIndex) const
  {
    return alphaD_[robotIndex];
  }

  int lambda(int contactIndex) const
  {
    return lambda_[contactIndex];
  }

  int alphaDBegin() const
  {
    return 0;
  }

  int alphaDBegin(int robotIndex) const
  {
    return alphaDBegin_[robotIndex];
  }

  int lambdaBegin() const
  {
    return totalAlphaD_;
  }

  int lambdaBegin(int contactIndex) const
  {
    return lambdaBegin_[contactIndex];
  }

  int nrUniLambda() const
  {
    return nrUniLambda_;
  }

  int nrBiLambda() const
  {
    return nrBiLambda_;
  }

  int unilateralBegin() const
  {
    return lambdaBegin();
  }

  int bilateralBegin() const
  {
    return unilateralBegin() + nrUniLambda_;
  }

  int nrContacts() const
  {
    return static_cast<int>(uniCont_.size() + biCont_.size());
  }

  const std::vector<UnilateralContact> & unilateralContacts() const
  {
    return uniCont_;
  }

  const std::vector<BilateralContact> & bilateralContacts() const
  {
    return biCont_;
  }

  const std::vector<BilateralContact> & allContacts() const
  {
    return allCont_;
  }

  void computeNormalAccB(const std::vector<rbd::MultiBody> & mbs, const std::vector<rbd::MultiBodyConfig> & mbcs);

  const std::vector<std::vector<sva::MotionVecd>> & normalAccB() const
  {
    return normalAccB_;
  }

  const std::vector<sva::MotionVecd> & normalAccB(int robotIndex) const
  {
    return normalAccB_[robotIndex];
  }

  const Eigen::VectorXd& lambdaVecPrev() const
  {
    return lambdaVecPrev_;
  }

// private:
protected:
  std::vector<int> alphaD_; //< each robot alphaD vector size
  std::vector<int> alphaDBegin_; //< each robot alphaD vector begin in x
  std::vector<int> lambda_; //< each contact lambda
  std::vector<int> lambdaBegin_; //< each contact lambda vector begin in x
  int totalAlphaD_, totalLambda_;
  int nrUniLambda_, nrBiLambda_;
  int nrVars_; //< total number of var

  std::vector<UnilateralContact> uniCont_;
  std::vector<BilateralContact> biCont_;
  std::vector<BilateralContact> allCont_;

  std::vector<int> mobileRobotIndex_; //< robot index with dof > 0
  /// normal acceleration of each body of each robot
  std::vector<std::vector<sva::MotionVecd>> normalAccB_;

  Eigen::VectorXd passiveTorque_;
  Eigen::VectorXd lambdaVecPrev_;
};

} // namespace qp

} // namespace tasks
