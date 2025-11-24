#include "Tasks/QPCoincidenceConstr.h"
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

#include <mc_rtc/logging.h>
#include <RBDyn/Jacobian.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include "Tasks/QPTasks.h"
#include "utils.h"
#include <Eigen/QR>
#include <fstream>
#include <iostream>
#include <set>

namespace tasks
{
namespace qp
{

const Eigen::MatrixXd & CoincidenceConstr::AEq() const
{
  return A_;
}

const Eigen::VectorXd & CoincidenceConstr::bEq() const
{
  return b_;
}

int CoincidenceConstr::maxEq() const
{
  return int(A_.rows());
}

void CoincidenceConstr::setJointSelector(const Eigen::VectorXd & selector)
{
  jointSelector_ = selector;
}

// =================== FixedCoincidenceConstr ===================

FixedCoincidenceConstr::FixedCoincidenceConstr(int robotIndex,
                                               const std::string & body1Name,
                                               const std::string & body2Name,
                                               const Eigen::Vector3d & point1,
                                               const Eigen::Vector3d & point2,
                                               const Eigen::VectorXd & jointSelector)
: CoincidenceConstr(robotIndex, body1Name, body2Name, jointSelector), point1_(point1), point2_(point2)
{
  A_.resize(0, 0);
  b_.resize(0);
  std::cout << "FixedCoincidenceConstr created: " << body1Name << " <-> " << body2Name << " (robot " << robotIndex
            << ")" << std::endl;
}

std::string FixedCoincidenceConstr::nameEq() const
{
  return "FixedCoincidenceConstr";
}

std::string FixedCoincidenceConstr::descEq(const std::vector<rbd::MultiBody> & mbs, int i)
{
  return "FixedCoincidenceConstr between " + body1Name_ + " and " + body2Name_ + ".";
}

void FixedCoincidenceConstr::updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data)
{
  const rbd::MultiBody & mb = mbs[robotIndex_];
  body1Index_ = mb.bodyIndexByName(body1Name_);
  body2Index_ = mb.bodyIndexByName(body2Name_);
  if(body1Index_ == -1 || body2Index_ == -1)
  {
    std::cerr << "Error: Invalid body indices in FixedCoincidenceConstr::update" << std::endl;
    return;
  }

  A_.setZero(6, data.nrVars());
  b_.setZero(6);
}

void FixedCoincidenceConstr::update(const std::vector<rbd::MultiBody> & mbs,
                                    const std::vector<rbd::MultiBodyConfig> & mbcs,
                                    const SolverData & data)
{
  using namespace Eigen;
  const rbd::MultiBody & mb = mbs[robotIndex_];
  const rbd::MultiBodyConfig & mbc = mbcs[robotIndex_];
  if(body1Index_ == -1 || body2Index_ == -1)
  {
    std::cerr << "Error: Invalid body indices in FixedCoincidenceConstr::update" << std::endl;
    return;
  }

  int alphaDBegin = data.alphaDBegin(robotIndex_);
  int nrDof = mb.nrDof();

  rbd::Jacobian jac1(mb, body1Name_);
  Eigen::MatrixXd jacMat1 = jac1.jacobian(mb, mbc);
  Eigen::MatrixXd fulljacobian1(6, nrDof);
  jac1.fullJacobian(mb, jacMat1.block(0, 0, 6, nrDof), fulljacobian1);

  rbd::Jacobian jac2(mb, body2Name_);
  Eigen::MatrixXd jacMat2 = jac2.jacobian(mb, mbc);
  Eigen::MatrixXd fulljacobian2(6, nrDof);
  jac2.fullJacobian(mb, jacMat2.block(0, 0, 6, nrDof), fulljacobian2);

  A_.block(0, alphaDBegin, 6, nrDof) = (fulljacobian1 - fulljacobian2) * jointSelector_.asDiagonal();

  const std::vector<sva::MotionVecd> & normalAccB = data.normalAccB(robotIndex_);
  Vector6d normalAcc1 = jac1.normalAcceleration(mb, mbc, normalAccB).vector();
  Vector6d normalAcc2 = jac2.normalAcceleration(mb, mbc, normalAccB).vector();

  double Kp = 1;
  double Kd = 2;
  Eigen::VectorXd errorp = mbc.bodyPosW[body1Index_].translation() - mbc.bodyPosW[body2Index_].translation();
  Eigen::VectorXd errord = mbc.bodyVelW[body1Index_].linear() - mbc.bodyVelW[body2Index_].linear();

  b_.segment(0, 6) = -(normalAcc1.tail(6) - normalAcc2.tail(6)) - Kp * errorp - Kd * errord;
}

// =============== RotationalCoincidenceConstr ===============

RotationalCoincidenceConstr::RotationalCoincidenceConstr(int robotIndex,
                                                         const std::string & body1Name,
                                                         const std::string & body2Name,
                                                         const Eigen::VectorXd & jointSelector)
: CoincidenceConstr(robotIndex, body1Name, body2Name, jointSelector)
{
  A_.resize(0, 0);
  b_.resize(0);
  std::cout << "RotationalCoincidenceConstr created: " << body1Name << " <-> " << body2Name << " (robot " << robotIndex
            << ")" << std::endl;
}

std::string RotationalCoincidenceConstr::nameEq() const
{
  return "RotationalCoincidenceConstr";
}

std::string RotationalCoincidenceConstr::descEq(const std::vector<rbd::MultiBody> & mbs, int i)
{
  return "RotationalCoincidenceConstr between " + body1Name_ + " and " + body2Name_ + ".";
}

void RotationalCoincidenceConstr::updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data)
{
  const rbd::MultiBody & mb = mbs[robotIndex_];
  body1Index_ = mb.bodyIndexByName(body1Name_);
  body2Index_ = mb.bodyIndexByName(body2Name_);
  if(body1Index_ == -1 || body2Index_ == -1)
  {
    std::cerr << "Error: Invalid body indices in RotationalCoincidenceConstr::update" << std::endl;
    return;
  }

  A_.setZero(3, data.nrVars());
  b_.setZero(3);
}

void RotationalCoincidenceConstr::update(const std::vector<rbd::MultiBody> & mbs,
                                         const std::vector<rbd::MultiBodyConfig> & mbcs,
                                         const SolverData & data)
{
  using namespace Eigen;
  const rbd::MultiBody & mb = mbs[robotIndex_];
  const rbd::MultiBodyConfig & mbc = mbcs[robotIndex_];
  if(body1Index_ == -1 || body2Index_ == -1)
  {
    std::cerr << "Error: Invalid body indices in RotationalCoincidenceConstr::update" << std::endl;
    return;
  }

  int alphaDBegin = data.alphaDBegin(robotIndex_);
  int nrDof = mb.nrDof();

  rbd::Jacobian jac1(mb, body1Name_);
  Eigen::MatrixXd jacMat1 = jac1.jacobian(mb, mbc);
  Eigen::MatrixXd fulljacobian1(3, nrDof);
  jac1.fullJacobian(mb, jacMat1.block(3, 0, 3, nrDof), fulljacobian1);

  rbd::Jacobian jac2(mb, body2Name_);
  Eigen::MatrixXd jacMat2 = jac2.jacobian(mb, mbc);
  Eigen::MatrixXd fulljacobian2(3, nrDof);
  jac2.fullJacobian(mb, jacMat2.block(3, 0, 3, nrDof), fulljacobian2);

  A_.block(0, alphaDBegin, 3, nrDof) = (fulljacobian1 - fulljacobian2) * jointSelector_.asDiagonal();

  const std::vector<sva::MotionVecd> & normalAccB = data.normalAccB(robotIndex_);
  Vector6d normalAcc1 = jac1.normalAcceleration(mb, mbc, normalAccB).vector();
  Vector6d normalAcc2 = jac2.normalAcceleration(mb, mbc, normalAccB).vector();

  double Kp = 1;
  double Kd = 2;
  Eigen::VectorXd errorp = mbc.bodyPosW[body1Index_].translation() - mbc.bodyPosW[body2Index_].translation();
  Eigen::VectorXd errord = mbc.bodyVelW[body1Index_].linear() - mbc.bodyVelW[body2Index_].linear();

  b_.segment(0, 3) = -(normalAcc1.tail(3) - normalAcc2.tail(3)) - Kp * errorp - Kd * errord;

  double distance = errorp.norm();
  double distancex = (mbc.bodyPosW[body1Index_].translation()).norm();
  std::ofstream file("close_loop_measurement.csv", std::ios::app);
  file << distance << "," << distancex << "\n";
  file.close();
}

} // namespace qp
} // namespace tasks
