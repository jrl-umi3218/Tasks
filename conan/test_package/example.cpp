/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// includes
// std
#include <fstream>
#include <tuple>

// Eigen
#include <unsupported/Eigen/Polynomials>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include <RBDyn/EulerIntegration.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/ID.h>
#include <RBDyn/Body.h>
#include <RBDyn/Joint.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>


// sch
#include <sch/CD/CD_Pair.h>
#include <sch/S_Object/S_Box.h>
#include <sch/S_Object/S_Sphere.h>

// Tasks
#include <Tasks/Bounds.h>
#include <Tasks/QPConstr.h>
#include <Tasks/QPContactConstr.h>
#include <Tasks/QPMotionConstr.h>
#include <Tasks/QPSolver.h>
#include <Tasks/QPTasks.h>

// arms.h
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeZXZArm(bool isFixed = true,
                                                            const sva::PTransformd X_base = sva::PTransformd::Identity())
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;

  MultiBodyGraph mbg;

  double mass = 1.;
  Matrix3d I = Matrix3d::Identity();
  Vector3d h = Vector3d::Zero();

  RBInertiad rbi(mass, h, I);

  Body b0(rbi, "b0");
  Body b1(rbi, "b1");
  Body b2(rbi, "b2");
  Body b3(rbi, "b3");

  mbg.addBody(b0);
  mbg.addBody(b1);
  mbg.addBody(b2);
  mbg.addBody(b3);

  Joint j0(Joint::RevZ, true, "j0");
  Joint j1(Joint::RevX, true, "j1");
  Joint j2(Joint::RevZ, true, "j2");

  mbg.addJoint(j0);
  mbg.addJoint(j1);
  mbg.addJoint(j2);

  //  Root     j0       j1     j2
  //  ---- b0 ---- b1 ---- b2 ----b3
  //  Fixed    Z       X       Z

  PTransformd to(Vector3d(0., 0.5, 0.));
  PTransformd from(Vector3d(0., 0., 0.));

  mbg.linkBodies("b0", PTransformd::Identity(), "b1", from, "j0");
  mbg.linkBodies("b1", to, "b2", from, "j1");
  mbg.linkBodies("b2", to, "b3", from, "j2");

  MultiBody mb = mbg.makeMultiBody("b0", isFixed, X_base);

  MultiBodyConfig mbc(mb);
  mbc.zero(mb);

  return std::make_tuple(mb, mbc);
}

int main()
{
  using namespace Eigen;
  using namespace sva;
  using namespace rbd;
  using namespace tasks;

  MultiBody mb;
  MultiBodyConfig mbcInit;

  std::tie(mb, mbcInit) = makeZXZArm();

  std::vector<MultiBody> mbs = {mb};
  std::vector<MultiBodyConfig> mbcs(1);

  forwardKinematics(mb, mbcInit);
  forwardVelocity(mb, mbcInit);

  qp::QPSolver solver;

  solver.nrVars(mbs, {}, {});

  solver.updateConstrSize();

  Vector3d posD = Vector3d(0.707106, 0.707106, 0.);
  qp::PositionTask posTask(mbs, 0, "b3", posD);
  qp::SetPointTask posTaskSp(mbs, 0, &posTask, 10., 1.);

  // Test addTask
  solver.addTask(&posTaskSp);

  // Test PositionTask
  mbcs[0] = mbcInit;
  for(int i = 0; i < 10000; ++i)
  {
    eulerIntegration(mb, mbcs[0], 0.001);

    forwardKinematics(mb, mbcs[0]);
    forwardVelocity(mb, mbcs[0]);
  }

  std::cout << "(C++) Final norm of position task: " << posTask.eval().norm() << "\n";
  return 0;
}
