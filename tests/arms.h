/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
// std
#include <tuple>

// RBDyn
#include <RBDyn/Body.h>
#include <RBDyn/Joint.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>

/// @return An simple ZXZ arm with Y as up axis.
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

/// @return A one body robot for the environnment.
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeEnv()
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

  mbg.addBody(b0);

  MultiBody mb = mbg.makeMultiBody("b0", true);

  MultiBodyConfig mbc(mb);
  mbc.zero(mb);

  return std::make_tuple(mb, mbc);
}
