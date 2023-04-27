/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// boost
#define BOOST_TEST_MODULE AllocationTest
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>

// Tasks
#include "Tasks/QPConstr.h"
#include "Tasks/QPTasks.h"

// Arms
#include "arms.h"

#pragma GCC diagnostic ignored "-Wunused-variable"

std::string class_name(const std::string & fn_name)
{
  auto eq_pos = fn_name.find_first_of('=') + 2;
  auto c_pos = fn_name.find_first_of(';');
  if(c_pos == std::string::npos) { c_pos = fn_name.find_first_of(','); }
  return fn_name.substr(eq_pos, c_pos - eq_pos);
}

template<typename T, typename... Args>
void test_raw_ptr_creation(Args &&... args)
{
#ifdef __GNUG__
  std::string T_name = class_name(__PRETTY_FUNCTION__);
  std::cout << "Testing raw pointer creation for " << T_name << std::endl;
#endif
  T * p = nullptr;
  for(size_t i = 0; i < 100; ++i) { new T(std::forward<Args>(args)...); }
}

template<typename T, typename... Args>
void test_shared_ptr_creation(Args &&... args)
{
#ifdef __GNUG__
  std::string T_name = class_name(__PRETTY_FUNCTION__);
  std::cout << "Testing shared_ptr creation for for " << T_name << std::endl;
#endif
  std::shared_ptr<T> p = nullptr;
  for(size_t i = 0; i < 100; ++i) { p = std::make_shared<T>(std::forward<Args>(args)...); }
}

BOOST_AUTO_TEST_CASE(AllocationTest)
{
  rbd::MultiBody mb1, mb2;
  rbd::MultiBodyConfig mbc1Init, mbc2Init;

  std::tie(mb1, mbc1Init) = makeZXZArm();
  std::tie(mb2, mbc2Init) = makeZXZArm();

  std::vector<rbd::MultiBody> mbs = {mb1, mb2};
  Eigen::Vector2d pt2d = Eigen::Vector2d::Zero();
  sva::PTransformd pt = sva::PTransformd::Identity();
  test_raw_ptr_creation<tasks::qp::ImageConstr>(mbs, 0, "b3", pt, 1.0);
  test_raw_ptr_creation<tasks::GazeTask>(mbs[0], "b3", pt2d, 1.0, pt);
  test_raw_ptr_creation<tasks::PositionBasedVisServoTask>(mbs[0], "b3", pt, pt);
  test_raw_ptr_creation<tasks::qp::GazeTask>(mbs, 0, "b3", pt2d, 1.0, pt);
  test_raw_ptr_creation<tasks::qp::PositionBasedVisServoTask>(mbs, 0, "b3", pt, pt);
  test_raw_ptr_creation<tasks::qp::MultiRobotTransformTask>(mbs, 0, 1, "b3", "b3", pt, pt, 1.0, 1.0);
  test_shared_ptr_creation<tasks::qp::ImageConstr>(mbs, 0, "b3", pt, 1.0);
  test_shared_ptr_creation<tasks::GazeTask>(mbs[0], "b3", pt2d, 1.0, pt);
  test_shared_ptr_creation<tasks::PositionBasedVisServoTask>(mbs[0], "b3", pt, pt);
  test_shared_ptr_creation<tasks::qp::GazeTask>(mbs, 0, "b3", pt2d, 1.0, pt);
  test_shared_ptr_creation<tasks::qp::PositionBasedVisServoTask>(mbs, 0, "b3", pt, pt);
  test_shared_ptr_creation<tasks::qp::MultiRobotTransformTask>(mbs, 0, 1, "b3", "b3", pt, pt, 1.0, 1.0);
}
