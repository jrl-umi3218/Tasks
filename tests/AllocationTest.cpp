// Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
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


// boost
#define BOOST_TEST_MODULE AllocationTest
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

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
  if(c_pos == std::string::npos)
  {
    c_pos = fn_name.find_first_of(',');
  }
  return fn_name.substr(eq_pos, c_pos - eq_pos);
}

template<typename T, typename ... Args>
void test_raw_ptr_creation(Args ... args)
{
#ifdef __GNUG__
  std::string T_name = class_name(__PRETTY_FUNCTION__);
  std::cout << "Testing raw pointer creation for " << T_name << std::endl;
#endif
  T * p = nullptr;
  for(size_t i = 0; i < 100; ++i)
  {
    new T(args...);
  }
}

template<typename T, typename ... Args>
void test_shared_ptr_creation(Args ... args)
{
#ifdef __GNUG__
  std::string T_name = class_name(__PRETTY_FUNCTION__);
  std::cout << "Testing shared_ptr creation for for " << T_name << std::endl;
#endif
  std::shared_ptr<T> p = nullptr;
  for(size_t i = 0; i < 100; ++i)
  {
    p = std::make_shared<T>(args...);
  }
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
