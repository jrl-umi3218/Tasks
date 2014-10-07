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

// includes
// std
#include <fstream>
#include <iostream>
#include <tuple>

// boost
#define BOOST_TEST_MODULE QPMultiRobotTest
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// RBDyn
#include <RBDyn/EulerIntegration.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/ID.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>

// Tasks
#include "Bounds.h"
#include "QPConstr.h"
#include "QPMotionConstr.h"
#include "QPSolver.h"
#include "QPTasks.h"

// Arms
#include "arms.h"


// Test contact between two robot.
// We set two identical robot at the same positio
// then we link the end effector and add a task
// to make it move on the second robot.
// The first robot must have the same motion.
BOOST_AUTO_TEST_CASE(TwoArmContactTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb1, mb2;
	MultiBodyConfig mbc1Init, mbc2Init;

	std::tie(mb1, mbc1Init) = makeZXZArm();
	std::tie(mb2, mbc2Init) = makeZXZArm();

	forwardKinematics(mb1, mbc1Init);
	forwardVelocity(mb1, mbc1Init);
	forwardKinematics(mb2, mbc2Init);
	forwardVelocity(mb2, mbc2Init);

	sva::PTransformd X_0_b1(mbc1Init.bodyPosW.back());
	sva::PTransformd X_0_b2(mbc2Init.bodyPosW.back());
	sva::PTransformd X_b1_b2(X_0_b2*X_0_b1.inv());

	std::vector<MultiBody> mbs = {mb1, mb2};
	std::vector<MultiBodyConfig> mbcs = {mbc1Init, mbc2Init};

	// Test ContactAccConstr constraint
	// Also test PositionTask on the second robot

	qp::QPSolver solver;

	std::vector<qp::UnilateralContact> contVec =
		{qp::UnilateralContact(0, 1, 3, 3,
			{Vector3d::Zero()}, RotX(cst::pi<double>()/2.), X_b1_b2,
			3, std::tan(cst::pi<double>()/4.))};

	Matrix3d oriD = RotZ(cst::pi<double>()/4.);
	Vector3d posD(oriD*mbc2Init.bodyPosW.back().translation());
	qp::PositionTask posTask(mbs, 1, 3, posD);
	qp::SetPointTask posTaskSp(mbs, 1, &posTask, 1000., 1.);

	qp::ContactAccConstr contCstrAcc;

	contCstrAcc.addToSolver(solver);
	solver.addTask(&posTaskSp);

	solver.nrVars(mbs, contVec, {});
	solver.updateConstrSize();

	// 3 dof + 3 dof + 3 lambda
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 3 + 3);

	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		for(std::size_t r = 0; r < mbs.size(); ++r)
		{
			eulerIntegration(mbs[r], mbcs[r], 0.001);

			forwardKinematics(mbs[r], mbcs[r]);
			forwardVelocity(mbs[r], mbcs[r]);
		}
		// check that the link hold
		sva::PTransformd X_0_b1_post(mbcs[0].bodyPosW.back());
		sva::PTransformd X_0_b2_post(mbcs[1].bodyPosW.back());
		sva::PTransformd X_b1_b2_post(X_0_b2*X_0_b1.inv());
		BOOST_CHECK_SMALL((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 1e-5);
	}
	// check that the task is well minimized
	BOOST_CHECK_SMALL(posTask.eval().norm(), 1e-5);

	contCstrAcc.removeFromSolver(solver);
	solver.removeTask(&posTaskSp);

	// Test ContactSpeedConstr constraint
	// Also test OrientationTask on the second robot

	mbcs = {mbc1Init, mbc2Init};
	qp::OrientationTask oriTask(mbs, 1, 3, oriD);
	qp::SetPointTask oriTaskSp(mbs, 1, &oriTask, 1000., 1.);

	qp::ContactSpeedConstr contCstrSpeed(0.001);

	contCstrSpeed.addToSolver(solver);
	solver.addTask(&oriTaskSp);

	solver.nrVars(mbs, contVec, {});
	solver.updateConstrSize();

	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		for(std::size_t r = 0; r < mbs.size(); ++r)
		{
			eulerIntegration(mbs[r], mbcs[r], 0.001);

			forwardKinematics(mbs[r], mbcs[r]);
			forwardVelocity(mbs[r], mbcs[r]);
		}
		// check that the link hold
		sva::PTransformd X_0_b1_post(mbcs[0].bodyPosW.back());
		sva::PTransformd X_0_b2_post(mbcs[1].bodyPosW.back());
		sva::PTransformd X_b1_b2_post(X_0_b2*X_0_b1.inv());
		BOOST_CHECK_SMALL((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 1e-5);
	}
	// check that the task is well minimized
	BOOST_CHECK_SMALL(oriTask.eval().norm(), 1e-5);
}
