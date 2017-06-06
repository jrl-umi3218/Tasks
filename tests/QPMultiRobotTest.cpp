// Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
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
#include "Tasks/Bounds.h"
#include "Tasks/QPConstr.h"
#include "Tasks/QPContactConstr.h"
#include "Tasks/QPMotionConstr.h"
#include "Tasks/QPSolver.h"
#include "Tasks/QPTasks.h"

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
		{qp::UnilateralContact(0, 1, "b3", "b3",
			{Vector3d::Zero()}, RotX(cst::pi<double>()/2.), X_b1_b2,
			3, std::tan(cst::pi<double>()/4.))};

	Matrix3d oriD = RotZ(cst::pi<double>()/4.);
	Vector3d posD(oriD*mbc2Init.bodyPosW.back().translation());
	qp::PositionTask posTask(mbs, 1, "b3", posD);
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
	qp::OrientationTask oriTask(mbs, 1, "b3", oriD);
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


// Test Motion constraint
// We setup two arm, one with a fixed base and the second
// with a freebase put on the body b3 of the first robot.
// First we launch an impossible motion to check the dynamics
// After we try with an unilateral contact
// Then we try with a bilateral contact.
BOOST_AUTO_TEST_CASE(TwoArmDDynamicContactTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb1, mb2;
	MultiBodyConfig mbc1Init, mbc2Init;

	std::tie(mb1, mbc1Init) = makeZXZArm();

	forwardKinematics(mb1, mbc1Init);
	forwardVelocity(mb1, mbc1Init);

	std::tie(mb2, mbc2Init) = makeZXZArm(false);
	Vector3d mb2InitPos = mbc1Init.bodyPosW.back().translation();
	Quaterniond mb2InitOri(RotY(cst::pi<double>()/2.));
	mbc2Init.q[0] = {mb2InitOri.w(), mb2InitOri.x(), mb2InitOri.y(), mb2InitOri.z(),
		mb2InitPos.x(), mb2InitPos.y()+ 1, mb2InitPos.z()};
	forwardKinematics(mb2, mbc2Init);
	forwardVelocity(mb2, mbc2Init);

	sva::PTransformd X_0_b1(mbc1Init.bodyPosW.back());
	sva::PTransformd X_0_b2(mbc2Init.bodyPosW.front());
	sva::PTransformd X_b1_b2(X_0_b2*X_0_b1.inv());

	std::vector<MultiBody> mbs = {mb1, mb2};
	std::vector<MultiBodyConfig> mbcs = {mbc1Init, mbc2Init};

	// Test ContactAccConstr constraint
	// Also test PositionTask on the second robot

	qp::QPSolver solver;

	std::vector<Eigen::Vector3d> points =
	{
		Vector3d(0.1, 0., 0.1),
		Vector3d(0.1, 0., -0.1),
		Vector3d(-0.1, 0., -0.1),
		Vector3d(-0.1, 0., 0.1),
	};

	std::vector<Eigen::Vector3d> biPoints =
	{
		Vector3d(0., 0., 0.),
		Vector3d(0., 0., 0.),
		Vector3d(0., 0., 0.),
		Vector3d(0., 0., 0.),
	};

	const int nrGen = 4;
	std::vector<Eigen::Matrix3d> biFrames =
	{
		RotX((1.*cst::pi<double>())/4.),
		RotX((3.*cst::pi<double>())/4.),
		Matrix3d(RotX((1.*cst::pi<double>())/4.)*RotY(cst::pi<double>()/2.)),
		Matrix3d(RotX((3.*cst::pi<double>())/4.)*RotY(cst::pi<double>()/2.)),
	};

	// The fixed robot can pull the other
	std::vector<qp::UnilateralContact> contVecFail =
		{qp::UnilateralContact(0, 1, "b3", "b0",
			points, RotX(-cst::pi<double>()/2.), X_b1_b2,
			nrGen, 0.7)};

	// The fixed robot can push the other
	std::vector<qp::UnilateralContact> contVec =
		{qp::UnilateralContact({0, 1, "b3", "b0"},
			points, RotX(cst::pi<double>()/2.), X_b1_b2,
			nrGen, 0.7)};

	// The fixed robot has non coplanar force apply on the other
	std::vector<qp::BilateralContact> contVecBi =
		{qp::BilateralContact({0, 1, "b3", "b0"},
			biPoints, biFrames, X_b1_b2,
			nrGen, 1.)};

	qp::PostureTask posture1Task(mbs, 0, mbc1Init.q, 2., 1.);
	qp::PostureTask posture2Task(mbs, 1, mbc2Init.q, 2., 1.);

	qp::ContactSpeedConstr contCstrSpeed(0.001);

	const double Inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double>> torqueMin1 = {{},{-Inf},{-Inf},{-Inf}};
	std::vector<std::vector<double>> torqueMax1 = {{},{Inf},{Inf},{Inf}};
	std::vector<std::vector<double>> torqueMin2 = {{0., 0., 0., 0., 0., 0.},
																							{-Inf},{-Inf},{-Inf}};
	std::vector<std::vector<double>> torqueMax2 = {{0., 0., 0., 0., 0., 0.},
																							{Inf},{Inf},{Inf}};
	qp::MotionConstr motion1(mbs, 0, {torqueMin1, torqueMax1});
	qp::MotionConstr motion2(mbs, 1, {torqueMin2, torqueMax2});
	qp::PositiveLambda plCstr;

	motion1.addToSolver(solver);
	motion2.addToSolver(solver);
	plCstr.addToSolver(solver);

	contCstrSpeed.addToSolver(solver);
	solver.addTask(&posture1Task);
	solver.addTask(&posture2Task);

	// check the impossible motion
	solver.nrVars(mbs, contVecFail, {});
	solver.updateConstrSize();

	// 3 dof + 9 dof + 4*nrGen lambda
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 9 + 4*nrGen);
	BOOST_REQUIRE(!solver.solve(mbs, mbcs));


	// check the unilateral motion
	mbcs = {mbc1Init, mbc2Init};
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
		sva::PTransformd X_0_b2_post(mbcs[1].bodyPosW.front());
		sva::PTransformd X_b1_b2_post(X_0_b2*X_0_b1.inv());
		BOOST_CHECK_SMALL((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 1e-5);

		// force in world frame must be the same
		auto f1 = contVec[0].force(solver.lambdaVec(0), contVec[0].r1Cone);
		auto f2 = contVec[0].force(solver.lambdaVec(0), contVec[0].r2Cone);
		BOOST_CHECK_SMALL((f1 + f2).norm(), 1e-5);
	}


	// check the bilateral motion
	mbcs = {mbc1Init, mbc2Init};
	solver.nrVars(mbs, {}, contVecBi);
	solver.updateConstrSize();
	// 3 dof + 9 dof + 4*nrGen lambda
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 9 + 4*nrGen);

	for(int i = 0; i < 1000; ++i)
	{
		//std::cout << i << std::endl;
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		for(std::size_t r = 0; r < mbs.size(); ++r)
		{
			eulerIntegration(mbs[r], mbcs[r], 0.001);

			forwardKinematics(mbs[r], mbcs[r]);
			forwardVelocity(mbs[r], mbcs[r]);
		}
		// check that the link hold
		sva::PTransformd X_0_b1_post(mbcs[0].bodyPosW.back());
		sva::PTransformd X_0_b2_post(mbcs[1].bodyPosW.front());
		sva::PTransformd X_b1_b2_post(X_0_b2*X_0_b1.inv());
		BOOST_CHECK_SMALL((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 1e-5);

		// force in world frame must be the same
		auto f1 = contVec[0].force(solver.lambdaVec(0), contVec[0].r1Cone);
		auto f2 = contVec[0].force(solver.lambdaVec(0), contVec[0].r2Cone);
		BOOST_CHECK_SMALL((f1 + f2).norm(), 1e-5);
	}

	plCstr.removeFromSolver(solver);
	motion2.removeFromSolver(solver);
	motion1.removeFromSolver(solver);
	contCstrSpeed.removeFromSolver(solver);

	solver.removeTask(&posture1Task);
	solver.removeTask(&posture2Task);
}


// Test the MultiCoMTask.
// We try to move the CoM of two arm at a specific position.
BOOST_AUTO_TEST_CASE(TwoArmMultiCoMTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb1, mb2;
	MultiBodyConfig mbc1Init, mbc2Init;

	std::tie(mb1, mbc1Init) = makeZXZArm(true,
		sva::PTransformd(Vector3d(-0.5, 0., 0.)));
	forwardKinematics(mb1, mbc1Init);
	forwardVelocity(mb1, mbc1Init);

	std::tie(mb2, mbc2Init) = makeZXZArm(true,
		sva::PTransformd(Vector3d(0.5, 0., 0.)));
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

	const int nrGen = 3;
	// The fixed robot can push the other
	std::vector<qp::UnilateralContact> contVec =
		{qp::UnilateralContact({0, 1, "b3", "b3"},
		 {Vector3d(0.,0.,0.)}, RotX(cst::pi<double>()/2.), X_b1_b2,
			nrGen, 0.7)};

	qp::PostureTask posture1Task(mbs, 0, mbc1Init.q, 2., 1.);
	qp::PostureTask posture2Task(mbs, 1, mbc2Init.q, 2., 1.);
	Vector3d comD(
		(rbd::computeCoM(mb1, mbc1Init) + rbd::computeCoM(mb2, mbc2Init))/2.
		 + Vector3d(0., 0., 0.5));
	qp::MultiCoMTask multiCoM(mbs, {0,1}, comD, 10., 500.);
	// call this method just for test coverage
	multiCoM.updateInertialParameters(mbs);

	qp::ContactSpeedConstr contCstrSpeed(0.001);

	solver.addTask(&posture1Task);
	solver.addTask(&posture2Task);

	solver.nrVars(mbs, contVec, {});

	// Add MultiCoMTask and ContactSpeedConstr after the nrVars call to test
	// addTask and addToSolver with nrVars init
	solver.addTask(mbs, &multiCoM);
	contCstrSpeed.addToSolver(mbs, solver);

	solver.updateConstrSize();

	// 3 dof + 3 dof + 1*nrGen lambda
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 3 + 1*nrGen);

	for(int i = 0; i < 2000; ++i)
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
	BOOST_CHECK_SMALL(multiCoM.speed().norm(), 1e-3);

	contCstrSpeed.removeFromSolver(solver);

	solver.removeTask(&posture1Task);
	solver.removeTask(&posture2Task);
	solver.removeTask(&multiCoM);
}


// Test the MultiRobotTransformTask
// We try to set he two end effector at the same frame.
BOOST_AUTO_TEST_CASE(MultiRobotTransformTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb1, mb2;
	MultiBodyConfig mbc1Init, mbc2Init;

	std::tie(mb1, mbc1Init) = makeZXZArm(true,
		sva::PTransformd(sva::RotZ(-cst::pi<double>()/4.), Vector3d(-0.5, 0., 0.)));
	forwardKinematics(mb1, mbc1Init);
	forwardVelocity(mb1, mbc1Init);

	std::tie(mb2, mbc2Init) = makeZXZArm(false,
		 sva::PTransformd(sva::RotZ(cst::pi<double>()/2.), Vector3d(0.5, 0., 0.)));
	forwardKinematics(mb2, mbc2Init);
	forwardVelocity(mb2, mbc2Init);

	std::vector<MultiBody> mbs = {mb1, mb2};
	std::vector<MultiBodyConfig> mbcs = {mbc1Init, mbc2Init};

	// Test ContactAccConstr constraint
	// Also test PositionTask on the second robot

	qp::QPSolver solver;

	qp::PostureTask posture1Task(mbs, 0, mbc1Init.q, 0.1, 10.);
	qp::PostureTask posture2Task(mbs, 1, mbc2Init.q, 0.1, 10.);
	qp::MultiRobotTransformTask mrtt(mbs, 0, 1, "b3", "b3",
		sva::PTransformd(sva::RotZ(-cst::pi<double>()/8.)),
		sva::PTransformd::Identity(), 100., 1000.);
	mrtt.dimWeight((Vector6d() << 0., 0., 1., 1., 1.,0.).finished());

	solver.addTask(&posture1Task);
	solver.addTask(&posture2Task);
	solver.addTask(&mrtt);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();
	// 3 dof + 9 dof
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 9);

	for(int i = 0; i < 2000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		for(std::size_t r = 0; r < mbs.size(); ++r)
		{
			eulerIntegration(mbs[r], mbcs[r], 0.005);

			forwardKinematics(mbs[r], mbcs[r]);
			forwardVelocity(mbs[r], mbcs[r]);
		}
	}
	BOOST_CHECK_SMALL(mrtt.eval().norm(), 1e-3);

	solver.removeTask(&posture1Task);
	solver.removeTask(&posture2Task);
	solver.removeTask(&mrtt);
}


// Test the TorqueTask
BOOST_AUTO_TEST_CASE(TorqueTaskTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb1, mb2;
	MultiBodyConfig mbc1Init, mbc2Init;

	std::tie(mb1, mbc1Init) = makeZXZArm(true,
		sva::PTransformd(sva::RotZ(-cst::pi<double>()/4.), Vector3d(-0.5, 0., 0.)));
	forwardKinematics(mb1, mbc1Init);
	forwardVelocity(mb1, mbc1Init);

	std::tie(mb2, mbc2Init) = makeZXZArm(false,
		 sva::PTransformd(sva::RotZ(cst::pi<double>()/2.), Vector3d(0.5, 0., 0.)));
	forwardKinematics(mb2, mbc2Init);
	forwardVelocity(mb2, mbc2Init);

	std::vector<MultiBody> mbs = {mb1, mb2};
	std::vector<MultiBodyConfig> mbcs = {mbc1Init, mbc2Init};

	// Test ContactAccConstr constraint
	// Also test PositionTask on the second robot

	qp::QPSolver solver;

        std::vector<std::vector<double>> lsup;
        std::vector<std::vector<double>> linf;
        std::vector<double> sup;
        std::vector<double> inf;

        for(const auto j: mb1.joints())
        {
          sup.resize(j.dof());
          inf.resize(j.dof());
          std::fill(sup.begin(), sup.end(), 1e4);
          std::fill(inf.begin(), inf.end(), -1e4);
          lsup.push_back(sup);
          linf.push_back(inf);
        }

        TorqueBound tb(lsup, linf);

	qp::PostureTask posture1Task(mbs, 0, mbc1Init.q, 0.1, 10.);
	qp::PostureTask posture2Task(mbs, 1, mbc2Init.q, 0.1, 10.);
	qp::TorqueTask tt(mbs, 0, tb, 1);

	solver.addTask(&posture1Task);
	solver.addTask(&posture2Task);
	solver.addTask(&tt);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();
	// 3 dof + 9 dof
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 9);

	for(int i = 0; i < 2000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		for(std::size_t r = 0; r < mbs.size(); ++r)
		{
			eulerIntegration(mbs[r], mbcs[r], 0.005);

			forwardKinematics(mbs[r], mbcs[r]);
			forwardVelocity(mbs[r], mbcs[r]);
		}
	}
	solver.removeTask(&posture1Task);
	solver.removeTask(&posture2Task);
	solver.removeTask(&tt);
}
