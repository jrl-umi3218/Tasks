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
#include <iostream>
#include <tuple>

// boost
#define BOOST_TEST_MODULE QPSolverTest
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// SpaceVecAlg
#include <SpaceVecAlg>

// RBDyn
#include <EulerIntegration.h>
#include <FK.h>
#include <FV.h>
#include <ID.h>
#include <MultiBody.h>
#include <MultiBodyConfig.h>
#include <MultiBodyGraph.h>

// Tasks
#include "QPConstr.h"
#include "QPSolver.h"
#include "QPTasks.h"

std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeZXZArm()
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;

	MultiBodyGraph mbg;

	double mass = 1.;
	Matrix3d I = Matrix3d::Identity();
	Vector3d h = Vector3d::Zero();

	RBInertia rbi(mass, h, I);

	Body b0(rbi, 0, "b0");
	Body b1(rbi, 1, "b1");
	Body b2(rbi, 2, "b2");
	Body b3(rbi, 3, "b3");

	mbg.addBody(b0);
	mbg.addBody(b1);
	mbg.addBody(b2);
	mbg.addBody(b3);

	Joint j0(Joint::RevZ, true, 0, "j0");
	Joint j1(Joint::RevX, true, 1, "j1");
	Joint j2(Joint::RevZ, true, 2, "j2");

	mbg.addJoint(j0);
	mbg.addJoint(j1);
	mbg.addJoint(j2);

	//  Root     j0       j1     j2
	//  ---- b0 ---- b1 ---- b2 ----b3
	//  Fixed    Z       X       Z


	PTransform to(Vector3d(0., 0.5, 0.));
	PTransform from(Vector3d(0., 0., 0.));


	mbg.linkBodies(0, PTransform::Identity(), 1, from, 0);
	mbg.linkBodies(1, to, 2, from, 1);
	mbg.linkBodies(2, to, 3, from, 2);

	MultiBody mb = mbg.makeMultiBody(0, true);

	MultiBodyConfig mbc(mb);

	mbc.q = {{}, {0.}, {0.}, {0.}};
	mbc.alpha = {{}, {0.}, {0.}, {0.}};
	mbc.alphaD = {{}, {0.}, {0.}, {0.}};
	mbc.jointTorque = {{}, {0.}, {0.}, {0.}};
	ForceVec f0(Vector6d::Zero());
	mbc.force = {f0, f0, f0, f0};

	return std::make_tuple(mb, mbc);
}



BOOST_AUTO_TEST_CASE(QPTaskTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);


	qp::QPSolver solver;

	solver.nrVars(mb, {});
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 3);

	solver.updateEqConstrSize();
	solver.updateInEqConstrSize();


	Vector3d posD = Vector3d(0.707106, 0.707106, 0.);
	qp::PositionTask posTask(mb, 3, posD);
	qp::SetPointTask posTaskSp(mb, &posTask, 10., 1.);

	// Test addTask
	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test PositionTask
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.update(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(posTask.eval().norm(), 0.00001);

	// test removeTask
	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);


	qp::OrientationTask oriTask(mb, 3, RotZ(cst::pi<double>()/2.));
	qp::SetPointTask oriTaskSp(mb, &oriTask, 10., 1.);

	// Test addTask
	solver.addTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test OrientatioTask
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.update(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(oriTask.eval().norm(), 0.00001);

	// test removeTask
	solver.removeTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);


	qp::PostureTask postureTask(mb, {{}, {0.2}, {0.4}, {-0.8}}, 10., 1.);
	solver.addTask(&postureTask);

	// Test PostureTask
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.update(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(postureTask.task().eval().norm(), 0.00001);

	solver.removeTask(&postureTask);


	Vector3d comD = RotZ(cst::pi<double>()/4.)*rbd::computeCoM(mb, mbcInit);
	qp::CoMTask comTask(mb, comD);
	qp::SetPointTask comTaskSp(mb, &comTask, 10., 1.);

	solver.addTask(&comTaskSp);

	// Test CoMTask
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.update(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(comTask.eval().norm(), 0.00001);

	solver.removeTask(&comTaskSp);
}



BOOST_AUTO_TEST_CASE(QPConstrTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);


	qp::QPSolver solver;


	std::vector<qp::Contact> contVec = {{3, {Vector3d::Zero()}, {Vector3d(0., -1., 0.)}}};

	Vector3d posD = Vector3d(0.707106, 0.707106, 0.);
	qp::PositionTask posTask(mb, 3, posD);
	qp::SetPointTask posTaskSp(mb, &posTask, 10., 1.);

	qp::OrientationTask oriTask(mb, 3, RotZ(cst::pi<double>()/2.));
	qp::SetPointTask oriTaskSp(mb, &oriTask, 10., 1.);


	qp::ContactAccConstr contCstrAcc(mb);

	// Test addEqualityConstraint
	solver.addEqualityConstraint(&contCstrAcc);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 1);
	solver.addConstraint(&contCstrAcc);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, contVec);
	solver.updateEqConstrSize();
	solver.updateInEqConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);
	solver.addTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 2);

	posTask.update(mb, mbcInit);
	oriTask.update(mb, mbcInit);
	Vector3d evalPos = posTask.eval();
	Vector3d evalOri = oriTask.eval();

	// Test ContactConstr
	mbcSolv = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.update(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL((posTask.eval() - evalPos).norm(), 0.00001);
	BOOST_CHECK_SMALL((oriTask.eval() - evalOri).norm(), 0.00001);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);
	solver.removeTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test removeEqualityConstraint
	solver.removeEqualityConstraint(&contCstrAcc);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 0);
	solver.removeConstraint(&contCstrAcc);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);



	contVec = {};
	qp::MotionConstr motionCstr(mb);

	solver.addEqualityConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 1);
	solver.addBoundConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	solver.addConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, contVec);
	solver.updateEqConstrSize();
	solver.updateInEqConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test MotionConstr
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.update(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(posTask.eval().norm(), 0.00001);

	MultiBodyConfig mbcTest(mbcSolv);
	mbcTest.jointTorque = {{}, {0.}, {0.}, {0.}};
	InverseDynamics id(mb);
	id.inverseDynamics(mb, mbcTest);

	BOOST_CHECK_SMALL(
		(dofToVector(mb, mbcSolv.jointTorque) -
			dofToVector(mb, mbcTest.jointTorque)).norm(), 0.00001);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test removeEqualityConstraint
	solver.removeEqualityConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 0);
	solver.removeBoundConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 0);
	solver.removeConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);
}
