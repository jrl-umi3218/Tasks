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
#include <tuple>

// boost
#define BOOST_TEST_MODULE QPSolverTest
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// Eigen
#include <unsupported/Eigen/Polynomials>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include <RBDyn/EulerIntegration.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/ID.h>

// sch
#include <sch/S_Object/S_Sphere.h>
#include <sch/S_Object/S_Box.h>
#include <sch/CD/CD_Pair.h>

// Tasks
#include "Tasks/Bounds.h"
#include "Tasks/QPConstr.h"
#include "Tasks/QPContactConstr.h"
#include "Tasks/QPMotionConstr.h"
#include "Tasks/QPSolver.h"
#include "Tasks/QPTasks.h"

// Arms
#include "arms.h"


BOOST_AUTO_TEST_CASE(FrictionConeTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	double angle = cst::pi<double>()/4.;
	double mu = std::tan(angle);

	qp::FrictionCone cone(Matrix3d::Identity(), 4, mu);
	for(Vector3d v: cone.generators)
	{
		// check cone equation x^2 + y^2 = z^2(tan(angle)^2)
		BOOST_CHECK_SMALL(std::pow(v.x(), 2) + std::pow(v.y(), 2) -
			std::pow(v.z(), 2)*std::pow(mu, 2), 0.00001);
	}

	Matrix3d rep;
	rep << 0., 1., 0.,
					0., 0., 1.,
					1., 0., 0.;
	qp::FrictionCone cone2(rep, 4, mu);
	for(Vector3d v: cone2.generators)
	{
		// check cone equation in rep frame z^2 + y^2 = x^2(tan(angle)^2)
		BOOST_CHECK_SMALL(std::pow(v.z(), 2) + std::pow(v.y(), 2) -
			std::pow(v.x(), 2)*std::pow(mu, 2), 0.00001);
	}
}


// TODO contacts Test



BOOST_AUTO_TEST_CASE(QPTaskTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	std::vector<MultiBody> mbs = {mb};
	std::vector<MultiBodyConfig> mbcs(1);

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);


	qp::QPSolver solver;

	solver.nrVars(mbs, {}, {});
	BOOST_CHECK_EQUAL(solver.nrVars(), 3);

	solver.updateConstrSize();

	Vector3d posD = Vector3d(0.707106, 0.707106, 0.);
	qp::PositionTask posTask(mbs, 0, "b3", posD);
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 10., 1.);

	// Test addTask
	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test PositionTask
	mbcs[0] = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mb, mbcs[0], 0.001);

		forwardKinematics(mb, mbcs[0]);
		forwardVelocity(mb, mbcs[0]);
	}

	BOOST_CHECK_SMALL(posTask.eval().norm(), 0.00001);

	// test removeTask
	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);


	qp::OrientationTask oriTask(mbs, 0, "b3", RotZ(cst::pi<double>()/2.));
	qp::SetPointTask oriTaskSp(mbs, 0, &oriTask, 10., 1.);

	// Test addTask
	solver.addTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test OrientatioTask
	mbcs[0] = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mb, mbcs[0], 0.001);

		forwardKinematics(mb, mbcs[0]);
		forwardVelocity(mb, mbcs[0]);
	}

	BOOST_CHECK_SMALL(oriTask.eval().norm(), 0.00001);

	// test removeTask
	solver.removeTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);


	qp::PostureTask postureTask(mbs, 0, {{}, {0.2}, {0.4}, {-0.8}}, 10., 1.);
	postureTask.jointsStiffness(mbs, {{"j2", 10.}});
	solver.addTask(&postureTask);

	// Test PostureTask
	mbcs[0] = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mb, mbcs[0], 0.001);

		forwardKinematics(mb, mbcs[0]);
		forwardVelocity(mb, mbcs[0]);
	}

	BOOST_CHECK_SMALL(postureTask.task().eval().norm(), 0.00001);

	solver.removeTask(&postureTask);


	Vector3d comD(RotZ(cst::pi<double>()/4.)*rbd::computeCoM(mb, mbcInit));
	qp::CoMTask comTask(mbs, 0, comD);
	qp::SetPointTask comTaskSp(mbs, 0, &comTask, 10., 1.);

	solver.addTask(&comTaskSp);

	// Test CoMTask
	mbcs[0] = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mb, mbcs[0], 0.001);

		forwardKinematics(mb, mbcs[0]);
		forwardVelocity(mb, mbcs[0]);
	}

	BOOST_CHECK_SMALL(comTask.eval().norm(), 0.00001);

	solver.removeTask(&comTaskSp);


	qp::LinVelocityTask linVelocityTask(mbs, 0, "b3", -Vector3d::UnitX()*0.005);
	qp::SetPointTask linVelocityTaskSp(mbs, 0, &linVelocityTask, 1000., 10000.);

	solver.addTask(&linVelocityTaskSp);

	// Test LinVelocityTask
	mbcs[0] = mbcInit;
	for(int i = 0; i < 4000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mb, mbcs[0], 0.001);

		forwardKinematics(mb, mbcs[0]);
		forwardVelocity(mb, mbcs[0]);
	}

	// error is huge, why ?
	BOOST_CHECK_SMALL(linVelocityTask.eval().norm(), 0.001);

	solver.removeTask(&linVelocityTaskSp);
}



BOOST_AUTO_TEST_CASE(QPConstrTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb, mbEnv;
	MultiBodyConfig mbcInit, mbcEnv;

	std::tie(mb, mbcInit) = makeZXZArm();
	std::tie(mbEnv, mbcEnv) = makeEnv();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);
	forwardKinematics(mbEnv, mbcEnv);
	forwardVelocity(mbEnv, mbcEnv);

	std::vector<MultiBody> mbs = {mb, mbEnv};
	std::vector<MultiBodyConfig> mbcs = {mbcInit, mbcEnv};


	qp::QPSolver solver;


	std::vector<qp::UnilateralContact> contVec =
		{qp::UnilateralContact(0, 1, "b3", "b0", {Vector3d::Zero()}, Matrix3d::Identity(),
			sva::PTransformd::Identity(), 3, std::tan(cst::pi<double>()/4.))};

	Vector3d posD = Vector3d(0.707106, 0.707106, 0.);
	qp::PositionTask posTask(mbs, 0, "b3", posD);
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 10., 1.);

	qp::OrientationTask oriTask(mbs, 0, "b3", RotZ(cst::pi<double>()/2.));
	qp::SetPointTask oriTaskSp(mbs, 0, &oriTask, 10., 1.);


	qp::ContactAccConstr contCstrAcc;

	// Test addEqualityConstraint
	mbcs[0] = mbcInit;
	contCstrAcc.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mbs, contVec, {});
	solver.updateConstrSize();

	// check vars number is 3 dof + 3 lambda
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 3);

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);
	solver.addTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 2);

	solver.data().computeNormalAccB(mbs, mbcs);
	posTask.update(mbs, mbcs, solver.data());
	oriTask.update(mbs, mbcs, solver.data());
	Vector3d evalPos = posTask.eval();
	Vector3d evalOri = oriTask.eval();

	// Test ContactAccConstr
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}

	BOOST_CHECK_SMALL((posTask.eval() - evalPos).norm(), 0.00001);
	BOOST_CHECK_SMALL((oriTask.eval() - evalOri).norm(), 0.00001);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);
	solver.removeTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test removeEqualityConstraint
	contCstrAcc.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);



	// Test ContactSpeedConstr
	mbcs[0] = mbcInit;
	qp::ContactSpeedConstr contCstrSpeed(0.001);

	// Test addEqualityConstraint
	contCstrSpeed.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mbs, contVec, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);
	solver.addTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 2);

	solver.data().computeNormalAccB(mbs, mbcs);
	posTask.update(mbs, mbcs, solver.data());
	oriTask.update(mbs, mbcs, solver.data());
	evalPos = posTask.eval();
	evalOri = oriTask.eval();

	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}

	BOOST_CHECK_SMALL((posTask.eval() - evalPos).norm(), 0.00001);
	BOOST_CHECK_SMALL((oriTask.eval() - evalOri).norm(), 0.00001);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);
	solver.removeTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test removeEqualityConstraint
	solver.removeEqualityConstraint(&contCstrSpeed);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 0);
	solver.removeConstraint(&contCstrSpeed);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);



	// MotionConstr test
	mbcs[0] = mbcInit;
	contVec = {};
	double Inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double>> torqueMin = {{},{-Inf},{-Inf},{-Inf}};
	std::vector<std::vector<double>> torqueMax = {{},{Inf},{Inf},{Inf}};
	qp::MotionConstr motionCstr(mbs, 0, {torqueMin, torqueMax});
	qp::PositiveLambda plCstr;

	motionCstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	plCstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 2);

	solver.nrVars(mbs, contVec, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test MotionConstr
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}
	motionCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
	motionCstr.torque(mbs, mbcs);

	BOOST_CHECK_SMALL(posTask.eval().norm(), 5e-05);

	MultiBodyConfig mbcTest(mbcs[0]);
	mbcTest.jointTorque = {{}, {0.}, {0.}, {0.}};
	InverseDynamics id(mb);
	id.inverseDynamics(mb, mbcTest);

	BOOST_CHECK_SMALL(
		(dofToVector(mb, mbcs[0].jointTorque) -
			dofToVector(mb, mbcTest.jointTorque)).norm(), 0.00001);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test removeEqualityConstraint
	motionCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);


	// MotionConstr test with contact
	mbcs[0] = mbcInit;
	contVec =
		{qp::UnilateralContact(0, 1, "b3", "b0", {Vector3d::Zero()}, Matrix3d::Identity(),
			sva::PTransformd::Identity(), 3, std::tan(cst::pi<double>()/4.))};

	solver.addGenInequalityConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 1);
	solver.addConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 2);

	solver.nrVars(mbs, contVec, {});
	solver.updateConstrSize();

	// check vars number is 3 dof + 3 lambda
	BOOST_CHECK_EQUAL(solver.nrVars(), 3 + 3);

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test MotionConstr
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}

	BOOST_CHECK_SMALL(posTask.eval().norm(), 5e-5);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test removeEqualityConstraint
	motionCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	plCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);
}



BOOST_AUTO_TEST_CASE(QPJointLimitsTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<rbd::MultiBody> mbs = {mb};
	std::vector<rbd::MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	int bodyI = mb.bodyIndexByName("b3");
	qp::PositionTask posTask(mbs, 0, "b3",
		RotZ(cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 10., 1.);

	double inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double> > lBound = {{}, {-cst::pi<double>()/4.}, {-inf}, {-inf}};
	std::vector<std::vector<double> > uBound = {{}, {cst::pi<double>()/4.}, {inf}, {inf}};

	qp::JointLimitsConstr jointConstr(mbs, 0, {lBound, uBound}, 0.001);

	// Test add*Constraint
	solver.addBoundConstraint(&jointConstr);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	solver.addConstraint(&jointConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	// Test JointLimitsConstr
	mbcs[0] = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		BOOST_REQUIRE_GT(mbcs[0].q[1][0], -cst::pi<double>()/4. - 0.01);
	}

	posTask.position(RotZ(-cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	mbcs[0] = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		BOOST_REQUIRE_LT(mbcs[0].q[1][0], cst::pi<double>()/4. + 0.01);
	}

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test remove*Constraint
	solver.removeBoundConstraint(&jointConstr);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 0);
	solver.removeConstraint(&jointConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);
}



BOOST_AUTO_TEST_CASE(QPDamperJointLimitsTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<rbd::MultiBody> mbs = {mb};
	std::vector<rbd::MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	int bodyI = mb.bodyIndexByName("b3");
	qp::PositionTask posTask(mbs, 0, "b3",
		RotZ(cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 10., 1.);

	double inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double> > lBound = {{}, {-cst::pi<double>()/4.}, {-inf}, {-inf}};
	std::vector<std::vector<double> > uBound = {{}, {cst::pi<double>()/4.}, {inf}, {inf}};
	std::vector<std::vector<double> > lVel = {{}, {-inf}, {-inf}, {-inf}};
	std::vector<std::vector<double> > uVel = {{}, {inf}, {inf}, {inf}};

	qp::DamperJointLimitsConstr dampJointConstr(mbs, 0, {lBound, uBound},
		{lVel, uVel}, 0.125, 0.025, 1., 0.001);

	// Test add*Constraint
	dampJointConstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	// Test JointLimitsConstr
	mbcs[0] = mbcInit;
	for(int i = 0; i < 2000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		BOOST_REQUIRE_GT(mbcs[0].q[1][0], -(cst::pi<double>()/4.)*(1. - 0.025) - 0.01);
	}

	posTask.position(RotZ(-cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	mbcs[0] = mbcInit;
	for(int i = 0; i < 2000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		BOOST_REQUIRE_LT(mbcs[0].q[1][0], (cst::pi<double>()/4.)*(1. - 0.025) + 0.01);
	}

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test remove*Constraint
	dampJointConstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);
}

BOOST_AUTO_TEST_CASE(QPMimicJointTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	/* Make a very simple gripper */
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
	Joint j1(Joint::RevZ, true, "j1");
	j1.makeMimic("j0", -1.0, 0.0);
	Joint j2(Joint::RevZ, true, "j2");
	j2.makeMimic("j0", 0.0, 0.5);

	mbg.addJoint(j0);
	mbg.addJoint(j1);
	mbg.addJoint(j2);

	PTransformd from = sva::PTransformd::Identity();

	mbg.linkBodies("b0", from, "b1", from, "j0");
	mbg.linkBodies("b1", from, "b2", from, "j1");
	mbg.linkBodies("b2", from, "b3", from, "j2");

	MultiBody mb = mbg.makeMultiBody("b0", true, from);

	MultiBodyConfig mbcInit(mb);
	mbcInit.zero(mb);
	BOOST_CHECK_EQUAL(mbcInit.q[3][0], 0.5);

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<rbd::MultiBody> mbs = {mb};
	std::vector<rbd::MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	tasks::qp::PostureTask pt(mbs, 0, mbcInit.q, 100.0, 10.0);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	solver.addTask(&pt);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	// Test mimic joint in mbc
	mbcs[0] = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		BOOST_REQUIRE(mbcs[0].q[1][0] == -mbcs[0].q[2][0]);
		BOOST_REQUIRE(mbcs[0].q[3][0] == 0.5);
	}
	BOOST_CHECK_SMALL(mbcs[0].q[1][0], 1e-4);

	mbcs[0] = mbcInit;
	pt.posture({{}, {0.2}, {-0.2}, {0.5}});
	for(int i = 0; i < 2000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		BOOST_REQUIRE(mbcs[0].q[1][0] == -mbcs[0].q[2][0]);
		BOOST_REQUIRE(mbcs[0].q[3][0] == 0.5);
	}
	BOOST_CHECK_SMALL(mbcs[0].q[1][0] - 0.2, 1e-4);

	solver.removeTask(&pt);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);
}


BOOST_AUTO_TEST_CASE(QPTorqueLimitsTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<rbd::MultiBody> mbs = {mb};
	std::vector<rbd::MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	int bodyI = mb.bodyIndexByName("b3");
	qp::PositionTask posTask(mbs, 0, "b3",
		RotZ(cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 10., 1.);

	std::vector<std::vector<double> > lBound = {{}, {-30.}, {-30.}, {-30.}};
	std::vector<std::vector<double> > uBound = {{}, {30.}, {30.}, {30.}};

	qp::MotionConstr motionCstr(mbs, 0, {lBound, uBound});
	qp::PositiveLambda plCstr;

	// Test add*Constraint

	solver.addGenInequalityConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 1);
	solver.addConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	plCstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 2);

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();


	// Test MotionConstr torque limits
	mbcs[0] = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		motionCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
		motionCstr.torque(mbs, mbcs);
		for(int i = 0; i < 3; ++i)
		{
			BOOST_REQUIRE_GT(mbcs[0].jointTorque[i+1][0], lBound[i+1][0] - 0.001);
			BOOST_REQUIRE_LT(mbcs[0].jointTorque[i+1][0], uBound[i+1][0] + 0.001);
		}
	}

	posTask.position(mbcInit.bodyPosW[bodyI].translation());
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		motionCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
		motionCstr.torque(mbs, mbcs);
		for(int i = 0; i < 3; ++i)
		{
			BOOST_REQUIRE_GT(mbcs[0].jointTorque[i+1][0], lBound[i+1][0] - 0.001);
			BOOST_REQUIRE_LT(mbcs[0].jointTorque[i+1][0], uBound[i+1][0] + 0.001);
		}
	}

	motionCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);


	Eigen::VectorXd lpoly(2), upoly(2);
	Eigen::VectorXd null;
	lpoly << -30, 1.;
	upoly << 30, 1.;
	std::vector<std::vector<Eigen::VectorXd> > lBoundPoly = {{null}, {lpoly}, {lpoly}, {lpoly}};
	std::vector<std::vector<Eigen::VectorXd> > uBoundPoly = {{null}, {upoly}, {upoly}, {upoly}};
	qp::MotionPolyConstr motionPolyCstr(mbs, 0, {lBoundPoly, uBoundPoly});

	motionPolyCstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 2);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	// Test MotionPolyConstr torque limits
	mbcs[0] = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		std::vector<std::vector<double>> oldQ = mbcs[0].q;
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		motionPolyCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
		motionPolyCstr.torque(mbs, mbcs);
		for(int i = 0; i < 3; ++i)
		{
			BOOST_REQUIRE_GT(mbcs[0].jointTorque[i+1][0], Eigen::poly_eval(lBoundPoly[i+1][0], oldQ[i+1][0]));
			BOOST_REQUIRE_LT(mbcs[0].jointTorque[i+1][0], Eigen::poly_eval(uBoundPoly[i+1][0], oldQ[i+1][0]));
		}
	}

	posTask.position(mbcInit.bodyPosW[bodyI].translation());
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		std::vector<std::vector<double>> oldQ = mbcs[0].q;
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
		motionPolyCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
		motionPolyCstr.torque(mbs, mbcs);
		for(int i = 0; i < 3; ++i)
		{
			BOOST_REQUIRE_GT(mbcs[0].jointTorque[i+1][0], Eigen::poly_eval(lBoundPoly[i+1][0], oldQ[i+1][0]));
			BOOST_REQUIRE_LT(mbcs[0].jointTorque[i+1][0], Eigen::poly_eval(uBoundPoly[i+1][0], oldQ[i+1][0]));
		}
	}

	motionPolyCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	plCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);
}



BOOST_AUTO_TEST_CASE(QPAutoCollTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<MultiBody> mbs = {mb};
	std::vector<MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	int bodyI = mb.bodyIndexByName("b3");
	qp::PositionTask posTask(mbs, 0, "b3", mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 50., 1.);

	sch::S_Sphere b0(0.25), b3(0.25);
	sch::CD_Pair pair(&b0, &b3);

	PTransformd I = PTransformd::Identity();
	qp::CollisionConstr autoCollConstr(mbs, 0.001);
	int collId1 = 10;
	autoCollConstr.addCollision(mbs, collId1,
		0, "b0", &b0, I,
		0, "b3", &b3, I,
		0.01, 0.005, 1.);
	BOOST_CHECK_EQUAL(autoCollConstr.nrCollisions(), 1);

	// Test addInequalityConstraint
	solver.addInequalityConstraint(&autoCollConstr);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 1);
	solver.addConstraint(&autoCollConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	mbcs[0] = mbcInit;
	std::ofstream distHard("selfDistHard.py");
	distHard << "dist = [";
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0. ? std::sqrt(dist) : -std::sqrt(-dist);
		distHard << dist << ", ";
		BOOST_REQUIRE_GT(dist, 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}
	distHard << "]" << std::endl;

	autoCollConstr.rmCollision(collId1);
	BOOST_CHECK_EQUAL(autoCollConstr.nrCollisions(), 0);


	// test automatic damping computation
	autoCollConstr.addCollision(mbs, collId1,
		0, "b0", &b0, I,
		0, "b3", &b3, I,
		0.1, 0.01, 0., 0.1);
	BOOST_CHECK_EQUAL(autoCollConstr.nrCollisions(), 1);
	posTask.position(mbcInit.bodyPosW[bodyI].translation());

	mbcs[0] = mbcInit;
	std::ofstream distSoft("selfDistSoft.py");
	distSoft << "dist = [";
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0. ? std::sqrt(dist) : -std::sqrt(-dist);
		distSoft << dist << ", ";
		BOOST_REQUIRE_GT(dist, 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}
	distSoft << "]" << std::endl;

	autoCollConstr.rmCollision(collId1);
	BOOST_CHECK_EQUAL(autoCollConstr.nrCollisions(), 0);


	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test remove*Constraint
	solver.removeInequalityConstraint(&autoCollConstr);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 0);
	solver.removeConstraint(&autoCollConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);
}


BOOST_AUTO_TEST_CASE(QPStaticEnvCollTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb, mbEnv;
	MultiBodyConfig mbcInit, mbcEnv;

	std::tie(mb, mbcInit) = makeZXZArm();
	std::tie(mbEnv, mbcEnv) = makeEnv();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);
	forwardKinematics(mbEnv, mbcEnv);
	forwardVelocity(mbEnv, mbcEnv);

	std::vector<MultiBody> mbs = {mb, mbEnv};
	std::vector<MultiBodyConfig> mbcs = {mbcInit, mbcEnv};


	qp::QPSolver solver;

	int bodyI = mb.bodyIndexByName("b3");
	qp::PositionTask posTask(mbs, 0, "b3", mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 50., 1.);

	sch::S_Sphere b0(0.25), b3(0.25);
	sch::CD_Pair pair(&b0, &b3);

	b0.setTransformation(qp::tosch(mbcInit.bodyPosW[0]));

	PTransformd I = PTransformd::Identity();
	qp::CollisionConstr seCollConstr(mbs, 0.001);
	int collId1 = 10;
	seCollConstr.addCollision(mbs, collId1,
		0, "b3", &b3, I,
		1, "b0", &b0, I,
		0.01, 0.005, 1.);
	BOOST_CHECK_EQUAL(seCollConstr.nrCollisions(), 1);

	// Test addInequalityConstraint
	solver.addInequalityConstraint(&seCollConstr);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 1);
	solver.addConstraint(&seCollConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	mbcs[0] = mbcInit;
	std::ofstream distHard("staticDistHard.py");
	distHard << "dist = [";
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
		distHard << dist << ", ";
		BOOST_REQUIRE_GT(dist, 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}
	distHard << "]" << std::endl;

	seCollConstr.rmCollision(collId1);
	BOOST_CHECK_EQUAL(seCollConstr.nrCollisions(), 0);

	// test damping computation
	seCollConstr.addCollision(mbs, collId1,
		0, "b3", &b3, I,
		1, "b0", &b0, I,
		0.1, 0.01, 0., 0.1);
	BOOST_CHECK_EQUAL(seCollConstr.nrCollisions(), 1);
	posTask.position(mbcInit.bodyPosW[bodyI].translation());

	mbcs[0] = mbcInit;
	std::ofstream distSoft("staticDistSoft.py");
	distSoft << "dist = [";
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
		distSoft << dist << ", ";
		BOOST_REQUIRE_GT(dist, 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}
	distSoft << "]" << std::endl;


	seCollConstr.rmCollision(collId1);
	BOOST_CHECK_EQUAL(seCollConstr.nrCollisions(), 0);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test remove*Constraint
	solver.removeInequalityConstraint(&seCollConstr);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 0);
	solver.removeConstraint(&seCollConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);
}


BOOST_AUTO_TEST_CASE(QPBilatContactTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb, mbEnv;
	MultiBodyConfig mbcInit, mbcEnv;

	std::tie(mb, mbcInit) = makeZXZArm(false);
	std::tie(mbEnv, mbcEnv) = makeEnv();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);
	forwardKinematics(mbEnv, mbcEnv);
	forwardVelocity(mbEnv, mbcEnv);

	std::vector<MultiBody> mbs = {mb, mbEnv};
	std::vector<MultiBodyConfig> mbcs = {mbcInit, mbcEnv};


	qp::QPSolver solver;

	double Inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double>> torqueMin = {{0., 0., 0., 0., 0., 0.},{-Inf},{-Inf},{-Inf}};
	std::vector<std::vector<double>> torqueMax = {{0., 0., 0., 0., 0., 0.},{Inf},{Inf},{Inf}};
	qp::MotionConstr motionCstr(mbs, 0, {torqueMin, torqueMax});
	qp::PositiveLambda plCstr;
	qp::ContactAccConstr contCstrAcc;

	solver.addGenInequalityConstraint(&motionCstr);
	solver.addConstraint(&motionCstr);

	solver.addEqualityConstraint(&contCstrAcc);
	solver.addConstraint(&contCstrAcc);

	plCstr.addToSolver(solver);

	std::vector<Eigen::Vector3d> points =
		{
			Vector3d(0.1, 0.1, 0.),
			Vector3d(-0.1, 0.1, 0.),
			Vector3d(-0.1, -0.1, 0.),
			Vector3d(0.1, -0.1, 0.)
		};

	std::vector<Eigen::Matrix3d> biFrames =
		{
			sva::RotY((0.*cst::pi<double>())/2.),
			sva::RotY((1.*cst::pi<double>())/2.),
			sva::RotY((2.*cst::pi<double>())/2.),
			sva::RotY((3.*cst::pi<double>())/2.),
		};

	std::vector<qp::UnilateralContact> uni =
		{qp::UnilateralContact(0, 1, "b0", "b0",
			points, Matrix3d::Identity(), sva::PTransformd::Identity(),
			3, 0.7)};
	std::vector<qp::BilateralContact> bi =
		{qp::BilateralContact(0, 1, "b0", "b0",
			points, biFrames, sva::PTransformd::Identity(),
			3, 0.7)};

	solver.nrVars(mbs, uni, {});
	solver.updateConstrSize();
	BOOST_CHECK_EQUAL(solver.nrVars(), 9 + 4*3);

	// This stance with unilateral contac is impossible so the solver must fail
	// The forces are apply on the Z axis and the gravity come frome the Y axis
	mbcs[0] = mbcInit;
	BOOST_REQUIRE(!solver.solve(mbs, mbcs));


	// We test it again with bilateral contact to check that the stance is now
	// valid
	// the forces are apply on the Z, Y, -Z and -Y axis
	solver.nrVars(mbs, {}, bi);
	solver.updateConstrSize();
	BOOST_CHECK_EQUAL(solver.nrVars(), 9 + 4*3);

	mbcs[0] = mbcInit;
	for(int i = 0; i < 10; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}


	plCstr.removeFromSolver(solver);

	solver.removeEqualityConstraint(&contCstrAcc);
	solver.removeConstraint(&contCstrAcc);

	solver.removeConstraint(&motionCstr);
	solver.removeGenInequalityConstraint(&motionCstr);
}


Eigen::Vector6d compute6dError(const sva::PTransformd& b1, const sva::PTransformd& b2)
{
	Eigen::Vector6d error;
	error.head<3>() = sva::rotationError(b1.rotation(), b2.rotation());
	error.tail<3>() = b1.translation() - b2.translation();
	return error;
}


Eigen::Vector6d compute6dErrorInB1(const sva::PTransformd& b1, const sva::PTransformd& b2)
{
	sva::MotionVecd error;
	error.angular() = sva::rotationError(b1.rotation(), b2.rotation());
	error.linear() = b1.translation() - b2.translation();
	return (sva::PTransformd(b1.rotation())*error).vector();
}


double computeDofError(const sva::PTransformd& b1, const sva::PTransformd& b2,
										const Eigen::MatrixXd& dof)
{
	return (dof*compute6dError(b1, b2)).norm();
}


double computeDofErrorInB1(const sva::PTransformd& b1, const sva::PTransformd& b2,
										const Eigen::MatrixXd& dof)
{
	return (dof*compute6dErrorInB1(b1, b2)).norm();
}


BOOST_AUTO_TEST_CASE(QPDofContactsTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb, mbEnv;
	MultiBodyConfig mbcInit, mbcEnv;

	std::tie(mb, mbcInit) = makeZXZArm(false);
	std::tie(mbEnv, mbcEnv) = makeEnv();

	Vector4d quat(Vector4d::Random().normalized());
	mbcInit.q[0] = {quat.w(), quat.x(), quat.y(), quat.z(), 0., 0., 0.};

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);
	forwardKinematics(mbEnv, mbcEnv);
	forwardVelocity(mbEnv, mbcEnv);

	std::vector<MultiBody> mbs = {mb, mbEnv};
	std::vector<MultiBodyConfig> mbcs = {mbcInit, mbcEnv};


	qp::QPSolver solver;

	sva::PTransformd X_b1_cf(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());

	qp::ContactPosConstr contCstrSpeed(0.005);
	// target in cf coordinate is transform into world frame
	qp::PositionTask posTask(mbs, 0, "b0",
		(X_b1_cf*mbcInit.bodyPosW[0]).rotation().transpose()*Vector3d(1., 1., -1.));
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 10., 1.);
	// Rotation in cf coordinate
	qp::OrientationTask oriTask(mbs, 0, "b0",
		(X_b1_cf*mbcInit.bodyPosW[0]).rotation().transpose()*sva::RotY(0.1)*sva::RotX(0.5));
	qp::SetPointTask oriTaskSp(mbs, 0, &oriTask, 10., 1.);

	contCstrSpeed.addToSolver(solver);
	solver.addTask(&posTaskSp);

	std::vector<Eigen::Vector3d> points =
		{
			Vector3d(0.1, 0.1, 0.),
			Vector3d(-0.1, 0.1, 0.),
			Vector3d(-0.1, -0.1, 0.),
			Vector3d(0.1, -0.1, 0.)
		};

	sva::PTransformd X_b1_b2(mbcEnv.bodyPosW[0]*mbcInit.bodyPosW[0].inv());
	std::vector<qp::UnilateralContact> uni =
		{qp::UnilateralContact(0, 1, "b0", "b0",
			points, Matrix3d::Identity(), X_b1_b2,
			3, 0.7, X_b1_cf)};

	// contactDof must be provide in r1BodyId frame
	MatrixXd contactDof(5, 6);
	contactDof.setZero();
	contactDof(0, 0) = 1.;
	contactDof(1, 1) = 1.;
	contactDof(2, 2) = 1.;
	contactDof(3, 3) = 1.;
	contactDof(4, 4) = 1.;

	// test Z free
	contCstrSpeed.addDofContact({0, 1, "b0", "b0"}, contactDof);
	solver.nrVars(mbs, uni, {});
	solver.updateConstrSize();

	mbcs[0] = mbcInit;
	for(int i = 0; i < 100; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.005);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}

	BOOST_CHECK_SMALL(computeDofErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0],
		X_b1_cf*mbcInit.bodyPosW[0], contactDof), 1e-6);
	BOOST_CHECK_GT(std::pow(
		compute6dErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0],
			X_b1_cf*mbcInit.bodyPosW[0])(5), 2), 0.1);


	// test Y Free and updateDofContacts
	contactDof.setZero();
	contactDof(0, 0) = 1.;
	contactDof(1, 1) = 1.;
	contactDof(2, 2) = 1.;
	contactDof(3, 3) = 1.;
	contactDof(4, 5) = 1.;
	contCstrSpeed.resetDofContacts();
	contCstrSpeed.addDofContact({0, 1, "b0", "b0"}, contactDof);
	contCstrSpeed.updateDofContacts();
	mbcs[0] = mbcInit;
	for(int i = 0; i < 100; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.005);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}

	BOOST_CHECK_SMALL(computeDofErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0],
		X_b1_cf*mbcInit.bodyPosW[0], contactDof), 1e-6);
	BOOST_CHECK_GT(std::pow(
		compute6dErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0],
			X_b1_cf*mbcInit.bodyPosW[0])(4), 2), 0.1);


	// test WX Free and updateDofContacts
	contactDof.setZero();
	contactDof(0, 1) = 1.;
	contactDof(1, 2) = 1.;
	contactDof(2, 3) = 1.;
	contactDof(3, 4) = 1.;
	contactDof(4, 5) = 1.;
	contCstrSpeed.resetDofContacts();
	contCstrSpeed.addDofContact({0, 1, "b0", "b0"}, contactDof);
	contCstrSpeed.updateDofContacts();

	// add the orientation task in cf coordinate
	solver.addTask(&oriTaskSp);
	solver.updateTasksNrVars(mbs);
	mbcs[0] = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.005);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}

	BOOST_CHECK_SMALL(computeDofErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0],
		X_b1_cf*mbcInit.bodyPosW[0], contactDof), 1e-4);
	BOOST_CHECK_GT(std::abs(
		compute6dErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0],
			X_b1_cf*mbcInit.bodyPosW[0])(0)), 0.1);
}


BOOST_AUTO_TEST_CASE(QPBoundedSpeedTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<MultiBody> mbs = {mb};
	std::vector<MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;
	solver.solver("QLD");

	std::string bodyName("b3");
	int bodyIndex = mb.bodyIndexByName(bodyName);
	sva::PTransformd bodyPoint(Vector3d(0., 0.1, 0.));

	qp::BoundedSpeedConstr constSpeed(mbs, 0, 0.005);
	qp::PostureTask postureTask(mbs, 0, {{}, {0.}, {0.}, {0.}}, 1., 0.01);
	qp::PositionTask posTask(mbs, 0, bodyName, Vector3d(1., -1., 1.), bodyPoint.translation());
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 20., 1.);
	MatrixXd dof(1, 6);
	VectorXd speed(1);

	// X body axis must have 0 velocity
	dof << 0.,0.,0.,1.,0.,0.;
	speed << 0.;
	constSpeed.addBoundedSpeed(mbs, bodyName, bodyPoint.translation(), dof, speed);
	BOOST_CHECK_EQUAL(constSpeed.nrBoundedSpeeds(), 1);

	constSpeed.addToSolver(solver);
	solver.addTask(&postureTask);
	solver.addTask(&posTaskSp);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	mbcs[0] = mbcInit;
	for(int i = 0; i < 100; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.005);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}

	sva::PTransformd initPos(bodyPoint*mbcInit.bodyPosW[bodyIndex]);
	sva::PTransformd finalPos(bodyPoint*mbcs[0].bodyPosW[bodyIndex]);
	BOOST_CHECK_SMALL(computeDofError(finalPos, initPos,
		dof), 1e-6);
	BOOST_CHECK_GT(std::pow(
		compute6dError(finalPos, initPos)(4), 2), 0.1);


	// same test but with Z axis
	BOOST_CHECK(constSpeed.removeBoundedSpeed(bodyName));
	constSpeed.updateBoundedSpeeds();
	BOOST_CHECK_EQUAL(constSpeed.nrBoundedSpeeds(), 0);
	BOOST_CHECK_EQUAL(constSpeed.maxGenInEq(), 0);

	// must resize constraint matrix since nrMaxIneq has changed
	solver.updateConstrSize();
	dof << 0.,0.,0.,0.,0.,1.;
	constSpeed.addBoundedSpeed(mbs, bodyName, bodyPoint.translation(), dof, speed);
	constSpeed.updateBoundedSpeeds();
	BOOST_CHECK_EQUAL(constSpeed.nrBoundedSpeeds(), 1);
	BOOST_CHECK_EQUAL(constSpeed.maxGenInEq(), 1);
	solver.updateConstrSize();

	mbcs[0] = mbcInit;
	for(int i = 0; i < 100; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.005);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);
	}

	initPos = bodyPoint*mbcInit.bodyPosW[bodyIndex];
	finalPos = bodyPoint*mbcs[0].bodyPosW[bodyIndex];
	BOOST_CHECK_SMALL(computeDofError(finalPos, initPos,
		dof), 1e-6);
	BOOST_CHECK_GT(std::pow(
		compute6dError(finalPos, initPos)(4), 2), 0.1);
}


BOOST_AUTO_TEST_CASE(MomentumTask)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<MultiBody> mbs = {mb};
	std::vector<MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	sva::ForceVecd momTarget(Vector3d(1., 1., 1.), Vector3d(0., 0., 0.));

	qp::MomentumTask momTask(mbs, 0, momTarget);
	qp::SetPointTask momTaskSp(mbs, 0, &momTask, 10., 1.);

	// Test addTask
	solver.addTask(&momTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test MomentumTask
	BOOST_REQUIRE(solver.solve(mbs, mbcs));
}


BOOST_AUTO_TEST_CASE(SurfaceOrientationTask)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<MultiBody> mbs = {mb};
	std::vector<MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	qp::SurfaceOrientationTask surfOriTask(mbs, 0, "b3",
		RotZ(cst::pi<double>()/2.), sva::PTransformd(Vector3d(0., 0., 0.)));
	qp::SetPointTask surfOriTaskSp(mbs, 0, &surfOriTask, 10., 1.);

	// Test addTask
	solver.addTask(&surfOriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test MomentumTask
	BOOST_REQUIRE(solver.solve(mbs, mbcs));
}


BOOST_AUTO_TEST_CASE(QPCoMPlaneTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<MultiBody> mbs = {mb};
	std::vector<MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	int bodyI = mb.bodyIndexByName("b3");
	Eigen::Vector3d initPos(mbcInit.bodyPosW[bodyI].translation());
	qp::PositionTask posTask(mbs, 0, "b3", initPos);
	qp::SetPointTask posTaskSp(mbs, 0, &posTask, 50., 1.);
	Eigen::Vector3d n1(0., 0., -1.);
	Eigen::Vector3d n2(0., 0., 1.);
	Eigen::Vector3d p1(0., 0., 0.1);
	Eigen::Vector3d p2(0., 0., -0.1);
	double offset1 = -n1.dot(p1);
	double offset2 = -n2.dot(p2);

	qp::CoMIncPlaneConstr comPlaneConstr(mbs, 0, 0.001);
	int planeId1 = 10;
	int planeId2 = 20;
	comPlaneConstr.addPlane(planeId1, n1, offset1, 0.01, 0.005, 0., 0.1);
	comPlaneConstr.addPlane(planeId2, n2, offset2, 0.01, 0.005, 0.1, 0.);
	BOOST_CHECK_EQUAL(comPlaneConstr.nrPlanes(), 2);

	comPlaneConstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	mbcs[0] = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);

		// check if the CoM is on the good side of each plane
		Eigen::Vector3d com = rbd::computeCoM(mbs[0], mbcs[0]);
		double dist1 = n1.dot(com) + offset1;
		double dist2 = n2.dot(com) + offset2;
		BOOST_REQUIRE_GT(dist1, 0.005);
		BOOST_REQUIRE_GT(dist2, 0.005);
	}

	// inverse rotation side
	mbcs[0] = mbcInit;
	posTask.position(initPos);
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(-0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mbs[0], mbcs[0], 0.001);

		forwardKinematics(mbs[0], mbcs[0]);
		forwardVelocity(mbs[0], mbcs[0]);

		// check if the CoM is on the good side of each plane
		Eigen::Vector3d com = rbd::computeCoM(mbs[0], mbcs[0]);
		double dist1 = n1.dot(com) + offset1;
		double dist2 = n2.dot(com) + offset2;
		BOOST_REQUIRE_GT(dist1, 0.005);
		BOOST_REQUIRE_GT(dist2, 0.005);
	}

	comPlaneConstr.rmPlane(planeId1);
	comPlaneConstr.rmPlane(planeId2);
	BOOST_CHECK_EQUAL(comPlaneConstr.nrPlanes(), 0);

	comPlaneConstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);
}


BOOST_AUTO_TEST_CASE(JointsSelector)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm(false);

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	std::vector<MultiBody> mbs = {mb};
	std::vector<MultiBodyConfig> mbcs = {mbcInit};

	qp::QPSolver solver;

	qp::PositionTask pt(mbs, 0, "b3", Vector3d::Zero());
	// construct two JointSelector that should have the same joints
	qp::JointsSelector js1(qp::JointsSelector::ActiveJoints(mbs, 0, &pt, {"j2", "j1"}));
	qp::JointsSelector js2(qp::JointsSelector::UnactiveJoints(mbs, 0, &pt, {"j0", "Root"}));
	qp::SetPointTask js1Sp(mbs, 0, &js1, 1., 1.);
	qp::SetPointTask js2Sp(mbs, 0, &js2, 1., 1.);

	BOOST_REQUIRE_EQUAL(js1.selectedJoints().size(), 2);
	BOOST_REQUIRE_EQUAL(js2.selectedJoints().size(), 2);
	// check they have the same joints selected
	for(std::size_t i = 0; i < 2; ++i)
	{
		BOOST_REQUIRE_EQUAL(js1.selectedJoints()[i].posInDof,
			js2.selectedJoints()[i].posInDof);
		BOOST_REQUIRE_EQUAL(js1.selectedJoints()[i].dof,
			js2.selectedJoints()[i].dof);
	}
	// check joints are sorted
	BOOST_REQUIRE_LT(js1.selectedJoints()[0].posInDof,
		js1.selectedJoints()[1].posInDof);

	solver.addTask(&js1Sp);
	solver.addTask(&js2Sp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 2);

	solver.nrVars(mbs, {}, {});
	solver.updateConstrSize();

	// Test JointsSelector
	BOOST_REQUIRE(solver.solve(mbs, mbcs));

	// jacobian first column should be zero
	BOOST_REQUIRE_EQUAL(js1.jac().block(0, 0, 3, 7), MatrixXd::Zero(3,7));

	// matrix should be equals
	BOOST_REQUIRE_EQUAL(js1.jac(), js2.jac());
	BOOST_REQUIRE_EQUAL(js1.eval(), js2.eval());
	BOOST_REQUIRE_EQUAL(js1.speed(), js2.speed());
	BOOST_REQUIRE_EQUAL(js1.normalAcc(), js2.normalAcc());
}


BOOST_AUTO_TEST_CASE(QPTransformTaskTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeZXZArm(false);

	std::vector<MultiBody> mbs = {mb};
	std::vector<MultiBodyConfig> mbcs(1);

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);


	qp::QPSolver solver;

	solver.nrVars(mbs, {}, {});
	BOOST_CHECK_EQUAL(solver.nrVars(), 9);

	solver.updateConstrSize();


	// Test TransformTask
	sva::PTransformd X_b_s(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());
	sva::PTransformd X_0_t(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());
	Quaterniond E_0_c(Vector4d::Random().normalized());

	qp::TransformTask transTask(mbs, 0, "b3", X_0_t, X_b_s, E_0_c.matrix());
	VectorXd dimW(6);
	dimW << 1., 1., 1., 1., 1., 1.;
	qp::SetPointTask transTaskSp(mbs, 0, &transTask, 10., dimW, 100.);
	// must add postureTask to avoid singularity
	qp::PostureTask postureTask(mbs, 0, {{}, {0.}, {0.}, {0.}}, .5, 10.);

	// addTask all the tasks
	solver.addTask(&transTaskSp);
	solver.addTask(&postureTask);
	solver.updateTasksNrVars(mbs);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 2);

	mbcs[0] = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mb, mbcs[0], 0.001);

		forwardKinematics(mb, mbcs[0]);
		forwardVelocity(mb, mbcs[0]);
	}

	BOOST_CHECK_SMALL(transTask.eval().norm(), 1e-5);

	// removeTask transTaskSp
	solver.removeTask(&transTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	// Test SurfaceTransformTask
	qp::SurfaceTransformTask surfTransTask(mbs, 0, "b3", X_0_t, X_b_s);
	qp::SetPointTask surfTransTaskSp(mbs, 0, &surfTransTask, 10., dimW, 100.);
	solver.addTask(&surfTransTaskSp);
	solver.updateTasksNrVars(mbs);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 2);

	mbcs[0] = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mbs, mbcs));
		eulerIntegration(mb, mbcs[0], 0.001);

		forwardKinematics(mb, mbcs[0]);
		forwardVelocity(mb, mbcs[0]);
	}

	BOOST_CHECK_SMALL(surfTransTask.eval().norm(), 1e-5);

	solver.removeTask(&surfTransTaskSp);
	solver.removeTask(&postureTask);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);
}
