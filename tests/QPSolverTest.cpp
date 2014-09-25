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
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>

// sch
#include <sch/S_Object/S_Sphere.h>
#include <sch/S_Object/S_Box.h>
#include <sch/CD/CD_Pair.h>

// Tasks
#include "QPConstr.h"
#include "QPMotionConstr.h"
#include "QPSolver.h"
#include "QPTasks.h"

/// @return An simple ZXZ arm with Y as up axis.
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeZXZArm(bool isFixed=true)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;

	MultiBodyGraph mbg;

	double mass = 1.;
	Matrix3d I = Matrix3d::Identity();
	Vector3d h = Vector3d::Zero();

	RBInertiad rbi(mass, h, I);

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


	PTransformd to(Vector3d(0., 0.5, 0.));
	PTransformd from(Vector3d(0., 0., 0.));


	mbg.linkBodies(0, PTransformd::Identity(), 1, from, 0);
	mbg.linkBodies(1, to, 2, from, 1);
	mbg.linkBodies(2, to, 3, from, 2);

	MultiBody mb = mbg.makeMultiBody(0, isFixed);

	MultiBodyConfig mbc(mb);
	mbc.zero(mb);

	return std::make_tuple(mb, mbc);
}



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



BOOST_AUTO_TEST_CASE(BilateralContactTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	const double radius = 0.1;
	double angle = cst::pi<double>()/4.;
	double mu = std::tan(angle);

	Vector3d center(1., 0.5, -0.2);
	center = Vector3d::Zero();

	qp::BilateralContact bi(0, center, radius, 4, Matrix3d::Identity(), 4, mu);
	for(const Vector3d& p: bi.points)
	{
		// check if the point belong to the circle of radius radius
		BOOST_CHECK_SMALL((p - center).squaredNorm() - std::pow(radius, 2), 1e-6);
	}

	std::vector<Vector3d> T = {Vector3d::UnitX(), Vector3d::UnitX(),
															Vector3d::UnitX(), Vector3d::UnitX()};
	std::vector<Vector3d> B = {Vector3d::UnitY(), Vector3d::UnitZ(),
															-Vector3d::UnitY(), -Vector3d::UnitZ()};
	std::vector<Vector3d> N = {Vector3d::UnitZ(), -Vector3d::UnitY(),
															-Vector3d::UnitZ(), Vector3d::UnitY()};

	for(int i = 0; i < 4; ++i)
	{
		// check cone equation T^2 + B^2 = N^2(tan(angle)^2)
		for(const Vector3d& g: bi.cones[i].generators)
		{
			BOOST_CHECK_SMALL(std::pow(T[i].dot(g), 2) + std::pow(B[i].dot(g), 2) -
				std::pow(N[i].dot(g), 2)*std::pow(mu, 2), 1e-6);
		}
	}
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

	solver.nrVars(mb, {}, {});
	BOOST_CHECK_EQUAL(solver.nrVars(), 3);

	solver.updateConstrSize();
	solver.updateConstrSize();


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
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
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
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(oriTask.eval().norm(), 0.00001);

	// test removeTask
	solver.removeTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);


	qp::PostureTask postureTask(mb, {{}, {0.2}, {0.4}, {-0.8}}, 10., 1.);
	postureTask.jointsStiffness(mb, {{2, 10.}});
	solver.addTask(&postureTask);

	// Test PostureTask
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(postureTask.task().eval().norm(), 0.00001);

	solver.removeTask(&postureTask);


	Vector3d comD(RotZ(cst::pi<double>()/4.)*rbd::computeCoM(mb, mbcInit));
	qp::CoMTask comTask(mb, comD);
	qp::SetPointTask comTaskSp(mb, &comTask, 10., 1.);

	solver.addTask(&comTaskSp);

	// Test CoMTask
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(comTask.eval().norm(), 0.00001);

	solver.removeTask(&comTaskSp);


	qp::LinVelocityTask linVelocityTask(mb, 3, -Vector3d::UnitX()*0.005);
	qp::SetPointTask linVelocityTaskSp(mb, &linVelocityTask, 1000., 10000.);

	solver.addTask(&linVelocityTaskSp);

	// Test LinVelocityTask
	mbcSolv = mbcInit;
	for(int i = 0; i < 4000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
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

	MultiBody mb;
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);


	qp::QPSolver solver;


	std::vector<qp::UnilateralContact> contVec =
		{qp::UnilateralContact(3, {Vector3d::Zero()}, Matrix3d::Identity(), 3,
														std::tan(cst::pi<double>()/4.))};

	Vector3d posD = Vector3d(0.707106, 0.707106, 0.);
	qp::PositionTask posTask(mb, 3, posD);
	qp::SetPointTask posTaskSp(mb, &posTask, 10., 1.);

	qp::OrientationTask oriTask(mb, 3, RotZ(cst::pi<double>()/2.));
	qp::SetPointTask oriTaskSp(mb, &oriTask, 10., 1.);


	qp::ContactAccConstr contCstrAcc(mb);

	// Test addEqualityConstraint
	contCstrAcc.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, contVec, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);
	solver.addTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 2);

	solver.data().computeNormalAccB(mb, mbcInit);
	posTask.update(mb, mbcInit, solver.data());
	oriTask.update(mb, mbcInit, solver.data());
	Vector3d evalPos = posTask.eval();
	Vector3d evalOri = oriTask.eval();

	// Test ContactAccConstr
	mbcSolv = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
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
	contCstrAcc.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);



	// Test ContactSpeedConstr
	qp::ContactSpeedConstr contCstrSpeed(mb, 0.001);

	// Test addEqualityConstraint
	contCstrSpeed.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, contVec, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);
	solver.addTask(&oriTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 2);

	solver.data().computeNormalAccB(mb, mbcInit);
	posTask.update(mb, mbcInit, solver.data());
	oriTask.update(mb, mbcInit, solver.data());
	evalPos = posTask.eval();
	evalOri = oriTask.eval();

	mbcSolv = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
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
	solver.removeEqualityConstraint(&contCstrSpeed);
	BOOST_CHECK_EQUAL(solver.nrEqualityConstraints(), 0);
	solver.removeConstraint(&contCstrSpeed);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);



	// MotionConstr test
	contVec = {};
	double Inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double>> torqueMin = {{},{-Inf},{-Inf},{-Inf}};
	std::vector<std::vector<double>> torqueMax = {{},{Inf},{Inf},{Inf}};
	qp::MotionConstr motionCstr(mb, torqueMin, torqueMax);

	motionCstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, contVec, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test MotionConstr
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}
	motionCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
	motionCstr.torque(mb, mbcSolv);

	BOOST_CHECK_SMALL(posTask.eval().norm(), 5e-05);

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
	motionCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);



	// MotionConstr test with contact
	contVec =
		{qp::UnilateralContact(3, {Vector3d::Zero()}, Matrix3d::Identity(), 3,
														std::tan(cst::pi<double>()/4.))};

	solver.addGenInequalityConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 1);
	solver.addBoundConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	solver.addConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, contVec, {});
	solver.updateConstrSize();
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test MotionConstr
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(posTask.eval().norm(), 5e-5);

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test removeEqualityConstraint
	motionCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 0);
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
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	qp::QPSolver solver;

	int bodyI = mb.bodyIndexById(3);
	qp::PositionTask posTask(mb, 3,
		RotZ(cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mb, &posTask, 10., 1.);

	double inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double> > lBound = {{}, {-cst::pi<double>()/4.}, {-inf}, {-inf}};
	std::vector<std::vector<double> > uBound = {{}, {cst::pi<double>()/4.}, {inf}, {inf}};

	qp::JointLimitsConstr jointConstr(mb, lBound, uBound, 0.001);

	// Test add*Constraint
	solver.addBoundConstraint(&jointConstr);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	solver.addConstraint(&jointConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	// Test JointLimitsConstr
	mbcSolv = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
		BOOST_REQUIRE_GT(mbcSolv.q[1][0], -cst::pi<double>()/4. - 0.01);
	}

	posTask.position(RotZ(-cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	mbcSolv = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
		BOOST_REQUIRE_LT(mbcSolv.q[1][0], cst::pi<double>()/4. + 0.01);
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
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	qp::QPSolver solver;

	int bodyI = mb.bodyIndexById(3);
	qp::PositionTask posTask(mb, 3,
		RotZ(cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mb, &posTask, 10., 1.);

	double inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double> > lBound = {{}, {-cst::pi<double>()/4.}, {-inf}, {-inf}};
	std::vector<std::vector<double> > uBound = {{}, {cst::pi<double>()/4.}, {inf}, {inf}};
	std::vector<std::vector<double> > lVel = {{}, {-inf}, {-inf}, {-inf}};
	std::vector<std::vector<double> > uVel = {{}, {inf}, {inf}, {inf}};

	qp::DamperJointLimitsConstr dampJointConstr(mb, lBound, uBound, lVel, uVel,
																						0.125, 0.025, 1., 0.001);

	// Test add*Constraint
	dampJointConstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	// Test JointLimitsConstr
	mbcSolv = mbcInit;
	for(int i = 0; i < 2000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
		BOOST_REQUIRE_GT(mbcSolv.q[1][0], -(cst::pi<double>()/4.)*(1. - 0.025) - 0.01);
	}

	posTask.position(RotZ(-cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	mbcSolv = mbcInit;
	for(int i = 0; i < 2000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
		BOOST_REQUIRE_LT(mbcSolv.q[1][0], (cst::pi<double>()/4.)*(1. - 0.025) + 0.01);
	}

	solver.removeTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 0);

	// Test remove*Constraint
	dampJointConstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);
}



BOOST_AUTO_TEST_CASE(QPTorqueLimitsTest)
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

	int bodyI = mb.bodyIndexById(3);
	qp::PositionTask posTask(mb, 3,
		RotZ(cst::pi<double>()/2.)*mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mb, &posTask, 10., 1.);

	std::vector<std::vector<double> > lBound = {{}, {-30.}, {-30.}, {-30.}};
	std::vector<std::vector<double> > uBound = {{}, {30.}, {30.}, {30.}};

	qp::MotionConstr motionCstr(mb, lBound, uBound);

	// Test add*Constraint

	solver.addGenInequalityConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 1);
	solver.addBoundConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	solver.addConstraint(&motionCstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();


	// Test MotionConstr torque limits
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
		motionCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
		motionCstr.torque(mb, mbcSolv);
		for(int i = 0; i < 3; ++i)
		{
			BOOST_REQUIRE_GT(mbcSolv.jointTorque[i+1][0], lBound[i+1][0] - 0.001);
			BOOST_REQUIRE_LT(mbcSolv.jointTorque[i+1][0], uBound[i+1][0] + 0.001);
		}
	}

	posTask.position(mbcInit.bodyPosW[bodyI].translation());
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
		motionCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
		motionCstr.torque(mb, mbcSolv);
		for(int i = 0; i < 3; ++i)
		{
			BOOST_REQUIRE_GT(mbcSolv.jointTorque[i+1][0], lBound[i+1][0] - 0.001);
			BOOST_REQUIRE_LT(mbcSolv.jointTorque[i+1][0], uBound[i+1][0] + 0.001);
		}
	}

	motionCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 0);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 0);



	Eigen::VectorXd lpoly(2), upoly(2);
	Eigen::VectorXd null;
	lpoly << -30, 1.;
	upoly << 30, 1.;
	std::vector<std::vector<Eigen::VectorXd> > lBoundPoly = {{null}, {lpoly}, {lpoly}, {lpoly}};
	std::vector<std::vector<Eigen::VectorXd> > uBoundPoly = {{null}, {upoly}, {upoly}, {upoly}};
	qp::MotionPolyConstr motionPolyCstr(mb, lBoundPoly, uBoundPoly);

	motionPolyCstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrBoundConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();

	// Test MotionPolyConstr torque limits
	mbcSolv = mbcInit;
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		std::vector<std::vector<double>> oldQ = mbcSolv.q;
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
		motionPolyCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
		motionPolyCstr.torque(mb, mbcSolv);
		for(int i = 0; i < 3; ++i)
		{
			BOOST_REQUIRE_GT(mbcSolv.jointTorque[i+1][0], Eigen::poly_eval(lBoundPoly[i+1][0], oldQ[i+1][0]));
			BOOST_REQUIRE_LT(mbcSolv.jointTorque[i+1][0], Eigen::poly_eval(uBoundPoly[i+1][0], oldQ[i+1][0]));
		}
	}

	posTask.position(mbcInit.bodyPosW[bodyI].translation());
	for(int i = 0; i < 10000; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		std::vector<std::vector<double>> oldQ = mbcSolv.q;
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
		motionPolyCstr.computeTorque(solver.alphaDVec(), solver.lambdaVec());
		motionPolyCstr.torque(mb, mbcSolv);
		for(int i = 0; i < 3; ++i)
		{
			BOOST_REQUIRE_GT(mbcSolv.jointTorque[i+1][0], Eigen::poly_eval(lBoundPoly[i+1][0], oldQ[i+1][0]));
			BOOST_REQUIRE_LT(mbcSolv.jointTorque[i+1][0], Eigen::poly_eval(uBoundPoly[i+1][0], oldQ[i+1][0]));
		}
	}

	motionPolyCstr.removeFromSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrGenInequalityConstraints(), 0);
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
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);


	qp::QPSolver solver;

	int bodyI = mb.bodyIndexById(3);
	qp::PositionTask posTask(mb, 3, mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mb, &posTask, 50., 1.);

	sch::S_Sphere b0(0.25), b3(0.25);
	sch::CD_Pair pair(&b0, &b3);

	PTransformd I = PTransformd::Identity();
	qp::SelfCollisionConstr autoCollConstr(mb, 0.001);
	int collId1 = 10;
	autoCollConstr.addCollision(mb, collId1, 0, &b0, I, 3, &b3, I, 0.01, 0.005, 1.);
	BOOST_CHECK_EQUAL(autoCollConstr.nrCollisions(), 1);

	// Test addInequalityConstraint
	solver.addInequalityConstraint(&autoCollConstr);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 1);
	solver.addConstraint(&autoCollConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	mbcSolv = mbcInit;
	std::ofstream distHard("selfDistHard.py");
	distHard << "dist = [";
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0. ? std::sqrt(dist) : -std::sqrt(-dist);
		distHard << dist << ", ";
		BOOST_REQUIRE_GT(dist, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}
	distHard << "]" << std::endl;

	autoCollConstr.rmCollision(collId1);
	BOOST_CHECK_EQUAL(autoCollConstr.nrCollisions(), 0);


	// test automatic damping computation
	autoCollConstr.addCollision(mb, collId1, 0, &b0, I, 3, &b3, I, 0.1, 0.01, 0., 0.1);
	BOOST_CHECK_EQUAL(autoCollConstr.nrCollisions(), 1);
	posTask.position(mbcInit.bodyPosW[bodyI].translation());

	mbcSolv = mbcInit;
	std::ofstream distSoft("selfDistSoft.py");
	distSoft << "dist = [";
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0. ? std::sqrt(dist) : -std::sqrt(-dist);
		distSoft << dist << ", ";
		BOOST_REQUIRE_GT(dist, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
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

	MultiBody mb;
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm();

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);


	qp::QPSolver solver;

	int bodyI = mb.bodyIndexById(3);
	qp::PositionTask posTask(mb, 3, mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mb, &posTask, 50., 1.);

	sch::S_Sphere b0(0.25), b3(0.25);
	sch::CD_Pair pair(&b0, &b3);

	b0.setTransformation(qp::tosch(mbcInit.bodyPosW[0]));

	PTransformd I = PTransformd::Identity();
	qp::StaticEnvCollisionConstr seCollConstr(mb, 0.001);
	int collId1 = 10;
	seCollConstr.addCollision(mb, collId1, 3, &b3, I, 0, &b0, 0.01, 0.005, 1.);
	BOOST_CHECK_EQUAL(seCollConstr.nrCollisions(), 1);

	// Test addInequalityConstraint
	solver.addInequalityConstraint(&seCollConstr);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 1);
	solver.addConstraint(&seCollConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	mbcSolv = mbcInit;
	std::ofstream distHard("staticDistHard.py");
	distHard << "dist = [";
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
		distHard << dist << ", ";
		BOOST_REQUIRE_GT(dist, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}
	distHard << "]" << std::endl;

	seCollConstr.rmCollision(collId1);
	BOOST_CHECK_EQUAL(seCollConstr.nrCollisions(), 0);

	// test damping computation
	seCollConstr.addCollision(mb, collId1, 3, &b3, I, 0, &b0, 0.1, 0.01, 0., 0.1);
	BOOST_CHECK_EQUAL(seCollConstr.nrCollisions(), 1);
	posTask.position(mbcInit.bodyPosW[bodyI].translation());

	mbcSolv = mbcInit;
	std::ofstream distSoft("staticDistSoft.py");
	distSoft << "dist = [";
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		dist = dist >= 0 ? std::sqrt(dist) : -std::sqrt(-dist);
		distSoft << dist << ", ";
		BOOST_REQUIRE_GT(dist, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
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

BOOST_AUTO_TEST_CASE(QPCoMCollTest)
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

	int bodyI = mb.bodyIndexById(3);
	qp::PositionTask posTask(mb, 3, mbcInit.bodyPosW[bodyI].translation());
	qp::SetPointTask posTaskSp(mb, &posTask, 50., 1.);
	sch::S_Box b0(0.25, 0.25, 0.25);
	sch::S_Sphere comSphere(0.001);
	sch::CD_Pair pair(&comSphere, &b0);

	Eigen::Vector3d com = rbd::computeCoM(mb, mbcInit);
	sva::PTransformd comT(com);

	b0.setTransformation(qp::tosch(comT));

	qp::CoMCollisionConstr comCollConstr(mb, 0.001);
	int collId1 = 10;
	comCollConstr.addCollision(mb, collId1, &b0, 0.01, 0.005, 1.);
	BOOST_CHECK_EQUAL(comCollConstr.nrCollisions(), 1);

	// Test addInequalityConstraint
	solver.addInequalityConstraint(&comCollConstr);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 1);
	solver.addConstraint(&comCollConstr);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);


	mbcSolv = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		//Update comSphere position because comCollConstr attributes are private
		com = rbd::computeCoM(mb, mbcInit);
		comT = sva::PTransformd(com);
		comSphere.setTransformation(qp::tosch(comT));
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		double norm = (com - Vector3d(pb2Tmp[0],pb2Tmp[1],pb2Tmp[2])).norm();
		dist = dist >= 0 ? -norm : norm;
		BOOST_REQUIRE_GT(dist, 0.005);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	comCollConstr.rmCollision(collId1);
	// try the heuristic to compute the damping
	comCollConstr.addCollision(mb, collId1, &b0, 0.01, 0.005, 0., 0.1);
	mbcSolv = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(-0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		sch::Point3 pb1Tmp, pb2Tmp;
		com = rbd::computeCoM(mb, mbcInit);
		comT = sva::PTransformd(com);
		comSphere.setTransformation(qp::tosch(comT));
		double dist = pair.getClosestPoints(pb1Tmp, pb2Tmp);
		double norm = (com - Vector3d(pb2Tmp[0],pb2Tmp[1],pb2Tmp[2])).norm();
		dist = dist >= 0 ? -norm : norm;
		BOOST_REQUIRE_GT(dist, 0.005);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	comCollConstr.rmCollision(collId1);
	BOOST_CHECK_EQUAL(comCollConstr.nrCollisions(), 0);
}


BOOST_AUTO_TEST_CASE(QPBilatContactTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm(false);

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);


	qp::QPSolver solver;

	double Inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double>> torqueMin = {{},{-Inf},{-Inf},{-Inf}};
	std::vector<std::vector<double>> torqueMax = {{},{Inf},{Inf},{Inf}};
	qp::MotionConstr motionCstr(mb, torqueMin, torqueMax);
	qp::ContactAccConstr contCstrAcc(mb);

	solver.addGenInequalityConstraint(&motionCstr);
	solver.addBoundConstraint(&motionCstr);
	solver.addConstraint(&motionCstr);

	solver.addEqualityConstraint(&contCstrAcc);
	solver.addConstraint(&contCstrAcc);

	std::vector<Eigen::Vector3d> points =
		{
			Vector3d(0.1, 0.1, 0.),
			 Vector3d(-0.1, 0.1, 0.),
			Vector3d(-0.1, -0.1, 0.),
			Vector3d(0.1, -0.1, 0.)
		};

	std::vector<qp::UnilateralContact> uni =
		{qp::UnilateralContact(0, points, Matrix3d::Identity(), 3, 0.7)};
	std::vector<qp::BilateralContact> bi =
		{qp::BilateralContact(0, Vector3d::Zero(), 0.1, 4, Matrix3d::Identity(), 3., 0.7)};

	solver.nrVars(mb, uni, {});
	solver.updateConstrSize();

	// This stance with unilateral contac is impossible so the solver must fail
	mbcSolv = mbcInit;
	BOOST_REQUIRE(!solver.solve(mb, mbcSolv));


	// We test it again with bilateral contact to check that the stance is now
	// valid
	solver.nrVars(mb, {}, bi);
	solver.updateConstrSize();

	mbcSolv = mbcInit;
	for(int i = 0; i < 10; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}


	solver.removeEqualityConstraint(&contCstrAcc);
	solver.removeConstraint(&contCstrAcc);

	solver.removeConstraint(&motionCstr);
	solver.removeBoundConstraint(&motionCstr);
	solver.removeGenInequalityConstraint(&motionCstr);
}


Eigen::Vector6d compute6dError(const sva::PTransformd& b1, const sva::PTransformd& b2)
{
	Eigen::Vector6d error;
	error.head<3>() = sva::rotationError(b1.rotation(), b2.rotation(), 1e-7);
	error.tail<3>() = b1.translation() - b2.translation();
	return error;
}


double computeDofError(const sva::PTransformd& b1, const sva::PTransformd& b2,
										const Eigen::MatrixXd& dof)
{
	return (dof*compute6dError(b1, b2)).norm();
}


BOOST_AUTO_TEST_CASE(QPDofContactsTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit, mbcSolv;

	std::tie(mb, mbcInit) = makeZXZArm(false);

	forwardKinematics(mb, mbcInit);
	forwardVelocity(mb, mbcInit);

	qp::QPSolver solver;

	qp::ContactSpeedConstr contCstrSpeed(mb, 0.005);
	qp::PositionTask posTask(mb, 0, Vector3d(1., 1., -1.));
	qp::SetPointTask posTaskSp(mb, &posTask, 10., 1.);

	contCstrSpeed.addToSolver(solver);
	solver.addTask(&posTaskSp);

	std::vector<Eigen::Vector3d> points =
		{
			Vector3d(0.1, 0.1, 0.),
			 Vector3d(-0.1, 0.1, 0.),
			Vector3d(-0.1, -0.1, 0.),
			Vector3d(0.1, -0.1, 0.)
		};

	std::vector<qp::UnilateralContact> uni =
		{qp::UnilateralContact(0, points, Matrix3d::Identity(), 3, 0.7)};

	MatrixXd contactDof(5, 6);
	contactDof.setZero();
	contactDof(0, 0) = 1.;
	contactDof(1, 1) = 1.;
	contactDof(2, 2) = 1.;
	contactDof(3, 3) = 1.;
	contactDof(4, 4) = 1.;

	// test Z free
	contCstrSpeed.addDofContact(0, contactDof);
	solver.nrVars(mb, uni, {});
	solver.updateConstrSize();

	mbcSolv = mbcInit;
	for(int i = 0; i < 100; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.005);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(computeDofError(mbcSolv.bodyPosW[0], mbcInit.bodyPosW[0],
		contactDof), 1e-6);
	BOOST_CHECK_GT(std::pow(
		compute6dError(mbcSolv.bodyPosW[0], mbcInit.bodyPosW[0])(5), 2), 0.1);


	// test Y Free and updateDofContacts
	contactDof.setZero();
	contactDof(0, 0) = 1.;
	contactDof(1, 1) = 1.;
	contactDof(2, 2) = 1.;
	contactDof(3, 3) = 1.;
	contactDof(4, 5) = 1.;
	contCstrSpeed.resetDofContacts();
	contCstrSpeed.addDofContact(0, contactDof);
	contCstrSpeed.updateDofContacts();
	mbcSolv = mbcInit;
	for(int i = 0; i < 100; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.005);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	BOOST_CHECK_SMALL(computeDofError(mbcSolv.bodyPosW[0], mbcInit.bodyPosW[0],
		contactDof), 1e-6);
	BOOST_CHECK_GT(std::pow(
		compute6dError(mbcSolv.bodyPosW[0], mbcInit.bodyPosW[0])(4), 2), 0.1);
}


BOOST_AUTO_TEST_CASE(QPConstantSpeedTest)
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
	solver.solver("QLD");

	int bodyId = 3;
	int bodyIndex = mb.bodyIndexById(bodyId);
	sva::PTransformd bodyPoint(Vector3d(0., 0.1, 0.));

	qp::ConstantSpeedConstr constSpeed(mb, 0.005);
	qp::PostureTask postureTask(mb, {{}, {0.}, {0.}, {0.}}, 1., 0.01);
	qp::PositionTask posTask(mb, bodyId, Vector3d(1., -1., 1.), bodyPoint.translation());
	qp::SetPointTask posTaskSp(mb, &posTask, 20., 1.);
	MatrixXd dof(1, 6);
	VectorXd speed(1);

	// X body axis must have 0 velocity
	dof << 0.,0.,0.,1.,0.,0.;
	speed << 0.;
	constSpeed.addConstantSpeed(mb, bodyId, bodyPoint.translation(), dof, speed);
	BOOST_CHECK_EQUAL(constSpeed.nrConstantSpeed(), 1);

	constSpeed.addToSolver(solver);
	solver.addTask(&postureTask);
	solver.addTask(&posTaskSp);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();

	mbcSolv = mbcInit;
	for(int i = 0; i < 100; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.005);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	sva::PTransformd initPos(bodyPoint*mbcInit.bodyPosW[bodyIndex]);
	sva::PTransformd finalPos(bodyPoint*mbcSolv.bodyPosW[bodyIndex]);
	BOOST_CHECK_SMALL(computeDofError(finalPos, initPos,
		dof), 1e-6);
	BOOST_CHECK_GT(std::pow(
		compute6dError(finalPos, initPos)(4), 2), 0.1);


	// same test but with Z axis
	BOOST_CHECK(constSpeed.removeConstantSpeed(bodyId));
	BOOST_CHECK_EQUAL(constSpeed.nrConstantSpeed(), 0);
	BOOST_CHECK_EQUAL(constSpeed.maxEq(), 0);
	// must resize constraint matrix since nrMaxIneq has changed
	solver.updateConstrSize();
	dof << 0.,0.,0.,0.,0.,1.;
	constSpeed.addConstantSpeed(mb, bodyId, bodyPoint.translation(), dof, speed);
	BOOST_CHECK_EQUAL(constSpeed.nrConstantSpeed(), 1);
	BOOST_CHECK_EQUAL(constSpeed.maxEq(), 1);
	solver.updateConstrSize();

	mbcSolv = mbcInit;
	for(int i = 0; i < 100; ++i)
	{
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.005);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);
	}

	initPos = bodyPoint*mbcInit.bodyPosW[bodyIndex];
	finalPos = bodyPoint*mbcSolv.bodyPosW[bodyIndex];
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
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	forwardKinematics(mb, mbc);
	forwardVelocity(mb, mbc);

	qp::QPSolver solver;

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();

	sva::ForceVecd momTarget(Vector3d(1., 1., 1.), Vector3d(0., 0., 0.));

	qp::MomentumTask momTask(mb, momTarget);
	qp::SetPointTask momTaskSp(mb, &momTask, 10., 1.);

	// Test addTask
	solver.addTask(&momTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	// Test MomentumTask
	BOOST_REQUIRE(solver.solve(mb, mbc));
}


BOOST_AUTO_TEST_CASE(QPCoMPlaneTest)
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

	int bodyI = mb.bodyIndexById(3);
	Eigen::Vector3d initPos(mbcInit.bodyPosW[bodyI].translation());
	qp::PositionTask posTask(mb, 3, initPos);
	qp::SetPointTask posTaskSp(mb, &posTask, 50., 1.);
	Eigen::Vector3d n1(0., 0., -1.);
	Eigen::Vector3d n2(0., 0., 1.);
	Eigen::Vector3d p1(0., 0., 0.1);
	Eigen::Vector3d p2(0., 0., -0.1);
	double offset1 = -n1.dot(p1);
	double offset2 = -n2.dot(p2);

	qp::CoMIncPlaneConstr comPlaneConstr(mb, 0.001);
	int planeId1 = 10;
	int planeId2 = 20;
	comPlaneConstr.addPlane(planeId1, n1, offset1, 0.01, 0.005, 0., 0.1);
	comPlaneConstr.addPlane(planeId2, n2, offset2, 0.01, 0.005, 0.1, 0.);
	BOOST_CHECK_EQUAL(comPlaneConstr.nrPlanes(), 2);

	comPlaneConstr.addToSolver(solver);
	BOOST_CHECK_EQUAL(solver.nrInequalityConstraints(), 1);
	BOOST_CHECK_EQUAL(solver.nrConstraints(), 1);

	solver.nrVars(mb, {}, {});
	solver.updateConstrSize();

	solver.addTask(&posTaskSp);
	BOOST_CHECK_EQUAL(solver.nrTasks(), 1);

	mbcSolv = mbcInit;
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);

		// check if the CoM is on the good side of each plane
		Eigen::Vector3d com = rbd::computeCoM(mb, mbcSolv);
		double dist1 = n1.dot(com) + offset1;
		double dist2 = n2.dot(com) + offset2;
		BOOST_REQUIRE_GT(dist1, 0.005);
		BOOST_REQUIRE_GT(dist2, 0.005);
	}

	// inverse rotation side
	mbcSolv = mbcInit;
	posTask.position(initPos);
	for(int i = 0; i < 1000; ++i)
	{
		posTask.position(RotX(-0.01)*posTask.position());
		BOOST_REQUIRE(solver.solve(mb, mbcSolv));
		eulerIntegration(mb, mbcSolv, 0.001);

		forwardKinematics(mb, mbcSolv);
		forwardVelocity(mb, mbcSolv);

		// check if the CoM is on the good side of each plane
		Eigen::Vector3d com = rbd::computeCoM(mb, mbcSolv);
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
