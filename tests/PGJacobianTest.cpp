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
#define BOOST_TEST_MODULE PGJacobianTest
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// eigen3
#include <Eigen/Core>
#include <Eigen/Geometry>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include <RBDyn/EulerIntegration.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/MultiBodyGraph.h>
#include <RBDyn/Jacobian.h>

// Tasks
#include "PGJacobian.h"


const double TOL = 1e-8;


std::vector<double> randomFree()
{
	using namespace Eigen;
	Quaterniond q(Vector4d::Random());
	q.normalize();
	Vector3d p(Vector3d::Random());

	return {q.w(), q.x(), q.y(), q.z(), p.x(), p.y(), p.z()};
}


std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeFreeXArm()
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

	mbg.addBody(b0);
	mbg.addBody(b1);

	Joint j0(Joint::RevX, true, 0, "j0");

	mbg.addJoint(j0);

	//  Root     j0
	//  ---- b0 ---- b1
	//  Free     X

	PTransformd to(Vector3d(0., 0.5, 0.));

	mbg.linkBodies(0, to, 1, PTransformd::Identity(), 0);

	MultiBody mb = mbg.makeMultiBody(0, false);

	MultiBodyConfig mbc(mb);

	mbc.q = {{1., 0., 0., 0., 0., 0., 0.}, {0.}};
	mbc.alpha = {{0., 0., 0., 0., 0., 0.}, {0.}};
	mbc.alphaD = {{0., 0., 0., 0., 0., 0.}, {0.}};
	mbc.jointTorque = {{0., 0., 0., 0., 0., 0.}, {0.}, {0.}, {0.}};
	ForceVecd f0(Vector6d::Zero());
	mbc.force = {f0, f0};

	return std::make_tuple(mb, mbc);
}



BOOST_AUTO_TEST_CASE(angularVelToQuatVelVSquatVelToAngularVel)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	for(int i = 0; i < 10; ++i)
	{
		std::vector<double> q = randomFree();
		BOOST_REQUIRE_EQUAL(pg::angularVelToQuatVel(q).transpose()*2.,
			pg::quatVelToAngularVel(q)*0.5);
		BOOST_REQUIRE_EQUAL(pg::angularVelToQuatVel(q)*2.,
			pg::quatVelToAngularVel(q).transpose()*0.5);
	}

	for(int i = 0; i < 10; ++i)
	{
		std::vector<double> q = randomFree();
		Vector3d av = Vector3d::Random();
		Vector4d qv = pg::angularVelToQuatVel(q)*av;
		Vector3d av2 = pg::quatVelToAngularVel(q)*qv;

		BOOST_REQUIRE_SMALL((av - av2).norm(), TOL);
	}

	for(int i = 0; i < 10; ++i)
	{
		std::vector<double> q = randomFree();
		Vector3d av = Vector3d::Random();
		Vector4d qv = pg::angularVelToQuatVel(q)*av;
		Vector3d av2 = pg::quatVelToAngularVel(q)*qv;
		Vector4d qv2 = pg::angularVelToQuatVel(q)*av2;

		BOOST_REQUIRE_SMALL((qv - qv2).norm(), TOL);
	}
}



BOOST_AUTO_TEST_CASE(RotPGJacobianTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeFreeXArm();
	for(int i = 0; i < 10; ++i)
	{
		mbc.q[0] = randomFree();
		mbc.q[1][0] = Matrix<double, 1, 1>::Random()(0);

		forwardKinematics(mb, mbc);
		forwardVelocity(mb, mbc);

		Jacobian jacStd(mb, 1);
		pg::PGJacobian jacPg(mb, jacStd);

		MatrixXd fullStd(6, mb.nrDof());
		MatrixXd fullPg(6, mb.nrParams());

		const MatrixXd& jacRefStd = jacStd.jacobian(mb, mbc);
		jacStd.fullJacobian(mb, jacRefStd, fullStd);

		const MatrixXd& jacRefPg = jacPg.jacobian(mb, mbc, jacRefStd);
		jacPg.fullJacobian(mb, jacStd, jacRefPg, fullPg);

		VectorXd alpha = VectorXd::Random(7);
		VectorXd alphaPg(8);
		alphaPg << pg::angularVelToQuatVel(mbc.q[0])*alpha.head<3>(),
							 alpha.tail<4>();

		VectorXd resStd = fullStd*alpha;
		VectorXd resPg = fullPg*alphaPg;

		BOOST_CHECK_SMALL((resStd.head<3>() - resPg.head<3>()).norm(), TOL);
	}
}



BOOST_AUTO_TEST_CASE(TransPGJacobianTest)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;

	MultiBody mb;
	MultiBodyConfig mbcInit;

	std::tie(mb, mbcInit) = makeFreeXArm();
	for(int i = 0; i < 5; ++i)
	{
		mbcInit.q[0] = randomFree();
		mbcInit.q[1][0] = Matrix<double, 1, 1>::Random()(0);

		forwardKinematics(mb, mbcInit);
		forwardVelocity(mb, mbcInit);

		MultiBodyConfig mbcStd(mbcInit);
		MultiBodyConfig mbcPg(mbcInit);

		Jacobian jacStd(mb, 1);
		pg::PGJacobian jacPg(mb, jacStd);

		MatrixXd fullStd(6, mb.nrDof());
		MatrixXd fullPg(6, mb.nrParams());

		const int nrIt = 5000;
		const double step = 0.005;
		const Vector3d obj(100., 50., 10.);

		for(int i = 0; i < nrIt; ++i)
		{
			Vector3d err = obj - mbcStd.bodyPosW[1].translation();

			const MatrixXd& jac = jacStd.jacobian(mb, mbcStd);
			jacStd.fullJacobian(mb, jac, fullStd);

			VectorXd res =
				fullStd.block(3, 0., 3, mb.nrDof()).jacobiSvd(ComputeThinU | ComputeThinV).solve(err);

			mbcStd.alpha = vectorToDof(mb, res);

			eulerIntegration(mb, mbcStd, step);

			forwardKinematics(mb, mbcStd);
			forwardVelocity(mb, mbcStd);
		}
		BOOST_CHECK_SMALL((mbcStd.bodyPosW[1].translation() - obj).norm(), TOL);

		for(int i = 0; i < nrIt; ++i)
		{
			Vector3d err = obj - mbcPg.bodyPosW[1].translation();

			const MatrixXd& jac = jacPg.jacobian(mb, mbcPg, jacStd.jacobian(mb, mbcPg));
			jacPg.fullJacobian(mb, jacStd, jac, fullPg);

			VectorXd res =
				fullPg.block(3, 0., 3, mb.nrParams()).jacobiSvd(ComputeThinU | ComputeThinV).solve(err);

			VectorXd q = paramToVector(mb, mbcPg.q);
			q += res*step;
			q.head<4>().normalize();
			mbcPg.q = vectorToParam(mb, q);

			forwardKinematics(mb, mbcPg);
			forwardVelocity(mb, mbcPg);
		}
		BOOST_CHECK_SMALL((mbcPg.bodyPosW[1].translation() - obj).norm(), TOL);
	}
}

