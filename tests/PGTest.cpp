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
#define BOOST_TEST_MODULE PGTest
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// eigen3
#include <Eigen/Core>
#include <Eigen/Geometry>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include <EulerIntegration.h>
#include <FK.h>
#include <FV.h>
#include <MultiBody.h>
#include <MultiBodyConfig.h>
#include <MultiBodyGraph.h>
#include <Jacobian.h>

// Tasks
#include "PostureGenerator.h"
#include "PGConstr.h"


const double TOL = 1e-8;


std::tuple<rbd::MultiBody, rbd::MultiBodyConfig> makeFreeLegs()
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;

	MultiBodyGraph mbg;

	double mass = 1.;
	Matrix3d I = Matrix3d::Identity();
	Vector3d h = Vector3d::Zero();

	RBInertia rbi(mass, h, I);

	Body b0(rbi, 0, "waist");
	Body b1(rbi, 1, "rleg");
	Body b2(rbi, 2, "rfoot");
	Body b3(rbi, 3, "lleg");
	Body b4(rbi, 4, "lfoot");

	mbg.addBody(b0);
	mbg.addBody(b1);
	mbg.addBody(b2);
	mbg.addBody(b3);
	mbg.addBody(b4);

	Joint j0(Joint::RevZ, true, 0, "waist_rleg");
	Joint j1(Joint::RevZ, true, 1, "rleg_rfoot");
	Joint j2(Joint::RevZ, true, 2, "waist_lleg");
	Joint j3(Joint::RevZ, true, 3, "llegt_lfoot");

	mbg.addJoint(j0);
	mbg.addJoint(j1);
	mbg.addJoint(j2);
	mbg.addJoint(j3);

	//  Root     j0
	//  ---- b0 ---- b1
	//  Free     X

	PTransform toRLeg(Vector3d(0., 0., 0.1));
	PTransform toLLeg(Vector3d(0., 0., -0.1));
	PTransform toFoot(Vector3d(0., -0.3, 0.));

	mbg.linkBodies(0, toRLeg, 1, PTransform::Identity(), 0);
	mbg.linkBodies(1, toFoot, 2, PTransform::Identity(), 1);
	mbg.linkBodies(0, toLLeg, 3, PTransform::Identity(), 2);
	mbg.linkBodies(3, toFoot, 4, PTransform::Identity(), 3);

	MultiBody mb = mbg.makeMultiBody(0, false);

	MultiBodyConfig mbc(mb);

	mbc.q = {{1., 0., 0., 0., 0., 0., 0.}, {0.}, {0.}, {0.}, {0.}};
	mbc.alpha = {{0., 0., 0., 0., 0., 0.}, {0.}, {0.}, {0.}, {0.}};
	mbc.alphaD = {{0., 0., 0., 0., 0., 0.}, {0.}, {0.}, {0.}, {0.}};
	mbc.jointTorque = {{0., 0., 0., 0., 0., 0.}, {0.}, {0.}, {0.}};
	ForceVec f0(Vector6d::Zero());
	mbc.force = {f0, f0, f0, f0, f0};

	return std::make_tuple(mb, mbc);
}



BOOST_AUTO_TEST_CASE(DummyContact)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;
	using namespace tasks;
	namespace cst = boost::math::constants;


	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeFreeLegs();

	rbd::forwardKinematics(mb, mbc);
	rbd::forwardVelocity(mb, mbc);

	pg::PostureGenerator pg(mb, mbc);

	pg::UnitQuaternion unitQConstr(mb);
	pg::DummyContact contactConstr(mb, 2, Vector3d(1., 1., 1.));

	pg.addConstraint(&unitQConstr);
	pg.addConstraint(&contactConstr);

	BOOST_CHECK(pg.solve());
}
