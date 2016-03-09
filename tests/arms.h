// Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
//
// This file is part of Tasks.
//
// RBDyn is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RBDyn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with RBDyn.  If not, see <http://www.gnu.org/licenses/>.

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
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig>
makeZXZArm(bool isFixed=true,
	const sva::PTransformd X_base=sva::PTransformd::Identity())
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

	MultiBody mb = mbg.makeMultiBody(0, isFixed, X_base);

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

	Body b0(rbi, 0, "b0");

	mbg.addBody(b0);

	MultiBody mb = mbg.makeMultiBody(0, true);

	MultiBodyConfig mbc(mb);
	mbc.zero(mb);

	return std::make_tuple(mb, mbc);
}
