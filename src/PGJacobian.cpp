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

// associated header
#include "PGJacobian.h"

// includes
// RBDyn
#include <Jacobian.h>
#include <MultiBody.h>
#include <MultiBodyConfig.h>

namespace tasks
{

namespace pg
{

Eigen::Matrix<double, 4, 3> angularVelToQuatVel(const std::vector<double>& q)
{
	return 0.5*(Eigen::Matrix<double, 4, 3>() <<
		-q[1], -q[2], -q[3],
		q[0], -q[3], q[2],
		q[3], q[0], -q[1],
		-q[2], q[1], q[0]).finished();
}


Eigen::Matrix<double, 3, 4> quatVelToAngularVel(const std::vector<double>& q)
{
	return 2.*(Eigen::Matrix<double, 3, 4>() <<
		-q[1], q[0], q[3], -q[2],
		-q[2], -q[3], q[0], q[1],
		-q[3], q[2], -q[1], q[0]).finished();
}


PGJacobian::PGJacobian(const rbd::MultiBody& mb, const rbd::Jacobian& jac):
  free_(),
  spherical_(),
  other_(),
  jac_()
{
	int curPos = 0;
	int curPosPG = 0;
	for(int i: jac.jointsPath())
	{
		const rbd::Joint& j = mb.joint(i);
		if(j.type() == rbd::Joint::Free)
		{
			free_.emplace_back(i, curPos, curPosPG);
		}
		else if(j.type() == rbd::Joint::Spherical)
		{
			spherical_.emplace_back(i, curPos, curPosPG);
		}
		else if(j.dof() != 0) // 0 dof joint (fixed) are not considered
		{
			other_.emplace_back(i, curPos, curPosPG);
		}
		curPos += j.dof();
		curPosPG += j.params();
	}

	jac_.resize(6, curPosPG);
}


const Eigen::MatrixXd& PGJacobian::jacobian(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc,
	const Eigen::MatrixXd& jacobian)
{
	for(Joint j: free_)
	{
		Eigen::Matrix3d E = rbd::QuatToE(mbc.q[j.index]);
		Eigen::Matrix<double, 3, 4> W = quatVelToAngularVel(mbc.q[j.index]);

		jac_.block(0, j.posPG, 3, 4) = jacobian.block(0, j.pos, 3, 3)*W;
		jac_.block(3, j.posPG, 3, 4) = jacobian.block(3, j.pos, 3, 3)*W;

		jac_.block(3, j.posPG + 4, 3, 3) =
			jacobian.block(3, j.pos + 3, 3, 3)*E.transpose();
	}

	for(Joint j: spherical_)
	{
		Eigen::Matrix<double, 3, 4> W = quatVelToAngularVel(mbc.q[j.index]);

		jac_.block(0, j.posPG, 3, 4) = jacobian.block(0, j.pos, 3, 3)*W;
		jac_.block(3, j.posPG, 3, 4) = jacobian.block(3, j.pos, 3, 3)*W;
	}

	for(Joint j: other_)
	{
		int nrP = mb.joint(j.index).params();

		jac_.block(0, j.posPG, 6, nrP) = jacobian.block(0, j.pos, 6, nrP);
	}

	return jac_;
}


void PGJacobian::fullJacobian(const rbd::MultiBody& mb,
	const rbd::Jacobian& jac,
	const Eigen::MatrixXd& jacobian, Eigen::MatrixXd& res)
{
	res.setZero();
	int jacPGPos = 0;
	for(int i: jac.jointsPath())
	{
		int params = mb.joint(i).params();
		res.block(0, mb.jointPosInParam(i), res.rows(), params) =
			jacobian.block(0, jacPGPos, res.rows(), params);
		jacPGPos += params;
	}
}


} // namespace pg

} // namespace tasks
