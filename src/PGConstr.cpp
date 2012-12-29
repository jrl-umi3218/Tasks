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
#include "PGConstr.h"


namespace tasks
{

namespace pg
{


UnitQuaternion::UnitQuaternion(const rbd::MultiBody& mb):
  quat_(),
  val_(),
  jac_()
{
	for(std::size_t i = 0; i < mb.nrJoints(); ++i)
	{
		if(mb.joint(i).type() == rbd::Joint::Spherical ||
			 mb.joint(i).type() == rbd::Joint::Free)
		{
			quat_.push_back(i);
		}
	}
	val_.resize(quat_.size());
	jac_.resize(quat_.size()*4);
	size(quat_.size());
}


void UnitQuaternion::update(const rbd::MultiBody& /* mb */,
	const rbd::MultiBodyConfig& mbc)
{
	for(std::size_t i = 0; i < quat_.size(); ++i)
	{
		int j = quat_[i];
		Eigen::Quaterniond q(mbc.q[j][0], mbc.q[j][1], mbc.q[j][2], mbc.q[j][3]);
		val_(i) = q.coeffs().squaredNorm();
		jac_.segment(i*4, 4) << q.w()*2., q.x()*2., q.y()*2., q.z()*2.;
	}
}


std::vector<std::pair<int, int> >
	UnitQuaternion::structure(const rbd::MultiBody& mb) const
{
	std::vector<std::pair<int, int> > s;
	for(std::size_t i = 0; i < quat_.size(); ++i)
	{
		int pos = mb.jointPosInParam(quat_[i]);
		for(int j = 0; j < 4; ++j)
		{
			s.emplace_back(i, pos + j);
		}
	}
	return s;
}


const Eigen::VectorXd& UnitQuaternion::value() const
{
	return val_;
}


const Eigen::VectorXd& UnitQuaternion::jac() const
{
	return jac_;
}


Eigen::VectorXd UnitQuaternion::lower() const
{
	Eigen::VectorXd l(quat_.size());
	l.fill(1.);
	return l;
}


Eigen::VectorXd UnitQuaternion::upper() const
{
	Eigen::VectorXd u(quat_.size());
	u.fill(1.);
	return u;
}

} // namespace pg

} // namespace tasks
