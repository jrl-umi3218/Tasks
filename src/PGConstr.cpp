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



/**
	*													UnitQuaternion
	*/



UnitQuaternion::UnitQuaternion(const rbd::MultiBody& mb):
  quatPos_(),
  val_(),
  jac_()
{
	for(std::size_t i = 0; i < mb.nrJoints(); ++i)
	{
		if(mb.joint(i).type() == rbd::Joint::Spherical ||
			 mb.joint(i).type() == rbd::Joint::Free)
		{
			quatPos_.push_back(mb.jointPosInParam(i));
		}
	}
	val_.resize(quatPos_.size());
	jac_.resize(quatPos_.size()*4);
	size(quatPos_.size());
}


void UnitQuaternion::update(const rbd::MultiBody& /* mb */,
	const rbd::MultiBodyConfig& /* mbc */,
	int n, const double* x)
{
	using namespace Eigen;
	Map<const VectorXd> qVec(x, n);
	for(std::size_t i = 0; i < quatPos_.size(); ++i)
	{
		int pos = quatPos_[i];
		Quaterniond q(qVec[pos + 0], qVec[pos + 1],
			qVec[pos + 2], qVec[pos + 3]);
		/*
		val_(i) = q.coeffs().squaredNorm();
		jac_.segment(i*4, 4) << q.w()*2., q.x()*2., q.y()*2., q.z()*2.;
		*/
		val_(i) = q.coeffs().norm();
		jac_.segment(i*4, 4) << q.w()/val_(i), q.x()/val_(i), q.y()/val_(i), q.z()/val_(i);
	}
}


std::vector<std::pair<int, int> >
	UnitQuaternion::structure(const rbd::MultiBody& mb) const
{
	std::vector<std::pair<int, int> > s;
	for(std::size_t i = 0; i < quatPos_.size(); ++i)
	{
		int pos = quatPos_[i];
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
	Eigen::VectorXd l(quatPos_.size());
	l.fill(1.);
	return l;
}


Eigen::VectorXd UnitQuaternion::upper() const
{
	Eigen::VectorXd u(quatPos_.size());
	u.fill(1.);
	return u;
}



/**
	*																DummyContact
	*/



DummyContact::DummyContact(const rbd::MultiBody& mb, int bodyId,
  const Eigen::Vector3d& obj):
  Constraint(3),
  obj_(obj),
  body_(mb.bodyIndexById(bodyId)),
  bodyJac_(mb, bodyId),
  bodyPgJac_(mb, bodyJac_),
  val_(3),
  jac_(bodyPgJac_.params()*3)
{
}


void DummyContact::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc,
	int /* n */, const double*  /* x */)
{
	val_ = mbc.bodyPosW[body_].translation() - obj_;
	const Eigen::MatrixXd& jac =
		bodyPgJac_.jacobian(mb, mbc, bodyJac_.jacobian(mb, mbc));
	jac_.segment(0, jac.cols()) = jac.row(3);
	jac_.segment(jac.cols(), jac.cols()) = jac.row(4);
	jac_.segment(jac.cols()*2, jac.cols()) = jac.row(5);
}


std::vector<std::pair<int, int> >
	DummyContact::structure(const rbd::MultiBody& mb) const
{
	std::vector<std::pair<int, int> > s;

	for(int i = 0; i < 3; ++i)
	{
		for(int j: bodyJac_.jointsPath())
		{
			int pos = mb.jointPosInParam(j);
			for(int d = 0; d < mb.joint(j).params(); ++d)
			{
				s.emplace_back(i, pos + d);
			}
		}
	}

	return s;
}


const Eigen::VectorXd& DummyContact::value() const
{
	return val_;
}


const Eigen::VectorXd& DummyContact::jac() const
{
	return jac_;
}


Eigen::VectorXd DummyContact::lower() const
{
	return Eigen::Vector3d::Zero();
}

Eigen::VectorXd DummyContact::upper() const
{
	return Eigen::Vector3d::Zero();
}


} // namespace pg

} // namespace tasks
