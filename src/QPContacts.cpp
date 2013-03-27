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
#include "QPContacts.h"

// includes
// std
#include <stdexcept>

// Eigen
#include <Eigen/Geometry>

//boost
#include <boost/math/constants/constants.hpp>


namespace tasks
{

namespace qp
{

/**
	*													FrictionCone
	*/



FrictionCone::FrictionCone(Eigen::Matrix3d frame, int nrGen, double mu):
	generators(nrGen)
{
	Eigen::Vector3d normal = frame.row(2);
	Eigen::Vector3d tan = frame.row(0);
	double angle = std::atan(mu);

	Eigen::Vector3d gen = Eigen::AngleAxisd(angle, tan)*normal;
	double step = (boost::math::constants::pi<double>()*2.)/nrGen;

	for(int i = 0; i < nrGen; ++i)
	{
		generators[i] = Eigen::AngleAxisd(step*i, normal)*gen;
	}
}



/**
	*													UnilateralContact
	*/



UnilateralContact::UnilateralContact(int bodyId,
	const std::vector<Eigen::Vector3d>& points,
	Eigen::Matrix3d frame, int nrGen, double mu):
	bodyId(bodyId),
	points(points),
	cone(frame, nrGen, mu)
{
}


Eigen::Vector3d UnilateralContact::force(const Eigen::VectorXd& lambda) const
{
	Eigen::Vector3d F(Eigen::Vector3d::Zero());

	for(std::size_t i = 0; i < cone.generators.size(); ++i)
	{
		F += cone.generators[i]*lambda(i);
	}

	return F;
}


Eigen::Vector3d UnilateralContact::sForce(const Eigen::VectorXd& lambda) const
{
	if(static_cast<std::size_t>(lambda.rows()) != cone.generators.size())
	{
		std::ostringstream str;
		str << "number of lambda and generator mismatch: expected ("
				<< cone.generators.size() << ") gived (" << lambda.rows() << ")";
		throw std::domain_error(str.str());
	}

	return force(lambda);
}



/**
	*													BilateralContact
	*/


BilateralContact::BilateralContact(int bId,
	const std::vector<Eigen::Vector3d>& pts,
	Eigen::Matrix3d f):
	bodyId(bId),
	points(pts),
	frame(f)
{
}


Eigen::Vector3d BilateralContact::force(const Eigen::Vector3d& lambda) const
{
	return (frame.array().rowwise()*lambda.transpose().array()).colwise().sum();
}


/**
	*													SlidingContact
	*/


SlidingContact::SlidingContact(int bId,
	const std::vector<Eigen::Vector3d>& pts,
	Eigen::Matrix3d f, double m):
	bodyId(bId),
	points(pts),
	frame(f),
	mu(m)
{
}


Eigen::Vector3d SlidingContact::force(const Eigen::Vector3d& lambda) const
{
	return (frame.array().rowwise()*lambda.transpose().array()).colwise().sum();
}


} // namespace qp

} // namespace tasks
