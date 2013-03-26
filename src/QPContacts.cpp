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



FrictionCone::FrictionCone(Eigen::Matrix3d frame, int nrGen, double angle):
	generators(nrGen)
{
	Eigen::Vector3d normal = frame.col(2);
	Eigen::Vector3d tan = frame.col(0);

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
	Eigen::Matrix3d frame, int nrGen, double angle):
	bodyId(bodyId),
	points(points),
	cone(frame, nrGen, angle)
{
}

} // namespace qp

} // namespace tasks
