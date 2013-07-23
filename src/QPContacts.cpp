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


/// @throw std::domain_error if points is not a valid points index.
void checkRange(int point, const std::vector<Eigen::Vector3d>& points)
{
	if(point < 0 || point >= static_cast<int>(points.size()))
	{
		std::ostringstream str;
		str << "invalid point index: must be in the range [0," << points.size()
				<< "[";
		throw std::domain_error(str.str());
	}
}



/**
	*													FrictionCone
	*/



FrictionCone::FrictionCone(const Eigen::Matrix3d& frame, int nrGen, double mu):
	generators(nrGen)
{
	Eigen::Vector3d normal(frame.row(2));
	Eigen::Vector3d tan(frame.row(0));
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
	const Eigen::Matrix3d& frame, int nrGen, double mu):
	bodyId(bodyId),
	points(points),
	cone(frame, nrGen, mu)
{
}


Eigen::Vector3d UnilateralContact::force(const Eigen::VectorXd& lambda,
	int /* point */) const
{
	Eigen::Vector3d F(Eigen::Vector3d::Zero());

	for(std::size_t i = 0; i < cone.generators.size(); ++i)
	{
		F += cone.generators[i]*lambda(i);
	}

	return F;
}


Eigen::Vector3d UnilateralContact::force(const Eigen::VectorXd& lambda) const
{
	Eigen::Vector3d F(Eigen::Vector3d::Zero());
	int pos = 0;

	for(int i = 0; i < int(points.size()); ++i)
	{
		F += force(lambda.segment(pos, nrLambda(i)), i);
		pos += nrLambda(i);
	}

	return F;
}


int UnilateralContact::nrLambda(int /* point */) const
{
	return static_cast<int>(cone.generators.size());
}


int UnilateralContact::nrLambda() const
{
	int totalLambda = 0;
	for(int i = 0; i < int(points.size()); ++i)
	{
		totalLambda += nrLambda(i);
	}
	return totalLambda;
}


Eigen::Vector3d UnilateralContact::sForce(const Eigen::VectorXd& lambda,
	int point) const
{
	checkRange(point, points);
	if(static_cast<int>(lambda.rows()) != nrLambda(point))
	{
		std::ostringstream str;
		str << "number of lambda and generator mismatch: expected ("
				<< nrLambda(point) << ") gived (" << lambda.rows() << ")";
		throw std::domain_error(str.str());
	}

	return force(lambda, point);
}


Eigen::Vector3d UnilateralContact::sForce(const Eigen::VectorXd& lambda) const
{
	int totalLambda = nrLambda();

	if(static_cast<int>(lambda.rows()) != totalLambda)
	{
		std::ostringstream str;
		str << "number of lambda and generator mismatch: expected ("
				<< totalLambda << ") gived (" << lambda.rows() << ")";
		throw std::domain_error(str.str());
	}

	return force(lambda);
}


int UnilateralContact::sNrLambda(int point) const
{
	checkRange(point, points);
	return nrLambda(point);
}


/**
	*													BilateralContact
	*/


BilateralContact::BilateralContact(int bId,
	const Eigen::Vector3d& center, double radius, int nrPoints,
	const Eigen::Matrix3d& frame, int nrGen, double mu):
	bodyId(bId),
	points(nrPoints),
	cones(nrPoints)
{
	Eigen::Vector3d normal(frame.row(2));
	Eigen::Vector3d tan(frame.row(0));

	double step = (boost::math::constants::pi<double>()*2.)/nrPoints;
	Eigen::Vector3d sPoint = normal*radius;

	for(int i = 0; i < nrPoints; ++i)
	{
		Eigen::AngleAxisd rot(step*i, tan);
		points[i] = center + rot*sPoint;
		// we must inverse rot because we are in anti trig frame
		cones[i] = FrictionCone(rot.inverse()*frame, nrGen, mu);
	}
}


BilateralContact::BilateralContact(int bId, std::vector<Eigen::Vector3d>& p,
	const std::vector<Eigen::Matrix3d>& frames,
	int nrGen, double mu):
	bodyId(bId),
	points(p),
	cones(p.size())
{
	for(std::size_t i = 0; i < frames.size(); ++i)
	{
		cones[i] = FrictionCone(frames[i], nrGen, mu);
	}
}


Eigen::Vector3d BilateralContact::force(const Eigen::VectorXd& lambda,
	int point) const
{
	Eigen::Vector3d F(Eigen::Vector3d::Zero());

	for(std::size_t i = 0; i < cones[point].generators.size(); ++i)
	{
		F += cones[point].generators[i]*lambda(i);
	}

	return F;
}


Eigen::Vector3d BilateralContact::force(const Eigen::VectorXd& lambda) const
{
	Eigen::Vector3d F(Eigen::Vector3d::Zero());
	int pos = 0;

	for(int i = 0; i < int(points.size()); ++i)
	{
		F += force(lambda.segment(pos, nrLambda(i)), i);
		pos += nrLambda(i);
	}

	return F;
}


int BilateralContact::nrLambda(int point) const
{
	return static_cast<int>(cones[point].generators.size());
}


int BilateralContact::nrLambda() const
{
	int totalLambda = 0;
	for(int i = 0; i < int(points.size()); ++i)
	{
		totalLambda += nrLambda(i);
	}
	return totalLambda;
}


Eigen::Vector3d BilateralContact::sForce(const Eigen::VectorXd& lambda,
	int point) const
{
	checkRange(point, points);
	if(static_cast<int>(lambda.rows()) != nrLambda(point))
	{
		std::ostringstream str;
		str << "number of lambda and generator mismatch: expected ("
				<< nrLambda(point) << ") gived (" << lambda.rows() << ")";
		throw std::domain_error(str.str());
	}

	return force(lambda, point);
}


Eigen::Vector3d BilateralContact::sForce(const Eigen::VectorXd& lambda) const
{
	int totalLambda = nrLambda();

	if(static_cast<int>(lambda.rows()) != totalLambda)
	{
		std::ostringstream str;
		str << "number of lambda and generator mismatch: expected ("
				<< totalLambda << ") gived (" << lambda.rows() << ")";
		throw std::domain_error(str.str());
	}

	return force(lambda);
}


int BilateralContact::sNrLambda(int point) const
{
	checkRange(point, points);
	return nrLambda(point);
}


} // namespace qp

} // namespace tasks
