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

#pragma once

// include
// std
#include <vector>

// Eigen
#include <Eigen/Core>


namespace tasks
{

namespace qp
{


struct FrictionCone
{
	FrictionCone(){}
	FrictionCone(const Eigen::Matrix3d& frame, int nrGen, double mu);

	std::vector<Eigen::Vector3d> generators;
};


struct UnilateralContact
{
	UnilateralContact(){}
	UnilateralContact(int bodyId, const std::vector<Eigen::Vector3d>& points,
		const Eigen::Matrix3d& frame, int nrGen, double mu);

	/// @return Force vector of the contact point.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda, int point) const;

	/// @return Number of lambda needed to compute the force vector the contact point.
	int nrLambda(int point) const;

	/**
		* Safe version of @see force.
		* @throw std::domain_error If lambda don't match the number of generator
		* or if point is not a valid index.
		*/
	Eigen::Vector3d sForce(const Eigen::VectorXd& lambda, int point) const;

	/**
		* Safe version of @see nrLambda.
		* @throw std::domain_error If point is not a valid index.
		*/
	int sNrLambda(int point) const;


	int bodyId;
	std::vector<Eigen::Vector3d> points;
	FrictionCone cone;
};


struct BilateralContact
{
	BilateralContact(){}
	BilateralContact(int bodyId, const Eigen::Vector3d& center,
		double radius, int nrPoints,
		const Eigen::Matrix3d& frame, int nrGen, double mu);
	BilateralContact(int bodyId, std::vector<Eigen::Vector3d>& points,
		const std::vector<Eigen::Matrix3d>& frames,
		int nrGen, double mu);

	/// @return Force vector of the contact point.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda, int point) const;

	/// @return Number of lambda needed to compute the force vector the contact point.
	int nrLambda(int point) const;

	/**
		* Safe version of @see force. Also the generic one for python bindingx
		* @throw std::domain_error If lambda don't match the number of generator
		* or if point it not a valid index.
		*/
	Eigen::Vector3d sForce(const Eigen::VectorXd& lambda, int point) const;

	/**
		* Safe version of @see nrLambda.
		* @throw std::domain_error If point is not a valid index.
		*/
	int sNrLambda(int point) const;


	int bodyId;
	std::vector<Eigen::Vector3d> points;
	std::vector<FrictionCone> cones;
};


} // namespace qp

} // namespace tasks
