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
	FrictionCone(Eigen::Matrix3d frame, int nrGen, double mu);

	std::vector<Eigen::Vector3d> generators;
};


struct UnilateralContact
{
	UnilateralContact(){}
	UnilateralContact(int bodyId, const std::vector<Eigen::Vector3d>& points,
		Eigen::Matrix3d frame, int nrGen, double mu);

	/// @return Force vector of the contact.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda) const;

	/// @return Number of lambda needed to compute the force vector
	int nrLambda() const;

	/**
		* Safe version of @see force.
		* @throw std::domain_error If lambda don't match the number of generator.
		*/
	Eigen::Vector3d sForce(const Eigen::VectorXd& lambda) const;

	int bodyId;
	std::vector<Eigen::Vector3d> points;
	FrictionCone cone;
};


struct BilateralContact
{
	BilateralContact(){}
	BilateralContact(int bodyId, const std::vector<Eigen::Vector3d>& points,
		Eigen::Matrix3d frame);

	/// @return Force vector of the contact.
	Eigen::Vector3d force(const Eigen::Vector3d& lambda) const;

	/// @return Number of lambda needed to compute the force vector
	int nrLambda() const;

	/**
		* Safe version of @see force. Also the generic one for python bindingx
		* @throw std::domain_error If lambda don't match the number of generator.
		*/
	Eigen::Vector3d sForce(const Eigen::VectorXd& lambda) const;

	int bodyId;
	std::vector<Eigen::Vector3d> points;
	Eigen::Matrix3d frame;
};


} // namespace qp

} // namespace tasks
