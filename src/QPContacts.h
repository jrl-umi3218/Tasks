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

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>


namespace tasks
{

namespace qp
{


struct FrictionCone
{
	FrictionCone(){}
	FrictionCone(const Eigen::Matrix3d& frame, int nrGen, double mu,
		double direction=1.);

	std::vector<Eigen::Vector3d> generators;
};


struct ContactId
{
	ContactId();
	ContactId(int r1Index, int r2Index, int r1BodyId, int r2BodyId);

	bool operator==(const ContactId& cId) const;
	bool operator!=(const ContactId& cId) const;

	bool operator<(const ContactId& cId) const;

	int r1Index, r2Index;
	int r1BodyId, r2BodyId;
};



struct UnilateralContact
{
	UnilateralContact(){}

	UnilateralContact(int r1Index, int r2Index, int r1BodyId, int r2BodyId,
		std::vector<Eigen::Vector3d> r1Points,
		const Eigen::Matrix3d& r1Frame,
		const sva::PTransformd& X_b1_b2,
        int nrGen, double mu,
        const sva::PTransformd& X_b1_s1 = sva::PTransformd::Identity());

	UnilateralContact(const ContactId& cId,
		std::vector<Eigen::Vector3d> r1Points,
		const Eigen::Matrix3d& r1Frame,
		const sva::PTransformd& X_b1_b2,
        int nrGen, double mu,
        const sva::PTransformd& X_b1_s1 = sva::PTransformd::Identity());

	/// @return Cone c, point p force vector in body coordinate.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda, int p,
		const FrictionCone& c) const;
	/// @return Cone c, force vector in body coordinate.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda,
		const FrictionCone& c) const;

	/// @return Number of lambda needed to compute the force vector of the contact point.
	int nrLambda(int point) const;
	/// @return Number of lambda needed to compute the force vector.
	int nrLambda() const;

	/**
		* Safe version of @see force.
		* @throw std::domain_error If lambda don't match the number of generator
		* or if point is not a valid index.
		*/
	Eigen::Vector3d sForce(const Eigen::VectorXd& lambda, int point,
		const FrictionCone& c) const;
	/**
		* Safe version of @see force.
		* @throw std::domain_error If lambda don't match the number of generator.
		*/
	Eigen::Vector3d sForce(const Eigen::VectorXd& lambda,
		const FrictionCone& c) const;

	/**
		* Safe version of @see nrLambda.
		* @throw std::domain_error If point is not a valid index.
		*/
	int sNrLambda(int point) const;


	ContactId contactId;
	std::vector<Eigen::Vector3d> r1Points, r2Points;
	FrictionCone r1Cone, r2Cone;
	sva::PTransformd X_b1_b2;
    sva::PTransformd X_b1_s1;

private:
	void construct(const Eigen::MatrixXd& r1Frame, int nrGen, double mu);
};


struct BilateralContact
{
	BilateralContact(){}

	BilateralContact(int r1Index, int r2Index,
		int r1BodyId, int r2BodyId,
		std::vector<Eigen::Vector3d> r1Points,
		const std::vector<Eigen::Matrix3d>& r1Frames,
		const sva::PTransformd& X_b1_b2,
        int nrGen, double mu,
        const sva::PTransformd& X_b1_s1 = sva::PTransformd::Identity());

	BilateralContact(const ContactId& cId,
		std::vector<Eigen::Vector3d> r1Points,
		const std::vector<Eigen::Matrix3d>& r1Frames,
		const sva::PTransformd& X_b1_b2,
        int nrGen, double mu,
        const sva::PTransformd& X_b1_s1 = sva::PTransformd::Identity());

	BilateralContact(const UnilateralContact& c);

	/// @return Cone c[point] force vector in body coordinate.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda, int point,
		const std::vector<FrictionCone>& c) const;
	/// @return Cones c force vector in body coordinate.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda,
		const std::vector<FrictionCone>& c) const;

	/// @return Number of lambda needed to compute the force vector of the contact point.
	int nrLambda(int point) const;
	/// @return Number of lambda needed to compute the force vector.
	int nrLambda() const;

	/**
		* Safe version of @see force. Also the generic one for python binding.
		* @throw std::domain_error If lambda don't match the number of generator
		* or if point it not a valid index.
		*/
	Eigen::Vector3d sForce(const Eigen::VectorXd& lambda, int point,
		const std::vector<FrictionCone>& c) const;
	/**
		* Safe version of @see force. Also the generic one for python binding.
		* @throw std::domain_error If lambda don't match the number of generator.
		*/
	Eigen::Vector3d sForce(const Eigen::VectorXd& lambda,
		const std::vector<FrictionCone>& c) const;

	/**
		* Safe version of @see nrLambda.
		* @throw std::domain_error If point is not a valid index.
		*/
	int sNrLambda(int point) const;


	ContactId contactId;
	std::vector<Eigen::Vector3d> r1Points, r2Points;
	std::vector<FrictionCone> r1Cones, r2Cones;
	sva::PTransformd X_b1_b2;
    sva::PTransformd X_b1_s1;

private:
	void construct(const std::vector<Eigen::Matrix3d>& r1Frames, int nrGen, double mu);
};


} // namespace qp

} // namespace tasks
