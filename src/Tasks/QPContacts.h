// Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
//
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

#include <tasks/config.hh>

namespace tasks
{

namespace qp
{


/**
	* Linearized friction cone generator.
	* Compute the vector that linearize the friction cone (generatrix).
	*/
struct TASKS_DLLAPI FrictionCone
{
	FrictionCone(){}

	/**
		* @param frame Friction cone frame. The friction cone is define along
		* the frame normal axis (last line).
		* @param nrGen Number of generatrix.
		* @param mu Coefficient of friction.
		* @param dir Cone direction.
		*/
	FrictionCone(const Eigen::Matrix3d& frame, int nrGen, double mu,
		double direction=1.);

	/// Vector of generatrix
	std::vector<Eigen::Vector3d> generators;
};



/**
	* Unique identifier for a contact.
	*/
struct TASKS_DLLAPI ContactId
{
	ContactId();
	/**
		* @param r1Index First robot imply in the contact.
		* @param r2Index Second robot imply in the contact.
		* @param r1BodyName Body name of robot r1Index imply in the contact.
		* @param r2BodyName Body name of robot r2Index imply in the contact.
		* @param ambiguityId If two or more contacts have the same r1Index, r2Index,
		* r1BodyId, r2BodyId the ambiguityId is used to dissociate them.
		*/
	ContactId(int r1Index, int r2Index,
		const std::string& r1BodyName, const std::string& r2BodyName,
		int ambiguityId=-1);

	bool operator==(const ContactId& cId) const;
	bool operator!=(const ContactId& cId) const;

	bool operator<(const ContactId& cId) const;

	int r1Index, r2Index;
	std::string r1BodyName, r2BodyName;
	int ambiguityId;
};



/**
	* Model of a planar contacts.
	* All friction cone are in the same direction.
	*/
struct TASKS_DLLAPI UnilateralContact
{
	UnilateralContact(){}

	/**
		* @param r1Index First robot imply in the contact.
		* @param r2Index Second robot imply in the contact.
		* @param r1BodyName Body name of robot r1Index imply in the contact.
		* @param r2BodyName Body name of robot r2Index imply in the contact.
		* @param r1Points Contact points in r1BodyId frame.
		* @param r1Frame Friction cone frame in rbBodyId frame.
		* @param X_b1_b2 Transformation between r1BodyId and r2BodyId frame.
		* @param nrGen Number of generatrix.
		* @param mu Coefficient of friction.
		* @param X_b1_cf Define the \f$ cf \f$ frame (common frame) used in
		* the contact constraints.
		* @see ContactConstrCommon::addDofContact
		*/
	UnilateralContact(int r1Index, int r2Index,
		const std::string& r1BodyName, const std::string& r2BodyName,
		std::vector<Eigen::Vector3d> r1Points,
		const Eigen::Matrix3d& r1Frame,
		const sva::PTransformd& X_b1_b2,
		int nrGen, double mu,
		const sva::PTransformd& X_b1_cf=sva::PTransformd::Identity());

	/**
		* @param r1Index First robot imply in the contact.
		* @param r2Index Second robot imply in the contact.
		* @param r1BodyName Body name of robot r1Index imply in the contact.
		* @param r2BodyName Body name of robot r2Index imply in the contact.
		* @param ambId ambiguityId.
		* @param r1Points Contact points in r1BodyId frame.
		* @param r1Frame Friction cone frame in rbBodyId frame.
		* @param X_b1_b2 Transformation between r1BodyId and r2BodyId frame.
		* @param nrGen Number of generatrix.
		* @param mu Coefficient of friction.
		* @param X_b1_cf Define the \f$ cf \f$ frame (common frame) used in
		* the contact constraints.
		* @see ContactConstrCommon::addDofContact
		* @see ContactId::ContactId
		*/
	UnilateralContact(int r1Index, int r2Index,
		const std::string& r1BodyName, const std::string& r2BodyName,
		int ambId, std::vector<Eigen::Vector3d> r1Points,
		const Eigen::Matrix3d& r1Frame,
		const sva::PTransformd& X_b1_b2,
		int nrGen, double mu,
		const sva::PTransformd& X_b1_cf=sva::PTransformd::Identity());

	/**
		* @param cId Contact id.
		* @param r1Points Contact points in r1BodyId frame.
		* @param r1Frame Friction cone frame in rbBodyId frame.
		* @param X_b1_b2 Transformation between r1BodyId and r2BodyId frame.
		* @param nrGen Number of generatrix.
		* @param mu Coefficient of friction.
		* @param X_b1_cf Define the \f$ cf \f$ frame (common frame) used in
		* the contact constraints.
		* @see ContactConstrCommon::addDofContact
		*/
	UnilateralContact(const ContactId& cId,
		std::vector<Eigen::Vector3d> r1Points,
		const Eigen::Matrix3d& r1Frame,
		const sva::PTransformd& X_b1_b2,
		int nrGen, double mu,
		const sva::PTransformd& X_b1_cf=sva::PTransformd::Identity());

	/// @return Cone c, point p force vector in body coordinate.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda, int p,
		const FrictionCone& c) const;
	/// @return Cone c, force vector in body coordinate.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda,
		const FrictionCone& c) const;

	/**
	 * Compute force applied on the body origin in the body frame.
	 * @param lambda Linearized contact forces.
	 * @param r_b_pi List of transformation r_b_pi (body origin to point i).
	 * @param c_pi_b Friction cone associated with each point in body frame.
	 * @return F_b, the 6D force applied on the body origin at the body frame.
	 */
	sva::ForceVecd force(const Eigen::VectorXd& lambda,
		const std::vector<Eigen::Vector3d>& r_b_pi,
		const FrictionCone& c_pi_b) const;

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
		* Safe version of @see force.
		* @throw std::domain_error If lambda don't match the number of generator.
		*/
	sva::ForceVecd sForce(const Eigen::VectorXd& lambda,
		const std::vector<Eigen::Vector3d>& r_b_pi,
		const FrictionCone& c_pi_b) const;

	/**
		* Safe version of @see nrLambda.
		* @throw std::domain_error If point is not a valid index.
		*/
	int sNrLambda(int point) const;


	ContactId contactId;
	std::vector<Eigen::Vector3d> r1Points, r2Points;
	FrictionCone r1Cone, r2Cone;
	sva::PTransformd X_b1_b2;
	sva::PTransformd X_b1_cf;

private:
	void construct(const Eigen::MatrixXd& r1Frame, int nrGen, double mu);
};



/**
	* Model of a generic contacts.
	* Friction cone frame can be different for each contact points.
	*/
struct TASKS_DLLAPI BilateralContact
{
	BilateralContact(){}

	/**
		* @param r1Index First robot imply in the contact.
		* @param r2Index Second robot imply in the contact.
		* @param r1BodyName Body name of robot r1Index imply in the contact.
		* @param r2BodyName Body name of robot r2Index imply in the contact.
		* @param r1Points Contact points in r1BodyId frame.
		* @param r1Frames Friction cone frame in r1Points order in rbBodyId frame.
		* @param X_b1_b2 Transformation between r1BodyId and r2BodyId frame.
		* @param nrGen Number of generatrix.
		* @param mu Coefficient of friction.
		* @param X_b1_cf Define the \f$ cf \f$ frame (common frame) used in
		* the contact constraints.
		* @see ContactConstrCommon::addDofContact
		*/
	BilateralContact(int r1Index, int r2Index,
		const std::string& r1BodyName, const std::string& r2BodyName,
		std::vector<Eigen::Vector3d> r1Points,
		const std::vector<Eigen::Matrix3d>& r1Frames,
		const sva::PTransformd& X_b1_b2,
		int nrGen, double mu,
		const sva::PTransformd& X_b1_cf=sva::PTransformd::Identity());

	/**
		* @param r1Index First robot imply in the contact.
		* @param r2Index Second robot imply in the contact.
		* @param r1BodyName Body id of robot r1Index imply in the contact.
		* @param r2BodyName Body id of robot r2Index imply in the contact.
		* @param ambId ambiguityId.
		* @param r1Points Contact points in r1BodyId frame.
		* @param r1Frames Friction cone frame in r1Points order in rbBodyId frame.
		* @param X_b1_b2 Transformation between r1BodyId and r2BodyId frame.
		* @param nrGen Number of generatrix.
		* @param mu Coefficient of friction.
		* @param X_b1_cf Define the \f$ cf \f$ frame (common frame) used in
		* the contact constraints.
		* @see ContactConstrCommon::addDofContact
		* @see ContactId::ContactId
		*/
	BilateralContact(int r1Index, int r2Index,
		const std::string& r1BodyName, const std::string& r2BodyName,
		int ambId, std::vector<Eigen::Vector3d> r1Points,
		const std::vector<Eigen::Matrix3d>& r1Frames,
		const sva::PTransformd& X_b1_b2,
		int nrGen, double mu,
		const sva::PTransformd& X_b1_cf=sva::PTransformd::Identity());

	/**
		* @param cId Contact id.
		* @param r1Points Contact points in r1BodyId frame.
		* @param r1Frames Friction cone frame in r1Points order in rbBodyId frame.
		* @param X_b1_b2 Transformation between r1BodyId and r2BodyId frame.
		* @param nrGen Number of generatrix.
		* @param mu Coefficient of friction.
		* @param X_b1_cf Define the \f$ cf \f$ frame (common frame) used in
		* the contact constraints.
		* @see ContactConstrCommon::addDofContact
		*/
	BilateralContact(const ContactId& cId,
		std::vector<Eigen::Vector3d> r1Points,
		const std::vector<Eigen::Matrix3d>& r1Frames,
		const sva::PTransformd& X_b1_b2,
		int nrGen, double mu,
		const sva::PTransformd& X_b1_cf=sva::PTransformd::Identity());

	/**
		* Construct a BilateralContact from an UnilateralContact.
		*/
	BilateralContact(const UnilateralContact& c);

	/// @return Cone c[point] force vector in body coordinate.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda, int point,
		const std::vector<FrictionCone>& c) const;
	/// @return Cones c force vector in body coordinate.
	Eigen::Vector3d force(const Eigen::VectorXd& lambda,
		const std::vector<FrictionCone>& c) const;

	/**
	 * Compute force applied on the body origin in the body frame.
	 * @param lambda Linearized contact forces.
	 * @param r_b_pi List of transformation r_b_pi (body origin to point i).
	 * @param c_pi_b Frictions cones associated with each point in body frame.
	 * @return F_b, the 6D force applied on the body origin at the body frame.
	 */
	sva::ForceVecd force(const Eigen::VectorXd& lambda,
		const std::vector<Eigen::Vector3d>& r_b_pi,
		const std::vector<FrictionCone>& c_pi_b) const;

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
		* Safe version of @see force.
		* @throw std::domain_error If lambda don't match the number of generator.
		*/
	sva::ForceVecd sForce(const Eigen::VectorXd& lambda,
		const std::vector<Eigen::Vector3d>& r_b_pi,
		const std::vector<FrictionCone>& c_pi_b) const;

	/**
		* Safe version of @see nrLambda.
		* @throw std::domain_error If point is not a valid index.
		*/
	int sNrLambda(int point) const;


	ContactId contactId;
	std::vector<Eigen::Vector3d> r1Points, r2Points;
	std::vector<FrictionCone> r1Cones, r2Cones;
	sva::PTransformd X_b1_b2;
	sva::PTransformd X_b1_cf;

private:
	void construct(const std::vector<Eigen::Matrix3d>& r1Frames, int nrGen, double mu);
};


} // namespace qp

} // namespace tasks
