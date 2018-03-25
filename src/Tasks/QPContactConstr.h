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

// includes
// std
#include <map>
#include <set>

// Eigen
#include <Eigen/Core>

// RBDyn
#include <RBDyn/Jacobian.h>

// Tasks
#include "QPSolver.h"

namespace tasks
{

namespace qp
{


/**
	* Manage modifier on contact constraint.
	*/
class TASKS_DLLAPI ContactConstrCommon
{
public:
	/**
		* Add a virtual contact.
		* A virtual contact don't have any motion constraint apply on it (move freely).
		* It can be seen like a Dof contact with \f$ S \in \mathbb{R}^{0 \times 6} \f$
		*/
	bool addVirtualContact(const ContactId& contactId);

	/// Remove a virtual contact.
	bool removeVirtualContact(const ContactId& contactId);

	/// Remove all virtual contact.
	void resetVirtualContacts();

	/**
		* Free some degree of freedom of a contact like in the BoundedSpeedConstr.
		* \f[ \bar{v} = S v \f]
		* with \f$ v \f$ the velocity in the UnilateralContact and BilateralContact
		* \f$ cf \f$ frame.
		* @param dof \f$ S \in \mathbb{R}^{n \times 6} \f$
		* @see BoundedSpeedConstr
		* @see ContactConstr::updateDofContacts
		*/
	bool addDofContact(const ContactId& contactId, const Eigen::MatrixXd& dof);

	/**
		* Remove a Dof contact.
		* @see ContactConstr::updateDofContacts
		*/
	bool removeDofContact(const ContactId& contactId);

	/**
		* Remove all Dof contact.
		* @see ContactConstr::updateDofContacts
		*/
	void resetDofContacts();

protected:
	struct ContactCommon
	{
		ContactId cId;
		sva::PTransformd X_b1_cf;
		sva::PTransformd X_b1_b2;

		bool operator==(const ContactCommon& cc) const;
		bool operator<(const ContactCommon& cc) const;
	};

protected:
	std::set<ContactCommon> contactCommonInContact(
		const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

protected:
	std::set<ContactId> virtualContacts_;
	std::map<ContactId, Eigen::MatrixXd> dofContacts_;
};


/**
	* Common contact constraint computation.
	*/
class TASKS_DLLAPI ContactConstr : public ConstraintFunction<Equality>,
	public ContactConstrCommon
{
public:
	ContactConstr();

	/**
		* Update \f$ S \f$ matrix based on Dof contact.
		* You must call this method after calling ContactConstrCommon::addDofContact,
		* ContactConstrCommon::removeDofContact
		* and ContactConstrCommon::resetDofContacts.
		*/
	void updateDofContacts();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data) override;

	virtual std::string descEq(const std::vector<rbd::MultiBody>& mbs, int line) override;

	// Inequality Constraint
	virtual int nrEq() const override;
	virtual int maxEq() const override;

	virtual const Eigen::MatrixXd& AEq() const override;
	virtual const Eigen::VectorXd& bEq() const override;

protected:
	struct ContactSideData
	{
		ContactSideData(int rI, int aDB, double s, const rbd::Jacobian& j,
			const sva::PTransformd& Xbp):
			robotIndex(rI), alphaDBegin(aDB), bodyIndex(j.jointsPath().back()),
			sign(s), jac(j), X_b_p(Xbp)
		{}

		int robotIndex, alphaDBegin, bodyIndex;
		double sign;
		rbd::Jacobian jac;
		sva::PTransformd X_b_p;
	};

	struct ContactData
	{
		ContactData(std::vector<ContactSideData> csds,
			const Eigen::MatrixXd& d, int b1, int b2, const sva::PTransformd& X_bb,
			const sva::PTransformd& X_bcf, const ContactId& cId):
			contacts(std::move(csds)),
			dof(d),
			b1Index(b1),
			b2Index(b2),
			X_b1_b2(X_bb),
			X_b1_cf(X_bcf),
			contactId(cId)
		{}

		std::vector<ContactSideData> contacts;
		Eigen::MatrixXd dof;
		int b1Index, b2Index;
		sva::PTransformd X_b1_b2;
		sva::PTransformd X_b1_cf;
		ContactId contactId;
	};

protected:
	void updateNrEq();

protected:
	std::vector<ContactData> cont_;

	Eigen::MatrixXd fullJac_, dofJac_;

	Eigen::MatrixXd A_;
	Eigen::VectorXd b_;

	int nrEq_, totalAlphaD_;
	double timeStep_;
};



/**
	* Contact constraint by targeting a null acceleration.
	* \f[ a = 0 \f]
	* The contact velocity must be null to make this constraint work.
	*
	* This constraint formulation is usually unstable, prefer
	* ContactSpeedConstr and ContactPosConstr.
	*/
class TASKS_DLLAPI ContactAccConstr : public ContactConstr
{
public:
	ContactAccConstr();

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual std::string nameEq() const override;
};



/**
	* Contact constraint by targeting a null velocity.
	* \f[ v + a \Delta_{dt} = 0 \f]
	* This constraint formulation is stable for robot/environment contact.
	* For robot/robot contact prefer ContactPosConstr.
	*/
class TASKS_DLLAPI ContactSpeedConstr : public ContactConstr
{
public:
	/// @param timeStep Time step in second.
	ContactSpeedConstr(double timeStep);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual std::string nameEq() const override;
private:
	double timeStep_;
};



/**
	* Contact constraint by targeting a constant frame.
	* \f[ v + a \Delta_{dt} = \frac{\epsilon(q)}{\Delta_{dt}} \f]
	* Where \f$ \epsilon \f$ is the error between the initial
	* and current contact frame.
	*
	* This constraint formulation is stable in all case.
	*/
class TASKS_DLLAPI ContactPosConstr : public ContactConstr
{
public:
	/// @param timeStep Time step in second.
	ContactPosConstr(double timeStep);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data) override;

	virtual std::string nameEq() const override;
private:
	double timeStep_;
};


} // qp

} // tasks
