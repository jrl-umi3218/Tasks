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

class ContactConstrCommon
{
public:
	bool addVirtualContact(const ContactId& contactId);
	bool removeVirtualContact(const ContactId& contactId);
	void resetVirtualContacts();

	bool addDofContact(const ContactId& contactId, const Eigen::MatrixXd& dof);
	bool removeDofContact(const ContactId& contactId);
	void resetDofContacts();

protected:
	struct ContactCommon
	{
		ContactId cId;
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



class ContactConstr : public ConstraintFunction<Equality>,
	public ContactConstrCommon
{
public:
	ContactConstr();

	void updateDofContacts();

	// Constraint
	virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
		const SolverData& data);

	virtual std::string descEq(const std::vector<rbd::MultiBody>& mbs, int line);

	// Inequality Constraint
	virtual int nrEq() const;
	virtual int maxEq() const;

	virtual const Eigen::MatrixXd& AEq() const;
	virtual const Eigen::VectorXd& bEq() const;

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
			const ContactId& cId):
			contacts(std::move(csds)),
			dof(d),
			b1Index(b1),
			b2Index(b2),
			X_b1_b2(X_bb),
			contactId(cId)
		{}

		std::vector<ContactSideData> contacts;
		Eigen::MatrixXd dof;
		int b1Index, b2Index;
		sva::PTransformd X_b1_b2;
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



class ContactAccConstr : public ContactConstr
{
public:
	ContactAccConstr();

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameEq() const;
};



class ContactSpeedConstr : public ContactConstr
{
public:
	ContactSpeedConstr(double timeStep);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameEq() const;
private:
	double timeStep_;
};



class ContactPosConstr : public ContactConstr
{
public:
	ContactPosConstr(double timeStep);

	virtual void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const SolverData& data);

	virtual std::string nameEq() const;
private:
	double timeStep_;
};

} // qp

} // tasks
