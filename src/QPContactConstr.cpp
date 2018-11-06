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

// associated header
#include "Tasks/QPContactConstr.h"

// includes
// std
#include <set>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// Tasks
#include "utils.h"

namespace tasks
{

namespace qp
{

/**
	*															ContactCommon
	*/


bool ContactConstrCommon::ContactCommon::operator==(const ContactCommon& cc) const
{
	return cId == cc.cId;
}


bool ContactConstrCommon::ContactCommon::operator<(const ContactCommon& cc) const
{
	return cId < cc.cId;
}


bool ContactConstrCommon::addVirtualContact(const ContactId& cId)
{
	return virtualContacts_.insert(cId).second;
}


bool ContactConstrCommon::removeVirtualContact(const ContactId& cId)
{
	return virtualContacts_.erase(cId) == 1;
}


void ContactConstrCommon::resetVirtualContacts()
{
	virtualContacts_.clear();
}


bool ContactConstrCommon::addDofContact(const ContactId& cId,
	const Eigen::MatrixXd& dof)
{
	return dofContacts_.insert({cId, dof}).second;
}


bool ContactConstrCommon::removeDofContact(const ContactId& cId)
{
	return dofContacts_.erase(cId) == 1;
}


void ContactConstrCommon::resetDofContacts()
{
	dofContacts_.clear();
}


std::set<ContactConstrCommon::ContactCommon>
ContactConstrCommon::contactCommonInContact(
	const std::vector<rbd::MultiBody>& mbs, const SolverData& data)
{
	std::set<ContactCommon> ret;
	auto isValid = [&mbs, this](const ContactId& contactId)
	{
		// if is virtualContacts we don't add it
		return (virtualContacts_.find(contactId) == virtualContacts_.end());
	};

	for(const BilateralContact& c: data.allContacts())
	{
		if(isValid(c.contactId))
		{
			ret.insert({c.contactId, c.X_b1_cf, c.X_b1_b2});
		}
	}

	return std::move(ret);
}


/**
	*															ContactConstr
	*/


ContactConstr::ContactConstr():
	cont_(),
	fullJac_(),
	dofJac_(),
	A_(),
	b_(),
	nrEq_(0),
	totalAlphaD_(0)
{}


void ContactConstr::updateDofContacts()
{
	for(ContactData& c: cont_)
	{
		auto it = dofContacts_.find(c.contactId);
		if(it != dofContacts_.end())
		{
			c.dof = it->second;
		}
		else
		{
			c.dof.setIdentity(6, 6);
		}
	}
	updateNrEq();
}


void ContactConstr::updateNrVars(const std::vector<rbd::MultiBody>& mbs,
	const SolverData& data)
{
	cont_.clear();
	totalAlphaD_ = data.totalAlphaD();

	int maxDof = std::max_element(mbs.begin(), mbs.end(), compareDof)->nrDof();
	fullJac_.resize(6, maxDof);
	dofJac_.resize(6, maxDof);

	std::set<ContactCommon> contactCSet = contactCommonInContact(mbs, data);
	for(const ContactCommon& cC: contactCSet)
	{
		Eigen::MatrixXd dof(Eigen::MatrixXd::Identity(6, 6));
		auto it = dofContacts_.find(cC.cId);
		if(it != dofContacts_.end())
		{
			dof = it->second;
		}
		std::vector<ContactSideData> contacts;
		auto addContact = [&mbs, &data, &contacts](int rIndex,
			const std::string& bName,
			double sign, const sva::PTransformd& point)
		{
			if(mbs[rIndex].nrDof() > 0)
			{
				contacts.emplace_back(rIndex, data.alphaDBegin(rIndex), sign,
							rbd::Jacobian(mbs[rIndex], bName), point);
			}
		};
		addContact(cC.cId.r1Index, cC.cId.r1BodyName, 1., cC.X_b1_cf);
		addContact(cC.cId.r2Index, cC.cId.r2BodyName, -1., cC.X_b1_cf*cC.X_b1_b2.inv());

		int b1Index = mbs[cC.cId.r1Index].bodyIndexByName(cC.cId.r1BodyName);
		int b2Index = mbs[cC.cId.r2Index].bodyIndexByName(cC.cId.r2BodyName);
		cont_.emplace_back(std::move(contacts), dof, b1Index, b2Index, cC.X_b1_b2,
			cC.X_b1_cf, cC.cId);
	}
	updateNrEq();

	A_.setZero(cont_.size()*6, data.nrVars());
	b_.setZero(cont_.size()*6);
}


int ContactConstr::nrEq() const
{
	return nrEq_;
}


std::string ContactConstr::descEq(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	std::ostringstream oss;
	int contact = line/6;
	for(const ContactSideData& csd: cont_[contact].contacts)
	{
		int body = csd.jac.jointsPath().back();
		oss << "Contact: " << mbs[csd.robotIndex].body(body).name() << std::endl;
	}
	return oss.str();
}


int ContactConstr::maxEq() const
{
	return int(A_.rows());
}


const Eigen::MatrixXd& ContactConstr::AEq() const
{
	return A_;
}


const Eigen::VectorXd& ContactConstr::bEq() const
{
	return b_;
}


void ContactConstr::updateNrEq()
{
	nrEq_ = 0;
	for(const ContactData& c: cont_)
	{
		nrEq_ += int(c.dof.rows());
	}
}


/**
	*															ContactAccConstr
	*/


ContactAccConstr::ContactAccConstr():
	ContactConstr()
{}


void ContactAccConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	using namespace Eigen;

	A_.setZero();
	b_.setZero();
	// J_i*alphaD + JD_i*alpha = 0

	int index = 0;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		ContactData& cd = cont_[i];
		int rows = int(cd.dof.rows());

		for(std::size_t j = 0; j < cd.contacts.size(); ++j)
		{
			ContactSideData& csd = cd.contacts[j];
			const rbd::MultiBody& mb = mbs[csd.robotIndex];
			const rbd::MultiBodyConfig& mbc = mbcs[csd.robotIndex];

			// AEq = J_i
			sva::PTransformd X_0_p = csd.X_b_p*mbc.bodyPosW[csd.bodyIndex];
			const MatrixXd& jacMat = csd.jac.jacobian(mb, mbc, X_0_p);
			dofJac_.block(0, 0, rows, csd.jac.dof()).noalias() =
				csd.sign*cd.dof*jacMat;
			csd.jac.fullJacobian(mb, dofJac_.block(0, 0, rows, csd.jac.dof()),
				fullJac_);
			A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() +=
					fullJac_.block(0, 0, rows, mb.nrDof());

			// BEq = -JD_i*alpha
			Vector6d normalAcc = csd.jac.normalAcceleration(
				mb, mbc, data.normalAccB(csd.robotIndex), csd.X_b_p,
				sva::MotionVecd(Vector6d::Zero())).vector();
			b_.segment(index, rows).noalias() -= csd.sign*cd.dof*normalAcc;
		}
		index += rows;
	}
}


std::string ContactAccConstr::nameEq() const
{
	return "ContactAccConstr";
}


/**
	*															ContactSpeedConstr
	*/


ContactSpeedConstr::ContactSpeedConstr(double timeStep):
	ContactConstr(),
	timeStep_(timeStep)
{}


void ContactSpeedConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	using namespace Eigen;

	A_.block(0, 0, nrEq_, totalAlphaD_).setZero();
	b_.head(nrEq_).setZero();
	// J_i*alphaD + JD_i*alpha = 0

	int index = 0;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		ContactData& cd = cont_[i];
		int rows = int(cd.dof.rows());

		for(std::size_t j = 0; j < cd.contacts.size(); ++j)
		{
			ContactSideData& csd = cd.contacts[j];
			const rbd::MultiBody& mb = mbs[csd.robotIndex];
			const rbd::MultiBodyConfig& mbc = mbcs[csd.robotIndex];

			// AEq = J_i
			sva::PTransformd X_0_p = csd.X_b_p*mbc.bodyPosW[csd.bodyIndex];
			const MatrixXd& jacMat = csd.jac.jacobian(mb, mbc, X_0_p);
			dofJac_.block(0, 0, rows, csd.jac.dof()).noalias() =
				csd.sign*cd.dof*jacMat;
			csd.jac.fullJacobian(mb, dofJac_.block(0, 0, rows, csd.jac.dof()),
				fullJac_);
			A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() +=
					fullJac_.block(0, 0, rows, mb.nrDof());

			// BEq = -JD_i*alpha
			Vector6d normalAcc = csd.jac.normalAcceleration(
				mb, mbc, data.normalAccB(csd.robotIndex), csd.X_b_p,
				sva::MotionVecd(Vector6d::Zero())).vector();
			Vector6d velocity = csd.jac.velocity(mb, mbc, csd.X_b_p).vector();
			b_.segment(index, rows).noalias() -= csd.sign*cd.dof*(normalAcc +
				velocity/timeStep_);
		}

		index += rows;
	}
}


std::string ContactSpeedConstr::nameEq() const
{
	return "ContactSpeedConstr";
}


/**
	*															ContactPosConstr
	*/


ContactPosConstr::ContactPosConstr(double timeStep):
	ContactConstr(),
	timeStep_(timeStep)
{}


void ContactPosConstr::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const SolverData& data)
{
	using namespace Eigen;

	A_.block(0, 0, nrEq_, totalAlphaD_).setZero();
	b_.head(nrEq_).setZero();
	// J_i*alphaD + JD_i*alpha = 0

	int index = 0;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		ContactData& cd = cont_[i];
		int rows = int(cd.dof.rows());

		for(std::size_t j = 0; j < cd.contacts.size(); ++j)
		{
			ContactSideData& csd = cd.contacts[j];
			const rbd::MultiBody& mb = mbs[csd.robotIndex];
			const rbd::MultiBodyConfig& mbc = mbcs[csd.robotIndex];

			// AEq = J_i
			sva::PTransformd X_0_p = csd.X_b_p*mbc.bodyPosW[csd.bodyIndex];
			const MatrixXd& jacMat = csd.jac.jacobian(mb, mbc, X_0_p);
			dofJac_.block(0, 0, rows, csd.jac.dof()).noalias() =
				csd.sign*cd.dof*jacMat;
			csd.jac.fullJacobian(mb, dofJac_.block(0, 0, rows, csd.jac.dof()),
				fullJac_);
			A_.block(index, csd.alphaDBegin, rows, mb.nrDof()).noalias() +=
					fullJac_.block(0, 0, rows, mb.nrDof());

			// BEq = -JD_i*alpha
			Vector6d normalAcc = csd.jac.normalAcceleration(
				mb, mbc, data.normalAccB(csd.robotIndex), csd.X_b_p,
				sva::MotionVecd(Vector6d::Zero())).vector();
			Vector6d velocity = csd.jac.velocity(mb, mbc, csd.X_b_p).vector();
			b_.segment(index, rows).noalias() -= csd.sign*cd.dof*(normalAcc +
				velocity/timeStep_);
		}

		// target the derivative of the position error
		sva::PTransformd X_0_b1cf =
			cd.X_b1_cf*mbcs[cd.contactId.r1Index].bodyPosW[cd.b1Index];
		sva::PTransformd X_0_b2cf =
			cd.X_b1_cf*cd.X_b1_b2.inv()*mbcs[cd.contactId.r2Index].bodyPosW[cd.b2Index];

		sva::PTransformd X_b1cf_b2cf = X_0_b2cf*X_0_b1cf.inv();
		Eigen::Vector6d error;
		error.head<3>() = sva::rotationVelocity(X_b1cf_b2cf.rotation());
		error.tail<3>() = X_b1cf_b2cf.translation();
		b_.segment(index, rows) += cd.dof*(error/timeStep_);

		index += rows;
	}
}


std::string ContactPosConstr::nameEq() const
{
	return "ContactPosConstr";
}

} // qp

} // tasks
