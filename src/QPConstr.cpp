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
#include "QPConstr.h"

// includes
// RBDyn
#include <MultiBody.h>
#include <MultiBodyConfig.h>


namespace tasks
{

namespace qp
{

MotionConstr::MotionConstr(const rbd::MultiBody& mb):
	fd_(mb),
	cont_(),
	fullJac_(6, mb.nrDof()),
	AEq_(),
	BEq_(),
	XL_(),
	XU_(),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0)
{
}


void MotionConstr::updateNrVars(const rbd::MultiBody& mb,
	int alphaD, int lambda, int torque, const std::vector<Contact>& cont)
{
	cont_.resize(cont.size());

	nrDof_ = alphaD;
	nrFor_ = lambda;
	nrTor_ = torque;

	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		cont_[i].jac = rbd::Jacobian(mb, cont[i].bodyId);
		cont_[i].transJac.resize(6, cont_[i].jac.dof());
		cont_[i].points = cont[i].points;
		cont_[i].normals = cont[i].normals;
	}

	AEq_.resize(nrDof_, nrDof_ + nrFor_ + nrTor_);
	BEq_.resize(nrDof_);

	AEq_.setZero();
	BEq_.setZero();

	XL_.resize(nrFor_);
	XU_.resize(nrFor_);

	XL_.fill(0.);
	XU_.fill(std::numeric_limits<double>::infinity());
}


void MotionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	fd_.computeH(mb, mbc);
	fd_.computeC(mb, mbc);

	// H*alphaD - tau - tau_c = -C

	// AEq
	//         nrDof      nrFor            nrTor
	// nrDof [   H      -Sum J_i^t*ni     [0 ... -1]

	AEq_.block(0, 0, nrDof_, nrDof_) = fd_.H();

	int contPos = nrDof_;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		const MatrixXd& jac = cont_[i].jac.jacobian(mb, mbc);

		for(std::size_t j = 0; j < cont_[i].points.size(); ++j)
		{
			cont_[i].jac.translateJacobian(jac, mbc, cont_[i].points[j], cont_[i].transJac);
			cont_[i].jac.fullJacobian(mb, cont_[i].transJac, fullJac_);

			AEq_.block(0, contPos, nrDof_, 1) =
				-fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*cont_[i].normals[j];

			contPos += 1;
		}
	}

	AEq_.block(mb.joint(0).dof(), contPos, nrTor_, nrTor_) =
		-MatrixXd::Identity(nrTor_, nrTor_);


	// BEq = -C
	BEq_ = -fd_.C();
}


int MotionConstr::nrEqLine()
{
	return nrDof_;
}


const Eigen::MatrixXd& MotionConstr::AEq() const
{
	return AEq_;
}


const Eigen::VectorXd& MotionConstr::BEq() const
{
	return BEq_;
}


int MotionConstr::beginVar()
{
	return nrDof_;
}


const Eigen::VectorXd& MotionConstr::Lower() const
{
	return XL_;
}


const Eigen::VectorXd& MotionConstr::Upper() const
{
	return XU_;
}



ContactAccConstr::ContactAccConstr(const rbd::MultiBody& mb):
	cont_(),
	fullJac_(6, mb.nrDof()),
	alphaVec_(mb.nrDof()),
	AEq_(),
	BEq_(),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0)
{}


void ContactAccConstr::updateNrVars(const rbd::MultiBody& mb,
	int alphaD, int lambda, int torque, const std::vector<Contact>& cont)
{
	cont_.resize(cont.size());

	nrDof_ = alphaD;
	nrFor_ = lambda;
	nrTor_ = torque;

	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		cont_[i].jac = rbd::Jacobian(mb, cont[i].bodyId);
	}

	AEq_.resize(cont_.size()*6 , nrDof_ + nrFor_ + nrTor_);
	BEq_.resize(cont_.size()*6);

	AEq_.setZero();
	BEq_.setZero();
}


void ContactAccConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	rbd::paramToVector(mbc.alpha, alphaVec_);

	// J_i*alphaD + JD_i*alpha = 0

	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		// AEq = J_i
		const MatrixXd& jac = cont_[i].jac.jacobian(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jac, fullJac_);
		AEq_.block(i*6, 0, 6, mb.nrDof()) = fullJac_;

		// BEq = -JD_i*alpha
		const MatrixXd& jacDot = cont_[i].jac.jacobianDot(mb, mbc);
		cont_[i].jac.fullJacobian(mb, jacDot, fullJac_);
		BEq_.segment(i*6, 6) = -fullJac_*alphaVec_;
	}
}


int ContactAccConstr::nrEqLine()
{
	return AEq_.rows();
}


const Eigen::MatrixXd& ContactAccConstr::AEq() const
{
	return AEq_;
}


const Eigen::VectorXd& ContactAccConstr::BEq() const
{
	return BEq_;
}

} // namespace qp

} // namespace tasks
