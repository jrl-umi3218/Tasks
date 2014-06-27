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
#include "QPMotionConstr.h"

// includes
// Eigen
#include <unsupported/Eigen/Polynomials>

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
	*															MotionConstrCommon
	*/


MotionConstrCommon::ContactData::ContactData(const rbd::MultiBody& mb, int b,
	std::vector<Eigen::Vector3d> p,
	const std::vector<FrictionCone>& cones):
	jac(mb, b),
	body(jac.jointsPath().back()),
	points(std::move(p)),
	generators(cones.size()),
	jacTrans(6, jac.dof()),
	generatorsComp(cones.size())
{
	for(std::size_t i = 0; i < cones.size(); ++i)
	{
		generators[i].resize(3, cones[i].generators.size());
		generatorsComp[i].resize(3, cones[i].generators.size());
		for(std::size_t j = 0; j < cones[i].generators.size(); ++j)
		{
			generators[i].col(j) = cones[i].generators[j];
		}
	}
}



MotionConstrCommon::MotionConstrCommon(const rbd::MultiBody& mb):
	fd_(mb),
	cont_(),
	fullJac_(6, mb.nrDof()),
	A_(),
	AL_(),
	AU_(),
	XL_(),
	XU_(),
	curTorque_(mb.nrDof()),
	nrDof_(0),
	nrFor_(0),
	nrTor_(0)
{
}


void MotionConstrCommon::computeTorque(const Eigen::VectorXd& alphaD,
															 const Eigen::VectorXd& lambda)
{
	curTorque_ = fd_.H()*alphaD;
	curTorque_ += fd_.C();
	curTorque_ += A_.block(0, nrDof_, A_.rows(), A_.cols() - nrDof_)*lambda;
}


const Eigen::VectorXd& MotionConstrCommon::torque() const
{
	return curTorque_;
}


void MotionConstrCommon::torque(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc) const
{
	int pos = mb.joint(0).dof();
	for(std::size_t i = 1; i < mbc.jointTorque.size(); ++i)
	{
		for(double& d: mbc.jointTorque[i])
		{
			d = curTorque_(pos);
			++pos;
		}
	}
}


void MotionConstrCommon::updateNrVars(const rbd::MultiBody& mb,
	const SolverData& data)
{
	const auto& uniCont = data.unilateralContacts();
	const auto& biCont = data.bilateralContacts();
	cont_.resize(data.nrContacts());

	nrDof_ = data.alphaD();
	nrFor_ = data.lambda();
	nrTor_ = data.torque();

	std::size_t iCont = 0;
	for(const UnilateralContact& c: uniCont)
	{
		cont_[iCont] = ContactData(mb, c.bodyId, c.points,
			std::vector<FrictionCone>(c.points.size(), c.cone));

		++iCont;
	}

	for(const BilateralContact& c: biCont)
	{
		cont_[iCont] = ContactData(mb, c.bodyId, c.points, c.cones);

		++iCont;
	}

	A_.setZero(nrDof_, data.nrVars());
	AL_.setZero(nrDof_);
	AU_.setZero(nrDof_);
	curTorque_.resize(nrDof_);

	XL_.resize(data.lambda());
	XU_.resize(data.lambda());

	XL_.fill(0.);
	XU_.fill(std::numeric_limits<double>::infinity());
}


void MotionConstrCommon::computeMatrix(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	fd_.computeH(mb, mbc);
	fd_.computeC(mb, mbc);

	// H*alphaD - tau - tau_c = -C

	// AEq
	//         nrDof      nrFor            nrTor
	// nrDof [   H      -Sum J_i^t*ni     [0 ... -1]

	A_.block(0, 0, nrDof_, nrDof_) = fd_.H();

	int contPos = nrDof_;
	for(std::size_t i = 0; i < cont_.size(); ++i)
	{
		const MatrixXd& jac = cont_[i].jac.jacobian(mb, mbc);

		// for each contact point we compute all the torques
		// due to each generator of the friction cone
		for(std::size_t j = 0; j < cont_[i].points.size(); ++j)
		{
			cont_[i].generatorsComp[j].noalias() =
				mbc.bodyPosW[cont_[i].body].rotation().transpose()*cont_[i].generators[j];

			cont_[i].jac.translateJacobian(jac, mbc,
				cont_[i].points[j], cont_[i].jacTrans);
			cont_[i].jac.fullJacobian(mb, cont_[i].jacTrans, fullJac_);

			A_.block(0, contPos, nrDof_, cont_[i].generatorsComp[j].cols()).noalias() =
				-fullJac_.block(3, 0, 3, fullJac_.cols()).transpose()*
					cont_[i].generatorsComp[j];

			contPos += int(cont_[i].generatorsComp[j].cols());
		}
	}

	// BEq = -C
	AL_ = -fd_.C();
	AU_ = -fd_.C();
}


int MotionConstrCommon::maxInEq() const
{
	return nrDof_;
}


const Eigen::MatrixXd& MotionConstrCommon::AInEq() const
{
	return A_;
}


const Eigen::VectorXd& MotionConstrCommon::LowerInEq() const
{
	return AL_;
}


const Eigen::VectorXd& MotionConstrCommon::UpperInEq() const
{
	return AU_;
}


int MotionConstrCommon::beginVar() const
{
	return nrDof_;
}


const Eigen::VectorXd& MotionConstrCommon::Lower() const
{
	return XL_;
}


const Eigen::VectorXd& MotionConstrCommon::Upper() const
{
	return XU_;
}


std::string MotionConstrCommon::nameInEq() const
{
	return "MotionConstr";
}


std::string MotionConstrCommon::descInEq(const rbd::MultiBody& mb, int line)
{
	int jIndex = findJointFromVector(mb, line, true);
	return std::string("Joint: ") + mb.joint(jIndex).name();
}


std::string MotionConstrCommon::nameBound() const
{
	return "MotionConstr";
}


std::string MotionConstrCommon::descBound(const rbd::MultiBody& mb, int line)
{
	int start = 0;
	int end = 0;
	for(const ContactData& cd: cont_)
	{
		end += int(cd.points.size());
		if(line >= start && line < end)
		{
			return std::string("Body: ") + mb.body(cd.body).name();
		}
		start = end;
	}
	return "";
}


/**
	*															MotionConstr
	*/


MotionConstr::MotionConstr(const rbd::MultiBody& mb,
													std::vector<std::vector<double>> lTorqueBounds,
													std::vector<std::vector<double>> uTorqueBounds):
	MotionConstrCommon(mb),
	torqueL_(),
	torqueU_()
{
	int vars = mb.nrDof() - mb.joint(0).dof();
	torqueL_.resize(vars);
	torqueU_.resize(vars);

	// remove the joint 0
	lTorqueBounds[0] = {};
	uTorqueBounds[0] = {};

	rbd::paramToVector(lTorqueBounds, torqueL_);
	rbd::paramToVector(uTorqueBounds, torqueU_);
}


void MotionConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	const SolverData& /* data */)
{
	computeMatrix(mb, mbc);

	AL_.segment(mb.joint(0).dof(), nrTor_) += torqueL_;
	AU_.segment(mb.joint(0).dof(), nrTor_) += torqueU_;
}


/**
	*															MotionSpringConstr
	*/


MotionSpringConstr::MotionSpringConstr(const rbd::MultiBody& mb,
																		 std::vector<std::vector<double>> lTorqueBounds,
																		 std::vector<std::vector<double>> uTorqueBounds,
																		 const std::vector<SpringJoint>& springs):
	MotionConstrCommon(mb),
	torqueL_(),
	torqueU_(),
	springs_()
{
	int vars = mb.nrDof() - mb.joint(0).dof();
	torqueL_.resize(vars);
	torqueU_.resize(vars);

	// remove the joint 0
	lTorqueBounds[0] = {};
	uTorqueBounds[0] = {};

	rbd::paramToVector(lTorqueBounds, torqueL_);
	rbd::paramToVector(uTorqueBounds, torqueU_);

	springs_.reserve(springs.size());
	for(const SpringJoint& sj: springs)
	{
		int index = mb.jointIndexById(sj.jointId);
		int posInDof = mb.jointPosInDof(index) - mb.joint(0).dof();
		springs_.push_back({index, posInDof, sj.K, sj.C, sj.O});
	}
}


void MotionSpringConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	const SolverData& /* data */)
{
	computeMatrix(mb, mbc);

	for(const SpringJointData& sj: springs_)
	{
		double spring = mbc.q[sj.index][0]*sj.K + mbc.alpha[sj.index][0]*sj.C + sj.O;
		torqueL_(sj.posInDof) = -spring;
		torqueU_(sj.posInDof) = -spring;
	}

	AL_.segment(mb.joint(0).dof(), nrTor_) += torqueL_;
	AU_.segment(mb.joint(0).dof(), nrTor_) += torqueU_;
}


/**
	*															MotionPolyConstr
	*/


MotionPolyConstr::MotionPolyConstr(const rbd::MultiBody& mb,
																 const std::vector<std::vector<Eigen::VectorXd>>& lTorqueBounds,
																 const std::vector<std::vector<Eigen::VectorXd>>& uTorqueBounds):
	MotionConstrCommon(mb),
	torqueL_(),
	torqueU_()
{
	// remove non managed joint
	for(int i = 0; i < mb.nrJoints(); ++i)
	{
		if(mb.joint(i).dof() == 1)
		{
			jointIndex_.push_back(i);
			torqueL_.push_back(lTorqueBounds[i][0]);
			torqueU_.push_back(uTorqueBounds[i][0]);
		}
	}
}


void MotionPolyConstr::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	const SolverData& /* data */)
{
	computeMatrix(mb, mbc);

	for(std::size_t i = 0; i < jointIndex_.size(); ++i)
	{
		int index = jointIndex_[i];
		int dofPos = mb.jointPosInDof(index);
		AL_(dofPos) += Eigen::poly_eval(torqueL_[i], mbc.q[index][0]);
		AU_(dofPos) += Eigen::poly_eval(torqueU_[i], mbc.q[index][0]);
	}
}


} // namespace qp

} // namespace tasks
