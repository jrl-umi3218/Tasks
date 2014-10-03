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


MotionConstrCommon::ContactData::ContactData(const rbd::MultiBody& mb,
	int bId, int lB,
	std::vector<Eigen::Vector3d> pts,
	const std::vector<FrictionCone>& cones):
	bodyIndex(),
	lambdaBegin(lB),
	jac(mb, bId),
	points(std::move(pts)),
	generators(cones.size()),
	jacTrans(6, jac.dof())//,
//	generatorsComp(cones.size())
{
	bodyIndex = jac.jointsPath().back();
	for(std::size_t i = 0; i < cones.size(); ++i)
	{
		generators[i].resize(3, cones[i].generators.size());
		//generatorsComp[i].resize(3, cones[i].generators.size());
		for(std::size_t j = 0; j < cones[i].generators.size(); ++j)
		{
			generators[i].col(j) = cones[i].generators[j];
		}
	}
}


MotionConstrCommon::RobotData::RobotData(const rbd::MultiBody& mb, int rI, int aDB,
	std::vector<ContactData> contacts):
	robotIndex(rI),
	alphaDBegin(aDB),
	nrDof(mb.nrDof()),
	fd(mb),
	cont(std::move(contacts)),
	fullJac(6, mb.nrDof()),
	curTorque(mb.nrDof())
{ }


MotionConstrCommon::MotionConstrCommon():
	robots_(),
	A_(),
	AL_(),
	AU_(),
	XL_(),
	XU_(),
	rIndexToRobot_()
{
}


void MotionConstrCommon::computeTorque(int robotIndex,
	const Eigen::VectorXd& alphaD, const Eigen::VectorXd& lambda)
{
	int r = rIndexToRobot_.at(robotIndex);
	RobotData& rd = robots_[r];

	rd.curTorque = rd.fd.H()*alphaD;
	rd.curTorque += rd.fd.C();
	// TODO FORCE !!!
	// curTorque_ += A_.block(rd.alphaDBegin, rd.nrDof, A_.rows(), A_.cols() - nrDof_)*lambda;
}


const Eigen::VectorXd& MotionConstrCommon::torque(int robotIndex) const
{
	int r = rIndexToRobot_.at(robotIndex);
	return robots_[r].curTorque;
}


void MotionConstrCommon::torque(const std::vector<rbd::MultiBody>& mbs,
	std::vector<rbd::MultiBodyConfig>& mbcs, int robotIndex) const
{
	int r = rIndexToRobot_.at(robotIndex);
	const RobotData& rd = robots_[r];
	const rbd::MultiBody& mb = mbs[robotIndex];
	rbd::MultiBodyConfig& mbc = mbcs[robotIndex];

	int pos = mb.joint(0).dof();
	for(std::size_t i = 1; i < mbc.jointTorque.size(); ++i)
	{
		for(double& d: mbc.jointTorque[i])
		{
			d = rd.curTorque(pos);
			++pos;
		}
	}
}


void MotionConstrCommon::updateNrVars(const std::vector<rbd::MultiBody>& mbs,
	const SolverData& data)
{
	lambdaBegin_ = data.lambdaBegin();

	const auto& cCont = data.allContacts();

	std::vector<std::vector<int>> rToC(mbs.size());
	for(std::size_t i = 0; i < cCont.size(); ++i)
	{
		const BilateralContact& c = cCont[i];

		// if it's a self contact we can add this constraint juste once
		if(c.contactId.r1Index == c.contactId.r2Index)
		{
			rToC[c.contactId.r1Index].push_back(int(i));
		}
		else
		{
			rToC[c.contactId.r1Index].push_back(int(i));
			rToC[c.contactId.r2Index].push_back(int(i));
		}
	}

	robots_.clear();
	rIndexToRobot_.clear();
	int totalDof = 0;
	for(int r = 0; r < int(mbs.size()); ++r)
	{
		const rbd::MultiBody& mb = mbs[r];

		if(mb.nrDof() > 0)
		{
			const std::vector<int>& cIndex = rToC[r];
			std::vector<ContactData> cd;
			cd.reserve(cIndex.size());

			for(int ci: cIndex)
			{
				const BilateralContact& c = cCont[ci];
				if(r == c.contactId.r1Index)
				{
					cd.emplace_back(mb, c.contactId.r1BodyId, data.lambdaBegin(ci),
						c.r1Points, c.r1Cones);
				}
				// we don't use else to manage self contact on the robot
				if(r == c.contactId.r2Index)
				{
					cd.emplace_back(mb, c.contactId.r2BodyId, data.lambdaBegin(ci),
						c.r2Points, c.r2Cones);
				}
			}

			rIndexToRobot_[r] = int(robots_.size());
			robots_.emplace_back(mb, r, data.alphaD(r), std::move(cd));
			totalDof += mb.nrDof();
		}
	}

	A_.setZero(totalDof, data.nrVars());
	AL_.setZero(totalDof);
	AU_.setZero(totalDof);

	XL_.resize(data.totalLambda());
	XU_.resize(data.totalLambda());

	XL_.fill(0.);
	XU_.fill(std::numeric_limits<double>::infinity());
}


void MotionConstrCommon::computeMatrix(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs)
{
	using namespace Eigen;

	int totalDof = 0;
	for(RobotData& rd: robots_)
	{
		const rbd::MultiBody& mb = mbs[rd.robotIndex];
		const rbd::MultiBodyConfig& mbc = mbcs[rd.robotIndex];

		rd.fd.computeH(mb, mbc);
		rd.fd.computeC(mb, mbc);

		// tauMin -C <= H*alphaD - J^t G lambda <= tauMax - C

		// fill inertia matrix part
		A_.block(totalDof, rd.alphaDBegin, rd.nrDof, rd.nrDof) = rd.fd.H();

		for(std::size_t i = 0; i < rd.cont.size(); ++i)
		{
			const MatrixXd& jac = rd.cont[i].jac.bodyJacobian(mb, mbc);

			ContactData& cd = rd.cont[i];
			for(std::size_t j = 0; j < cd.points.size(); ++j)
			{
				/*
				cd.generatorsComp[j].noalias() =
					mbc.bodyPosW[cd.bodyIndex].rotation().transpose()*cd.generators[j];
					*/

				cd.jac.translateBodyJacobian(jac, mbc, cd.points[j], cd.jacTrans);
				cd.jac.fullJacobian(mb, cd.jacTrans, rd.fullJac);

				A_.block(totalDof, cd.lambdaBegin, rd.nrDof, cd.generators[j].cols()).noalias() =
					-rd.fullJac.block(3, 0, 3, rd.fullJac.cols()).transpose()*
						cd.generators[j];
			}
		}

		// BEq = -C
		AL_ = -rd.fd.C();
		AU_ = -rd.fd.C();

		totalDof += rd.nrDof;
	}
}


int MotionConstrCommon::maxGenInEq() const
{
	return int(A_.rows());
}


const Eigen::MatrixXd& MotionConstrCommon::AGenInEq() const
{
	return A_;
}


const Eigen::VectorXd& MotionConstrCommon::LowerGenInEq() const
{
	return AL_;
}


const Eigen::VectorXd& MotionConstrCommon::UpperGenInEq() const
{
	return AU_;
}


int MotionConstrCommon::beginVar() const
{
	return lambdaBegin_;
}


const Eigen::VectorXd& MotionConstrCommon::Lower() const
{
	return XL_;
}


const Eigen::VectorXd& MotionConstrCommon::Upper() const
{
	return XU_;
}


std::string MotionConstrCommon::nameGenInEq() const
{
	return "MotionConstr";
}


std::string MotionConstrCommon::descGenInEq(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	int totalDof = 0;
	for(const rbd::MultiBody& mb: mbs)
	{
		totalDof += mb.nrDof();
		if(line < totalDof)
		{
			int jIndex = findJointFromVector(mb, line, true);
			return std::string("Joint: ") + mb.joint(jIndex).name();
		}
	}
	return "";
}


std::string MotionConstrCommon::nameBound() const
{
	return "MotionConstr";
}


std::string MotionConstrCommon::descBound(const std::vector<rbd::MultiBody>& mbs,
	int line)
{
	for(const RobotData& rd: robots_)
	{
		for(const ContactData& cd: rd.cont)
		{
			int begin = cd.lambdaBegin - lambdaBegin_;
			int nrLambda = std::accumulate(cd.generators.begin(), cd.generators.end(),
				0, [](int acc, const Eigen::Matrix<double, 3, Eigen::Dynamic>& g)
				{return acc + g.cols();});

			int end = begin + nrLambda;
			if(line >= begin && line < end)
			{
				return std::string("Body: ") +
					mbs[rd.robotIndex].body(cd.bodyIndex).name();
			}
		}
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
