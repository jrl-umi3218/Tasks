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
#include "QPSolver.h"

// includes
// std
#include <limits>
#include <cmath>

//boost
#include <boost/math/constants/constants.hpp>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// Tasks
#include "QL.h"


namespace tasks
{

namespace qp
{


/**
	*													FrictionCone
	*/



FrictionCone::FrictionCone(Eigen::Matrix3d frame, int nrGen, double angle):
	generators(nrGen)
{
	Eigen::Vector3d normal = frame.col(2);
	Eigen::Vector3d tan = frame.col(0);

	Eigen::Vector3d gen = Eigen::AngleAxisd(angle, tan)*normal;
	double step = (boost::math::constants::pi<double>()*2.)/nrGen;

	for(int i = 0; i < nrGen; ++i)
	{
		generators[i] = Eigen::AngleAxisd(step*i, normal)*gen;
	}
}



/**
	*													Contact
	*/



Contact::Contact(int bodyId, const std::vector<Eigen::Vector3d>& points,
  Eigen::Matrix3d frame, int nrGen, double angle):
  bodyId(bodyId),
  point(Eigen::Vector3d::Zero()),
  cone(frame, nrGen, angle)
{
	for(const Eigen::Vector3d& p: points)
		point += p;
	point /= points.size();
}



/**
	*													QPSolver
	*/



QPSolver::QPSolver():
  constr_(),
  eqConstr_(),
  inEqConstr_(),
  boundConstr_(),
  tasks_(),
  alphaD_(0),
  lambda_(0),
  torque_(0),
  nrVars_(0),
  cont_(),
  A1_(),B1_(),A2_(),B2_(),
  XL_(),XU_(),
  Q_(),C_(),
  res_()
{ }


bool QPSolver::update(const rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc)
{
	for(std::size_t i = 0; i < constr_.size(); ++i)
	{
		constr_[i]->update(mb, mbc);
	}

	for(std::size_t i = 0; i < tasks_.size(); ++i)
	{
		tasks_[i]->update(mb, mbc);
	}

	A1_.setZero();
	B1_.setZero();
	A2_.setZero();
	B2_.setZero();
	XL_.fill(-std::numeric_limits<double>::infinity());
	XU_.fill(std::numeric_limits<double>::infinity());
	Q_.setZero();
	C_.setZero();

	int count = 0;
	for(std::size_t i = 0; i < eqConstr_.size(); ++i)
	{
		const Eigen::MatrixXd& A1 = eqConstr_[i]->AEq();
		const Eigen::VectorXd& B1 = eqConstr_[i]->BEq();

		A1_.block(count, 0, A1.rows(), nrVars_) = A1;
		B1_.segment(count, A1.rows()) = B1;

		count += A1.rows();
	}

	count = 0;
	for(std::size_t i = 0; i < inEqConstr_.size(); ++i)
	{
		const Eigen::MatrixXd& A2 = inEqConstr_[i]->AInEq();
		const Eigen::VectorXd& B2 = inEqConstr_[i]->BInEq();

		A2_.block(count, 0, A2.rows(), nrVars_) = A2;
		B2_.segment(count, A2.rows()) = B2;

		count += A2.rows();
	}

	for(std::size_t i = 0; i < boundConstr_.size(); ++i)
	{
		const Eigen::VectorXd& XL = boundConstr_[i]->Lower();
		const Eigen::VectorXd& XU = boundConstr_[i]->Upper();
		int bv = boundConstr_[i]->beginVar();

		XL_.segment(bv, XL.size()) = XL;
		XU_.segment(bv, XU.size()) = XU;
	}

	for(std::size_t i = 0; i < tasks_.size(); ++i)
	{
		const Eigen::MatrixXd& Q = tasks_[i]->Q();
		const Eigen::VectorXd& C = tasks_[i]->C();
		std::pair<int, int> b = tasks_[i]->begin();

		int r = Q.rows();
		int c = Q.cols();

		Q_.block(b.first, b.second, r, c) += tasks_[i]->weight()*Q;
		C_.segment(b.first, r) += tasks_[i]->weight()*C;
	}

	/*
	int tPos = alphaD_ + lambda_;
	Q_.block(tPos, tPos, torque_, torque_).setIdentity();
	*/

	res_.setZero();
	bool success = false;
	double iter = 1e-8;
	while(!success && iter < 1e-3)
	{
		success = solveQP(A1_.cols(), A1_.rows(), A2_.rows(),
			Q_, C_, A1_, B1_, A2_, B2_, XL_, XU_, res_, iter);
		iter *= 10.;
	}

	if(success)
	{
		rbd::vectorToParam(res_.head(alphaD_), mbc.alphaD);
		rbd::vectorToParam(res_.tail(torque_), mbc.jointTorque);

		// don't write contact force to the structure since there are
		// to compute C vector.

		/// @todo Change QPSolver api to allow to write contact force to mbc.
		/*
		int lambdaPos = alphaD_;
		for(std::size_t i = 0; i < cont_.size(); ++i)
		{
			sva::ForceVec GF(Eigen::Vector6d::Zero());

			for(std::size_t j = 0; j < cont_[i].points.size(); ++j)
			{
				double lambdaCoef = res_(lambdaPos);

				sva::PTransform T(cont_[i].points[j]);
				sva::ForceVec F(Eigen::Vector3d::Zero(), cont_[i].normals[j]*lambdaCoef);

				GF = GF + T.transMul(F);
				++lambdaPos;
			}

			mbc.force[mb.bodyIndexById(cont_[i].bodyId)] = GF;
		}
		*/
	}

	return success;
}


void QPSolver::updateEqConstrSize()
{
	int nbEq = 0;
	for(std::size_t i = 0; i < eqConstr_.size(); ++i)
	{
		nbEq += eqConstr_[i]->nrEqLine();
	}

	A1_.resize(nbEq, nrVars_);
	B1_.resize(nbEq);
}


void QPSolver::updateInEqConstrSize()
{
	int nbInEq = 0;
	for(std::size_t i = 0; i < inEqConstr_.size(); ++i)
	{
		nbInEq += inEqConstr_[i]->nrInEqLine();
	}

	A2_.resize(nbInEq, nrVars_);
	B2_.resize(nbInEq);
}


void QPSolver::nrVars(const rbd::MultiBody& mb, std::vector<Contact> cont)
{
	alphaD_ = mb.nrDof();
	lambda_ = 0;
	torque_ = (mb.nrDof() - mb.joint(0).dof());
	cont_ = cont;

	for(const Contact& c: cont_)
	{
		lambda_ += c.cone.generators.size();
	}

	nrVars_ = alphaD_ + lambda_ + torque_;

	if(XL_.rows() != nrVars_)
	{
		XL_.resize(nrVars_);
		XU_.resize(nrVars_);

		Q_.resize(nrVars_, nrVars_);
		C_.resize(nrVars_);

		res_.resize(nrVars_);
	}

	for(Task* t: tasks_)
	{
		t->updateNrVars(mb, alphaD_, lambda_, torque_, cont_);
	}

	for(Constraint* c: constr_)
	{
		c->updateNrVars(mb, alphaD_, lambda_, torque_, cont_);
	}
}


int QPSolver::nrVars() const
{
	return nrVars_;
}


void QPSolver::addEqualityConstraint(EqualityConstraint* co)
{
	eqConstr_.push_back(co);
}


void QPSolver::removeEqualityConstraint(EqualityConstraint* co)
{
	eqConstr_.erase(std::find(eqConstr_.begin(), eqConstr_.end(), co));
}


int QPSolver::nrEqualityConstraints() const
{
	return eqConstr_.size();
}


void QPSolver::addInequalityConstraint(InequalityConstraint* co)
{
	inEqConstr_.push_back(co);
}


void QPSolver::removeInequalityConstraint(InequalityConstraint* co)
{
	inEqConstr_.erase(std::find(inEqConstr_.begin(), inEqConstr_.end(), co));
}


int QPSolver::nrInequalityConstraints() const
{
	return inEqConstr_.size();
}


void QPSolver::addBoundConstraint(BoundConstraint* co)
{
	boundConstr_.push_back(co);
}


void QPSolver::removeBoundConstraint(BoundConstraint* co)
{
	boundConstr_.erase(std::find(boundConstr_.begin(), boundConstr_.end(), co));
}


int QPSolver::nrBoundConstraints() const
{
	return boundConstr_.size();
}


void QPSolver::addConstraint(Constraint* co)
{
	if(std::find(constr_.begin(), constr_.end(), co) == constr_.end())
	{
		constr_.push_back(co);
	}
}


void QPSolver::removeConstraint(Constraint* co)
{
	auto it = std::find(constr_.begin(), constr_.end(), co);
	if(it != constr_.end())
	{
		constr_.erase(it);
	}
}

int QPSolver::nrConstraints() const
{
	return constr_.size();
}


void QPSolver::addTask(Task* task)
{
	tasks_.push_back(task);
}


void QPSolver::removeTask(Task* task)
{
	tasks_.erase(std::find(tasks_.begin(), tasks_.end(), task));
}


int QPSolver::nrTasks() const
{
	return tasks_.size();
}


void QPSolver::resetTasks()
{
	tasks_.clear();
}


const Eigen::VectorXd& QPSolver::result() const
{
	return res_;
}


Eigen::VectorXd QPSolver::alphaDVec() const
{
	return res_.head(alphaD_);
}


Eigen::VectorXd QPSolver::lambdaVec() const
{
	return res_.segment(alphaD_, lambda_);
}


Eigen::VectorXd QPSolver::torqueVec() const
{
	return res_.tail(torque_);
}


} // namespace qp

} // namespace tasks
