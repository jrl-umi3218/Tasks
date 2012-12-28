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
#include <vector>
#include <utility>

// Eigen
#include <Eigen/Core>

// Ipopt
// Remove unused parameter warning from IpTNLP header
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <IpTNLP.hpp>
#pragma GCC diagnostic pop

// RBDyn
#include <MultiBody.h>
#include <MultiBodyConfig.h>



namespace tasks
{

namespace pg
{

class Objective
{
public:
	Objective(double weight):
		weight_(weight)
	{}

	double weight() const
	{
		return weight_;
	}

	void weight(double w)
	{
		weight_ = w;
	}

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc) = 0;

	virtual double value() const = 0;
	virtual const Eigen::VectorXd& gradient() const = 0;

private:
	double weight_;
};



class Constraint
{
public:
	Constraint(int size):
		size_(size)
	{}

	int size() const
	{
		return size_;
	}

	virtual void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc) = 0;

	virtual std::vector<std::pair<int, int>> structure() const = 0;
	virtual const Eigen::VectorXd& value() const = 0;
	virtual const Eigen::VectorXd& jac() const = 0;

	virtual Eigen::VectorXd lower() const = 0;
	virtual Eigen::VectorXd upper() const = 0;

private:
	int size_;
};



class PostureGenerator : Ipopt::TNLP
{
public:
	PostureGenerator(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	void addObjective(Objective* obj);
	void removeObjective(Objective* obj);
	int nrObjectives() const;

	void addConstraint(Constraint* co);
	void removeConstraint(Constraint* co);
	int nrConstraints() const;

	// TNLP overloaded functions
	virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
		Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
		IndexStyleEnum& index_style);

	virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l,
		Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

	virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
		bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
		Ipopt::Index m, bool init_lambda,
		Ipopt::Number* lambda);

	virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x,
		bool new_x, Ipopt::Number& obj_value);
	virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x,
		bool new_x, Ipopt::Number* grad_f);

	virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
		Ipopt::Index m, Ipopt::Number* g);
	virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
		Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
		Ipopt::Index *jCol, Ipopt::Number* values);

private:
	rbd::MultiBody mb_;
	rbd::MultiBodyConfig mbc_;

	std::vector<Objective*> obj_;
	std::vector<Constraint*> constr_;
	std::vector<std::pair<int, int>> struct_;
};

} // namespace pg

} // namespace tasks
