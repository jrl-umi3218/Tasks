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



class PostureGenerator
{
public:
	PostureGenerator(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	void addObjective(Objective* obj);
	void removeObjective(Objective* obj);
	int nrObjectives() const;

	void addConstraint(Constraint* co);
	void removeConstraint(Constraint* co);
	int nrConstraints() const;

	// TNLP functions
	bool get_nlp_info(int& n, int& m,
		int& nnz_jac_g, int& nnz_h_lag);

	bool get_bounds_info(int n, double* x_l,
		double* x_u, int m, double* g_l, double* g_u);

	bool get_starting_point(int n, bool init_x, double* x,
		bool init_z, double* z_L, double* z_U,
		int m, bool init_lambda,
		double* lambda);

	bool eval_f(int n, const double* x,
		bool new_x, double& obj_value);
	bool eval_grad_f(int n, const double* x,
		bool new_x, double* grad_f);

	bool eval_g(int n, const double* x, bool new_x,
		int m, double* g);
	bool eval_jac_g(int n, const double* x, bool new_x,
		int m, int nele_jac, int* iRow,
		int *jCol, double* values);

private:
	rbd::MultiBody mb_;
	rbd::MultiBodyConfig mbc_;

	std::vector<Objective*> obj_;
	std::vector<Constraint*> constr_;
	std::vector<std::pair<int, int>> struct_;
};

} // namespace pg

} // namespace tasks
