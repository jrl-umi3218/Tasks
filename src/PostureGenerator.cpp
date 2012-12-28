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
#include "PostureGenerator.h"

// includes
// RBDyn
#include <FK.h>
#include <FV.h>

namespace tasks
{

namespace pg
{

PostureGenerator::PostureGenerator(const rbd::MultiBody& mb,
  const rbd::MultiBodyConfig& mbc):
  mb_(mb),
  mbc_(mbc)
{}


void PostureGenerator::addObjective(Objective* obj)
{
  obj_.push_back(obj);
}


void PostureGenerator::removeObjective(Objective* obj)
{
  obj_.erase(std::find(obj_.begin(), obj_.end(), obj));
}


int PostureGenerator::nrObjectives() const
{
  return obj_.size();
}


void PostureGenerator::addConstraint(Constraint* co)
{
  constr_.push_back(co);
}


void PostureGenerator::removeConstraint(Constraint* co)
{
  constr_.erase(std::find(constr_.begin(), constr_.end(), co));
}


int PostureGenerator::nrConstraints() const
{
	return constr_.size();
}


bool PostureGenerator::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
	Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	n = mb_.nrParams();

	int size = 0;
	int nbNNZ = 0;

	for(Constraint* c: constr_)
	{
		size += c->size();
		std::vector<std::pair<int, int> > v = c->structure();
		nbNNZ += v.size();
		struct_.insert(struct_.end(), v.begin(), v.end());
	}

	m = size;
	nnz_jac_g = nbNNZ;

	nnz_h_lag = 0;
	index_style = C_STYLE;

	return true;
}


bool PostureGenerator::get_bounds_info(Ipopt::Index n,
	Ipopt::Number* x_l, Ipopt::Number* x_u,
	Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
	using namespace Eigen;
	Map<VectorXd>(x_l, n).fill(-std::numeric_limits<double>::infinity());
	Map<VectorXd>(x_u, n).fill(std::numeric_limits<double>::infinity());

	Map<VectorXd>(g_l, m).fill(0);
	Map<VectorXd>(g_u, m).fill(0);

	return true;
}


bool PostureGenerator::get_starting_point(Ipopt::Index n,
	bool /* init_x */, Ipopt::Number* x,
	bool /* init_z */, Ipopt::Number* /* z_L */, Ipopt::Number* /* z_U */,
	Ipopt::Index /* m */, bool /* init_lambda */,
	Ipopt::Number* /* lambda */)
{
	using namespace Eigen;
	Map<VectorXd>(x, n) = rbd::paramToVector(mb_, mbc_.q);

	return true;
}


bool PostureGenerator::eval_f(Ipopt::Index n, const Ipopt::Number* x,
	bool new_x, Ipopt::Number& obj_value)
{
	using namespace Eigen;
	if(new_x)
	{
		rbd::vectorToParam(VectorXd(Map<const VectorXd>(x, n)), mbc_.q);
		rbd::forwardKinematics(mb_, mbc_);
		rbd::forwardVelocity(mb_, mbc_);
	}

	obj_value = 0.;
	for(Objective* o: obj_)
	{
		obj_value += o->value()*o->weight();
	}

	return true;
}


bool PostureGenerator::eval_grad_f(Ipopt::Index n,
	const Ipopt::Number* /* x */, bool /* new_x */, Ipopt::Number* grad_f)
{
	using namespace Eigen;
	Map<VectorXd> grad(grad_f, n);

	grad.setZero();
	for(Objective* o: obj_)
	{
		grad += o->gradient()*o->weight();
	}

	return true;
}


bool PostureGenerator::eval_g(
	Ipopt::Index /* n */, const Ipopt::Number* /* x */, bool /* new_x */,
	Ipopt::Index m, Ipopt::Number* g)
{
	using namespace Eigen;
	Map<VectorXd> gVal(g, m);

	int pos = 0;
	for(Constraint* c: constr_)
	{
		int s = c->size();
		gVal.segment(pos, s) = c->value();
		pos += s;
	}
	return true;
}


bool PostureGenerator::eval_jac_g(
	Ipopt::Index /* n */, const Ipopt::Number* /* x */, bool /* new_x */,
	Ipopt::Index /* m */, Ipopt::Index nele_jac,
	Ipopt::Index* iRow, Ipopt::Index *jCol, Ipopt::Number* values)
{
	using namespace Eigen;

	if(values == NULL)
	{
		for(int i = 0; i < nele_jac; ++i)
		{
			const std::pair<int, int>& p = struct_[i];
			iRow[i] = p.first;
			jCol[i] = p.second;
		}
	}
	else
	{
		Map<VectorXd> jac(values, nele_jac);

		int pos = 0;
		for(Constraint* c: constr_)
		{
			const VectorXd& cJac = c->jac();
			int s = cJac.rows();
			jac.segment(pos, s) = cJac;
			pos += s;
		}
	}

	return true;
}

} // namespace pg

} // namespace tasks
