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

// Ipopt
// Remove unused parameter warning from IpTNLP header
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>
#pragma GCC diagnostic pop

// RBDyn
#include <FK.h>
#include <FV.h>

namespace tasks
{

namespace pg
{



/**
	*													PostureGeneratorWrapper
	*/



class PostureGeneratorWrapper : public Ipopt::TNLP
{
public:
	PostureGeneratorWrapper(PostureGenerator* pg):
	pg_(pg)
	{}

	virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
		Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
		IndexStyleEnum& index_style)
	{
		index_style = C_STYLE;
		return pg_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag);
	}

	virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l,
		Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
	{
		return pg_->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
	}

	virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
		bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
		Ipopt::Index m, bool init_lambda,
		Ipopt::Number* lambda)
	{
		return pg_->get_starting_point(n, init_x, x, init_z, z_L, z_U, m,
			init_lambda, lambda);
	}

	virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x,
		bool new_x, Ipopt::Number& obj_value)
	{
		return pg_->eval_f(n, x, new_x, obj_value);
	}

	virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x,
		bool new_x, Ipopt::Number* grad_f)
	{
		return pg_->eval_grad_f(n, x, new_x, grad_f);
	}

	virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
		Ipopt::Index m, Ipopt::Number* g)
	{
		return pg_->eval_g(n, x, new_x, m, g);
	}

	virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
		Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
		Ipopt::Index *jCol, Ipopt::Number* values)
	{
		return pg_->eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
	}

	virtual void finalize_solution(Ipopt::SolverReturn status, int n,
		const Ipopt::Number* x, const Ipopt::Number* /* z_L */,
		const Ipopt::Number* /* z_U */, Ipopt::Index m, const Ipopt::Number* g,
		const Ipopt::Number* /* lambda */, Ipopt::Number obj_value,
		const Ipopt::IpoptData* /* ip_data */,
		Ipopt::IpoptCalculatedQuantities* /* ip_cq */)
	{
		pg_->finalize_solution(status, n, x, m, g, obj_value);
	}

private:
	PostureGenerator* pg_;
};



/**
	*														PostureGenerator
	*/



PostureGenerator::PostureGenerator(const rbd::MultiBody& mb,
  const rbd::MultiBodyConfig& mbc):
  mb_(mb),
  mbc_(mbc)
{}


bool PostureGenerator::solve()
{
	using namespace Ipopt;
	SmartPtr<TNLP> pg = new PostureGeneratorWrapper(this);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->Options()->SetNumericValue("tol", 1e-9);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

	ApplicationReturnStatus status;
	status = app->Initialize();
	if(status != Solve_Succeeded)
		return false;

	status = app->OptimizeTNLP(pg);

	return status == Solve_Succeeded;
}


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


bool PostureGenerator::get_nlp_info(int& n, int& m,
	int& nnz_jac_g, int& nnz_h_lag)
{
	n = mb_.nrParams();

	int size = 0;
	int nbNNZ = 0;

	for(Constraint* c: constr_)
	{
		std::vector<std::pair<int, int> > v = c->structure(mb_);

		size += c->size();
		nbNNZ += v.size();
		struct_.insert(struct_.end(), v.begin(), v.end());
	}

	m = size;
	nnz_jac_g = nbNNZ;

	nnz_h_lag = 0;

	return true;
}


bool PostureGenerator::get_bounds_info(int n,
	double* x_l, double* x_u,
	int m, double* g_l, double* g_u)
{
	using namespace Eigen;
	Map<VectorXd>(x_l, n).fill(-std::numeric_limits<double>::infinity());
	Map<VectorXd>(x_u, n).fill(std::numeric_limits<double>::infinity());

	Map<VectorXd> gLow(g_l, m);
	Map<VectorXd> gUpp(g_u, m);

	int pos = 0;
	for(Constraint* c: constr_)
	{
		int s = c->size();
		gLow.segment(pos, s) = c->lower();
		gUpp.segment(pos, s) = c->upper();
		pos += s;
	}

	return true;
}


bool PostureGenerator::get_starting_point(int n,
	bool /* init_x */, double* x,
	bool /* init_z */, double* /* z_L */, double* /* z_U */,
	int /* m */, bool /* init_lambda */,
	double* /* lambda */)
{
	using namespace Eigen;
	Map<VectorXd>(x, n) = rbd::paramToVector(mb_, mbc_.q);

	return true;
}


bool PostureGenerator::eval_f(int n, const double* x,
	bool new_x, double& obj_value)
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


bool PostureGenerator::eval_grad_f(int n,
	const double* /* x */, bool /* new_x */, double* grad_f)
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
	int /* n */, const double* /* x */, bool /* new_x */,
	int m, double* g)
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
	int /* n */, const double* /* x */, bool /* new_x */,
	int /* m */, int nele_jac,
	int* iRow, int *jCol, double* values)
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


void PostureGenerator::finalize_solution(int status, int n,
	const double* x, int  /* m */, const double* /* g */, double /* obj_value */)
{
	using namespace Eigen;

	if(status == Success)
	{
		rbd::vectorToParam(VectorXd(Map<const VectorXd>(x, n)), mbc_.q);
		rbd::forwardKinematics(mb_, mbc_);
		rbd::forwardVelocity(mb_, mbc_);
	}
}


} // namespace pg

} // namespace tasks
