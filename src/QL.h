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
#include <Eigen/Core>

namespace tasks
{

namespace qp
{

extern "C" int ql_(int *m,int *me,int *mmax,int *n,int *nmax,int *mnn,
			 double *c,double *d,double *a,double *b,double *xl,
			 double *xu,double *x,double *u,double *eps, int* mode,int *iout,int *ifail,
			 int *iprint,double *war,int *lwar,int *iwar,int *liwar
	);


bool solveQP(int N, int ME, int mi,
		/*const*/ Eigen::MatrixXd& Q, /* const */ Eigen::VectorXd& c,
		const Eigen::MatrixXd& A1, const Eigen::VectorXd& B1,
		const Eigen::MatrixXd& A2, const Eigen::VectorXd& B2,
		/*const*/ Eigen::VectorXd& XLower,
		/*const*/ Eigen::VectorXd& XUpper,
		Eigen::VectorXd& X,
		double epsilon);

} // namespace tasks

} // namespace qp

