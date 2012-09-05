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
#include "QL.h"

// includes
// std
#include <iostream>

namespace tasks
{

namespace qp
{

bool solveQP(int N, int ME, int mi,
		/*const*/ Eigen::MatrixXd& Q, /* const */ Eigen::VectorXd& c,
		const Eigen::MatrixXd& A1, const Eigen::VectorXd& B1,
		const Eigen::MatrixXd& A2, const Eigen::VectorXd& B2,
		/*const*/ Eigen::VectorXd& XLower,
		/*const*/ Eigen::VectorXd& XUpper,
		Eigen::VectorXd& X,
		double epsilon)
{
	int M = ME+mi;
	int MMAX=M+1;
	int NMAX=N;
	int MNN=M+N+N;

	// Note: the fusion of A1 and A2 could be realized with
	//  Eigen::MatrixXd A (MMAX,N);
	//  A << -A1, -A2, Eigen::MatrixXd::Zero(1,N);
	// but the computation time is higher.
	double * A=new double[MMAX*N];
	for(int j=0; j<N; ++j)
		for(int i=0; i<ME; ++i)
			A[j*(MMAX)+i]=-A1(i,j);
	for(int j=0; j<N; ++j)
		for(int i=0; i<mi; ++i)
			A[j*(MMAX)+i+ME]=-A2(i,j);
	for(int j=0; j<N; ++j)
		A[j*(MMAX)+mi+ME]=0;

	// But this is ok.
	Eigen::VectorXd B(M);
	B << B1, B2;

	double *U=new double[(MNN)];
	int IOUT=6;
	int IFAIL;
	int IPRINT=1;

	int LWAR= std::ceil(3.0*(NMAX)*(NMAX)/2.0 + 10*(NMAX) + 2*(MMAX))+1;
	int LIWAR= N;

	double * WAR=new double[LWAR];
	int * IWAR=new int[LIWAR];
	IWAR[0]=1;

	int MODE = 1;
	ql_(&M, &ME, &MMAX, &N, &NMAX, &MNN,
		 Q.data(), c.data(), A, B.data(),
		 XLower.data(), XUpper.data(),
		 X.data(), U, &epsilon, &MODE, &IOUT, &IFAIL,
		 &IPRINT, WAR, &LWAR, IWAR, &LIWAR);

	// cleaning
	delete[] A;
	delete[] U;
	delete[] WAR;
	delete[] IWAR;


	if(IFAIL!=0)
		std::cerr<<"QP resolution failed. IFAIL: "<<IFAIL<<std::endl;

	return (IFAIL==0);
}

} // namespace qp

} // namespace tasks
