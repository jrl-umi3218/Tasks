// Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
//
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
#include "Tasks/QPSolverData.h"

// includes
// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>


namespace tasks
{

namespace qp
{

SolverData::SolverData():
	alphaD_(),
	alphaDBegin_(),
	lambda_(),
	totalAlphaD_(0),
	totalLambda_(0),
	nrUniLambda_(0),
	nrBiLambda_(0),
	nrVars_(0),
	uniCont_(),
	biCont_(),
	allCont_(),
	mobileRobotIndex_(),
	normalAccB_()
{}


void SolverData::computeNormalAccB(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs)
{
	// we just need to update mobile robot normal acceleration
	for(int r: mobileRobotIndex_)
	{
		const rbd::MultiBody& mb = mbs[r];
		const rbd::MultiBodyConfig& mbc = mbcs[r];
		std::vector<sva::MotionVecd>& normalAccBr = normalAccB_[r];

		const std::vector<int>& pred = mb.predecessors();
		const std::vector<int>& succ = mb.successors();

		for(int i = 0; i < mb.nrJoints(); ++i)
		{
			const sva::PTransformd& X_p_i = mbc.parentToSon[i];
			const sva::MotionVecd& vj_i = mbc.jointVelocity[i];
			const sva::MotionVecd& vb_i = mbc.bodyVelB[i];

			if(pred[i] != -1)
				normalAccBr[succ[i]] = X_p_i*normalAccBr[pred[i]] + vb_i.cross(vj_i);
			else
				normalAccBr[succ[i]] = vb_i.cross(vj_i);
		}
	}
}

} // namespace qp

} // namespace tasks
