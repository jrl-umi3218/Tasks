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
#include "QPSolverData.h"

// includes
// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>


namespace tasks
{

namespace qp
{

SolverData::SolverData():
	alphaD_(0),
	lambda_(0),
	torque_(0),
	nrVars_(0),
	uniCont_(),
	biCont_()
{}


void SolverData::computeNormalAccB(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc)
{
	const std::vector<int>& pred = mb.predecessors();
	const std::vector<int>& succ = mb.successors();

	for(int i = 0; i < mb.nrJoints(); ++i)
	{
		const sva::PTransformd& X_p_i = mbc.parentToSon[i];
		const sva::MotionVecd& vj_i = mbc.jointVelocity[i];
		const sva::MotionVecd& vb_i = mbc.bodyVelB[i];

		if(pred[i] != -1)
			normalAccB_[succ[i]] = X_p_i*normalAccB_[pred[i]] + vb_i.cross(vj_i);
		else
			normalAccB_[succ[i]] = vb_i.cross(vj_i);
	}
}

} // namespace qp

} // namespace tasks
