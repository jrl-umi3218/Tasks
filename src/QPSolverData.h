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
// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// Tasks
#include "QPContacts.h"


// forward declaration
namespace rbd
{
class MultiBody;
class MultiBodyConfig;
}

namespace tasks
{

namespace qp
{


class SolverData
{
public:
	friend class QPSolver;

	SolverData();

	int nrVars() const
	{
		return nrVars_;
	}

	int alphaD() const
	{
		return alphaD_;
	}

	int lambda() const
	{
		return lambda_;
	}

	int torque() const
	{
		return torque_;
	}


	int alphaDBegin() const
	{
		return 0;
	}

	int lambdaBegin() const
	{
		return alphaD_;
	}

	int torqueBegin() const
	{
		return alphaD_ + lambda_;
	}


	int nrContacts() const
	{
		return static_cast<int>(uniCont_.size() + biCont_.size());
	}

	int nrUnilateralLambda() const
	{
		return lambdaUni_;
	}

	int nrBilateralLambda() const
	{
		return lambdaBi_;
	}

	int unilateralBegin() const
	{
		return alphaD_;
	}

	int bilateralBegin() const
	{
		return alphaD_ + lambdaUni_;
	}


	const std::vector<UnilateralContact>& unilateralContacts() const
	{
		return uniCont_;
	}

	const std::vector<BilateralContact>& bilateralContacts() const
	{
		return biCont_;
	}


	void computeNormalAccB(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

private:
	int alphaD_;
	int lambda_;
	int torque_;
	int nrVars_;

	int lambdaUni_;
	int lambdaBi_;

	std::vector<UnilateralContact> uniCont_;
	std::vector<BilateralContact> biCont_;

	std::vector<sva::MotionVecd> normalAccB_;
};


} // namespace qp

} // namespace tasks

