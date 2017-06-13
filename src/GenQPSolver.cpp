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
#include "Tasks/GenQPSolver.h"

// includes
// std
#include <map>

// Tasks
#include "QLDQPSolver.h"

#ifdef LSSOL_SOLVER_FOUND
	#include "LSSOLQPSolver.h"
#endif


namespace tasks
{

namespace qp
{


#ifdef LSSOL_SOLVER_FOUND
	const std::string GenQPSolver::default_qp_solver("LSSOL");
#else
	const std::string GenQPSolver::default_qp_solver("QLD");
#endif


template <typename T>
T* allocateQP()
{
	return new T;
}


static const
std::map<std::string, std::function<GenQPSolver*(void)> > qpFactory = {
#ifdef LSSOL_SOLVER_FOUND
	{"LSSOL", allocateQP<LSSOLQPSolver>},
#endif
	{"QLD", allocateQP<QLDQPSolver>}
};


GenQPSolver* createQPSolver(const std::string& name)
{
	return qpFactory.at(name)();
}

void GenQPSolver::setDependencies(int nrVars, std::vector<std::tuple<int, int, double>> dependencies)
{
	dependencies_ = dependencies;
	fullToReduced_.resize(nrVars, -1);
	reducedToFull_.resize(nrVars - static_cast<int>(dependencies_.size()), -1);
	/* Retrieve the variables which are removed due to the dependencies */
	std::vector<int> removedVars; removedVars.reserve(dependencies.size() + 1);
	for(const auto & d : dependencies)
	{
		removedVars.push_back(std::get<1>(d));
	}
	/* Prevent issue once we have gone past the last removed variable */
	removedVars.push_back(nrVars);
	std::sort(removedVars.begin(), removedVars.end());
	/* Number of removed variables encountered so far */
	size_t shift = 0;
	for(size_t i = 0; i < fullToReduced_.size(); ++i)
	{
		if(i == static_cast<size_t>(removedVars[shift])) { shift++; continue; }
		fullToReduced_[i] = static_cast<int>(i - shift);
		reducedToFull_[i-shift] = static_cast<int>(i);
	}
}

} // namespace qp

} // namespace tasks
