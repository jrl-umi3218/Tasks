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
#include "GenQPSolver.h"

// includes
// std
#include <map>

// Tasks
#include "QLDQPSolver.h"
#include "LSSOLQPSolver.h"


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


} // namespace qp

} // namespace tasks
