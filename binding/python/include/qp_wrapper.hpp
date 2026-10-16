/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <Tasks/QPSolver.h>
#include <Tasks/QPTasks.h>

#include <chrono>

namespace tasks
{
namespace qp
{

JointsSelector * ActiveJoints2Ptr(const std::vector<rbd::MultiBody> & mbs,
                                  int robotIndex,
                                  HighLevelTask * hl,
                                  const std::vector<std::string> & activeJointsNames)
{ return new JointsSelector(JointsSelector::ActiveJoints(mbs, robotIndex, hl, activeJointsNames)); }

JointsSelector * UnactiveJoints2Ptr(const std::vector<rbd::MultiBody> & mbs,
                                    int robotIndex,
                                    HighLevelTask * hl,
                                    const std::vector<std::string> & unactiveJointsNames)
{ return new JointsSelector(JointsSelector::UnactiveJoints(mbs, robotIndex, hl, unactiveJointsNames)); }

double QPSolverSolveTime(const QPSolver & solver)
{ return std::chrono::duration<double>(solver.solveTime()).count(); }

double QPSolverSolveAndBuildTime(const QPSolver & solver)
{ return std::chrono::duration<double>(solver.solveAndBuildTime()).count(); }

} // namespace qp
} // namespace tasks
