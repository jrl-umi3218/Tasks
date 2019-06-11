/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <Tasks/QPTasks.h>

namespace tasks
{
namespace qp
{

JointsSelector* ActiveJoints2Ptr(const std::vector<rbd::MultiBody>& mbs, int robotIndex, HighLevelTask* hl, const std::vector<std::string> & activeJointsNames)
{
  return new JointsSelector(JointsSelector::ActiveJoints(mbs, robotIndex, hl, activeJointsNames));
}

JointsSelector* UnactiveJoints2Ptr(const std::vector<rbd::MultiBody>& mbs, int robotIndex, HighLevelTask* hl, const std::vector<std::string> & unactiveJointsNames)
{
  return new JointsSelector(JointsSelector::UnactiveJoints(mbs, robotIndex, hl, unactiveJointsNames));
}

}
}
