// Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
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
