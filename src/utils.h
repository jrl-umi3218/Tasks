/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

namespace tasks
{

namespace qp
{

inline int findJointFromVector(const rbd::MultiBody & mb, int line, bool withBase)
{
  int start = withBase ? 0 : 1;
  for(int j = start; j < int(mb.nrJoints()); ++j)
  {
    if(line >= start && line <= (start + mb.joint(j).dof()))
    {
      return j;
    }
    start += mb.joint(j).dof();
  }
  return -1;
}

inline bool compareDof(const rbd::MultiBody & mb1, const rbd::MultiBody & mb2)
{
  return mb1.nrDof() < mb2.nrDof();
}

} // namespace qp

} // namespace tasks
