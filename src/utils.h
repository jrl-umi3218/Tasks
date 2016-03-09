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

#pragma once

namespace tasks
{

namespace qp
{

inline int findJointFromVector(const rbd::MultiBody& mb, int line, bool withBase)
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

inline bool compareDof(const rbd::MultiBody& mb1, const rbd::MultiBody& mb2)
{
	return mb1.nrDof() < mb2.nrDof();
}

} // namespace qp

} // namespace tasks
