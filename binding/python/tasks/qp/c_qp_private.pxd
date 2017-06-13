# Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
#
# This file is part of Tasks.
#
# Tasks is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Tasks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Tasks.  If not, see <http://www.gnu.org/licenses/>.

from rbdyn.c_rbdyn cimport *
from c_qp cimport *
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "qp_wrapper.hpp" namespace "tasks::qp":
    JointsSelector* ActiveJoints2Ptr(const vector[MultiBody]&, int, HighLevelTask*, const vector[string])
    JointsSelector* UnactiveJoints2Ptr(const vector[MultiBody]&, int, HighLevelTask*, const vector[string])
