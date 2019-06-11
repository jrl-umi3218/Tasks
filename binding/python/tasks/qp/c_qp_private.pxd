#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from rbdyn.c_rbdyn cimport *
from c_qp cimport *
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "qp_wrapper.hpp" namespace "tasks::qp":
    JointsSelector* ActiveJoints2Ptr(const vector[MultiBody]&, int, HighLevelTask*, const vector[string])
    JointsSelector* UnactiveJoints2Ptr(const vector[MultiBody]&, int, HighLevelTask*, const vector[string])
