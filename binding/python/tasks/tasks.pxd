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

cimport c_tasks
from libcpp.vector cimport vector

cdef class cpu_times(object):
  cdef c_tasks.cpu_times impl

cdef cpu_times cpu_timesFromC(const c_tasks.cpu_times &)

cdef class PositionTask(object):
  cdef c_tasks.PositionTask * impl

cdef class OrientationTask(object):
  cdef c_tasks.OrientationTask * impl

cdef class SurfaceOrientationTask(object):
  cdef c_tasks.SurfaceOrientationTask * impl

cdef class GazeTask(object):
  cdef c_tasks.GazeTask * impl

cdef class PositionBasedVisServoTask(object):
  cdef c_tasks.PositionBasedVisServoTask * impl

cdef class PostureTask(object):
  cdef c_tasks.PostureTask * impl

cdef class CoMTask(object):
  cdef c_tasks.CoMTask * impl

cdef class MultiCoMTask(object):
  cdef c_tasks.MultiCoMTask * impl

cdef class MultiRobotTransformTask(object):
  cdef c_tasks.MultiRobotTransformTask * impl

cdef class MomentumTask(object):
  cdef c_tasks.MomentumTask * impl

cdef class LinVelocityTask(object):
  cdef c_tasks.LinVelocityTask * impl

cdef class OrientationTrackingTask(object):
  cdef c_tasks.OrientationTrackingTask * impl

cdef class TransformTask(object):
  cdef c_tasks.TransformTask * impl

cdef class SurfaceTransformTask(object):
  cdef c_tasks.SurfaceTransformTask * impl

cdef class QBound(object):
  cdef c_tasks.QBound impl

cdef class AlphaBound(object):
  cdef c_tasks.AlphaBound impl

cdef class TorqueBound(object):
  cdef c_tasks.TorqueBound impl

cdef class PolyTorqueBound(object):
  cdef c_tasks.PolyTorqueBound impl
