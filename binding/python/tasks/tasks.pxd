#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

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
