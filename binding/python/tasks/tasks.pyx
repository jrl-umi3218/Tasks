#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

# distutils: language = c++
cimport c_tasks
cimport eigen.c_eigen as c_eigen
cimport eigen.eigen as eigen
cimport rbdyn.rbdyn as rbdyn
cimport sva.sva as sva

from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef class cpu_times(object):
  property wall:
    def __get__(self):
      return self.impl.wall
    def __set__(self, value):
      self.impl.wall = value
  property user:
    def __get__(self):
      return self.impl.user
    def __set__(self, value):
      self.impl.user = value
  property system:
    def __get__(self):
      return self.impl.system
    def __set__(self, value):
      self.impl.system = value

cdef cpu_times cpu_timesFromC(const c_tasks.cpu_times & arg):
  cdef cpu_times ret = cpu_times()
  ret.impl = arg
  return ret

cdef class PositionTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, eigen.Vector3d pos, eigen.Vector3d bodyPoint = eigen.Vector3d.Zero()):
    self.impl = new c_tasks.PositionTask(deref(mb.impl), bodyName, pos.impl, bodyPoint.impl)
  def position(self, eigen.Vector3d pos = None):
    if pos is None:
      return eigen.Vector3dFromC(self.impl.position())
    else:
      self.impl.position(pos.impl)
  def bodyPoint(self, eigen.Vector3d point = None):
    if point is None:
      return eigen.Vector3dFromC(self.impl.bodyPoint())
    else:
      self.impl.bodyPoint(point.impl)
  # Common to *Task
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.update(deref(mb.impl), deref(mbc.impl))
  def updateDot(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.updateDot(deref(mb.impl), deref(mbc.impl))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class OrientationTask(object):
  def __dealloc__(self):
    del self.impl
  def __m3ctor__(self, rbdyn.MultiBody mb, bodyName, eigen.Matrix3d ori):
    self.impl = new c_tasks.OrientationTask(deref(mb.impl), <string>bodyName, ori.impl)
  def __qdctor__(self, rbdyn.MultiBody mb, bodyName, eigen.Quaterniond ori):
    self.impl = new c_tasks.OrientationTask(deref(mb.impl), <string>bodyName, ori.impl)
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, ori):
    if isinstance(ori, eigen.Matrix3d):
      self.__m3ctor__(mb, bodyName, ori)
    elif isinstance(ori, eigen.Quaterniond):
      self.__qdctor__(mb, bodyName, ori)
    else:
      raise TypeError("Wrong arguments passed to OrientationTask ctor")
  def __orientationm3(self, eigen.Matrix3d ori):
    self.impl.orientation(ori.impl)
  def __orientationqd(self, eigen.Quaterniond ori):
    self.impl.orientation(ori.impl)
  def orientation(self, ori = None):
    if ori is None:
      return eigen.Matrix3dFromC(self.impl.orientation())
    elif isinstance(ori, eigen.Matrix3d):
      return self.__orientationm3(ori)
    elif isinstance(ori, eigen.Quaterniond):
      return self.__orientationqd(ori)
    else:
      raise TypeError("Wrong argument passed to orientation")
  # Common to *Task
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.update(deref(mb.impl), deref(mbc.impl))
  def updateDot(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.updateDot(deref(mb.impl), deref(mbc.impl))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class SurfaceOrientationTask(object):
  def __dealloc__(self):
    del self.impl
  def __m3ctor__(self, rbdyn.MultiBody mb, bodyName, eigen.Matrix3d ori, sva.PTransformd X_b_s):
    self.impl = new c_tasks.SurfaceOrientationTask(deref(mb.impl), <string>bodyName, ori.impl, deref(X_b_s.impl))
  def __qdctor__(self, rbdyn.MultiBody mb, bodyName, eigen.Quaterniond ori, sva.PTransformd X_b_s):
    self.impl = new c_tasks.SurfaceOrientationTask(deref(mb.impl), <string>bodyName, ori.impl, deref(X_b_s.impl))
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, ori, sva.PTransformd X_b_s):
    if isinstance(ori, eigen.Matrix3d):
      self.__m3ctor__(mb, bodyName, ori, X_b_s)
    elif isinstance(ori, eigen.Quaterniond):
      self.__qdctor__(mb, bodyName, ori, X_b_s)
    else:
      raise TypeError("Wrong arguments passed to SurfaceOrientationTask ctor")
  def __orientationm3(self, eigen.Matrix3d ori):
    self.impl.orientation(ori.impl)
  def __orientationqd(self, eigen.Quaterniond ori):
    self.impl.orientation(ori.impl)
  def orientation(self, ori = None):
    if ori is None:
      return eigen.Matrix3dFromC(self.impl.orientation())
    elif isinstance(ori, eigen.Matrix3d):
      return self.__orientationm3(ori)
    elif isinstance(ori, eigen.Quaterniond):
      return self.__orientationqd(ori)
    else:
      raise TypeError("Wrong argument passed to orientation")
  # Common to *Task
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.update(deref(mb.impl), deref(mbc.impl))
  def updateDot(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.updateDot(deref(mb.impl), deref(mbc.impl))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class GazeTask(object):
  def __dealloc__(self):
    del self.impl
  def __v2ctor__(self, rbdyn.MultiBody mb, bodyName, eigen.Vector2d point2d, double depthEstimate, sva.PTransformd X_b_gaze, eigen.Vector2d point2d_ref):
    self.impl = new c_tasks.GazeTask(deref(mb.impl), bodyName, point2d.impl, depthEstimate, deref(X_b_gaze.impl), point2d_ref.impl)
  def __v3ctor__(self, rbdyn.MultiBody mb, bodyName, eigen.Vector3d point3d, sva.PTransformd X_b_gaze, eigen.Vector2d point2d_ref):
    self.impl = new c_tasks.GazeTask(deref(mb.impl), bodyName, point3d.impl, deref(X_b_gaze.impl), point2d_ref.impl)
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, point, *args):
    if isinstance(point, eigen.Vector3d):
      if len(args) > 0 and isinstance(args[0], sva.PTransformd):
        if len(args) == 2 and isinstance(args[1], eigen.Vector2d):
          self.__v3ctor__(mb, bodyName, point, args[0], args[1])
        elif len(args) == 1:
          self.__v3ctor__(mb, bodyName, point, args[0], eigen.Vector2d.Zero())
    elif isinstance(point, eigen.Vector2d):
      if len(args) >= 2 and isinstance(args[1], sva.PTransformd):
        if len(args) == 2:
          self.__v2ctor__(mb, bodyName, point, args[0], args[1], eigen.Vector2d.Zero())
        elif len(args) == 3 and isinstance(args[2], eigen.Vector2d):
          self.__v2ctor__(mb, bodyName, point, args[0], args[1], args[2])
    else:
      raise TypeError("Wrong arguments passed to GazeTask ctor")
  def __v3error(self, eigen.Vector3d point3d, eigen.Vector2d point2d_ref):
    self.impl.error(point3d.impl, point2d_ref.impl)
  def __v2error(self, eigen.Vector2d point2d, eigen.Vector2d point2d_ref):
    self.impl.error(point2d.impl, point2d_ref.impl)
  def error(self, point, eigen.Vector2d point2d_ref = eigen.Vector2d.Zero()):
    if isinstance(point, eigen.Vector3d):
      return self.__v3error(point, point2d_ref)
    else:
      return self.__v2error(point, point2d_ref)
  def speed(self):
    return eigen.VectorXdFromC(self.impl.speed())
  def normalAcc(self):
    return eigen.VectorXdFromC(self.impl.normalAcc())
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc, normalAccB):
    self.impl.update(deref(mb.impl), deref(mbc.impl), sva.MotionVecdVector(normalAccB).v)
  # Common to *Task
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class PositionBasedVisServoTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, sva.PTransformd X_t_s, sva.PTransformd X_b_s):
    self.impl = new c_tasks.PositionBasedVisServoTask(deref(mb.impl), bodyName, deref(X_t_s.impl), deref(X_b_s.impl))
  def error(self, sva.PTransformd X_t_s):
    self.impl.error(deref(X_t_s.impl))
  def speed(self):
    return eigen.VectorXdFromC(self.impl.speed())
  def normalAcc(self):
    return eigen.VectorXdFromC(self.impl.normalAcc())
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc, normalAccB):
    self.impl.update(deref(mb.impl), deref(mbc.impl), sva.MotionVecdVector(normalAccB).v)
  # Common to *Task
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class PostureTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, q):
    self.impl = new c_tasks.PostureTask(deref(mb.impl), q)
  def posture(self, q = None):
    if q is None:
      return self.impl.posture()
    else:
      self.impl.posture(q)
  # Common to *Task
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.update(deref(mb.impl), deref(mbc.impl))
  def updateDot(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.updateDot(deref(mb.impl), deref(mbc.impl))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class CoMTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, eigen.Vector3d com, weight = None):
    if weight is None:
      self.impl = new c_tasks.CoMTask(deref(mb.impl), com.impl)
    else:
      self.impl = new c_tasks.CoMTask(deref(mb.impl), com.impl, weight)
  def com(self, eigen.Vector3d com = None):
    if com is None:
      return eigen.Vector3dFromC(self.impl.com())
    else:
      self.impl.com(com.impl)
  def updateInertialParameters(self, rbdyn.MultiBody mb):
    self.impl.updateInertialParameters(deref(mb.impl))
  # Common to *Task
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.update(deref(mb.impl), deref(mbc.impl))
  def updateDot(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.updateDot(deref(mb.impl), deref(mbc.impl))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class MultiCoMTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBodyVector mbs, robotIndexes, eigen.Vector3d com):
    self.impl = new c_tasks.MultiCoMTask(deref(mbs.v), robotIndexes, com.impl)
  def com(self, eigen.Vector3d com = None):
    if com is None:
      return eigen.Vector3dFromC(self.impl.com())
    else:
      self.impl.com(com.impl)
  def updateInertialParameters(self, rbdyn.MultiBodyVector mbs):
    self.impl.updateInertialParameters(deref(mbs.v))
  def update(self, rbdyn.MultiBodyVector mbs, rbdyn.MultiBodyConfigVector mbcs):
    self.impl.update(deref(mbs.v), deref(mbcs.v))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def speed(self):
    return eigen.VectorXdFromC(self.impl.speed())
  def normalAcc(self):
    return eigen.VectorXdFromC(self.impl.normalAcc())
  def jac(self, int index):
    return eigen.MatrixXdFromC(self.impl.jac(index))

cdef class MultiRobotTransformTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBodyVector mbs, int r1Index, int r2Index, r1BodyName, r2BodyName, sva.PTransformd X_r1b_r1s, sva.PTransformd X_r2b_r2s):
    self.impl = new c_tasks.MultiRobotTransformTask(deref(mbs.v), r1Index, r2Index, r1BodyName, r2BodyName, deref(X_r1b_r1s.impl), deref(X_r2b_r2s.impl))
  def r1Index(self):
    return self.impl.r1Index()
  def X_r1b_r1s(self, sva.PTransformd X_r1b_r1s = None):
    if X_r1b_r1s is None:
      return sva.PTransformdFromC(self.impl.X_r1b_r1s())
    else:
      self.impl.X_r1b_r1s(deref(X_r1b_r1s.impl))
  def r2Index(self):
    return self.impl.r2Index()
  def X_r2b_r2s(self, sva.PTransformd X_r2b_r2s = None):
    if X_r2b_r2s is None:
      return sva.PTransformdFromC(self.impl.X_r2b_r2s())
    else:
      self.impl.X_r2b_r2s(deref(X_r2b_r2s.impl))
  def update(self, rbdyn.MultiBodyVector mbs, rbdyn.MultiBodyConfigVector mbcs, normalAccB):
    self.impl.update(deref(mbs.v), deref(mbcs.v), sva.MotionVecdVectorVector(normalAccB).v)
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def speed(self):
    return eigen.VectorXdFromC(self.impl.speed())
  def normalAcc(self):
    return eigen.VectorXdFromC(self.impl.normalAcc())
  def jac(self, int index):
    return eigen.MatrixXdFromC(self.impl.jac(index))

cdef class MomentumTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, sva.ForceVecd mom):
    self.impl = new c_tasks.MomentumTask(deref(mb.impl), deref(mom.impl))
  def momentum(self, sva.ForceVecd mom = None):
    if mom is None:
      return sva.ForceVecdFromC(self.impl.momentum())
    else:
      self.impl.momentum(deref(mom.impl))
  # Common to *Task
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.update(deref(mb.impl), deref(mbc.impl))
  def updateDot(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.updateDot(deref(mb.impl), deref(mbc.impl))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class LinVelocityTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, eigen.Vector3d pos, eigen.Vector3d bodyPoint = eigen.Vector3d.Zero()):
    self.impl = new c_tasks.LinVelocityTask(deref(mb.impl), bodyName, pos.impl, bodyPoint.impl)
  def velocity(self, eigen.Vector3d pos = None):
    if pos is None:
      return eigen.Vector3dFromC(self.impl.velocity())
    else:
      self.impl.velocity(pos.impl)
  def bodyPoint(self, eigen.Vector3d point = None):
    if point is None:
      return eigen.Vector3dFromC(self.impl.bodyPoint())
    else:
      self.impl.bodyPoint(point.impl)
  # Common to *Task
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.update(deref(mb.impl), deref(mbc.impl))
  def updateDot(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.updateDot(deref(mb.impl), deref(mbc.impl))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class OrientationTrackingTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, eigen.Vector3d bodyPoint, eigen.Vector3d bodyAxis, trackingJointsNames, eigen.Vector3d trackedPoint):
    self.impl = new c_tasks.OrientationTrackingTask(deref(mb.impl), bodyName, bodyPoint.impl, bodyAxis.impl, trackingJointsNames, trackedPoint.impl)
  def trackedPoint(self, eigen.Vector3d point = None):
    if point is None:
      return eigen.Vector3dFromC(self.impl.trackedPoint())
    else:
      self.impl.trackedPoint(point.impl)
  def bodyPoint(self, eigen.Vector3d point = None):
    if point is None:
      return eigen.Vector3dFromC(self.impl.bodyPoint())
    else:
      self.impl.bodyPoint(point.impl)
  def bodyAxis(self, eigen.Vector3d axis = None):
    if axis is None:
      return eigen.Vector3dFromC(self.impl.bodyAxis())
    else:
      self.impl.bodyAxis(axis.impl)
  # Common to *Task
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.update(deref(mb.impl), deref(mbc.impl))
  def updateDot(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc):
    self.impl.updateDot(deref(mb.impl), deref(mbc.impl))
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def jacDot(self):
    return eigen.MatrixXdFromC(self.impl.jacDot())

cdef class TransformTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, sva.PTransformd X_0_t, sva.PTransformd X_b_p = sva.PTransformd.Identity(), eigen.Matrix3d E_0_c = eigen.Matrix3d.Identity()):
    self.impl = new c_tasks.TransformTask(deref(mb.impl), bodyName, deref(X_0_t.impl), deref(X_b_p.impl), E_0_c.impl)
  def E_0_c(self, eigen.Matrix3d X_0_t = None):
    if X_0_t is None:
      return eigen.Matrix3dFromC(self.impl.E_0_c())
    else:
      self.impl.E_0_c(X_0_t.impl)
  # Common to *TransformTask
  def target(self, sva.PTransformd X_0_t = None):
    if X_0_t is None:
      return sva.PTransformdFromC(self.impl.target())
    else:
      self.impl.target(deref(X_0_t.impl))
  def X_b_p(self, sva.PTransformd X_b_p = None):
    if X_b_p is None:
      return sva.PTransformdFromC(self.impl.X_b_p())
    else:
      self.impl.X_b_p(deref(X_b_p.impl))
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc, mbcs):
    self.impl.update(deref(mb.impl), deref(mbc.impl), sva.MotionVecdVector(mbcs).v)
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def speed(self):
    return eigen.VectorXdFromC(self.impl.speed())
  def normalAcc(self):
    return eigen.VectorXdFromC(self.impl.normalAcc())

cdef class SurfaceTransformTask(object):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, rbdyn.MultiBody mb, bodyName, sva.PTransformd X_0_t, sva.PTransformd X_b_p = sva.PTransformd.Identity()):
    self.impl = new c_tasks.SurfaceTransformTask(deref(mb.impl), bodyName, deref(X_0_t.impl), deref(X_b_p.impl))
  # Common to *TransformTask
  def target(self, sva.PTransformd X_0_t = None):
    if X_0_t is None:
      return sva.PTransformdFromC(self.impl.target())
    else:
      self.impl.target(deref(X_0_t.impl))
  def X_b_p(self, sva.PTransformd X_b_p = None):
    if X_b_p is None:
      return sva.PTransformdFromC(self.impl.X_b_p())
    else:
      self.impl.X_b_p(deref(X_b_p.impl))
  def update(self, rbdyn.MultiBody mb, rbdyn.MultiBodyConfig mbc, mbcs):
    self.impl.update(deref(mb.impl), deref(mbc.impl), sva.MotionVecdVector(mbcs).v)
  def eval(self):
    return eigen.VectorXdFromC(self.impl.eval())
  def jac(self):
    return eigen.MatrixXdFromC(self.impl.jac())
  def speed(self):
    return eigen.VectorXdFromC(self.impl.speed())
  def normalAcc(self):
    return eigen.VectorXdFromC(self.impl.normalAcc())

cdef class QBound(object):
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_tasks.QBound()
    elif len(args) == 2:
      self.impl = c_tasks.QBound(args[0], args[1])
    else:
      raise TypeError("Wrong arguments passed to QBound ctor")
  property lQBound:
    def __get__(self):
      return self.impl.lQBound
    def __set__(self, value):
      self.impl.lQBound = value
  property uQBound:
    def __get__(self):
      return self.impl.uQBound
    def __set__(self, value):
      self.impl.uQBound = value

cdef class AlphaBound(object):
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_tasks.AlphaBound()
    elif len(args) == 2:
      self.impl = c_tasks.AlphaBound(args[0], args[1])
    else:
      raise TypeError("Wrong arguments passed to AlphaBound ctor")
  property lAlphaBound:
    def __get__(self):
      return self.impl.lAlphaBound
    def __set__(self, value):
      self.impl.lAlphaBound = value
  property uAlphaBound:
    def __get__(self):
      return self.impl.uAlphaBound
    def __set__(self, value):
      self.impl.uAlphaBound = value

cdef class TorqueBound(object):
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_tasks.TorqueBound()
    elif len(args) == 2:
      self.impl = c_tasks.TorqueBound(args[0], args[1])
    else:
      raise TypeError("Wrong arguments passed to TorqueBound ctor")
  property lTorqueBound:
    def __get__(self):
      return self.impl.lTorqueBound
    def __set__(self, value):
      self.impl.lTorqueBound = value
  property uTorqueBound:
    def __get__(self):
      return self.impl.uTorqueBound
    def __set__(self, value):
      self.impl.uTorqueBound = value

cdef class PolyTorqueBound(object):
  def __vctor__(self, lPTB, uPTB):
    cdef vector[vector[c_eigen.VectorXd]] vlPTB
    cdef vector[vector[c_eigen.VectorXd]] vuPTB
    for vv in lPTB:
      vlPTB.push_back(eigen.VectorXdVector(vv).v)
    for vv in uPTB:
      vuPTB.push_back(eigen.VectorXdVector(vv).v)
    self.impl = c_tasks.PolyTorqueBound(vlPTB, vuPTB)
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_tasks.PolyTorqueBound()
    elif len(args) == 2:
      self.__vctor__(args[0], args[1])
    else:
      raise TypeError("Wrong arguments passed to PolyTorqueBound ctor")
  property lPolyTorqueBound:
    def __get__(self):
      ret = []
      for vv in self.impl.lPolyTorqueBound:
        ret.append([])
        for v in vv:
          ret[-1].append(eigen.VectorXdFromC(v))
      return ret
    def __set__(self, value):
      cdef vector[vector[c_eigen.VectorXd]] newV
      for vv in value:
        newV.push_back(eigen.VectorXdVector(vv).v)
      self.impl.lPolyTorqueBound = newV
  property uPolyTorqueBound:
    def __get__(self):
      ret = []
      for vv in self.impl.lPolyTorqueBound:
        ret.append([])
        for v in vv:
          ret[-1].append(eigen.VectorXdFromC(v))
      return ret
    def __set__(self, value):
      cdef vector[vector[c_eigen.VectorXd]] newV
      for vv in value:
        newV.push_back(eigen.VectorXdVector(vv).v)
      self.impl.lPolyTorqueBound = newV
