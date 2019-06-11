# distutils: language = c++

#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cimport sva.c_sva
cimport c_qp
cimport c_qp_private
from eigen.eigen cimport *
from rbdyn.rbdyn cimport *
from sva.sva cimport *
from sch.sch cimport *
cimport tasks.tasks as tasks

from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp.vector cimport vector

def check_args(argList, typeList):
  if len(argList) != len(typeList):
    return False
  for arg,type in zip(argList,typeList):
    if type is not None and not isinstance(arg, type):
      return False
  return True

cdef class FrictionCone(object):
  def __ctor__(self, Matrix3d frame, int nrGen, double mu, double direction = 1):
    self.impl = c_qp.FrictionCone(frame.impl, nrGen, mu, direction)
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_qp.FrictionCone()
    elif check_args(args, [Matrix3d, None, None]) or check_args(args, [Matrix3d, None, None, None]):
      self.__ctor__(*args)
    else:
      raise TypeError("Wrong arguments passed to FrictionCone ctor")
  property generators:
    def __get__(self):
      ret = []
      for v in self.impl.generators:
        ret.append(Vector3dFromC(v))
      return ret
    def __set__(self, value):
      self.impl.generators = Vector3dVector(value).v

cdef class FrictionConeVector(object):
  def __addFC(self, FrictionCone fc):
    self.v.push_back(fc.impl)
  def __cinit__(self, *args):
    if len(args) == 1 and isinstance(args[0], list):
      args = args[0]
    for fc in args:
      self.__addFC(fc)


cdef FrictionConeFromC(c_qp.FrictionCone fc):
  cdef FrictionCone ret = FrictionCone()
  ret.impl = fc
  return ret

cdef class ContactId(object):
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_qp.ContactId()
    elif len(args) == 4 or len(args) == 5:
      r1BodyName = args[2]
      r2BodyName = args[3]
      if isinstance(r1BodyName, unicode):
        r1BodyName = r1BodyName.encode(u'ascii')
      if isinstance(r2BodyName, unicode):
        r2BodyName = r2BodyName.encode(u'ascii')
      ambId = -1
      if len(args) == 5:
        ambId = args[4]
      self.impl = c_qp.ContactId(args[0], args[1], r1BodyName, r2BodyName, ambId)
    else:
      raise TypeError("Wrong arguments passed to ContactId ctor")
  property r1Index:
    def __get__(self):
      return self.impl.r1Index
    def __set__(self, value):
      self.impl.r1Index = value
  property r1BodyName:
    def __get__(self):
      return self.impl.r1BodyName
    def __set__(self, value):
      self.impl.r1BodyName = value
  property r2Index:
    def __get__(self):
      return self.impl.r2Index
    def __set__(self, value):
      self.impl.r2Index = value
  property r2BodyName:
    def __get__(self):
      return self.impl.r2BodyName
    def __set__(self, value):
      self.impl.r2BodyName = value
  property ambiguityId:
    def __get__(self):
      return self.impl.ambiguityId
    def __set__(self, value):
      self.impl.ambiguityId = value
  def __richcmp__(ContactId self, ContactId other, int op):
    if op == 0:
      return self.impl < other.impl
    elif op == 2:
      return self.impl == other.impl
    elif op == 3:
      return self.impl != other.impl
    else:
      raise NotImplementedError("This comparison is not supported")

cdef ContactId ContactIdFromC(const c_qp.ContactId & ci):
  cdef ContactId ret = ContactId()
  ret.impl = ci
  return ret

cdef class UnilateralContact(object):
  def __noambidctor__(self, int r1Index, int r2Index, r1BodyName, r2BodyName, r1Points, Matrix3d r1Frame, PTransformd X_b1_b2, int nrGen, double mu, PTransformd X_b1_cf = PTransformd.Identity()):
    if isinstance(r1BodyName, unicode):
      r1BodyName = r1BodyName.encode(u'ascii')
    if isinstance(r2BodyName, unicode):
      r2BodyName = r2BodyName.encode(u'ascii')
    self.impl = c_qp.UnilateralContact(r1Index, r2Index, r1BodyName, r2BodyName, Vector3dVector(r1Points).v, r1Frame.impl, deref(X_b1_b2.impl), nrGen, mu, deref(X_b1_cf.impl))
  def __ambidctor__(self, int r1Index, int r2Index, r1BodyName, r2BodyName, int ambiguityId, r1Points, Matrix3d r1Frame, PTransformd X_b1_b2, int nrGen, double mu, PTransformd X_b1_cf = PTransformd.Identity()):
    if isinstance(r1BodyName, unicode):
      r1BodyName = r1BodyName.encode(u'ascii')
    if isinstance(r2BodyName, unicode):
      r2BodyName = r2BodyName.encode(u'ascii')
    self.impl = c_qp.UnilateralContact(r1Index, r2Index, r1BodyName, r2BodyName, ambiguityId, Vector3dVector(r1Points).v, r1Frame.impl, deref(X_b1_b2.impl), nrGen, mu, deref(X_b1_cf.impl))
  def __cidctor__(self, ContactId contactId, r1Points, Matrix3d r1Frame, PTransformd X_b1_b2, int nrGen, double mu, PTransformd X_b1_cf = PTransformd.Identity()):
    self.impl = c_qp.UnilateralContact(contactId.impl, Vector3dVector(r1Points).v, r1Frame.impl, deref(X_b1_b2.impl), nrGen, mu, deref(X_b1_cf.impl))
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_qp.UnilateralContact()
    elif check_args(args, [None, None, None, None, list, Matrix3d, PTransformd, None, None]) or check_args(args, [None, None, None, None, list, Matrix3d, PTransformd, None, None, PTransformd]):
      self.__noambidctor__(*args)
    elif check_args(args, [None, None, None, None, None, list, Matrix3d, PTransformd, None, None]) or check_args(args, [None, None, None, None, None, list, Matrix3d, PTransformd, None, None, PTransformd]):
      self.__ambidctor__(*args)
    elif check_args(args, [ContactId, list, Matrix3d, PTransformd, None, None]) or check_args(args, [ContactId, list, Matrix3d, PTransformd, None, None, PTransformd]):
      self.__cidctor__(*args)
    else:
      raise TypeError("Wrong arguments passed to UnilateralContact ctor")
  property contactId:
    def __get__(self):
      return ContactIdFromC(self.impl.contactId)
    def __set__(self, ContactId ci):
      self.impl.contactId = ci.impl
  property r1Points:
    def __get__(self):
      ret = []
      for v in self.impl.r1Points:
        ret.append(Vector3dFromC(v))
      return ret
    def __set__(self, value):
      self.impl.r1Points = Vector3dVector(value).v
  property r2Points:
    def __get__(self):
      ret = []
      for v in self.impl.r2Points:
        ret.append(Vector3dFromC(v))
      return ret
    def __set__(self, value):
      self.impl.r2Points = Vector3dVector(value).v
  property r1Cone:
    def __get__(self):
      return FrictionConeFromC(self.impl.r1Cone)
    def __set__(self, FrictionCone value):
      self.impl.r1Cone = value.impl
  property r2Cone:
    def __get__(self):
      return FrictionConeFromC(self.impl.r2Cone)
    def __set__(self, FrictionCone value):
      self.impl.r2Cone = value.impl
  property X_b1_b2:
    def __get__(self):
      return PTransformdFromC(self.impl.X_b1_b2, copy = False)
    def __set__(self, PTransformd value):
      self.impl.X_b1_b2 = deref(value.impl)
  property X_b1_cf:
    def __get__(self):
      return PTransformdFromC(self.impl.X_b1_cf, copy = False)
    def __set__(self, PTransformd value):
      self.impl.X_b1_cf = deref(value.impl)
  def __forceFC(self, VectorXd _lambda, FrictionCone c):
    return Vector3dFromC(self.impl.sForce(_lambda.impl, c.impl))
  def __forceIFC(self, VectorXd _lambda, int point, FrictionCone c):
    return Vector3dFromC(self.impl.sForce(_lambda.impl, point, c.impl))
  def __forceVFC(self, VectorXd _lambda, r_b_pi, FrictionCone c):
    return ForceVecdFromC(self.impl.sForce(_lambda.impl, Vector3dVector(r_b_pi).v, c.impl))
  def force(self, VectorXd _lambda, *args):
    if check_args(args, [FrictionCone]):
      return self.__forceFC(_lambda, *args)
    elif check_args(args, [None, FrictionCone]):
      return self.__forceIFC(_lambda, *args)
    elif check_args(args, [list, FrictionCone]):
      return self.__forceVFC(_lambda, *args)
    else:
      raise TypeError("Wrong type passed for UnilateralContact.force")
  def nrLambda(self, point = None):
    if point is None:
      return self.impl.nrLambda()
    else:
      return self.impl.sNrLambda(point)

cdef UnilateralContactFromC(c_qp.UnilateralContact c):
  cdef UnilateralContact ret = UnilateralContact()
  ret.impl = c
  return ret

cdef class UnilateralContactVector(object):
  def __addUC(self, UnilateralContact uc):
    self.v.push_back(uc.impl)
  def __cinit__(self, *args):
    if len(args) == 1 and isinstance(args[0], list):
      args = args[0]
    for uc in args:
      self.__addUC(uc)

cdef class BilateralContact(object):
  def __noambidctor__(self, int r1Index, int r2Index, r1BodyName, r2BodyName, r1Points, r1Frames, PTransformd X_b1_b2, int nrGen, double mu, PTransformd X_b1_cf = PTransformd.Identity()):
    if isinstance(r1BodyName, unicode):
      r1BodyName = r1BodyName.encode(u'ascii')
    if isinstance(r2BodyName, unicode):
      r2BodyName = r2BodyName.encode(u'ascii')
    self.impl = c_qp.BilateralContact(r1Index, r2Index, r1BodyName, r2BodyName, Vector3dVector(r1Points).v, Matrix3dVector(r1Frames).v, deref(X_b1_b2.impl), nrGen, mu, deref(X_b1_cf.impl))
  def __ambidctor__(self, int r1Index, int r2Index, r1BodyName, r2BodyName, int ambiguityId, r1Points, r1Frames, PTransformd X_b1_b2, int nrGen, double mu, PTransformd X_b1_cf = PTransformd.Identity()):
    if isinstance(r1BodyName, unicode):
      r1BodyName = r1BodyName.encode(u'ascii')
    if isinstance(r2BodyName, unicode):
      r2BodyName = r2BodyName.encode(u'ascii')
    self.impl = c_qp.BilateralContact(r1Index, r2Index, r1BodyName, r2BodyName, ambiguityId, Vector3dVector(r1Points).v, Matrix3dVector(r1Frames).v, deref(X_b1_b2.impl), nrGen, mu, deref(X_b1_cf.impl))
  def __cidctor__(self, ContactId contactId, r1Points, r1Frames, PTransformd X_b1_b2, int nrGen, double mu, PTransformd X_b1_cf = PTransformd.Identity()):
    self.impl = c_qp.BilateralContact(contactId.impl, Vector3dVector(r1Points).v, Matrix3dVector(r1Frames).v, deref(X_b1_b2.impl), nrGen, mu, deref(X_b1_cf.impl))
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_qp.BilateralContact()
    elif check_args(args, [None, None, None, None, list, list, PTransformd, None, None]) or check_args(args, [None, None, None, None, list, list, PTransformd, None, None, PTransformd]):
      self.__noambidctor__(*args)
    elif check_args(args, [None, None, None, None, None, list, list, PTransformd, None, None]) or check_args(args, [None, None, None, None, list, list, PTransformd, None, None, PTransformd]):
      self.__ambidctor__(*args)
    elif check_args(args, [ContactId, list, list, PTransformd, None, None]) or check_args(args, [ContactId, list, list, PTransformd, None, None, PTransformd]):
      self.__cidctor__(*args)
  property contactId:
    def __get__(self):
      return ContactIdFromC(self.impl.contactId)
    def __set__(self, ContactId ci):
      self.impl.contactId = ci.impl
  property r1Points:
    def __get__(self):
      ret = []
      for v in self.impl.r1Points:
        ret.append(Vector3dFromC(v))
      return ret
    def __set__(self, value):
      self.impl.r1Points = Vector3dVector(value).v
  property r2Points:
    def __get__(self):
      ret = []
      for v in self.impl.r2Points:
        ret.append(Vector3dFromC(v))
      return ret
    def __set__(self, value):
      self.impl.r2Points = Vector3dVector(value).v
  property r1Cones:
    def __get__(self):
      ret = []
      for fc in self.impl.r1Cones:
        ret.append(FrictionConeFromC(fc))
      return ret
    def __set__(self, value):
      self.impl.r1Cones = FrictionConeVector(value).v
  property r2Cones:
    def __get__(self):
      ret = []
      for fc in self.impl.r2Cones:
        ret.append(FrictionConeFromC(fc))
      return ret
    def __set__(self, value):
      self.impl.r2Cones = FrictionConeVector(value).v
  property X_b1_b2:
    def __get__(self):
      return PTransformdFromC(self.impl.X_b1_b2, copy = False)
    def __set__(self, PTransformd value):
      self.impl.X_b1_b2 = deref(value.impl)
  property X_b1_cf:
    def __get__(self):
      return PTransformdFromC(self.impl.X_b1_cf, copy = False)
    def __set__(self, PTransformd value):
      self.impl.X_b1_cf = deref(value.impl)
  def __forceFC(self, VectorXd _lambda, c):
    return Vector3dFromC(self.impl.sForce(_lambda.impl, FrictionConeVector(c).v))
  def __forceIFC(self, VectorXd _lambda, int point, c):
    return Vector3dFromC(self.impl.sForce(_lambda.impl, point, FrictionConeVector(c).v))
  def __forceVFC(self, VectorXd _lambda, r_b_pi, c):
    return ForceVecdFromC(self.impl.sForce(_lambda.impl, Vector3dVector(r_b_pi).v, FrictionConeVector(c).v))
  def force(self, VectorXd _lambda, *args):
    if check_args(args, [list]):
      return self.__forceFC(_lambda, *args)
    elif check_args(args, [list, list]):
      return self.__forceVFC(_lambda, *args)
    elif check_args(args, [None, list]):
      return self.__forceIFC(_lambda, *args)
    else:
      raise TypeError("Wrong type passed for BilateralContact.force")
  def nrLambda(self, point = None):
    if point is None:
      return self.impl.nrLambda()
    else:
      return self.impl.sNrLambda(point)

cdef BilateralContactFromC(c_qp.BilateralContact c):
  cdef BilateralContact ret = BilateralContact()
  ret.impl = c
  return ret

cdef class BilateralContactVector(object):
  def __addUC(self, BilateralContact uc):
    self.v.push_back(uc.impl)
  def __cinit__(self, *args):
    if len(args) == 1 and isinstance(args[0], list):
      args = args[0]
    for uc in args:
      self.__addUC(uc)

cdef class SolverData(object):
  def nrVars(self):
    return self.impl.nrVars()
  def totalAlphaD(self):
    return self.impl.totalAlphaD()
  def totalLambda(self):
    return self.impl.totalLambda()
  def alphaD(self, int idx):
    return self.impl.alphaD(idx)
  def Lambda(self, int idx):
    return self.impl._lambda(idx)
  def alphaDBegin(self):
    return self.impl.alphaDBegin()
  def alphaDBegin(self, int idx):
    return self.impl.alphaDBegin(idx)
  def lambdaBegin(self):
    return self.impl.lambdaBegin()
  def lambdaBegin(self, int idx):
    return self.impl.lambdaBegin(idx)
  def nrUniLambda(self):
    return self.impl.nrUniLambda()
  def nrBiLambda(self):
    return self.impl.nrBiLambda()
  def unilateralBegin(self):
    return self.impl.unilateralBegin()
  def bilateralBegin(self):
    return self.impl.bilateralBegin()
  def nrContacts(self):
    return self.impl.nrContacts()
  def unilateralContacts(self):
    #FIXME Makes a copy of the vector because const iterator are not supported by cython 2.0?
    cdef vector[c_qp.UnilateralContact] v = self.impl.unilateralContacts()
    ret = []
    for c in v:
      ret.append(UnilateralContactFromC(c))
    return ret
  def bilateralContacts(self):
    #FIXME Makes a copy of the vector because const iterator are not supported by cython 2.0?
    cdef vector[c_qp.BilateralContact] v = self.impl.bilateralContacts()
    ret = []
    for c in v:
      ret.append(BilateralContactFromC(c))
    return ret
  def allContacts(self):
    #FIXME Makes a copy of the vector because const iterator are not supported by cython 2.0?
    cdef vector[c_qp.BilateralContact] v = self.impl.allContacts()
    ret = []
    for c in v:
      ret.append(BilateralContactFromC(c))
    return ret
  def computeNormalAccB(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs):
    self.impl.computeNormalAccB(deref(mbs.v), deref(mbcs.v))
  def normalAccB(self, int robotIndex):
    #FIXME Makes a copy of the vector because const iterator are not supported by cython 2.0?
    cdef vector[c_sva.MotionVecd] v = self.impl.normalAccB(robotIndex)
    ret = []
    for mv in v:
      ret.append(MotionVecdFromC(mv))
    return ret
cdef SolverData SolverDataFromC(const c_qp.SolverData& sd):
  cdef SolverData ret = SolverData()
  ret.impl = sd
  return ret

cdef class JointStiffness(object):
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_qp.JointStiffness()
    elif len(args) == 2:
      if isinstance(args[0], unicode):
        jName = args[0].encode(u'ascii')
      else:
        jName = args[0]
      self.impl = c_qp.JointStiffness(jName, args[1])
    else:
      raise TypeError("Wrong arguments passed to JointStiffness ctor")
  property jointName:
    def __get__(self):
      return self.impl.jointName
    def __set__(self, value):
      self.impl.jointName = value
  property stiffness:
    def __get__(self):
      return self.impl.stiffness
    def __set__(self, value):
      self.impl.stiffness = value

cdef class JointStiffnessVector(object):
  def __addJS(self, JointStiffness js):
    self.v.push_back(js.impl)
  def __cinit__(self, *args):
    if len(args) == 1 and isinstance(args[0], list):
      args = args[0]
    for js in args:
      self.__addJS(js)

cdef class JointGains(object):
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_qp.JointGains()
    elif len(args) == 2:
      self.impl = c_qp.JointGains(args[0], args[1])
    elif len(args) == 3:
      self.impl = c_qp.JointGains(args[0], args[1], args[2])
    else:
      raise TypeError("Wrong arguments passed to JointGains ctor")
  property jointName:
    def __get__(self):
      return self.impl.jointName
    def __set__(self, value):
      self.impl.jointName = value
  property stiffness:
    def __get__(self):
      return self.impl.stiffness
    def __set__(self, value):
      self.impl.stiffness = value
  property damping:
    def __get__(self):
      return self.impl.damping
    def __set__(self, value):
      self.impl.damping = value

cdef class JointGainsVector(object):
  def __addJG(self, JointGains jg):
    self.v.push_back(jg.impl)
  def __cinit__(self, *args):
    if len(args) == 1 and isinstance(args[0], list):
      args = args[0]
    for jg in args:
      self.__addJG(jg)

cdef class SpringJoint(object):
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = c_qp.SpringJoint()
    elif len(args) == 4:
      self.impl = c_qp.SpringJoint(args[0], args[1], args[2], args[3])
    else:
      raise TypeError("Wrong arguments passed to SpringJoint ctor")
  property jointName:
    def __get__(self):
      return self.impl.jointName
    def __set__(self, value):
      self.impl.jointName = value
  property K:
    def __get__(self):
      return self.impl.K
    def __set__(self, value):
      self.impl.K = value
  property C:
    def __get__(self):
      return self.impl.C
    def __set__(self, value):
      self.impl.C = value
  property O:
    def __get__(self):
      return self.impl.O
    def __set__(self, value):
      self.impl.O = value

cdef class SpringJointVector(object):
  def __addSJ(self, SpringJoint sj):
    self.v.push_back(sj.impl)
  def __cinit__(self, *args):
    if len(args) == 1 and isinstance(args[0], list):
      args = args[0]
    for sj in args:
      self.__addSJ(sj)

cdef class Constraint(object):
  def updateNrVars(self, MultiBodyVector mbs, SolverData sd):
    self.constraint_base.updateNrVars(deref(mbs.v), sd.impl)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData sd):
    self.constraint_base.update(deref(mbs.v), deref(mbcs.v), sd.impl)

cdef class Equality(Constraint):
  def maxEq(self):
    return self.eq_base.maxEq()
  def nrEq(self):
    return self.eq_base.nrEq()
  def AEq(self):
    return MatrixXdFromC(self.eq_base.AEq())
  def bEq(self):
    return VectorXdFromC(self.eq_base.bEq())

cdef class Inequality(Constraint):
  def maxInEq(self):
    return self.ineq_base.maxInEq()
  def nrInEq(self):
    return self.ineq_base.nrInEq()
  def AInEq(self):
    return MatrixXdFromC(self.ineq_base.AInEq())
  def bInEq(self):
    return VectorXdFromC(self.ineq_base.bInEq())

cdef class GenInequality(Constraint):
  def maxGenInEq(self):
    return self.genineq_base.maxGenInEq()
  def nrGenInEq(self):
    return self.genineq_base.nrGenInEq()
  def AGenInEq(self):
    return MatrixXdFromC(self.genineq_base.AGenInEq())
  def LowerGenInEq(self):
    return VectorXdFromC(self.genineq_base.LowerGenInEq())
  def UpperGenInEq(self):
    return VectorXdFromC(self.genineq_base.UpperGenInEq())

cdef class Bound(Constraint):
  def beginVar(self):
    return self.bound_base.beginVar()
  def Lower(self):
    return MatrixXdFromC(self.bound_base.Lower())
  def Upper(self):
    return VectorXdFromC(self.bound_base.Upper())

cdef class Task(object):
  def weight(self, weight = None):
    if weight is None:
      return self.base.weight()
    else:
      self.base.weight(weight)

cdef class HighLevelTask(object):
  def dim(self):
    return self.base.dim()
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData sd):
    self.base.update(deref(mbs.v), deref(mbcs.v), sd.impl)
  def jac(self):
    return MatrixXdFromC(self.base.jac())
  def eval(self):
    return VectorXdFromC(self.base.eval())
  def speed(self):
    return VectorXdFromC(self.base.speed())
  def normalAcc(self):
    return VectorXdFromC(self.base.normalAcc())

cdef class PositionTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, Vector3d pos, Vector3d bodyPoint = Vector3d.Zero()):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    self.impl = new c_qp.PositionTask(deref(mbs.v), robotIndex, bodyName, pos.impl, bodyPoint.impl)
    self.base = self.impl
  def position(self, Vector3d pos = None):
    if pos is None:
      return Vector3dFromC(self.impl.position())
    else:
      self.impl.position(pos.impl)
  def bodyPoint(self, Vector3d point = None):
    if point is None:
      return Vector3dFromC(self.impl.bodyPoint())
    else:
      self.impl.bodyPoint(point.impl)

cdef class OrientationTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __m3ctor__(self, MultiBodyVector mbs, int robotIndex, bodyName, Matrix3d ori):
    self.impl = new c_qp.OrientationTask(deref(mbs.v), robotIndex, <string>bodyName, ori.impl)
  def __qctor__(self, MultiBodyVector mbs, int robotIndex, bodyName, Quaterniond ori):
    self.impl = new c_qp.OrientationTask(deref(mbs.v), robotIndex, <string>bodyName, ori.impl)
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, ori):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    if isinstance(ori, Matrix3d):
      self.__m3ctor__(mbs, robotIndex, bodyName, ori)
    elif isinstance(ori, Quaterniond):
      self.__qctor__(mbs, robotIndex, bodyName, ori)
    else:
      raise TypeError("Wrong arguments passed to OrientationTask ctor")
    self.base = self.impl
  def __m3orientation__(self, Matrix3d ori):
    self.impl.orientation(ori.impl)
  def __qorientation__(self, Quaterniond ori):
    self.impl.orientation(ori.impl)
  def orientation(self, ori = None):
    if ori is None:
      return Matrix3dFromC(self.impl.orientation())
    else:
      if isinstance(ori, Matrix3d):
        self.__m3orientation__(ori)
      elif isinstance(ori, Quaterniond):
        self.__qorientation__(ori)
      else:
        raise TypeError("Wrong arguments passed to OrientationTask.orientation")

cdef class TransformTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, PTransformd X_0_t, PTransformd X_b_p = PTransformd.Identity(), Matrix3d X_0_c = Matrix3d.Identity()):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    self.impl = new c_qp.TransformTask(deref(mbs.v), robotIndex, bodyName, deref(X_0_t.impl), deref(X_b_p.impl), X_0_c.impl)
    self.base = self.impl
  def E_0_c(self, Matrix3d E_0_c = None):
    if E_0_c is None:
      return Matrix3dFromC(self.impl.E_0_c())
    else:
      self.impl.E_0_c(E_0_c.impl)
  # Common to *TransformTask
  def target(self, PTransformd target = None):
    if target is None:
      return PTransformdFromC(self.impl.target())
    else:
      self.impl.target(deref(target.impl))
  def X_b_p(self, PTransformd X_b_p = None):
    if X_b_p is None:
      return PTransformdFromC(self.impl.X_b_p())
    else:
      self.impl.X_b_p(deref(X_b_p.impl))

cdef class SurfaceTransformTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, PTransformd X_0_t, PTransformd X_b_p = PTransformd.Identity()):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    self.impl = new c_qp.SurfaceTransformTask(deref(mbs.v), robotIndex, bodyName, deref(X_0_t.impl), deref(X_b_p.impl))
    self.base = self.impl
  # Common to *TransformTask
  def target(self, PTransformd target = None):
    if target is None:
      return PTransformdFromC(self.impl.target())
    else:
      self.impl.target(deref(target.impl))
  def X_b_p(self, PTransformd X_b_p = None):
    if X_b_p is None:
      return PTransformdFromC(self.impl.X_b_p())
    else:
      self.impl.X_b_p(deref(X_b_p.impl))

cdef class SurfaceOrientationTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __qctor__(self, MultiBodyVector mbs, int robotIndex, bodyName, Quaterniond ori, PTransformd X_b_s):
    self.impl = new c_qp.SurfaceOrientationTask(deref(mbs.v), robotIndex, <string>bodyName, ori.impl, deref(X_b_s.impl))
  def __m3ctor__(self, MultiBodyVector mbs, int robotIndex, bodyName, Matrix3d ori, PTransformd X_b_s):
    self.impl = new c_qp.SurfaceOrientationTask(deref(mbs.v), robotIndex, <string>bodyName, ori.impl, deref(X_b_s.impl))
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, ori, PTransformd X_b_s):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    if isinstance(ori, Quaterniond):
      self.__qctor__(mbs, robotIndex, bodyName, ori, X_b_s)
    elif isinstance(ori, Matrix3d):
      self.__m3ctor__(mbs, robotIndex, bodyName, ori, X_b_s)
    else:
      raise TypeError("Wrong arguments given to SurfaceOrientationTask ctor")
    self.base = self.impl
  def __qorientation__(self, Quaterniond ori):
    self.impl.orientation(ori.impl)
  def __m3orientation__(self, Matrix3d ori):
    self.impl.orientation(ori.impl)
  def orientation(self, ori = None):
    if ori is None:
      return Matrix3dFromC(self.impl.orientation())
    else:
      if isinstance(ori, Quaterniond):
        self.__qorientation__(ori)
      elif isinstance(ori, Matrix3d):
        self.__m3orientation__(ori)
      else:
        raise TypeError("Wrong arguments passed to SurfaceOrientationTask.orientation")

cdef class GazeTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __v2ctor__(self, MultiBodyVector mbs, int rIndex, bName, Vector2d point2d, double depthEstimate, PTransformd X_b_gaze, Vector2d point2d_ref = Vector2d.Zero()):
    self.impl = new c_qp.GazeTask(deref(mbs.v), rIndex, bName, point2d.impl, depthEstimate, deref(X_b_gaze.impl), point2d_ref.impl)
  def __v3ctor__(self, MultiBodyVector mbs, int rIndex, bName, Vector3d point3d, PTransformd X_b_gaze, Vector2d point2d_ref = Vector2d.Zero()):
    self.impl = new c_qp.GazeTask(deref(mbs.v), rIndex, bName, point3d.impl, deref(X_b_gaze.impl), point2d_ref.impl)
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, *args):
    if check_args(args, [Vector2d, None, PTransformd]) or check_args(args, [Vector2d, None, PTransformd, Vector2d]):
      self.__v2ctor__(mbs, robotIndex, bodyName, *args)
    elif check_args(args, [Vector3d, PTransformd]) or check_args(args, [Vector3d, PTransformd, Vector3d]):
      self.__v3ctor__(mbs, robotIndex, bodyName, *args)
    else:
      raise TypeError("Wrong arguments passed to Gaze ctor")
    self.base = self.impl
  def __v2error__(self, Vector2d point, Vector2d point2d_ref):
    self.impl.error(point.impl, point2d_ref.impl)
  def __v3error__(self, Vector3d point, Vector2d point2d_ref):
    self.impl.error(point.impl, point2d_ref.impl)
  def error(self, point, Vector2d point2d_ref = Vector2d.Zero()):
    if isinstance(point, Vector2d):
      self.__v2error__(point, point2d_ref)
    elif isinstance(point, Vector3d):
      self.__v3error__(point, point2d_ref)
    else:
      raise TypeError("Wrong arguments passed to Gaze.error")

cdef class PositionBasedVisServoTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, PTransformd X_t_s, PTransformd X_b_s = PTransformd.Identity()):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    self.impl = new c_qp.PositionBasedVisServoTask(deref(mbs.v), robotIndex, bodyName, deref(X_t_s.impl), deref(X_b_s.impl))
    self.base = self.impl
  def error(self, PTransformd X_t_s):
    self.impl.error(deref(X_t_s.impl))

cdef class PostureTask(Task):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, q, double stiffness, double weight, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.PostureTask(deref(mbs.v), robotIndex, q, stiffness, weight)
      self.base = self.impl
  def stiffness(self, s = None):
    if s is None:
      return self.impl.stiffness()
    else:
      self.impl.stiffness(s)
  def damping(self):
    return self.impl.damping()
  def gains(self, double stiffness, damping = None):
    if damping is None:
      self.impl.gains(stiffness)
    else:
      self.impl.gains(stiffness, damping)
  def posture(self, q = None):
    if q is None:
      return self.impl.posture()
    else:
      self.impl.posture(q)
  def jointsStiffness(self, MultiBodyVector mbs, jss):
    self.impl.jointsStiffness(deref(mbs.v), JointStiffnessVector(jss).v)
  def jointsGains(self, MultiBodyVector mbs, jgs):
    self.impl.jointsGains(deref(mbs.v), JointGainsVector(jgs).v)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData data):
    self.impl.ptUpdate(deref(mbs.v), deref(mbcs.v), data.impl)
  def eval(self):
    return VectorXdFromC(self.impl.ptEval())

cdef PostureTask PostureTaskFromPtr(c_qp.PostureTask * p):
    cdef PostureTask ret = PostureTask(None, 0, None, 0, 0, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = ret.base = p
    return ret

cdef class CoMTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, Vector3d com, weight = None):
    if weight is None:
      self.impl = new c_qp.CoMTask(deref(mbs.v), robotIndex, com.impl)
    else:
      self.impl = new c_qp.CoMTask(deref(mbs.v), robotIndex, com.impl, weight)
    self.base = self.impl
  def com(self, Vector3d com = None):
    if com is None:
      return Vector3dFromC(self.impl.com())
    else:
      self.impl.com(com.impl)
  def updateInertialParameters(self, MultiBodyVector mbs):
    self.impl.updateInertialParameters(deref(mbs.v))

cdef class MomentumTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, ForceVecd mom):
    self.impl = new c_qp.MomentumTask(deref(mbs.v), robotIndex, deref(mom.impl))
    self.base = self.impl
  def momentum(self, ForceVecd mom = None):
    if mom is None:
      return ForceVecdFromC(self.impl.momentum())
    else:
      self.impl.momentum(deref(mom.impl))

cdef class LinVelocityTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, Vector3d pos, Vector3d bodyPoint = Vector3d.Zero()):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    self.impl = new c_qp.LinVelocityTask(deref(mbs.v), robotIndex, bodyName, pos.impl, bodyPoint.impl)
    self.base = self.impl
  def velocity(self, Vector3d pos = None):
    if pos is None:
      return Vector3dFromC(self.impl.velocity())
    else:
      self.impl.velocity(pos.impl)
  def bodyPoint(self, Vector3d point = None):
    if point is None:
      return Vector3dFromC(self.impl.bodyPoint())
    else:
      self.impl.bodyPoint(point.impl)

cdef class OrientationTrackingTask(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bodyName, Vector3d bodyPoint, Vector3d bodyAxis, trackingJointName, Vector3d trackedPoint):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    self.impl = new c_qp.OrientationTrackingTask(deref(mbs.v), robotIndex, bodyName, bodyPoint.impl, bodyAxis.impl, trackingJointName, trackedPoint.impl)
    self.base = self.impl
  def trackedPoint(self, Vector3d pt = None):
    if pt is None:
      return Vector3dFromC(self.impl.trackedPoint())
    else:
      self.impl.trackedPoint(pt.impl)
  def bodyPoint(self, Vector3d pt = None):
    if pt is None:
      return Vector3dFromC(self.impl.bodyPoint())
    else:
      self.impl.bodyPoint(pt.impl)
  def bodyAxis(self, Vector3d pt = None):
    if pt is None:
      return Vector3dFromC(self.impl.bodyAxis())
    else:
      self.impl.bodyAxis(pt.impl)

class JointsSelector_SelectedData(object):
  def __init__(self, posInDof, dof):
    self.posInDof = posInDof
    self.dof = dof

cdef class JointsSelector(HighLevelTask):
  def __dealloc__(self):
    del self.impl
  def __ctor__(self, MultiBodyVector mbs, int robotIndex, HighLevelTask hl, selectedJointsId):
    self.impl = new c_qp.JointsSelector(deref(mbs.v), robotIndex, hl.base, selectedJointsId)
  def __cinit__(self, *args):
    if len(args) == 0:
      self.impl = NULL
    elif check_args(args, [None, None, HighLevelTask, list]):
      self.__ctor__(*args)
    else:
      raise TypeError("Wrong arguments for JointsSelector ctor")
    self.base = self.impl
  def selectedJoints(self):
    cdef vector[c_qp.JointsSelector_SelectedData] sdv = self.impl.selectedJoints()
    ret = []
    for sd in sdv:
      ret.append(JointsSelector_SelectedData(sd.posInDof, sd.dof))
    return ret
  @staticmethod
  def ActiveJoints(MultiBodyVector mbs, int robotIndex, HighLevelTask hl, activeJointsId):
    cdef JointsSelector ret = JointsSelector()
    for i, n in enumerate(activeJointsId):
      if isinstance(n, unicode):
        activeJointsId[i] = n.encode(u'ascii')
    ret.impl = c_qp_private.ActiveJoints2Ptr(deref(mbs.v), robotIndex, hl.base, activeJointsId)
    ret.base = ret.impl
    return ret
  @staticmethod
  def UnactiveJoints(MultiBodyVector mbs, int robotIndex, HighLevelTask hl, unactiveJointsId):
    cdef JointsSelector ret = JointsSelector()
    for i, n in enumerate(unactiveJointsId):
      if isinstance(n, unicode):
        unactiveJointsId[i] = n.encode(u'ascii')
    ret.impl = c_qp_private.UnactiveJoints2Ptr(deref(mbs.v), robotIndex, hl.base, unactiveJointsId)
    ret.base = ret.impl
    return ret

cdef class SetPointTask(Task):
  def __dealloc__(self):
    del self.impl
  def __PositionTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionTask pt, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, pt.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, pt.impl, stiffness, dimWeight.impl, weight)
  def __OrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __SurfaceOrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceOrientationTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __GazeTaskctor__(self, MultiBodyVector mbs, int rIdx, GazeTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __PositionBasedVisServoTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionBasedVisServoTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __CoMTaskctor__(self, MultiBodyVector mbs, int rIdx, CoMTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __LinVelocityTaskctor__(self, MultiBodyVector mbs, int rIdx, LinVelocityTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __OrientationTrackingTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTrackingTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __MomentumTaskctor__(self, MultiBodyVector mbs, int rIdx, MomentumTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __JointsSelectorctor__(self, MultiBodyVector mbs, int rIdx, JointsSelector task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __TransformTaskctor__(self, MultiBodyVector mbs, int rIdx, TransformTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __SurfaceTransformTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceTransformTask task, double stiffness, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, weight)
    else:
      self.impl = new c_qp.SetPointTask(deref(mbs.v), rIdx, task.impl, stiffness, dimWeight.impl, weight)
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, task, double stiffness, *args):
    cdef VectorXd dimWeight = None
    cdef double weight = 0
    if check_args(args, [VectorXd, None]):
      dimWeight = args[0]
      weight = args[1]
    elif check_args(args, [None]):
      weight = args[0]
    else:
      raise TypeError("Wrong types passed to SetPointTask ctor")
    if isinstance(task, PositionTask):
      self.__PositionTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, OrientationTask):
      self.__OrientationTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, SurfaceOrientationTask):
      self.__SurfaceOrientationTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, GazeTask):
      self.__GazeTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, PositionBasedVisServoTask):
      self.__PositionBasedVisServoTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, CoMTask):
      self.__CoMTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, LinVelocityTask):
      self.__LinVelocityTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, OrientationTrackingTask):
      self.__OrientationTrackingTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, MomentumTask):
      self.__MomentumTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, JointsSelector):
      self.__JointsSelectorctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, TransformTask):
      self.__TransformTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    elif isinstance(task, SurfaceTransformTask):
      self.__SurfaceTransformTaskctor__(mbs, robotIndex, task, stiffness, weight, dimWeight)
    else:
      raise TypeError("Failed to create SetPointTask, cannot handle this kind of Task")
    self.base = self.impl
  def stiffness(self, s = None):
    if s is None:
      return self.impl.stiffness()
    else:
      self.impl.stiffness(s)
  # SetPointTaskCommon
  def dimWeight(self, VectorXd v = None):
    if v is None:
      return VectorXdFromC(self.impl.dimWeight())
    else:
      self.impl.dimWeight(v.impl)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData data):
    self.impl.update(deref(mbs.v), deref(mbcs.v), data.impl)
  def Q(self):
    return MatrixXdFromC(self.impl.Q())
  def C(self):
    return VectorXdFromC(self.impl.C())

cdef class TrackingTask(Task):
  def __dealloc__(self):
    del self.impl
  def __PositionTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionTask pt, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, pt.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, pt.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __OrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __SurfaceOrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceOrientationTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __GazeTaskctor__(self, MultiBodyVector mbs, int rIdx, GazeTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __PositionBasedVisServoTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionBasedVisServoTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __CoMTaskctor__(self, MultiBodyVector mbs, int rIdx, CoMTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __LinVelocityTaskctor__(self, MultiBodyVector mbs, int rIdx, LinVelocityTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __OrientationTrackingTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTrackingTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __MomentumTaskctor__(self, MultiBodyVector mbs, int rIdx, MomentumTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __JointsSelectorctor__(self, MultiBodyVector mbs, int rIdx, JointsSelector task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __TransformTaskctor__(self, MultiBodyVector mbs, int rIdx, TransformTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __SurfaceTransformTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceTransformTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrackingTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, task, double gainPos, double gainVel, *args):
    cdef VectorXd dimWeight = None
    cdef double weight = 0
    if check_args(args, [VectorXd, None]):
      dimWeight = args[0]
      weight = args[1]
    elif check_args(args, [None]):
      weight = args[0]
    else:
      raise TypeError("Wrong types passed to TrackingTask ctor")
    if isinstance(task, PositionTask):
      self.__PositionTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, OrientationTask):
      self.__OrientationTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, SurfaceOrientationTask):
      self.__SurfaceOrientationTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, GazeTask):
      self.__GazeTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, PositionBasedVisServoTask):
      self.__PositionBasedVisServoTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, CoMTask):
      self.__CoMTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, LinVelocityTask):
      self.__LinVelocityTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, OrientationTrackingTask):
      self.__OrientationTrackingTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, MomentumTask):
      self.__MomentumTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, JointsSelector):
      self.__JointsSelectorctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, TransformTask):
      self.__TransformTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, SurfaceTransformTask):
      self.__SurfaceTransformTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    else:
      raise TypeError("Failed to create TrackingTask, cannot handle this kind of Task")
    self.base = self.impl
  # SetPointTaskCommon
  def dimWeight(self, VectorXd v = None):
    if v is None:
      return VectorXdFromC(self.impl.dimWeight())
    else:
      self.impl.dimWeight(v.impl)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData data):
    self.impl.update(deref(mbs.v), deref(mbcs.v), data.impl)
  def Q(self):
    return MatrixXdFromC(self.impl.Q())
  def C(self):
    return VectorXdFromC(self.impl.C())
  def setGains(self, double gainPos, double gainVel):
    self.impl.setGains(gainPos, gainVel)
  def errorPos(self, VectorXd err):
    self.impl.errorPos(err.impl)
  def errorVel(self, VectorXd err):
    self.impl.errorVel(err.impl)
  def refAccel(self, VectorXd acc):
    self.impl.refAccel(acc.impl)

cdef class TrajectoryTask(Task):
  def __dealloc__(self):
    del self.impl
  def __PositionTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionTask pt, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, pt.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, pt.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __OrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __SurfaceOrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceOrientationTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __GazeTaskctor__(self, MultiBodyVector mbs, int rIdx, GazeTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __PositionBasedVisServoTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionBasedVisServoTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __CoMTaskctor__(self, MultiBodyVector mbs, int rIdx, CoMTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __LinVelocityTaskctor__(self, MultiBodyVector mbs, int rIdx, LinVelocityTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __OrientationTrackingTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTrackingTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __MomentumTaskctor__(self, MultiBodyVector mbs, int rIdx, MomentumTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __JointsSelectorctor__(self, MultiBodyVector mbs, int rIdx, JointsSelector task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __TransformTaskctor__(self, MultiBodyVector mbs, int rIdx, TransformTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __SurfaceTransformTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceTransformTask task, double gainPos, double gainVel, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, weight)
    else:
      self.impl = new c_qp.TrajectoryTask(deref(mbs.v), rIdx, task.impl, gainPos, gainVel, dimWeight.impl, weight)
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, task, double gainPos, double gainVel, *args):
    cdef VectorXd dimWeight = None
    cdef double weight = 0
    if check_args(args, [VectorXd, None]):
      dimWeight = args[0]
      weight = args[1]
    elif check_args(args, [None]):
      weight = args[0]
    else:
      raise TypeError("Wrong types passed to TrajectoryTask ctor")
    if isinstance(task, PositionTask):
      self.__PositionTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, OrientationTask):
      self.__OrientationTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, SurfaceOrientationTask):
      self.__SurfaceOrientationTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, GazeTask):
      self.__GazeTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, PositionBasedVisServoTask):
      self.__PositionBasedVisServoTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, CoMTask):
      self.__CoMTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, LinVelocityTask):
      self.__LinVelocityTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, OrientationTrackingTask):
      self.__OrientationTrackingTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, MomentumTask):
      self.__MomentumTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, JointsSelector):
      self.__JointsSelectorctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, TransformTask):
      self.__TransformTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    elif isinstance(task, SurfaceTransformTask):
      self.__SurfaceTransformTaskctor__(mbs, robotIndex, task, gainPos, gainVel, weight, dimWeight)
    else:
      raise TypeError("Failed to create TrajectoryTask, cannot handle this kind of Task")
    self.base = self.impl
  # SetPointTaskCommon
  def dimWeight(self, VectorXd v = None):
    if v is None:
      return VectorXdFromC(self.impl.dimWeight())
    else:
      self.impl.dimWeight(v.impl)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData data):
    self.impl.update(deref(mbs.v), deref(mbcs.v), data.impl)
  def Q(self):
    return MatrixXdFromC(self.impl.Q())
  def C(self):
    return VectorXdFromC(self.impl.C())
  def setGains(self, double gainPos, double gainVel):
    self.impl.setGains(gainPos, gainVel)
  def refVel(self, VectorXd vel):
    self.impl.refVel(vel.impl)
  def refAccel(self, VectorXd acc):
    self.impl.refAccel(acc.impl)

cdef class TargetObjectiveTask(Task):
  def __dealloc__(self):
    del self.impl
  def __PositionTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionTask pt, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, pt.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, pt.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __OrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __SurfaceOrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceOrientationTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __GazeTaskctor__(self, MultiBodyVector mbs, int rIdx, GazeTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __PositionBasedVisServoTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionBasedVisServoTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __CoMTaskctor__(self, MultiBodyVector mbs, int rIdx, CoMTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __LinVelocityTaskctor__(self, MultiBodyVector mbs, int rIdx, LinVelocityTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __OrientationTrackingTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTrackingTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __MomentumTaskctor__(self, MultiBodyVector mbs, int rIdx, MomentumTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __JointsSelectorctor__(self, MultiBodyVector mbs, int rIdx, JointsSelector task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __TransformTaskctor__(self, MultiBodyVector mbs, int rIdx, TransformTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __SurfaceTransformTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceTransformTask task, double timeStep, double duration, VectorXd objDot, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, weight)
    else:
      self.impl = new c_qp.TargetObjectiveTask(deref(mbs.v), rIdx, task.impl, timeStep, duration, objDot.impl, dimWeight.impl, weight)
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, task, double timeStep, double duration, VectorXd objDot, *args):
    cdef VectorXd dimWeight = None
    cdef double weight = 0
    if check_args(args, [VectorXd, None]):
      dimWeight = args[0]
      weight = args[1]
    elif check_args(args, [None]):
      weight = args[0]
    else:
      raise TypeError("Wrong types passed to TargetObjectiveTask ctor")
    if isinstance(task, PositionTask):
      self.__PositionTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, OrientationTask):
      self.__OrientationTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, SurfaceOrientationTask):
      self.__SurfaceOrientationTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, GazeTask):
      self.__GazeTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, PositionBasedVisServoTask):
      self.__PositionBasedVisServoTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, CoMTask):
      self.__CoMTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, LinVelocityTask):
      self.__LinVelocityTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, OrientationTrackingTask):
      self.__OrientationTrackingTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, MomentumTask):
      self.__MomentumTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, JointsSelector):
      self.__JointsSelectorctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, TransformTask):
      self.__TransformTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    elif isinstance(task, SurfaceTransformTask):
      self.__SurfaceTransformTaskctor__(mbs, robotIndex, task, timeStep, duration, objDot, weight, dimWeight)
    else:
      raise TypeError("Failed to create TargetObjectiveTask, cannot handle this kind of Task")
    self.base = self.impl
  # SetPointTaskCommon
  def dimWeight(self, VectorXd v = None):
    if v is None:
      return VectorXdFromC(self.impl.dimWeight())
    else:
      self.impl.dimWeight(v.impl)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData data):
    self.impl.update(deref(mbs.v), deref(mbcs.v), data.impl)
  def Q(self):
    return MatrixXdFromC(self.impl.Q())
  def C(self):
    return VectorXdFromC(self.impl.C())
  def duration(self, d = None):
    if d is None:
      return self.impl.duration()
    else:
      self.impl.duration(d)
  def iter(self, d = None):
    if d is None:
      return self.impl.iter()
    else:
      self.impl.iter(d)
  def nrIter(self, d = None):
    if d is None:
      return self.impl.nrIter()
    else:
      self.impl.nrIter(d)
  def objDot(self, VectorXd objDot = None):
    if objDot is None:
      return VectorXdFromC(self.impl.objDot())
    else:
      self.impl.objDot(objDot.impl)
  def phi(self):
    return VectorXdFromC(self.impl.phi())
  def psi(self):
    return VectorXdFromC(self.impl.psi())

cdef class PIDTask(Task):
  def __dealloc__(self):
    del self.impl
  def __PositionTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionTask pt, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, pt.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, pt.impl, P, I, D, dimWeight.impl, weight)
  def __OrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __SurfaceOrientationTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceOrientationTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __GazeTaskctor__(self, MultiBodyVector mbs, int rIdx, GazeTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __PositionBasedVisServoTaskctor__(self, MultiBodyVector mbs, int rIdx, PositionBasedVisServoTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __CoMTaskctor__(self, MultiBodyVector mbs, int rIdx, CoMTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __LinVelocityTaskctor__(self, MultiBodyVector mbs, int rIdx, LinVelocityTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __OrientationTrackingTaskctor__(self, MultiBodyVector mbs, int rIdx, OrientationTrackingTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __MomentumTaskctor__(self, MultiBodyVector mbs, int rIdx, MomentumTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __JointsSelectorctor__(self, MultiBodyVector mbs, int rIdx, JointsSelector task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __TransformTaskctor__(self, MultiBodyVector mbs, int rIdx, TransformTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __SurfaceTransformTaskctor__(self, MultiBodyVector mbs, int rIdx, SurfaceTransformTask task, double P, double I, double D, double weight, VectorXd dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, weight)
    else:
      self.impl = new c_qp.PIDTask(deref(mbs.v), rIdx, task.impl, P, I, D, dimWeight.impl, weight)
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, task, double P, double I, double D, *args):
    cdef VectorXd dimWeight = None
    cdef double weight = 0
    if check_args(args, [VectorXd, None]):
      dimWeight = args[0]
      weight = args[1]
    elif check_args(args, [None]):
      weight = args[0]
    else:
      raise TypeError("Wrong types passed to PIDTask ctor")
    if isinstance(task, PositionTask):
      self.__PositionTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, OrientationTask):
      self.__OrientationTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, SurfaceOrientationTask):
      self.__SurfaceOrientationTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, GazeTask):
      self.__GazeTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, PositionBasedVisServoTask):
      self.__PositionBasedVisServoTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, CoMTask):
      self.__CoMTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, LinVelocityTask):
      self.__LinVelocityTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, OrientationTrackingTask):
      self.__OrientationTrackingTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, MomentumTask):
      self.__MomentumTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, JointsSelector):
      self.__JointsSelectorctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, TransformTask):
      self.__TransformTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    elif isinstance(task, SurfaceTransformTask):
      self.__SurfaceTransformTaskctor__(mbs, robotIndex, task, P, I, D, weight, dimWeight)
    else:
      raise TypeError("Failed to create PIDTask, cannot handle this kind of Task")
    self.base = self.impl
  # SetPointTaskCommon
  def dimWeight(self, VectorXd v = None):
    if v is None:
      return VectorXdFromC(self.impl.dimWeight())
    else:
      self.impl.dimWeight(v.impl)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData data):
    self.impl.update(deref(mbs.v), deref(mbcs.v), data.impl)
  def Q(self):
    return MatrixXdFromC(self.impl.Q())
  def C(self):
    return VectorXdFromC(self.impl.C())
  def P(self, P = None):
    if P is None:
      return self.impl.P()
    else:
      self.impl.P(P)
  def I(self, I = None):
    if I is None:
      return self.impl.I()
    else:
      self.impl.I(I)
  def D(self, D = None):
    if D is None:
      return self.impl.D()
    else:
      self.impl.D(D)
  def error(self, VectorXd error):
    self.impl.error(error.impl)
  def errorI(self, VectorXd errorI):
    self.impl.errorI(errorI.impl)
  def errorD(self, VectorXd errorD):
    self.impl.errorD(errorD.impl)

cdef class MultiCoMTask(Task):
  def __dealloc__(self):
    del self.impl
  def __ctor__(self, MultiBodyVector mbs, robotIndexes, Vector3d com, double stiffness, double weight, Vector3d dimWeight = None):
    if dimWeight is None:
      self.impl = new c_qp.MultiCoMTask(deref(mbs.v), robotIndexes, com.impl, stiffness, weight)
    else:
      self.impl = new c_qp.MultiCoMTask(deref(mbs.v), robotIndexes, com.impl, stiffness, dimWeight.impl, weight)
  def __cinit__(self, MultiBodyVector mbs, robotIndexes, Vector3d com, double stiffness, *args):
    cdef VectorXd dimWeight = None
    cdef double weight = 0
    if check_args(args, [VectorXd, None]):
      dimWeight = args[0]
      weight = args[1]
    elif check_args(args, [None]):
      weight = args[0]
    else:
      raise TypeError("Wrong types passed to MultiCoMTask ctor")
    self.__ctor__(mbs, robotIndexes, com, stiffness, weight, dimWeight)
    self.base = self.impl
  def com(self, Vector3d com = None):
    if com is None:
      return Vector3dFromC(self.impl.com())
    else:
      self.impl.com(com.impl)
  def updateInertialParameters(self, MultiBodyVector mbs):
    self.impl.updateInertialParameters(deref(mbs.v))
  def stiffness(self, s = None):
    if s is None:
      return self.impl.stiffness()
    else:
      self.impl.stiffness(s)
  def eval(self):
    return VectorXdFromC(self.impl.eval())
  def speed(self):
    return VectorXdFromC(self.impl.speed())
  def dimWeight(self, VectorXd v = None):
    if v is None:
      return VectorXdFromC(self.impl.dimWeight())
    else:
      self.impl.dimWeight(v.impl)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData data):
    self.impl.update(deref(mbs.v), deref(mbcs.v), data.impl)

cdef class MultiRobotTransformTask(Task):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, r1Index, r2Index, r1BodyIndex, r2BodyIndex, PTransformd X_r1b_r1s, PTransformd X_r2b_r2s, double stiffness, double weight):
    if isinstance(r1BodyIndex, unicode):
      r1BodyIndex = r1BodyIndex.encode(u'ascii')
    if isinstance(r2BodyIndex, unicode):
      r2BodyIndex = r2BodyIndex.encode(u'ascii')
    self.impl = new c_qp.MultiRobotTransformTask(deref(mbs.v), r1Index, r2Index, r1BodyIndex, r2BodyIndex, deref(X_r1b_r1s.impl), deref(X_r2b_r2s.impl), stiffness, weight)
    self.base = self.impl
  def X_r1b_r1s(self, PTransformd X_r1b_r1s = None):
    if X_r1b_r1s is None:
      return PTransformdFromC(self.impl.X_r1b_r1s())
    else:
      self.impl.X_r1b_r1s(deref(X_r1b_r1s.impl))
  def X_r2b_r2s(self, PTransformd X_r2b_r2s = None):
    if X_r2b_r2s is None:
      return PTransformdFromC(self.impl.X_r2b_r2s())
    else:
      self.impl.X_r2b_r2s(deref(X_r2b_r2s.impl))
  def stiffness(self, s = None):
    if s is None:
      return self.impl.stiffness()
    else:
      self.impl.stiffness(s)
  def eval(self):
    return VectorXdFromC(self.impl.eval())
  def speed(self):
    return VectorXdFromC(self.impl.speed())
  def dimWeight(self, VectorXd v = None):
    if v is None:
      return VectorXdFromC(self.impl.dimWeight())
    else:
      self.impl.dimWeight(v.impl)
  def update(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs, SolverData data):
    self.impl.update(deref(mbs.v), deref(mbcs.v), data.impl)

cdef class ContactTask(Task):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, ContactId cId, double stiffness, double weight):
    self.impl = new c_qp.ContactTask(cId.impl, stiffness, weight)
    self.base = self.impl
  def error(self, Vector3d err):
    self.impl.error(err.impl)
  def errorD(self, Vector3d errD):
    self.impl.errorD(errD.impl)

cdef class GripperTorqueTask(Task):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, ContactId cId, Vector3d origin, Vector3d axis, double weight):
    self.impl = new c_qp.GripperTorqueTask(cId.impl, origin.impl, axis.impl, weight)
    self.base = self.impl

cdef class MotionConstr(GenInequality):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, tasks.TorqueBound tb, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.MotionConstr(deref(mbs.v), robotIndex, tb.impl)
      self.cf_base = self.impl
      self.genineq_base = self.impl
      self.constraint_base = self.impl
  # Motion default
  def computeTorque(self, VectorXd alphaD, VectorXd _lambda):
    self.impl.computeTorque(alphaD.impl, _lambda.impl)
  def torque(self, MultiBodyVector mbs = None, MultiBodyConfigVector mbcs = None):
    if mbs is None and mbcs is None:
      VectorXdFromC(self.impl.torque())
    elif mbs is not None and mbcs is not None:
      self.impl.torque(deref(mbs.v), deref(mbcs.v))
    else:
      raise TypeError("Wrong arguments passed to MotionConstr.torque")
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef MotionConstr MotionConstrFromPtr(c_qp.MotionConstr * p):
    cdef MotionConstr ret = MotionConstr(None, 0, None, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = ret.genineq_base = ret.constraint_base = p
    return ret

cdef class MotionPolyConstr(GenInequality):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, tasks.PolyTorqueBound tb, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.MotionPolyConstr(deref(mbs.v), robotIndex, tb.impl)
      self.cf_base = self.impl
      self.genineq_base = self.impl
      self.constraint_base = self.impl
  # Motion default
  def computeTorque(self, VectorXd alphaD, VectorXd _lambda):
    self.impl.computeTorque(alphaD.impl, _lambda.impl)
  def torque(self, MultiBodyVector mbs = None, MultiBodyConfigVector mbcs = None):
    if mbs is None and mbcs is None:
      VectorXdFromC(self.impl.torque())
    elif mbs is not None and mbcs is not None:
      self.impl.torque(deref(mbs.v), deref(mbcs.v))
    else:
      raise TypeError("Wrong arguments passed to MotionPolyConstr.torque")
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef MotionPolyConstr MotionPolyConstrFromPtr(c_qp.MotionPolyConstr * p):
    cdef MotionPolyConstr ret = MotionPolyConstr(None, None, None, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = ret.genineq_base = ret.constraint_base = p
    return ret

cdef class MotionSpringConstr(GenInequality):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, tasks.TorqueBound tb, sjs, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.MotionSpringConstr(deref(mbs.v), robotIndex, tb.impl, SpringJointVector(sjs).v)
      self.cf_base = self.impl
      self.genineq_base = self.impl
      self.constraint_base = self.impl
  # Motion default
  def computeTorque(self, VectorXd alphaD, VectorXd _lambda):
    self.impl.computeTorque(alphaD.impl, _lambda.impl)
  def torque(self, MultiBodyVector mbs = None, MultiBodyConfigVector mbcs = None):
    if mbs is None and mbcs is None:
      VectorXdFromC(self.impl.torque())
    elif mbs is not None and mbcs is not None:
      self.impl.torque(deref(mbs.v), deref(mbcs.v))
    else:
      raise TypeError("Wrong arguments passed to MotionSpringConstr.torque")
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef MotionSpringConstr MotionSpringConstrFromPtr(c_qp.MotionSpringConstr * p):
    cdef MotionSpringConstr ret = MotionSpringConstr(None, None, None, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = ret.genineq_base = ret.constraint_base = p
    return ret

cdef class PositiveLambda(Bound):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
     self.impl = self.bound_base = self.cf_base = self.constraint_base = new c_qp.PositiveLambda()
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef PositiveLambda PositiveLambdaFromPtr(c_qp.PositiveLambda * p):
    cdef PositiveLambda ret = PositiveLambda(skip_alloc = True)
    ret.__own_impl = False
    ret.impl = ret.bound_base = ret.constraint_base = p
    return ret

cdef class ContactConstrCommon(Equality):
  def addVirtualContact(self, ContactId cid):
    return self.contactconstr_base.addVirtualContact(cid.impl)
  def removeVirtualContact(self, ContactId cid):
    return self.contactconstr_base.removeVirtualContact(cid.impl)
  def resetVirualContacts(self):
    self.contactconstr_base.resetVirtualContacts()
  def addDofContact(self, ContactId cid, MatrixXd dof):
    return self.contactconstr_base.addDofContact(cid.impl, dof.impl)
  def removeDofContact(self, ContactId cid):
    return self.contactconstr_base.removeDofContact(cid.impl)
  def resetDofContacts(self):
    self.contactconstr_base.resetDofContacts()

cdef class ContactAccConstr(ContactConstrCommon):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.ContactAccConstr()
      self.contactconstr_base = self.impl
      self.cf_base = self.impl
      self.eq_base = self.impl
      self.constraint_base = self.impl
  def updateDofContacts(self):
    self.impl.updateDofContacts()
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef ContactAccConstr ContactAccConstrFromPtr(c_qp.ContactAccConstr * p):
    cdef ContactAccConstr ret = ContactAccConstr(skip_alloc = True)
    ret.__own_impl = False
    ret.impl = p
    ret.contactconstr_base = ret.impl
    ret.eq_base = ret.impl
    ret.constraint_base = ret.impl
    return ret

cdef class ContactSpeedConstr(ContactConstrCommon):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, double timeStep, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.ContactSpeedConstr(timeStep)
      self.contactconstr_base = self.impl
      self.cf_base = self.impl
      self.eq_base = self.impl
      self.constraint_base = self.impl
  def updateDofContacts(self):
    self.impl.updateDofContacts()
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef ContactSpeedConstr ContactSpeedConstrFromPtr(c_qp.ContactSpeedConstr * p):
    cdef ContactSpeedConstr ret = ContactSpeedConstr(0.0, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = p
    ret.contactconstr_base = ret.impl
    ret.eq_base = ret.impl
    ret.constraint_base = ret.impl
    return ret

cdef class ContactPosConstr(ContactConstrCommon):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, double timeStep, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.ContactPosConstr(timeStep)
      self.contactconstr_base = self.impl
      self.cf_base = self.impl
      self.eq_base = self.impl
      self.constraint_base = self.impl
  def updateDofContacts(self):
    self.impl.updateDofContacts()
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef ContactPosConstr ContactPosConstrFromPtr(c_qp.ContactPosConstr * p):
    cdef ContactPosConstr ret = ContactPosConstr(0.0, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = p
    ret.contactconstr_base = ret.impl
    ret.eq_base = ret.impl
    ret.constraint_base = ret.impl
    return ret

cdef class CollisionConstr(Inequality):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, MultiBodyVector mbs, double step, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.CollisionConstr(deref(mbs.v), step)
      self.cf_base = self.impl
      self.ineq_base = self.impl
      self.constraint_base = self.impl
  def addCollision(self, MultiBodyVector mbs, int collId, int r1Index, r1BodyName, S_Object body1, PTransformd X_op1_o1, int r2Index, r2BodyName, S_Object body2, PTransformd X_op2_o2, double di, double ds, double damping, double dampingOff = 0):
    if isinstance(r1BodyName, unicode):
      r1BodyName = r1BodyName.encode(u'ascii')
    if isinstance(r2BodyName, unicode):
      r2BodyName = r2BodyName.encode(u'ascii')
    self.impl.addCollision(deref(mbs.v), collId, r1Index, r1BodyName, body1.impl, deref(X_op1_o1.impl), r2Index, r2BodyName, body2.impl, deref(X_op2_o2.impl), di, ds, damping, dampingOff)
  def rmCollision(self, int collId):
    return self.impl.rmCollision(collId)
  def nrCollisions(self):
    return self.impl.nrCollisions()
  def reset(self):
    self.impl.reset()
  def updateNrCollisions(self):
    self.impl.updateNrCollisions()
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef CollisionConstr CollisionConstrFromPtr(c_qp.CollisionConstr * p):
    cdef CollisionConstr ret = CollisionConstr(None, 0, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = ret.ineq_base = ret.constraint_base = p
    return ret

cdef class CoMIncPlaneConstr(Inequality):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, double step):
    self.impl = new c_qp.CoMIncPlaneConstr(deref(mbs.v), robotIndex, step)
    self.cf_base = self.impl
    self.ineq_base = self.impl
    self.constraint_base = self.impl
  def addPlane(self, int planeId, Vector3d normal, double offset, double di, double ds, double damping, double dampingOff = 0):
    self.impl.addPlane(planeId, normal.impl, offset, di, ds, damping, dampingOff)
  def rmPlane(self, int planeId):
    return self.impl.rmPlane(planeId)
  def nrPlanes(self):
    return self.impl.nrPlanes()
  def reset(self):
    self.impl.reset()
  def updateNrPlanes(self):
    self.impl.updateNrPlanes()
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef class JointLimitsConstr(Bound):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, tasks.QBound qb, double step, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.JointLimitsConstr(deref(mbs.v), robotIndex, qb.impl, step)
      self.cf_base = self.impl
      self.bound_base = self.impl
      self.constraint_base = self.impl
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef JointLimitsConstr JointLimitsConstrFromPtr(c_qp.JointLimitsConstr * p):
    cdef JointLimitsConstr ret = JointLimitsConstr(None, 0, None, 0, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = ret.bound_base = ret.constraint_base = p
    return ret

cdef class DamperJointLimitsConstr(Bound):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, tasks.QBound qb, tasks.AlphaBound ab, double interPercent, double securityPercent, double damperOffset, double step, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.DamperJointLimitsConstr(deref(mbs.v), robotIndex, qb.impl, ab.impl, interPercent, securityPercent, damperOffset, step)
      self.cf_base = self.impl
      self.bound_base = self.impl
      self.constraint_base = self.impl
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef DamperJointLimitsConstr DamperJointLimitsConstrFromPtr(c_qp.DamperJointLimitsConstr * p):
    cdef DamperJointLimitsConstr ret = DamperJointLimitsConstr(None, 0, None, None, 0, 0, 0, 0, skip_alloc = True)
    ret.__own_impl = False
    ret.impl = ret.bound_base = ret.constraint_base = p
    return ret

cdef class GripperTorqueConstr(Inequality):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self):
    self.impl = new c_qp.GripperTorqueConstr()
    self.cf_base = self.impl
    self.ineq_base = self.impl
    self.constraint_base = self.impl
  def addGripper(self, ContactId cid, double torqueLimit, Vector3d origin, Vector3d axis):
    self.impl.addGripper(cid.impl, torqueLimit, origin.impl, axis.impl)
  def rmGripper(self, ContactId cid):
    return self.impl.rmGripper(cid.impl)
  def reset(self):
    self.impl.reset()
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef class BoundedSpeedConstr(GenInequality):
  def __dealloc__(self):
    del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, double timeStep):
    self.impl = new c_qp.BoundedSpeedConstr(deref(mbs.v), robotIndex, timeStep)
    self.cf_base = self.impl
    self.genineq_base = self.impl
    self.constraint_base = self.impl
  def addBoundedSpeed(self, MultiBodyVector mbs, bodyName, Vector3d bodyPoint, MatrixXd dof, VectorXd speed):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    self.impl.addBoundedSpeed(deref(mbs.v), bodyName, bodyPoint.impl, dof.impl, speed.impl)
  def removeBoundedSpeed(self, bodyName):
    if isinstance(bodyName, unicode):
      bodyName = bodyName.encode(u'ascii')
    return self.impl.removeBoundedSpeed(bodyName)
  def resetBoundedSpeeds(self):
    self.impl.resetBoundedSpeeds()
  def nrBoundedSpeeds(self):
    return self.impl.nrBoundedSpeeds()
  def updateBoundedSpeeds(self):
    self.impl.updateBoundedSpeeds()
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef class ImageConstr(Inequality):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, MultiBodyVector mbs, int robotIndex, bName, PTransformd X_b_gaze, double step, double constrDirection=1.):
    self.impl = new c_qp.ImageConstr(deref(mbs.v), robotIndex, bName, deref(X_b_gaze.impl), step, constrDirection)
    self.ineq_base = self.impl
    self.constraint_base = self.impl
    self.cf_base = self.impl
  def setLimits(self, Vector2d min, Vector2d max, double iPercent, double sPercent, double damping, double dampingOffsetPercent):
    self.impl.setLimits(min.impl, max.impl, iPercent, sPercent, damping, dampingOffsetPercent)
  def __v2addpt__(self, Vector2d point2d, double depthEstimate):
    return self.impl.addPoint(point2d.impl, depthEstimate)
  def __v3addpt__(self, Vector3d point3d):
    return self.impl.addPoint(point3d.impl)
  def __mbsaddpt__(self, MultiBodyVector mbs, bName, PTransformd X_b_p=PTransformd.Identity()):
    self.impl.addPoint(deref(mbs.v), bName, deref(X_b_p.impl))
  def addPoint(self, *args):
    if check_args(args, [Vector2d, None]):
      return self.__v2addpt__(*args)
    elif check_args(args, [Vector3d]):
      return self.__v3addpt__(*args)
    elif check_args(args, [MultiBodyVector, None, PTransformd]):
      self.__mbsaddpt__(*args)
    else:
      raise TypeError("Wrong arguments passed to ImageConstr addPoint")
  def reset(self):
    self.impl.reset()
  def __v2updpt__(self, int pointId, Vector2d point2d):
    self.impl.updatePoint(pointId, point2d.impl)
  def __v2dupdpt__(self, int pointId, Vector2d point2d, double depthEstimate):
    self.impl.updatePoint(pointId, point2d.impl, depthEstimate)
  def __v3updpt__(self, int pointId, Vector3d point3d):
    self.impl.updatePoint(pointId, point3d.impl)
  def updatePoint(self, int pointId, *args):
    if check_args(args, [Vector2d]):
      self.__v2updpt__(pointId, *args)
    elif check_args(args, [Vector2d, None]):
      self.__v2dupdpt__(pointId, *args)
    elif check_args(args, [Vector3d]):
      self.__v3updpt__(pointId, *args)
    else:
      raise TypeError("Wrong arguments passed to ImageConstr updatePoint")
  # no binding on this function, shouldn't be used in the general case
#  def computeComponents(...):
#    self.impl.computeComponents(...)
  def __addToSolver(self, QPSolver solver):
    self.impl.addToSolver(deref(solver.impl))
  def __addToSolverMBS(self, MultiBodyVector mbs, QPSolver solver):
    self.impl.addToSolver(deref(mbs.v), deref(solver.impl))
  def addToSolver(self, *args):
    if check_args(args, [MultiBodyVector, QPSolver]):
      self.__addToSolverMBS(*args)
    else:
      self.__addToSolver(args[0])
  def removeFromSolver(self, QPSolver solver):
    self.impl.removeFromSolver(deref(solver.impl))

cdef class QPSolver(object):
  def __dealloc__(self):
    if self.__own_impl:
      del self.impl
  def __cinit__(self, skip_alloc = False):
    self.__own_impl = True
    if not skip_alloc:
      self.impl = new c_qp.QPSolver()
  def solve(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs):
    return self.impl.solve(deref(mbs.v), deref(mbcs.v))
  def solveNoMbcUpdate(self, MultiBodyVector mbs, MultiBodyConfigVector mbcs):
    return self.impl.solveNoMbcUpdate(deref(mbs.v), deref(mbcs.v))
  def updateMbc(self, MultiBodyConfig mbc, int robotIndex):
    self.impl.updateMbc(deref(mbc.impl), robotIndex)
  def updateConstrSize(self):
    self.impl.updateConstrSize()
  def nrVars(self, MultiBodyVector mbs = None, uni = None, bi = None):
    if mbs is None and uni is None and bi is None:
      return self.impl.nrVars()
    elif mbs is not None and uni is not None and bi is not None:
      self.impl.nrVars(deref(mbs.v), UnilateralContactVector(uni).v, BilateralContactVector(bi).v)
    else:
      raise TypeError("Wrong arguments passed to QPSolver.nrVars")
  def updateTasksNrVars(self, MultiBodyVector mbs):
    self.impl.updateTasksNrVars(deref(mbs.v))
  def updateConstrsNrVars(self, MultiBodyVector mbs):
    self.impl.updateConstrsNrVars(deref(mbs.v))
  def updateNrVars(self, MultiBodyVector mbs):
    self.impl.updateNrVars(deref(mbs.v))
  def addEqualityConstraint(self, Equality eq):
    self.impl.addEqualityConstraint(eq.eq_base)
  def removeEqualityConstraint(self, Equality eq):
    self.impl.removeEqualityConstraint(eq.eq_base)
  def nrEqualityConstraints(self):
    return self.impl.nrEqualityConstraints()
  def addInequalityConstraint(self, Inequality ineq):
    self.impl.addInequalityConstraint(ineq.ineq_base)
  def removeInequalityConstraint(self, Inequality ineq):
    self.impl.removeInequalityConstraint(ineq.ineq_base)
  def nrInequalityConstraints(self):
    return self.impl.nrInequalityConstraints()
  def addGenInequalityConstraint(self, GenInequality genineq):
    self.impl.addGenInequalityConstraint(genineq.genineq_base)
  def removeGenInequalityConstraint(self, GenInequality genineq):
    self.impl.removeGenInequalityConstraint(genineq.genineq_base)
  def nrGenInequalityConstraints(self):
    return self.impl.nrGenInequalityConstraints()
  def addBoundConstraint(self, Bound bound):
    self.impl.addBoundConstraint(bound.bound_base)
  def removeBoundConstraint(self, Bound bound):
    self.impl.removeBoundConstraint(bound.bound_base)
  def nrBoundConstraints(self):
    return self.impl.nrBoundConstraints()
  def __addConstraintC(self, Constraint c):
    self.impl.addConstraint(c.constraint_base)
  def __addConstraintMB(self, MultiBodyVector mbs, Constraint c):
    self.impl.addConstraint(deref(mbs.v), c.constraint_base)
  def addConstraint(self, *args):
    if len(args) == 1 and isinstance(args[0], Constraint):
      self.__addConstraintC(args[0])
    elif check_args(args, [None, Constraint]):
      self.__addConstraintMB(*args)
    else:
      raise TypeError("Wrong type passed to QPSolver.addConstraint")
  def removeConstraint(self, Constraint constraint):
    self.impl.removeConstraint(constraint.constraint_base)
  def nrConstraints(self):
    return self.impl.nrConstraints()
  def __addTaskT(self, Task t):
    self.impl.addTask(t.base)
  def __addTaskMB(self, MultiBodyVector mbs, Task t):
    self.impl.addTask(deref(mbs.v), t.base)
  def addTask(self, *args):
    if len(args) == 1 and isinstance(args[0], Task):
      self.__addTaskT(args[0])
    elif check_args(args, [None, Task]):
      self.__addTaskMB(*args)
    else:
      raise TypeError("Wrong type passed to QPSolver.addTask")
  def removeTask(self, Task t):
    self.impl.removeTask(t.base)
  def nrTasks(self):
    return self.impl.nrTasks()
  def resetTasks(self):
    self.impl.resetTasks()
  def solver(self, name):
    if isinstance(name, unicode):
      name = name.encode(u'ascii')
    self.impl.solver(name)
  def result(self):
    return VectorXdFromC(self.impl.result())
  def alphaDVec(self, robotIndex = None):
    if robotIndex is None:
      return VectorXdFromC(self.impl.alphaDVec())
    else:
      return VectorXdFromC(self.impl.alphaDVec(robotIndex))
  def lambdaVec(self, contactIndex = None):
    if contactIndex is None:
      return VectorXdFromC(self.impl.lambdaVec())
    else:
      return VectorXdFromC(self.impl.lambdaVec(contactIndex))
  def contactLambdaPosition(self, ContactId cid):
    return self.impl.contactLambdaPosition(cid.impl)
  def data(self):
    return SolverDataFromC(self.impl.data())
  def solveTime(self):
    return tasks.cpu_timesFromC(self.impl.solveTime())
  def solveAndBuildTime(self):
    return tasks.cpu_timesFromC(self.impl.solveAndBuildTime())

cdef QPSolver QPSolverFromPtr(c_qp.QPSolver * p):
    cdef QPSolver ret = QPSolver(skip_alloc = True)
    ret.__own_impl = False
    ret.impl = p
    return ret
