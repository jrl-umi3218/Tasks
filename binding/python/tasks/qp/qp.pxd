#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cimport c_qp

from libcpp.vector cimport vector
from libcpp cimport bool as cppbool

cdef class FrictionCone(object):
  cdef c_qp.FrictionCone impl

cdef class FrictionConeVector(object):
  cdef vector[c_qp.FrictionCone] v

cdef class ContactId(object):
  cdef c_qp.ContactId impl

cdef ContactId ContactIdFromC(const c_qp.ContactId&)

cdef class UnilateralContact(object):
  cdef c_qp.UnilateralContact impl

cdef class UnilateralContactVector(object):
  cdef vector[c_qp.UnilateralContact] v

cdef class BilateralContact(object):
  cdef c_qp.BilateralContact impl

cdef class BilateralContactVector(object):
  cdef vector[c_qp.BilateralContact] v

cdef class JointStiffness(object):
  cdef c_qp.JointStiffness impl

cdef class JointStiffnessVector(object):
  cdef vector[c_qp.JointStiffness] v

cdef class JointGains(object):
  cdef c_qp.JointGains impl

cdef class JointGainsVector(object):
  cdef vector[c_qp.JointGains] v

cdef class SpringJoint(object):
  cdef c_qp.SpringJoint impl

cdef class SpringJointVector(object):
  cdef vector[c_qp.SpringJoint] v

cdef class SolverData(object):
  cdef c_qp.SolverData impl

cdef SolverData SolverDataFromC(const c_qp.SolverData&)

cdef class Constraint(object):
  cdef c_qp.Constraint * constraint_base

#XXX Due to the lack of support for multiple inheritance in Cython, we need to
#    deviate from the C++ class hierarchy a bit and thus
#    Equality/Inequality/GenInequality/Bound will inherit from Constraint

cdef class Equality(Constraint):
  cdef c_qp.Equality * eq_base
  cdef c_qp.ConstraintFunction[c_qp.Equality] * cf_base

cdef class Inequality(Constraint):
  cdef c_qp.Inequality * ineq_base
  cdef c_qp.ConstraintFunction[c_qp.Inequality] * cf_base

cdef class GenInequality(Constraint):
  cdef c_qp.GenInequality * genineq_base
  cdef c_qp.ConstraintFunction[c_qp.GenInequality] * cf_base

cdef class Bound(Constraint):
  cdef c_qp.Bound * bound_base
  cdef c_qp.ConstraintFunction[c_qp.Bound] * cf_base

cdef class Task(object):
  cdef c_qp.Task * base

cdef class HighLevelTask(object):
  cdef c_qp.HighLevelTask * base

cdef class PositionTask(HighLevelTask):
  cdef c_qp.PositionTask * impl

cdef class OrientationTask(HighLevelTask):
  cdef c_qp.OrientationTask * impl

cdef class SurfaceOrientationTask(HighLevelTask):
  cdef c_qp.SurfaceOrientationTask * impl

cdef class GazeTask(HighLevelTask):
  cdef c_qp.GazeTask * impl

cdef class PositionBasedVisServoTask(HighLevelTask):
  cdef c_qp.PositionBasedVisServoTask * impl

cdef class PostureTask(Task):
  cdef c_qp.PostureTask * impl
  cdef cppbool __own_impl

cdef PostureTask PostureTaskFromPtr(c_qp.PostureTask*)

cdef class CoMTask(HighLevelTask):
  cdef c_qp.CoMTask * impl

cdef class MomentumTask(HighLevelTask):
  cdef c_qp.MomentumTask * impl

cdef class LinVelocityTask(HighLevelTask):
  cdef c_qp.LinVelocityTask * impl

cdef class OrientationTrackingTask(HighLevelTask):
  cdef c_qp.OrientationTrackingTask * impl

cdef class TransformTask(HighLevelTask):
  cdef c_qp.TransformTask * impl

cdef class SurfaceTransformTask(HighLevelTask):
  cdef c_qp.SurfaceTransformTask * impl

cdef class JointsSelector(HighLevelTask):
  cdef c_qp.JointsSelector * impl

cdef class SetPointTask(Task):
  cdef c_qp.SetPointTask * impl

cdef class TrackingTask(Task):
  cdef c_qp.TrackingTask * impl

cdef class TrajectoryTask(Task):
  cdef c_qp.TrajectoryTask * impl

cdef class TargetObjectiveTask(Task):
  cdef c_qp.TargetObjectiveTask * impl

cdef class PIDTask(Task):
  cdef c_qp.PIDTask * impl

cdef class MultiCoMTask(Task):
  cdef c_qp.MultiCoMTask * impl

cdef class MultiRobotTransformTask(Task):
  cdef c_qp.MultiRobotTransformTask * impl

cdef class ContactTask(Task):
  cdef c_qp.ContactTask * impl

cdef class GripperTorqueTask(Task):
  cdef c_qp.GripperTorqueTask * impl

cdef class MotionConstr(GenInequality):
  cdef c_qp.MotionConstr * impl
  cdef cppbool __own_impl

cdef MotionConstr MotionConstrFromPtr(c_qp.MotionConstr *)

cdef class MotionPolyConstr(GenInequality):
  cdef c_qp.MotionPolyConstr * impl
  cdef cppbool __own_impl

cdef MotionPolyConstr MotionPolyConstrFromPtr(c_qp.MotionPolyConstr *)

cdef class MotionSpringConstr(GenInequality):
  cdef c_qp.MotionSpringConstr * impl
  cdef cppbool __own_impl

cdef MotionSpringConstr MotionSpringConstrFromPtr(c_qp.MotionSpringConstr *)

cdef class PositiveLambda(Bound):
  cdef c_qp.PositiveLambda * impl
  cdef cppbool __own_impl

cdef PositiveLambda PositiveLambdaFromPtr(c_qp.PositiveLambda *)

cdef class ContactConstrCommon(Equality):
  cdef c_qp.ContactConstrCommon * contactconstr_base

cdef class ContactAccConstr(ContactConstrCommon):
  cdef c_qp.ContactAccConstr * impl
  cdef cppbool __own_impl

cdef ContactAccConstr ContactAccConstrFromPtr(c_qp.ContactAccConstr *)

cdef class ContactSpeedConstr(ContactConstrCommon):
  cdef c_qp.ContactSpeedConstr * impl
  cdef cppbool __own_impl

cdef ContactSpeedConstr ContactSpeedConstrFromPtr(c_qp.ContactSpeedConstr *)

cdef class ContactPosConstr(ContactConstrCommon):
  cdef c_qp.ContactPosConstr * impl
  cdef cppbool __own_impl

cdef ContactPosConstr ContactPosConstrFromPtr(c_qp.ContactPosConstr *)

cdef class CollisionConstr(Inequality):
  cdef c_qp.CollisionConstr * impl
  cdef cppbool __own_impl

cdef CollisionConstr CollisionConstrFromPtr(c_qp.CollisionConstr*)

cdef class CoMIncPlaneConstr(Inequality):
  cdef c_qp.CoMIncPlaneConstr * impl

cdef class JointLimitsConstr(Bound):
  cdef c_qp.JointLimitsConstr * impl
  cdef cppbool __own_impl

cdef JointLimitsConstr JointLimitsConstrFromPtr(c_qp.JointLimitsConstr *)

cdef class DamperJointLimitsConstr(Bound):
  cdef c_qp.DamperJointLimitsConstr * impl
  cdef cppbool __own_impl

cdef DamperJointLimitsConstr DamperJointLimitsConstrFromPtr(c_qp.DamperJointLimitsConstr *)

cdef class GripperTorqueConstr(Inequality):
  cdef c_qp.GripperTorqueConstr * impl

cdef class BoundedSpeedConstr(GenInequality):
  cdef c_qp.BoundedSpeedConstr * impl

cdef class ImageConstr(Inequality):
  cdef c_qp.ImageConstr * impl

cdef class QPSolver(object):
  cdef c_qp.QPSolver * impl
  cdef cppbool __own_impl

cdef QPSolver QPSolverFromPtr(c_qp.QPSolver *)
