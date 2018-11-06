#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

LEGACY = False

import unittest
from utils import expected_failure

import math
try:
  import eigen
  import sva
except ImportError:
  LEGACY = True
  import eigen3 as eigen
  import spacevecalg as sva
  from mc_rbdyn import rbdList
import rbdyn
import tasks
import sch

import arms

class TestFrictionCone(unittest.TestCase):
  def test_cone1(self):
    angle = math.pi/4
    mu = math.tan(angle)
    cone = tasks.qp.FrictionCone(eigen.Matrix3d.Identity(), 4, mu)
    for v in cone.generators:
      # check cone equation x^2 + y^2 = z^2(tan(angle)^2)
      self.assertAlmostEqual(v.x()**2 + v.y()**2, v.z()**2*mu**2, delta = 0.00001)

  def test_cone2(self):
    angle = math.pi/4
    mu = math.tan(angle)
    rep = eigen.Matrix3d.Zero()
    rep.coeff(0, 1, 1)
    rep.coeff(1, 2, 1)
    rep.coeff(2, 0, 1)
    cone = tasks.qp.FrictionCone(rep, 4, mu)
    for v in cone.generators:
      # check cone equation x^2 + y^2 = z^2(tan(angle)^2)
      # check cone equation in rep frame z^2 + y^2 = x^2(tan(angle)^2)
      self.assertAlmostEqual(v.z()**2 + v.y()**2, v.x()**2*mu**2, delta = 0.00001)

class TestQPTask(unittest.TestCase):
  def setUp(self):
    self.nrIter = 10000
    mb, self.mbcInit = arms.makeZXZArm()

    rbdyn.forwardKinematics(mb, self.mbcInit)
    rbdyn.forwardVelocity(mb, self.mbcInit)
    if not LEGACY:
      self.mbs = rbdyn.MultiBodyVector([mb])
      self.mbcs = rbdyn.MultiBodyConfigVector([self.mbcInit])
    else:
      self.mbs = [mb]
      self.mbcs = [rbdyn.MultiBodyConfig(self.mbcInit)]

    self.solver = tasks.qp.QPSolver()

    self.solver.nrVars(self.mbs, [], [])
    self.assertEqual(self.solver.nrVars(), 3)

    self.solver.updateConstrSize()

  def tearDown(self):
    pass

  def run_solver(self):
    for i in range(self.nrIter):
      if not LEGACY:
        self.assertTrue(self.solver.solve(self.mbs, self.mbcs))
      else:
        self.assertTrue(self.solver.solveNoMbcUpdate(self.mbs, self.mbcs))
        self.solver.updateMbc(self.mbcs[0], 0)
      rbdyn.eulerIntegration(self.mbs[0], self.mbcs[0], 0.001)
      rbdyn.forwardKinematics(self.mbs[0], self.mbcs[0])
      rbdyn.forwardVelocity(self.mbs[0], self.mbcs[0])

  def test_position_task(self):
    posD = eigen.Vector3d(0.707106, 0.707106, 0.)
    posTask = tasks.qp.PositionTask(self.mbs, 0, "b3", posD)
    posTaskSp = tasks.qp.SetPointTask(self.mbs, 0, posTask, 10, 1)

    self.solver.addTask(posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)

    # Test PositionTask
    self.run_solver()

    self.assertAlmostEqual(posTask.eval().norm(), 0, delta = 0.00001)

    self.solver.removeTask(posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

  def test_orientation_task(self):
    oriTask = tasks.qp.OrientationTask(self.mbs, 0, "b3", sva.RotZ(math.pi/2))
    oriTaskSp = tasks.qp.SetPointTask(self.mbs, 0, oriTask, 10, 1)

    self.solver.addTask(oriTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)

    # Test OrientationTask
    self.run_solver()

    self.assertAlmostEqual(oriTask.eval().norm(), 0, delta = 0.00001)

    self.solver.removeTask(oriTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

  def test_posture_task(self):
    postureTask = tasks.qp.PostureTask(self.mbs, 0, [[], [0.2], [0.4], [-0.8]], 10, 1)
    postureTask.jointsStiffness(self.mbs, [tasks.qp.JointStiffness("j2",10)])
    self.solver.addTask(postureTask)

    self.run_solver()

    self.assertAlmostEqual(postureTask.eval().norm(), 0, delta = 0.00001)

    self.solver.removeTask(postureTask)
    self.assertEqual(self.solver.nrTasks(), 0)

  def test_com_task(self):
    comD = eigen.Vector3d(sva.RotZ(math.pi/4)*rbdyn.computeCoM(self.mbs[0], self.mbcInit))
    comTask = tasks.qp.CoMTask(self.mbs, 0, comD)
    comTaskSp = tasks.qp.SetPointTask(self.mbs, 0, comTask, 10, 1)

    self.solver.addTask(comTaskSp)

    self.run_solver()

    self.assertAlmostEqual(comTask.eval().norm(), 0, delta = 0.00001)

    self.solver.removeTask(comTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

  def test_lin_velocity_task(self):
    linVelocityTask = tasks.qp.LinVelocityTask(self.mbs, 0, "b3", -eigen.Vector3d.UnitX()*0.005)
    linVelocityTaskSp = tasks.qp.SetPointTask(self.mbs, 0, linVelocityTask, 1000, 1000)

    self.solver.addTask(linVelocityTaskSp)

    self.run_solver()

    # XXX Huge error (same as cpp test)
    self.assertAlmostEqual(linVelocityTask.eval().norm(), 0, delta = 0.001)

    self.solver.removeTask(linVelocityTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

#
#
class TestQPConstr(unittest.TestCase):
  def setUp(self):
    self.nrIter = 1000
    mb, self.mbcInit = arms.makeZXZArm()
    mbEnv, mbcEnv = arms.makeEnv()
    if not LEGACY:
      self.mbs = rbdyn.MultiBodyVector([mb, mbEnv])
    else:
      self.mbs = [mb, mbEnv]

    rbdyn.forwardKinematics(mb, self.mbcInit)
    rbdyn.forwardVelocity(mb, self.mbcInit)
    rbdyn.forwardKinematics(mbEnv, mbcEnv)
    rbdyn.forwardVelocity(mbEnv, mbcEnv)

    if not LEGACY:
      self.mbcs = rbdyn.MultiBodyConfigVector([self.mbcInit, mbcEnv])
    else:
      self.mbcs = [rbdyn.MultiBodyConfig(self.mbcInit), mbcEnv]

    self.solver = tasks.qp.QPSolver()

    self.contVec = []

    posD = eigen.Vector3d(0.707106, 0.707106, 0.)
    self.posTask = tasks.qp.PositionTask(self.mbs, 0, "b3", posD)
    self.posTaskSp = tasks.qp.SetPointTask(self.mbs, 0, self.posTask, 10, 1)

    self.oriTask = tasks.qp.OrientationTask(self.mbs, 0, "b3", sva.RotZ(math.pi/2))
    self.oriTaskSp = tasks.qp.SetPointTask(self.mbs, 0, self.oriTask, 10, 1)

  def tearDown(self):
    pass

  def run_solver(self):
    for i in range(self.nrIter):
      if not LEGACY:
        self.assertTrue(self.solver.solve(self.mbs, self.mbcs))
      else:
        self.assertTrue(self.solver.solveNoMbcUpdate(self.mbs, self.mbcs))
        self.solver.updateMbc(self.mbcs[0], 0)
      rbdyn.eulerIntegration(self.mbs[0], self.mbcs[0], 0.001)
      rbdyn.forwardKinematics(self.mbs[0], self.mbcs[0])
      rbdyn.forwardVelocity(self.mbs[0], self.mbcs[0])

  def check_equality_constr(self, ConstrClass, *args):
    self.contVec = [tasks.qp.UnilateralContact(0, 1, "b3", "b0", [eigen.Vector3d.Zero()], eigen.Matrix3d.Identity(), sva.PTransformd.Identity(), 3, math.tan(math.pi/4))]
    constr = ConstrClass(*args)
    constr.addToSolver(self.solver)
    self.assertEqual(self.solver.nrEqualityConstraints(), 1)
    self.assertEqual(self.solver.nrConstraints(), 1)

    self.solver.nrVars(self.mbs, self.contVec, [])
    self.solver.updateConstrSize()

    # Check vars number is 3 dof + 3 lambda
    self.assertEqual(self.solver.nrVars(), 3+3)

    self.solver.addTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)
    self.solver.addTask(self.oriTaskSp)
    self.assertEqual(self.solver.nrTasks(), 2)

    self.solver.data().computeNormalAccB(self.mbs, self.mbcs)
    self.posTask.update(self.mbs, self.mbcs, self.solver.data())
    self.oriTask.update(self.mbs, self.mbcs, self.solver.data())
    evalPos = self.posTask.eval()
    evalOri = self.oriTask.eval()

    self.run_solver()

    if not LEGACY:
        self.assertAlmostEqual((self.posTask.eval() - evalPos).norm(), 0, delta = 0.00001)
        self.assertAlmostEqual((self.oriTask.eval() - evalOri).norm(), 0, delta = 0.00001)

    self.solver.removeTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)
    self.solver.removeTask(self.oriTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

    constr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrEqualityConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 0)

  def test_contact_acc_constr(self):
    self.check_equality_constr(tasks.qp.ContactAccConstr)

  def test_contact_speed_constr(self):
    self.check_equality_constr(tasks.qp.ContactSpeedConstr, 0.001)

  @expected_failure
  def test_contact_pos_constr(self):
    self.check_equality_constr(tasks.qp.ContactPosConstr, 0.001)

  def test_motion_constr(self):
    Inf = float("inf")
    torqueMin = [[], [-Inf], [-Inf], [-Inf]]
    torqueMax = [[], [Inf], [Inf], [Inf]]
    motionCstr = tasks.qp.MotionConstr(self.mbs, 0, tasks.TorqueBound(torqueMin, torqueMax))
    plCstr = tasks.qp.PositiveLambda()

    motionCstr.addToSolver(self.solver)
    self.assertEqual(self.solver.nrGenInequalityConstraints(), 1)
    self.assertEqual(self.solver.nrConstraints(), 1)

    plCstr.addToSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 1)
    self.assertEqual(self.solver.nrConstraints(), 2)

    self.solver.nrVars(self.mbs, self.contVec, [])
    self.solver.updateConstrSize()

    self.solver.addTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)

    self.nrIter = 10000
    self.run_solver()

    motionCstr.computeTorque(self.solver.alphaDVec(), self.solver.lambdaVec())
    motionCstr.torque(self.mbs, self.mbcs)

    self.assertAlmostEqual(self.posTask.eval().norm(), 0, delta = 5e-5)

    mbcTest = rbdyn.MultiBodyConfig(self.mbcs[0])
    mbcTest.jointTorque = [[], [0], [0], [0]]
    _id = rbdyn.InverseDynamics(self.mbs[0])
    _id.inverseDynamics(self.mbs[0], mbcTest)

    # This assert only pass on the cython bindings
    if not LEGACY:
      self.assertAlmostEqual((rbdyn.dofToVector(self.mbs[0], self.mbcs[0].jointTorque) - rbdyn.dofToVector(self.mbs[0], mbcTest.jointTorque)).norm(), 0, delta = 0.00001)

    motionCstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrGenInequalityConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 1)

    plCstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 0)

  def test_motion_constr_w_contact(self):
    Inf = float("inf")
    torqueMin = [[], [-Inf], [-Inf], [-Inf]]
    torqueMax = [[], [Inf], [Inf], [Inf]]
    motionCstr = tasks.qp.MotionConstr(self.mbs, 0, tasks.TorqueBound(torqueMin, torqueMax))
    plCstr = tasks.qp.PositiveLambda()

    motionCstr.addToSolver(self.solver)
    self.assertEqual(self.solver.nrGenInequalityConstraints(), 1)
    self.assertEqual(self.solver.nrConstraints(), 1)

    plCstr.addToSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 1)
    self.assertEqual(self.solver.nrConstraints(), 2)

    self.contVec = [tasks.qp.UnilateralContact(0, 1, "b3", "b0", [eigen.Vector3d.Zero()], eigen.Matrix3d.Identity(), sva.PTransformd.Identity(), 3, math.tan(math.pi/4))]
    self.solver.nrVars(self.mbs, self.contVec, [])
    self.solver.updateConstrSize()
    self.assertEqual(self.solver.nrVars(), 3+3)

    self.solver.addTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)

    self.nrIter = 10000
    self.run_solver()

    self.assertAlmostEqual(self.posTask.eval().norm(), 0, delta = 5e-5)

    self.solver.removeTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

    motionCstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrGenInequalityConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 1)

    plCstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 0)

class TestQPJointLimits(unittest.TestCase):
  def setUp(self):
    self.nrIter = 1000
    mb, self.mbcInit = arms.makeZXZArm()
    if not LEGACY:
      self.mbs = rbdyn.MultiBodyVector([mb])
    else:
      self.mbs = [mb]
    self.mb = self.mbs[0]

    rbdyn.forwardKinematics(mb, self.mbcInit)
    rbdyn.forwardVelocity(mb, self.mbcInit)

    if not LEGACY:
      self.mbcs = rbdyn.MultiBodyConfigVector([self.mbcInit])
    else:
      self.mbcs = [rbdyn.MultiBodyConfig(self.mbcInit)]

    self.solver = tasks.qp.QPSolver()

    self.bodyI = self.mb.bodyIndexByName("b3")
    if not LEGACY:
      self.posTask = tasks.qp.PositionTask(self.mbs, 0, "b3", sva.RotZ(math.pi/2)*self.mbcInit.bodyPosW[self.bodyI].translation())
    else:
      self.posTask = tasks.qp.PositionTask(self.mbs, 0, "b3", sva.RotZ(math.pi/2)*list(self.mbcInit.bodyPosW)[self.bodyI].translation())
    self.posTaskSp = tasks.qp.SetPointTask(self.mbs, 0, self.posTask, 10, 1)

  def tearDown(self):
    pass

  def run_solver(self):
    for i in range(self.nrIter):
      if not LEGACY:
        self.assertTrue(self.solver.solve(self.mbs, self.mbcs))
      else:
        self.assertTrue(self.solver.solveNoMbcUpdate(self.mbs, self.mbcs))
        self.solver.updateMbc(self.mbcs[0], 0)
      rbdyn.eulerIntegration(self.mbs[0], self.mbcs[0], 0.001)
      rbdyn.forwardKinematics(self.mbs[0], self.mbcs[0])
      rbdyn.forwardVelocity(self.mbs[0], self.mbcs[0])
      if not LEGACY:
        self.assertGreater(self.mbcs[0].q[1][0], -math.pi/4 - 0.01)
        self.assertLess(self.mbcs[0].q[1][0], math.pi/4 + 0.01)
      else:
        self.assertGreater(rbdList(self.mbcs[0].q)[1][0], -math.pi/4 - 0.01)
        self.assertLess(rbdList(self.mbcs[0].q)[1][0], math.pi/4 + 0.01)

  def test_joint_limits(self):
    inf = float("inf")
    lBound = [[], [-math.pi/4], [-inf], [-inf]]
    uBound = [[], [math.pi/4], [inf], [inf]]

    jointConstr = tasks.qp.JointLimitsConstr(self.mbs, 0, tasks.QBound(lBound, uBound), 0.001)

    self.solver.addBoundConstraint(jointConstr)
    self.assertEqual(self.solver.nrBoundConstraints(), 1)
    self.solver.addConstraint(jointConstr)
    self.assertEqual(self.solver.nrConstraints(), 1)

    self.solver.nrVars(self.mbs, [], [])
    self.solver.updateConstrSize()

    self.solver.addTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)

    self.run_solver()

    if not LEGACY:
      self.posTask.position(sva.RotZ(-math.pi/2)*self.mbcInit.bodyPosW[self.bodyI].translation())
    else:
      self.posTask.position(sva.RotZ(-math.pi/2)*list(self.mbcInit.bodyPosW)[self.bodyI].translation())
    if not LEGACY:
      self.mbcs[0] = self.mbcInit
    else:
      self.mbcs = [rbdyn.MultiBodyConfig(self.mbcInit)]

    self.run_solver()

    self.solver.removeTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

    self.solver.removeBoundConstraint(jointConstr)
    self.assertEqual(self.solver.nrBoundConstraints(), 0)
    self.solver.removeConstraint(jointConstr)
    self.assertEqual(self.solver.nrConstraints(), 0)

  def test_damper_joint_limits(self):
    inf = float("inf")
    lBound = [[], [-math.pi/4], [-inf], [-inf]]
    uBound = [[], [math.pi/4], [inf], [inf]]
    lVel = [[], [-inf], [-inf], [-inf]]
    uVel = [[], [inf], [inf], [inf]]

    dampJointConstr = tasks.qp.DamperJointLimitsConstr(self.mbs, 0, tasks.QBound(lBound, uBound), tasks.AlphaBound(lVel, uVel), 0.125, 0.025, 1.0, 0.001)

    dampJointConstr.addToSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 1)
    self.assertEqual(self.solver.nrConstraints(), 1)

    self.solver.nrVars(self.mbs, [], [])
    self.solver.updateConstrSize()

    self.solver.addTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)

    self.nrIter = 2000
    self.run_solver()

    if not LEGACY:
      self.posTask.position(sva.RotZ(-math.pi/2)*self.mbcInit.bodyPosW[self.bodyI].translation())
    else:
      self.posTask.position(sva.RotZ(-math.pi/2)*list(self.mbcInit.bodyPosW)[self.bodyI].translation())
    if not LEGACY:
      self.mbcs[0] = self.mbcInit
    else:
      self.mbcs = [rbdyn.MultiBodyConfig(self.mbcInit)]
    self.run_solver()

    self.solver.removeTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

    dampJointConstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 0)

class TestQPTorqueLimits(unittest.TestCase):
  def setUp(self):
    self.nrIter = 10000
    mb, self.mbcInit = arms.makeZXZArm()
    if not LEGACY:
      self.mbs = rbdyn.MultiBodyVector([mb])
    else:
      self.mbs = [mb]
    self.mb = self.mbs[0]

    rbdyn.forwardKinematics(mb, self.mbcInit)
    rbdyn.forwardVelocity(mb, self.mbcInit)

    if not LEGACY:
      self.mbcs = rbdyn.MultiBodyConfigVector([self.mbcInit])
    else:
      self.mbcs = [rbdyn.MultiBodyConfig(self.mbcInit)]

    self.solver = tasks.qp.QPSolver()

    self.bodyI = self.mb.bodyIndexByName("b3")
    if not LEGACY:
      self.posTask = tasks.qp.PositionTask(self.mbs, 0, "b3", sva.RotZ(math.pi/2)*self.mbcInit.bodyPosW[self.bodyI].translation())
    else:
      self.posTask = tasks.qp.PositionTask(self.mbs, 0, "b3", sva.RotZ(math.pi/2)*list(self.mbcInit.bodyPosW)[self.bodyI].translation())
    self.posTaskSp = tasks.qp.SetPointTask(self.mbs, 0, self.posTask, 10, 1)

    self.plCstr = tasks.qp.PositiveLambda()
    self.plCstr.addToSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 1)
    self.assertEqual(self.solver.nrConstraints(), 1)

  def tearDown(self):
    pass

  def test_motion_constr(self):
    lBound = [[], [-30], [-30], [-30]]
    uBound = [[], [30], [30], [30]]
    motionCstr = tasks.qp.MotionConstr(self.mbs, 0, tasks.TorqueBound(lBound, uBound))

    self.solver.addGenInequalityConstraint(motionCstr)
    self.assertEqual(self.solver.nrGenInequalityConstraints(), 1)
    self.solver.addConstraint(motionCstr)
    self.assertEqual(self.solver.nrConstraints(), 2)


    self.solver.addTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 1)

    self.solver.nrVars(self.mbs, [], [])
    self.solver.updateConstrSize()

    for i in range(self.nrIter):
      if not LEGACY:
        self.assertTrue(self.solver.solve(self.mbs, self.mbcs))
      else:
        self.assertTrue(self.solver.solveNoMbcUpdate(self.mbs, self.mbcs))
        self.solver.updateMbc(self.mbcs[0], 0)
      rbdyn.eulerIntegration(self.mbs[0], self.mbcs[0], 0.001)
      rbdyn.forwardKinematics(self.mbs[0], self.mbcs[0])
      rbdyn.forwardVelocity(self.mbs[0], self.mbcs[0])
      motionCstr.computeTorque(self.solver.alphaDVec(), self.solver.lambdaVec())
      motionCstr.torque(self.mbs, self.mbcs)
      if not LEGACY:
        for i in range(3):
          self.assertGreater(self.mbcs[0].jointTorque[i+1][0], lBound[i+1][0] - 0.001)
          self.assertLess(self.mbcs[0].jointTorque[i+1][0], uBound[i+1][0] + 0.001)
      else:
        pass #This test will fail for the legacy bindings

    if not LEGACY:
      self.posTask.position(self.mbcInit.bodyPosW[self.bodyI].translation())
    else:
      self.posTask.position(list(self.mbcInit.bodyPosW)[self.bodyI].translation())

    for i in range(self.nrIter):
      if not LEGACY:
        self.assertTrue(self.solver.solve(self.mbs, self.mbcs))
      else:
        self.assertTrue(self.solver.solveNoMbcUpdate(self.mbs, self.mbcs))
        self.solver.updateMbc(self.mbcs[0], 0)
      rbdyn.eulerIntegration(self.mbs[0], self.mbcs[0], 0.001)
      rbdyn.forwardKinematics(self.mbs[0], self.mbcs[0])
      rbdyn.forwardVelocity(self.mbs[0], self.mbcs[0])
      motionCstr.computeTorque(self.solver.alphaDVec(), self.solver.lambdaVec())
      motionCstr.torque(self.mbs, self.mbcs)
      if not LEGACY:
        for i in range(3):
          self.assertGreater(self.mbcs[0].jointTorque[i+1][0], lBound[i+1][0] - 0.001)
          self.assertLess(self.mbcs[0].jointTorque[i+1][0], uBound[i+1][0] + 0.001)
      else:
        pass #This test will fail for the legacy bindings

    motionCstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrGenInequalityConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 1)

    self.plCstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 0)

    self.solver.removeTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

  def test_motion_poly_constr(self):
    if not LEGACY:
      lpoly = eigen.VectorXd(-30, 1)
      upoly = eigen.VectorXd(30, 1)
    else:
      lpoly = eigen.VectorXd(2)
      lpoly[0] = -30
      lpoly[1] = 1
      upoly = eigen.VectorXd(2)
      upoly[0] = 30
      upoly[1] = 1
    null = eigen.VectorXd()
    lBoundPoly = [[null], [lpoly], [lpoly], [lpoly]]
    uBoundPoly = [[null], [upoly], [upoly], [upoly]]
    motionPolyCstr = tasks.qp.MotionPolyConstr(self.mbs, 0, tasks.PolyTorqueBound(lBoundPoly, uBoundPoly))

    motionPolyCstr.addToSolver(self.solver)
    self.assertEqual(self.solver.nrGenInequalityConstraints(), 1)
    self.assertEqual(self.solver.nrBoundConstraints(), 1)
    self.assertEqual(self.solver.nrConstraints(), 2)

    self.solver.nrVars(self.mbs, [], [])
    self.solver.updateConstrSize()

    for i in range(self.nrIter):
      if not LEGACY:
        self.assertTrue(self.solver.solve(self.mbs, self.mbcs))
      else:
        self.assertTrue(self.solver.solveNoMbcUpdate(self.mbs, self.mbcs))
        self.solver.updateMbc(self.mbcs[0], 0)
      oldQ = self.mbcs[0].q
      rbdyn.eulerIntegration(self.mbs[0], self.mbcs[0], 0.001)
      rbdyn.forwardKinematics(self.mbs[0], self.mbcs[0])
      rbdyn.forwardVelocity(self.mbs[0], self.mbcs[0])
      motionPolyCstr.computeTorque(self.solver.alphaDVec(), self.solver.lambdaVec())
      motionPolyCstr.torque(self.mbs, self.mbcs)
      if not LEGACY:
        for i in range(3):
          self.assertGreater(self.mbcs[0].jointTorque[i+1][0], eigen.poly_eval(lBoundPoly[i+1][0], oldQ[i+1][0]))
          self.assertLess(self.mbcs[0].jointTorque[i+1][0], eigen.poly_eval(uBoundPoly[i+1][0], oldQ[i+1][0]))
      else:
        pass #This test will fail for the legacy bindings

    motionPolyCstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrGenInequalityConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 1)

    self.plCstr.removeFromSolver(self.solver)
    self.assertEqual(self.solver.nrBoundConstraints(), 0)
    self.assertEqual(self.solver.nrConstraints(), 0)

    self.solver.removeTask(self.posTaskSp)
    self.assertEqual(self.solver.nrTasks(), 0)

class TestQPAutoColl(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm()

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit])
    else:
      mbs = [mb]
      mbcs = [rbdyn.MultiBodyConfig(mbcInit)]

    solver = tasks.qp.QPSolver()

    bodyI = mb.bodyIndexByName("b3")
    if not LEGACY:
      posTask = tasks.qp.PositionTask(mbs, 0, "b3", mbcInit.bodyPosW[bodyI].translation())
    else:
      posTask = tasks.qp.PositionTask(mbs, 0, "b3", list(mbcInit.bodyPosW)[bodyI].translation())
    posTaskSp = tasks.qp.SetPointTask(mbs, 0, posTask, 50, 1)

    b0 = sch.Sphere(0.25)
    b3 = sch.Sphere(0.25)
    pair = sch.CD_Pair(b0, b3)

    I = sva.PTransformd.Identity()
    autoCollConstr = tasks.qp.CollisionConstr(mbs, 0.001)
    collId1 = 10
    autoCollConstr.addCollision(mbs, collId1,
      0, "b0", b0, I,
      0, "b3", b3, I,
      0.01, 0.005, 1)
    self.assertEqual(autoCollConstr.nrCollisions(), 1)

    solver.addInequalityConstraint(autoCollConstr)
    self.assertEqual(solver.nrInequalityConstraints(), 1)
    solver.addConstraint(autoCollConstr)
    self.assertEqual(solver.nrConstraints(), 1)

    solver.nrVars(mbs, [], [])
    solver.updateConstrSize()

    solver.addTask(posTaskSp)
    self.assertEqual(solver.nrTasks(), 1)

    mbcs[0] = mbcInit
    for i in range(1000):
      posTask.position(sva.RotX(0.01)*posTask.position())
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.001)

      self.assertGreater(math.sqrt(pair.distance()), 0.001)

      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    autoCollConstr.rmCollision(collId1)
    self.assertEqual(autoCollConstr.nrCollisions(), 0)

class TestQPStaticEnvColl(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm()
    mbEnv, mbcEnv = arms.makeEnv()

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)
    rbdyn.forwardKinematics(mbEnv, mbcEnv)
    rbdyn.forwardVelocity(mbEnv, mbcEnv)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb, mbEnv])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit, mbcEnv])
    else:
      mbs = [mb, mbEnv]
      mbcs = [rbdyn.MultiBodyConfig(mbcInit), mbcEnv]

    solver = tasks.qp.QPSolver()

    bodyI = mb.bodyIndexByName("b3")
    if not LEGACY:
      posTask = tasks.qp.PositionTask(mbs, 0, "b3", mbcInit.bodyPosW[bodyI].translation())
    else:
      posTask = tasks.qp.PositionTask(mbs, 0, "b3", list(mbcInit.bodyPosW)[bodyI].translation())
    posTaskSp = tasks.qp.SetPointTask(mbs, 0, posTask, 50, 1)

    b0 = sch.Sphere(0.25)
    b3 = sch.Sphere(0.25)
    pair = sch.CD_Pair(b0, b3)

    if not LEGACY:
      b0.transform(mbcInit.bodyPosW[0])
    else:
      b0.transform(list(mbcInit.bodyPosW)[0])

    I = sva.PTransformd.Identity()
    seCollConstr = tasks.qp.CollisionConstr(mbs, 0.001)
    collId1 = 10
    seCollConstr.addCollision(mbs, collId1,
      0, "b3", b3, I,
      1, "b0", b0, I,
      0.01, 0.005, 1)
    self.assertEqual(seCollConstr.nrCollisions(), 1)

    solver.addInequalityConstraint(seCollConstr)
    self.assertEqual(solver.nrInequalityConstraints(), 1)
    solver.addConstraint(seCollConstr)
    self.assertEqual(solver.nrConstraints(), 1)

    solver.nrVars(mbs, [], [])
    solver.updateConstrSize()

    solver.addTask(posTaskSp)
    self.assertEqual(solver.nrTasks(), 1)

    mbcs[0] = mbcInit
    for i in range(1000):
      posTask.position(sva.RotX(0.01)*posTask.position())
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.001)

      self.assertGreater(math.sqrt(pair.distance()), 0.001)

      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    seCollConstr.rmCollision(collId1)
    self.assertEqual(seCollConstr.nrCollisions(), 0)

    # Test damping computation
    seCollConstr.addCollision(mbs, collId1,
      0, "b3", b3, I,
      1, "b0", b0, I,
      0.1, 0.01, 0, 0.1)
    self.assertEqual(seCollConstr.nrCollisions(), 1)
    if not LEGACY:
      posTask.position(mbcInit.bodyPosW[bodyI].translation())
    else:
      posTask.position(list(mbcInit.bodyPosW)[bodyI].translation())

    mbcs[0] = mbcInit
    for i in range(1000):
      posTask.position(sva.RotX(0.01)*posTask.position())
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.001)

      self.assertGreater(math.sqrt(pair.distance()), 0.001)

      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    seCollConstr.rmCollision(collId1)
    self.assertEqual(seCollConstr.nrCollisions(), 0)

    solver.removeTask(posTaskSp)
    self.assertEqual(solver.nrTasks(), 0)

    solver.removeInequalityConstraint(seCollConstr)
    self.assertEqual(solver.nrInequalityConstraints(), 0)
    solver.removeConstraint(seCollConstr)
    self.assertEqual(solver.nrConstraints(), 0)

class TestQPBilatContact(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm(False)
    mbEnv, mbcEnv = arms.makeEnv()

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)
    rbdyn.forwardKinematics(mbEnv, mbcEnv)
    rbdyn.forwardVelocity(mbEnv, mbcEnv)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb, mbEnv])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit, mbcEnv])
    else:
      mbs = [mb, mbEnv]
      mbcs = [rbdyn.MultiBodyConfig(mbcInit), mbcEnv]

    solver = tasks.qp.QPSolver()

    Inf = float("inf")
    torqueMin = [[0, 0, 0, 0, 0, 0], [-Inf], [-Inf], [-Inf]]
    torqueMax = [[0, 0, 0, 0, 0, 0], [Inf], [Inf], [Inf]]
    motionCstr = tasks.qp.MotionConstr(mbs, 0, tasks.TorqueBound(torqueMin, torqueMax))
    plCstr = tasks.qp.PositiveLambda()
    contCstrAcc = tasks.qp.ContactAccConstr()

    motionCstr.addToSolver(solver)
    contCstrAcc.addToSolver(solver)
    plCstr.addToSolver(solver)

    points = [ eigen.Vector3d(0.1, 0.1, 0),
               eigen.Vector3d(-0.1, 0.1, 0),
               eigen.Vector3d(-0.1, -0.1, 0),
               eigen.Vector3d(0.1, -0.1, 0) ]

    biFrames = [ sva.RotY(0*math.pi/2),
                 sva.RotY(1*math.pi/2),
                 sva.RotY(2*math.pi/2),
                 sva.RotY(3*math.pi/2) ]

    uni = [ tasks.qp.UnilateralContact(0, 1, "b0", "b0", points, eigen.Matrix3d.Identity(), sva.PTransformd.Identity(), 3, 0.7) ]
    bi = [ tasks.qp.BilateralContact(0, 1, "b0", "b0", points, biFrames, sva.PTransformd.Identity(), 3, 0.7) ]

    solver.nrVars(mbs, uni, [])
    solver.updateConstrSize()
    self.assertEqual(solver.nrVars(), 9 + 4*3)

    mbcs[0] = mbcInit
    self.assertFalse(solver.solve(mbs, mbcs))

    solver.nrVars(mbs, [], bi)
    solver.updateConstrSize()
    self.assertEqual(solver.nrVars(), 9 + 4*3)

    mbcs[0] = mbcInit
    for i in range(10):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.001)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    plCstr.removeFromSolver(solver)
    contCstrAcc.removeFromSolver(solver)
    motionCstr.removeFromSolver(solver)

def compute6dError(b1, b2):
  assert(isinstance(b1, sva.PTransformd) and isinstance(b2, sva.PTransformd))
  error = eigen.Vector6d.Zero()
  rError = sva.rotationError(b1.rotation(), b2.rotation())
  tError = b1.translation() - b2.translation()
  for i in range(3):
    error[i] = rError[i]
    error[i+3] = tError[i]
  return error

def compute6dErrorInB1(b1, b2):
  assert(isinstance(b1, sva.PTransformd) and isinstance(b2, sva.PTransformd))
  error = sva.MotionVecd(sva.rotationError(b1.rotation(), b2.rotation()), b1.translation() - b2.translation())
  return (sva.PTransformd(b1.rotation())*error).vector()

def computeDofError(b1, b2, dof):
  assert(isinstance(b1, sva.PTransformd) and isinstance(b2, sva.PTransformd) and isinstance(dof, eigen.MatrixXd))
  return (dof*compute6dError(b1,b2)).norm()

def computeDofErrorInB1(b1, b2, dof):
  assert(isinstance(b1, sva.PTransformd) and isinstance(b2, sva.PTransformd) and isinstance(dof, eigen.MatrixXd))
  error = compute6dErrorInB1(b1,b2)
  return (dof*error).norm()

class TestQPDofContacts(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm(False)
    mbEnv, mbcEnv = arms.makeEnv()

    quat = eigen.Vector4d.Random().normalized()
    if not LEGACY:
      mbcInit.q[0] = [quat.w(), quat.x(), quat.y(), quat.z(), 0, 0, 0]
    else:
      mbcInit.q = [[quat.w(), quat.x(), quat.y(), quat.z(), 0, 0, 0]] + list(mbcInit.q)[1:]

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)
    rbdyn.forwardKinematics(mbEnv, mbcEnv)
    rbdyn.forwardVelocity(mbEnv, mbcEnv)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb, mbEnv])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit, mbcEnv])
    else:
      mbs = [mb, mbEnv]
      mbcs = [rbdyn.MultiBodyConfig(mbcInit), mbcEnv]

    solver = tasks.qp.QPSolver()
    X_b1_cf = sva.PTransformd(eigen.Quaterniond(eigen.Vector4d.Random().normalized()), eigen.Vector3d.Random())

    contCstrSpeed = tasks.qp.ContactPosConstr(0.005)
    # Target in cf coordinate is transform into world frame
    if not LEGACY:
      posTask = tasks.qp.PositionTask(mbs, 0, "b0", (X_b1_cf*mbcInit.bodyPosW[0]).rotation().transpose()*eigen.Vector3d(1, 1, -1))
    else:
      posTask = tasks.qp.PositionTask(mbs, 0, "b0", (X_b1_cf*list(mbcInit.bodyPosW)[0]).rotation().transpose()*eigen.Vector3d(1, 1, -1))
    posTaskSp = tasks.qp.SetPointTask(mbs, 0, posTask, 10, 1)
    # Rotation in cf coordinate
    if not LEGACY:
      oriTask = tasks.qp.OrientationTask(mbs, 0, "b0", (X_b1_cf*mbcInit.bodyPosW[0]).rotation().transpose()*sva.RotY(0.1)*sva.RotX(0.5))
    else:
      oriTask = tasks.qp.OrientationTask(mbs, 0, "b0", (X_b1_cf*list(mbcInit.bodyPosW)[0]).rotation().transpose()*sva.RotY(0.1)*sva.RotX(0.5))
    oriTaskSp = tasks.qp.SetPointTask(mbs, 0, oriTask, 10, 1)

    contCstrSpeed.addToSolver(solver)
    solver.addTask(posTaskSp)

    points = [
      eigen.Vector3d(0.1, 0.1, 0),
      eigen.Vector3d(-0.1, 0.1, 0),
      eigen.Vector3d(-0.1, -0.1, 0),
      eigen.Vector3d(0.1, -0.1, 0)
    ]

    if not LEGACY:
      X_b1_b2 = sva.PTransformd(mbcEnv.bodyPosW[0]*mbcInit.bodyPosW[0].inv())
    else:
      X_b1_b2 = sva.PTransformd(list(mbcEnv.bodyPosW)[0]*list(mbcInit.bodyPosW)[0].inv())
    uni = [ tasks.qp.UnilateralContact(0, 1, "b0", "b0", points, eigen.Matrix3d.Identity(), X_b1_b2, 3, 0.7, X_b1_cf) ]

    # contactDof must be provided in r1BodyId frame
    contactDof = eigen.MatrixXd.Zero(5,6)
    if not LEGACY:
      for i in range(5):
        contactDof[i,i] = 1
    else:
      for i in range(5):
        contactDof.coeff(i,i,1)

    # test Z free
    contCstrSpeed.addDofContact(tasks.qp.ContactId(0, 1, "b0", "b0"), contactDof)
    solver.nrVars(mbs, uni, [])
    solver.updateConstrSize()

    #mbcs[0] = mbcInit
    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    for i in range(100):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.005)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    if not LEGACY:
      self.assertAlmostEqual(computeDofErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0], X_b1_cf*mbcInit.bodyPosW[0], contactDof), 0, delta = 1e-6)
      self.assertGreater(compute6dErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0], X_b1_cf*mbcInit.bodyPosW[0])[5]**2, 0.1)
    else:
      pass # These tests cannot be done on the legacy bindings as MatrixXd*Vector6d is not possible

    # Test Y free and updateDofContacts
    if not LEGACY:
      contactDof.setZero()
      for i in range(4):
        contactDof[i,i] = 1
      contactDof[4,5] = 1
    else:
      contactDof = eigen.MatrixXd.Zero(5,6)
      for i in range(4):
        contactDof.coeff(i, i, 1)
      contactDof.coeff(4, 5, 1)
    contCstrSpeed.resetDofContacts()
    contCstrSpeed.addDofContact(tasks.qp.ContactId(0, 1, "b0", "b0"), contactDof)
    contCstrSpeed.updateDofContacts()

    #mbcs[0] = mbcInit
    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    for i in range(100):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.005)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    if not LEGACY:
      self.assertAlmostEqual(computeDofErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0], X_b1_cf*mbcInit.bodyPosW[0], contactDof), 0, delta = 1e-6)
      self.assertGreater(compute6dErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0], X_b1_cf*mbcInit.bodyPosW[0])[4]**2, 0.1)
    else:
      pass # These tests cannot be done on the legacy bindings as MatrixXd*Vector6d is not possible

    # Test WX Free and updateDofContacts
    if not LEGACY:
      contactDof.setZero()
      for i in range(5):
        contactDof[i,i+1] = 1
    else:
      contactDof = eigen.MatrixXd.Zero(5,6)
      for i in range(5):
        contactDof.coeff(i, i+1, 1)
    contCstrSpeed.resetDofContacts()
    contCstrSpeed.addDofContact(tasks.qp.ContactId(0, 1, "b0", "b0"), contactDof)
    contCstrSpeed.updateDofContacts()

    # Add the orientation task in cf coordinate
    solver.addTask(oriTaskSp)
    solver.updateTasksNrVars(mbs)
    #mbcs[0] = mbcInit
    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    for i in range(1000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.005)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    if not LEGACY:
      self.assertAlmostEqual(computeDofErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0], X_b1_cf*mbcInit.bodyPosW[0], contactDof), 0, delta = 1e-4)
      self.assertGreater(abs(compute6dErrorInB1(X_b1_cf*mbcs[0].bodyPosW[0], X_b1_cf*mbcInit.bodyPosW[0])[0]), 0.1)
    else:
      pass # These tests cannot be done on the legacy bindings as MatrixXd*Vector6d is not possible

class TestQPBoundedSpeed(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm()

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit])
    else:
      mbs = [mb]
      mbcs = [mbcInit]

    solver = tasks.qp.QPSolver()
    solver.solver("QLD")

    bodyName = "b3"
    bodyIndex = mb.bodyIndexByName(bodyName)
    bodyPoint= sva.PTransformd(eigen.Vector3d(0, 0.1, 0))

    constSpeed = tasks.qp.BoundedSpeedConstr(mbs, 0, 0.005)
    postureTask = tasks.qp.PostureTask(mbs, 0, [[], [0], [0], [0]], 1, 0.01)
    posTask = tasks.qp.PositionTask(mbs, 0, bodyName, eigen.Vector3d(1, -1, 1), bodyPoint.translation())
    posTaskSp = tasks.qp.SetPointTask(mbs, 0, posTask, 20, 1)
    dof  = eigen.MatrixXd.Zero(1,6)
    speed = eigen.VectorXd.Zero(1)

    if not LEGACY:
      dof[0,3] = 1
    else:
      dof.coeff(0, 3, 1)
    constSpeed.addBoundedSpeed(mbs, bodyName, bodyPoint.translation(), dof, speed)
    self.assertEqual(constSpeed.nrBoundedSpeeds(), 1)

    constSpeed.addToSolver(solver)
    solver.addTask(postureTask)
    solver.addTask(posTaskSp)

    solver.nrVars(mbs, [], [])
    solver.updateConstrSize()

    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    for i in range(100):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.005)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    if not LEGACY:
      initPos = sva.PTransformd(bodyPoint*mbcInit.bodyPosW[bodyIndex])
      finalPos = sva.PTransformd(bodyPoint*mbcs[0].bodyPosW[bodyIndex])
      self.assertAlmostEqual(computeDofError(finalPos, initPos, dof), 0, delta = 1e-6)
      self.assertGreater(compute6dError(finalPos, initPos)[4]**2, 0.1)
    else:
      pass # See previous test

    # Same test but with Z axis
    self.assertTrue(constSpeed.removeBoundedSpeed(bodyName))
    constSpeed.updateBoundedSpeeds()
    self.assertEqual(constSpeed.nrBoundedSpeeds(), 0)
    self.assertEqual(constSpeed.maxGenInEq(), 0)

    # Must resize constraint matrix since nrMaxInEq has changed
    solver.updateConstrSize()
    if not LEGACY:
      dof.setZero()
      dof[0, 5] = 1
    else:
      dof = eigen.MatrixXd.Zero(1,6)
      dof.coeff(0, 5, 1)
    constSpeed.addBoundedSpeed(mbs, bodyName, bodyPoint.translation(), dof, speed)
    constSpeed.updateBoundedSpeeds()
    self.assertEqual(constSpeed.nrBoundedSpeeds(), 1)
    self.assertEqual(constSpeed.maxGenInEq(), 1)
    solver.updateConstrSize()

    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    for i in range(100):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.005)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    if not LEGACY:
      initPos = sva.PTransformd(bodyPoint*mbcInit.bodyPosW[bodyIndex])
      finalPos = sva.PTransformd(bodyPoint*mbcs[0].bodyPosW[bodyIndex])
      self.assertAlmostEqual(computeDofError(finalPos, initPos, dof), 0, delta = 1e-6)
      self.assertGreater(compute6dError(finalPos, initPos)[4]**2, 0.1)
    else:
      pass # See previous test

class TestMomentumTask(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm()

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit])
    else:
      mbs = [mb]
      mbcs = [mbcInit]

    solver = tasks.qp.QPSolver()

    solver.nrVars(mbs, [], [])
    solver.updateConstrSize()

    momTarget = sva.ForceVecd(eigen.Vector3d(1,1,1), eigen.Vector3d(0,0,0))

    momTask = tasks.qp.MomentumTask(mbs, 0, momTarget)
    momTaskSp = tasks.qp.SetPointTask(mbs, 0, momTask, 10, 1)

    solver.addTask(momTaskSp)
    self.assertEqual(solver.nrTasks(), 1)

    self.assertTrue(solver.solve(mbs, mbcs))

class TestSurfaceOrientationTask(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm()

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit])
    else:
      mbs = [mb]
      mbcs = [mbcInit]

    solver = tasks.qp.QPSolver()

    solver.nrVars(mbs, [], [])
    solver.updateConstrSize()

    surfOriTask = tasks.qp.SurfaceOrientationTask(mbs, 0, u'b3', sva.RotZ(math.pi/2), sva.PTransformd(eigen.Vector3d(0,0,0)))
    surfOriTaskSp = tasks.qp.SetPointTask(mbs, 0, surfOriTask, 10, 1)

    solver.addTask(surfOriTaskSp)
    self.assertEqual(solver.nrTasks(), 1)

    self.assertTrue(solver.solve(mbs, mbcs))

class TestQPCoMPlane(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm()

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit])
    else:
      mbs = [mb]
      mbcs = [mbcInit]

    solver = tasks.qp.QPSolver()

    bodyI = mb.bodyIndexByName("b3")
    if not LEGACY:
      initPos = eigen.Vector3d(mbcInit.bodyPosW[bodyI].translation())
    else:
      initPos = eigen.Vector3d(list(mbcInit.bodyPosW)[bodyI].translation())
    posTask = tasks.qp.PositionTask(mbs, 0, "b3", initPos)
    posTaskSp = tasks.qp.SetPointTask(mbs, 0, posTask, 50, 1)
    n1 = eigen.Vector3d(0, 0, -1)
    n2 = eigen.Vector3d(0, 0, 1)
    p1 = eigen.Vector3d(0, 0, 0.1)
    p2 = eigen.Vector3d(0, 0, -0.1)
    offset1 = -n1.dot(p1)
    offset2 = -n2.dot(p2)

    comPlaneConstr = tasks.qp.CoMIncPlaneConstr(mbs, 0, 0.001)
    planeId1 = 10
    planeId2 = 20
    comPlaneConstr.addPlane(planeId1, n1, offset1, 0.01, 0.005, 0, 0.1)
    comPlaneConstr.addPlane(planeId2, n2, offset2, 0.01, 0.005, 0.1, 0)
    self.assertEqual(comPlaneConstr.nrPlanes(), 2)

    comPlaneConstr.addToSolver(solver)
    self.assertEqual(solver.nrInequalityConstraints(), 1)
    self.assertEqual(solver.nrConstraints(), 1)

    solver.nrVars(mbs, [], [])
    solver.updateConstrSize()

    solver.addTask(posTaskSp)
    self.assertEqual(solver.nrTasks(), 1)

    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    for i in range(1000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.001)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

      # Check if the CoM is on the good side of each plane
      com = rbdyn.computeCoM(mbs[0], mbcs[0])
      dist1 = n1.dot(com) + offset1
      dist2 = n2.dot(com) + offset2
      self.assertGreater(dist1, 0.005)
      self.assertGreater(dist2, 0.005)

    # Inverse rotation side
    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    posTask.position(initPos)
    for i in range(1000):
      posTask.position(sva.RotX(-0.01)*posTask.position())
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.001)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

      # Check if the CoM is on the good side of each plane
      com = rbdyn.computeCoM(mbs[0], mbcs[0])
      dist1 = n1.dot(com) + offset1
      dist2 = n2.dot(com) + offset2
      self.assertGreater(dist1, 0.005)
      self.assertGreater(dist2, 0.005)

    comPlaneConstr.rmPlane(planeId1)
    comPlaneConstr.rmPlane(planeId2)
    self.assertEqual(comPlaneConstr.nrPlanes(), 0)

    comPlaneConstr.removeFromSolver(solver)
    self.assertEqual(solver.nrInequalityConstraints(), 0)
    self.assertEqual(solver.nrConstraints(), 0)

    solver.removeTask(posTaskSp)
    self.assertEqual(solver.nrTasks(), 0)

class TestJointsSelector(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm(False)

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit])
    else:
      mbs = [mb]
      mbcs = [mbcInit]

    solver = tasks.qp.QPSolver()

    pt = tasks.qp.PositionTask(mbs, 0, "b3", eigen.Vector3d.Zero())
    # Construct two JointSelector that should have the same joints
    js1 = tasks.qp.JointsSelector.ActiveJoints(mbs, 0, pt, ["j2", "j1"])
    js2 = tasks.qp.JointsSelector.UnactiveJoints(mbs, 0, pt, ["j0", "Root"])
    js1Sp = tasks.qp.SetPointTask(mbs, 0, js1, 1, 1)
    js2Sp = tasks.qp.SetPointTask(mbs, 0, js2, 1, 1)

    if not LEGACY:
      self.assertEqual(len(js1.selectedJoints()), 2)
      self.assertEqual(len(js2.selectedJoints()), 2)
      # Check that they have the same joints selected
      for i in range(2):
        self.assertEqual(js1.selectedJoints()[i].posInDof, js2.selectedJoints()[i].posInDof)
        self.assertEqual(js1.selectedJoints()[i].dof, js2.selectedJoints()[i].dof)
      # Check joints are sorted
      self.assertLess(js1.selectedJoints()[0].posInDof, js1.selectedJoints()[1].posInDof)
    else:
      pass # selectedJoints not implemented in legacy bindings

    solver.addTask(js1Sp)
    solver.addTask(js2Sp)
    self.assertEqual(solver.nrTasks(), 2)

    solver.nrVars(mbs, [], [])
    solver.updateConstrSize()

    self.assertTrue(solver.solve(mbs, mbcs))

    if not LEGACY:
      # Jacobian first column should be zero
      self.assertTrue(js1.jac().block(0, 0, 3, 7), eigen.MatrixXd.Zero(3,7))
      # Matrix should be equals
      self.assertEqual(js1.jac(), js2.jac())
      self.assertEqual(js1.eval(), js2.eval())
      self.assertEqual(js1.speed(), js2.speed())
      self.assertEqual(js1.normalAcc(), js2.normalAcc())

class TestQPTransformTask(unittest.TestCase):
  def test(self):
    mb, mbcInit = arms.makeZXZArm(False)

    rbdyn.forwardKinematics(mb, mbcInit)
    rbdyn.forwardVelocity(mb, mbcInit)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb])
      mbcs = rbdyn.MultiBodyConfigVector([mbcInit])
    else:
      mbs = [mb]
      mbcs = [mbcInit]

    solver = tasks.qp.QPSolver()

    solver.nrVars(mbs, [], [])
    self.assertEqual(solver.nrVars(), 9)

    solver.updateConstrSize()

    # Test TransformTask
    X_b_s = sva.PTransformd(eigen.Quaterniond(eigen.Vector4d.Random().normalized()), eigen.Vector3d.Random())
    X_0_t = sva.PTransformd(eigen.Quaterniond(eigen.Vector4d.Random().normalized()), eigen.Vector3d.Random())
    E_0_c = eigen.Quaterniond(eigen.Vector4d.Random().normalized())

    transTask = tasks.qp.TransformTask(mbs, 0, "b3", X_0_t, X_b_s, E_0_c.matrix())
    dimW = eigen.VectorXd.Zero(6)
    for i in range(6):
      dimW[i] = 1
    transTaskSp = tasks.qp.SetPointTask(mbs, 0, transTask, 10, dimW, 100)
    postureTask = tasks.qp.PostureTask(mbs, 0, [[],[0],[0],[0]], 0.5, 10)

    solver.addTask(transTaskSp)
    solver.addTask(postureTask)
    solver.updateTasksNrVars(mbs)
    self.assertEqual(solver.nrTasks(), 2)

    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    for i in range(10000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.001)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    self.assertAlmostEqual(transTask.eval().norm(), 0, delta = 1e-5)

    solver.removeTask(transTaskSp)
    self.assertEqual(solver.nrTasks(), 1)

    # Test SurfaceTransformTask
    surfTransTask = tasks.qp.SurfaceTransformTask(mbs, 0, "b3", X_0_t, X_b_s)
    surfTransTaskSp = tasks.qp.SetPointTask(mbs, 0, surfTransTask, 10, dimW, 100)
    solver.addTask(surfTransTaskSp)
    solver.updateTasksNrVars(mbs)
    self.assertEqual(solver.nrTasks(), 2)

    mbcs[0] = rbdyn.MultiBodyConfig(mbcInit)
    for i in range(10000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
      rbdyn.eulerIntegration(mbs[0], mbcs[0], 0.001)
      rbdyn.forwardKinematics(mbs[0], mbcs[0])
      rbdyn.forwardVelocity(mbs[0], mbcs[0])

    self.assertAlmostEqual(surfTransTask.eval().norm(), 0, delta = 1e-5)

    solver.removeTask(surfTransTaskSp)
    solver.removeTask(postureTask)
    self.assertEqual(solver.nrTasks(), 0)

if __name__ == "__main__":
  if not LEGACY:
    print("Running tests with new bindings")
  else:
    print("Running tests with legacy bindings")
  suite = unittest.TestSuite()
  suite.addTest(TestFrictionCone('test_cone1'))
  suite.addTest(TestFrictionCone('test_cone2'))
  suite.addTest(TestQPTask('test_position_task'))
  suite.addTest(TestQPTask('test_orientation_task'))
  suite.addTest(TestQPTask('test_posture_task'))
  suite.addTest(TestQPTask('test_com_task'))
  suite.addTest(TestQPTask('test_lin_velocity_task'))
  suite.addTest(TestQPConstr('test_contact_acc_constr'))
  suite.addTest(TestQPConstr('test_contact_speed_constr'))
  suite.addTest(TestQPConstr('test_motion_constr'))
  suite.addTest(TestQPConstr('test_motion_constr_w_contact'))
  suite.addTest(TestQPJointLimits('test_joint_limits'))
  suite.addTest(TestQPJointLimits('test_damper_joint_limits'))
  suite.addTest(TestQPTorqueLimits('test_motion_constr'))
  suite.addTest(TestQPTorqueLimits('test_motion_poly_constr'))
  suite.addTest(TestQPAutoColl('test'))
  suite.addTest(TestQPStaticEnvColl('test'))
  suite.addTest(TestQPBilatContact('test'))
  suite.addTest(TestQPDofContacts('test'))
  suite.addTest(TestQPBoundedSpeed('test'))
  suite.addTest(TestMomentumTask('test'))
  suite.addTest(TestSurfaceOrientationTask('test'))
  suite.addTest(TestQPCoMPlane('test'))
  suite.addTest(TestJointsSelector('test'))
  suite.addTest(TestQPTransformTask('test'))
  unittest.TextTestRunner(verbosity=2).run(suite)
