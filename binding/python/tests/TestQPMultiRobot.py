#!/usr/bin/env python

#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from __future__ import print_function

LEGACY = False

import unittest

import math
import platform
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

class TestTwoArmContact(unittest.TestCase):
  def test(self):
    mb1, mbc1Init = arms.makeZXZArm()
    mb2, mbc2Init = arms.makeZXZArm()

    rbdyn.forwardKinematics(mb1, mbc1Init)
    rbdyn.forwardVelocity(mb1, mbc1Init)
    rbdyn.forwardKinematics(mb2, mbc2Init)
    rbdyn.forwardVelocity(mb2, mbc2Init)

    if not LEGACY:
      X_0_b1 = sva.PTransformd(mbc1Init.bodyPosW[-1])
      X_0_b2 = sva.PTransformd(mbc2Init.bodyPosW[-1])
    else:
      X_0_b1 = sva.PTransformd(list(mbc1Init.bodyPosW)[-1])
      X_0_b2 = sva.PTransformd(list(mbc2Init.bodyPosW)[-1])
    X_b1_b2 = X_0_b2*X_0_b1.inv()

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb1, mb2])
      mbcs = rbdyn.MultiBodyConfigVector([mbc1Init, mbc2Init])
    else:
      mbs = [mb1, mb2]
      mbcs = [rbdyn.MultiBodyConfig(mbc1Init), rbdyn.MultiBodyConfig(mbc2Init)]

    # Test ContactAccConstr contraint and test PositionTask on the second robot
    solver = tasks.qp.QPSolver()

    contVec = [ tasks.qp.UnilateralContact(0, 1, "b3", "b3", [eigen.Vector3d.Zero()], sva.RotX(math.pi/2), X_b1_b2, 3, math.tan(math.pi/4)) ]

    oriD = sva.RotZ(math.pi/4)
    if platform.architecture()[0] == "32bit" and solver.solver() == b"QLD":
      oriD = sva.RotZ(0)
    if not LEGACY:
      posD = oriD*mbc2Init.bodyPosW[-1].translation()
    else:
      posD = oriD*list(mbc2Init.bodyPosW)[-1].translation()
    posTask = tasks.qp.PositionTask(mbs, 1, "b3", posD)
    posTaskSp = tasks.qp.SetPointTask(mbs, 1, posTask, 1000, 1)

    contCstrAcc = tasks.qp.ContactAccConstr()

    contCstrAcc.addToSolver(solver)
    solver.addTask(posTaskSp)

    solver.nrVars(mbs, contVec, [])
    solver.updateConstrSize()

    self.assertEqual(solver.nrVars(), 3 + 3 + 3)

    for i in range(1000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
        solver.updateMbc(mbcs[1], 1)
      for i in range(2):
        rbdyn.eulerIntegration(mbs[i], mbcs[i], 0.001)
        rbdyn.forwardKinematics(mbs[i], mbcs[i])
        rbdyn.forwardVelocity(mbs[i], mbcs[i])

      # Check that the link hold
      if not LEGACY:
        X_0_b1_post = mbcs[0].bodyPosW[-1]
        X_0_b2_post = mbcs[1].bodyPosW[-1]
      else:
        X_0_b1_post = list(mbcs[0].bodyPosW)[-1]
        X_0_b2_post = list(mbcs[1].bodyPosW)[-1]
      X_b1_b2_post = X_0_b2*X_0_b1.inv()
      self.assertAlmostEqual((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 0, delta = 1e-5)

    self.assertAlmostEqual(posTask.eval().norm(), 0, delta = 1e-5)

    contCstrAcc.removeFromSolver(solver)
    solver.removeTask(posTaskSp)

    # Test ContactSpeedConstr constraint and OrientationTask on the second robot
    if not LEGACY:
      mbcs = rbdyn.MultiBodyConfigVector([mbc1Init, mbc2Init])
    else:
      mbcs = [rbdyn.MultiBodyConfig(mbc1Init), rbdyn.MultiBodyConfig(mbc2Init)]
    oriTask = tasks.qp.OrientationTask(mbs, 1, "b3", oriD)
    oriTaskSp = tasks.qp.SetPointTask(mbs, 1, oriTask, 1000, 1)

    contCstrSpeed = tasks.qp.ContactSpeedConstr(0.001)
    contCstrSpeed.addToSolver(solver)
    solver.addTask(oriTaskSp)

    solver.nrVars(mbs, contVec, [])
    solver.updateConstrSize()
    for i in range(1000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
        solver.updateMbc(mbcs[1], 1)
      for i in range(2):
        rbdyn.eulerIntegration(mbs[i], mbcs[i], 0.001)
        rbdyn.forwardKinematics(mbs[i], mbcs[i])
        rbdyn.forwardVelocity(mbs[i], mbcs[i])
      # Check that the link hold
      if not LEGACY:
        X_0_b1_post = mbcs[0].bodyPosW[-1]
        X_0_b2_post = mbcs[1].bodyPosW[-1]
      else:
        X_0_b1_post = list(mbcs[0].bodyPosW)[-1]
        X_0_b2_post = list(mbcs[1].bodyPosW)[-1]
      X_b1_b2_post = X_0_b2*X_0_b1.inv()
      self.assertAlmostEqual((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 0, delta = 1e-5)

    self.assertAlmostEqual(oriTask.eval().norm(), 0, delta = 1e-5)


# Test Motion constraint
# We setup two arm, one with a fixed base and the second
# with a freebase put on the body b3 of the first robot.
# First we launch an impossible motion to check the dynamics
# After we try with an unilateral contact
# Then we try with a bilateral contact.
class TestTwoArmDDynamicContact(unittest.TestCase):
  def test(self):
    mb1, mbc1Init = arms.makeZXZArm()
    rbdyn.forwardKinematics(mb1, mbc1Init)
    rbdyn.forwardVelocity(mb1, mbc1Init)

    mb2, mbc2Init = arms.makeZXZArm(False)
    if not LEGACY:
      mb2InitPos = mbc1Init.bodyPosW[-1].translation()
    else:
      mb2InitPos = list(mbc1Init.bodyPosW)[-1].translation()
    mb2InitOri = eigen.Quaterniond(sva.RotY(math.pi/2))
    if not LEGACY:
      mbc2Init.q[0] = [mb2InitOri.w(), mb2InitOri.x(), mb2InitOri.y(), mb2InitOri.z(), mb2InitPos.x(), mb2InitPos.y() + 1, mb2InitPos.z()]
      mbc2Init.q[0] = [mb2InitOri.w(), mb2InitOri.x(), mb2InitOri.y(), mb2InitOri.z(), mb2InitPos.x(), mb2InitPos.y() + 1, mb2InitPos.z()]
    rbdyn.forwardKinematics(mb2, mbc2Init)
    rbdyn.forwardVelocity(mb2, mbc2Init)

    if not LEGACY:
      X_0_b1 = sva.PTransformd(mbc1Init.bodyPosW[-1])
      X_0_b2 = sva.PTransformd(mbc2Init.bodyPosW[-1])
    else:
      X_0_b1 = sva.PTransformd(list(mbc1Init.bodyPosW)[-1])
      X_0_b2 = sva.PTransformd(list(mbc2Init.bodyPosW)[-1])
    X_b1_b2 = X_0_b2*X_0_b1.inv()

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb1, mb2])
      mbcs = rbdyn.MultiBodyConfigVector([mbc1Init, mbc2Init])
    else:
      mbs = [mb1, mb2]
      mbcs = [rbdyn.MultiBodyConfig(mbc1Init), rbdyn.MultiBodyConfig(mbc2Init)]

    # Test ContactAccConstr constraint and PositionTask on the second robot
    solver = tasks.qp.QPSolver()

    points = [
      eigen.Vector3d(0.1, 0, 0.1),
      eigen.Vector3d(0.1, 0, -0.1),
      eigen.Vector3d(-0.1, 0, -0.1),
      eigen.Vector3d(-0.1, 0, 0.1),
    ]

    biPoints = [
      eigen.Vector3d.Zero(),
      eigen.Vector3d.Zero(),
      eigen.Vector3d.Zero(),
      eigen.Vector3d.Zero(),
    ]

    nrGen = 4
    biFrames = [
      sva.RotX(math.pi/4),
      sva.RotX(3*math.pi/4),
      sva.RotX(math.pi/4)*sva.RotY(math.pi/2),
      sva.RotX(3*math.pi/4)*sva.RotY(math.pi/2),
    ]

    # The fixed robot can pull the other
    contVecFail = [ tasks.qp.UnilateralContact(0, 1, "b3", "b0", points, sva.RotX(-math.pi/2), X_b1_b2, nrGen, 0.7) ]

    # The fixed robot can push the other
    contVec = [ tasks.qp.UnilateralContact(0, 1, "b3", "b0", points, sva.RotX(math.pi/2), X_b1_b2, nrGen, 0.7) ]

    # The fixed robot has non coplanar force apply on the other
    contVecBi = [ tasks.qp.BilateralContact(tasks.qp.ContactId(0, 1, "b3", "b0"), biPoints, biFrames, X_b1_b2, nrGen, 1) ]

    if not LEGACY:
      posture1Task = tasks.qp.PostureTask(mbs, 0, mbc1Init.q, 2, 1)
      posture2Task = tasks.qp.PostureTask(mbs, 1, mbc2Init.q, 2, 1)
    else:
      posture1Task = tasks.qp.PostureTask(mbs, 0, rbdList(mbc1Init.q), 2, 1)
      posture2Task = tasks.qp.PostureTask(mbs, 1, rbdList(mbc2Init.q), 2, 1)

    contCstrSpeed = tasks.qp.ContactSpeedConstr(0.001)

    Inf = float("inf")
    torqueMin1 = [[], [-Inf], [-Inf], [-Inf]]
    torqueMax1 = [[], [Inf], [Inf], [Inf]]
    torqueMin2 = [[0,0,0,0,0,0], [-Inf], [-Inf], [-Inf]]
    torqueMax2 = [[0,0,0,0,0,0], [Inf], [Inf], [Inf]]
    motion1 = tasks.qp.MotionConstr(mbs, 0, tasks.TorqueBound(torqueMin1, torqueMax1))
    motion2 = tasks.qp.MotionConstr(mbs, 1, tasks.TorqueBound(torqueMin2, torqueMax2))
    plCstr = tasks.qp.PositiveLambda()

    motion1.addToSolver(solver)
    motion2.addToSolver(solver)
    plCstr.addToSolver(solver)

    contCstrSpeed.addToSolver(solver)
    solver.addTask(posture1Task)
    solver.addTask(posture2Task)

    # Check the impossible motion
    solver.nrVars(mbs, contVecFail, [])
    solver.updateConstrSize()
    self.assertEqual(solver.nrVars(), 3 + 9 + 4*nrGen)
    self.assertFalse(solver.solve(mbs, mbcs))

    # Check the unilateral motion
    if not LEGACY:
      mbcs = rbdyn.MultiBodyConfigVector([mbc1Init, mbc2Init])
    else:
      mbcs = [rbdyn.MultiBodyConfig(mbc1Init), rbdyn.MultiBodyConfig(mbc2Init)]
    solver.nrVars(mbs, contVec, [])
    solver.updateConstrSize()
    for i in range(1000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
        solver.updateMbc(mbcs[1], 1)
      for i in range(2):
        rbdyn.eulerIntegration(mbs[i], mbcs[i], 0.001)
        rbdyn.forwardKinematics(mbs[i], mbcs[i])
        rbdyn.forwardVelocity(mbs[i], mbcs[i])

      # Check that the link hold
      if not LEGACY:
        X_0_b1_post = mbcs[0].bodyPosW[-1]
        X_0_b2_post = mbcs[1].bodyPosW[-1]
      else:
        X_0_b1_post = list(mbcs[0].bodyPosW)[-1]
        X_0_b2_post = list(mbcs[1].bodyPosW)[-1]
      X_b1_b2_post = X_0_b2*X_0_b1.inv()
      self.assertAlmostEqual((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 0, delta = 1e-5)

      # Force in the world frame must be the same
      f1 = contVec[0].force(solver.lambdaVec(0), contVec[0].r1Cone)
      f2 = contVec[0].force(solver.lambdaVec(0), contVec[0].r2Cone)
      self.assertAlmostEqual((f1+f2).norm(), 0, delta = 1e-5)

    # Check the bilateral motion
    if not LEGACY:
      mbcs = rbdyn.MultiBodyConfigVector([mbc1Init, mbc2Init])
    else:
      mbcs = [rbdyn.MultiBodyConfig(mbc1Init), rbdyn.MultiBodyConfig(mbc2Init)]
    solver.nrVars(mbs, contVec, [])
    solver.updateConstrSize()
    self.assertEqual(solver.nrVars(), 3 + 9 + 4*nrGen)
    for i in range(1000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
        solver.updateMbc(mbcs[1], 1)
      for i in range(2):
        rbdyn.eulerIntegration(mbs[i], mbcs[i], 0.001)
        rbdyn.forwardKinematics(mbs[i], mbcs[i])
        rbdyn.forwardVelocity(mbs[i], mbcs[i])

      # Check that the link hold
      if not LEGACY:
        X_0_b1_post = mbcs[0].bodyPosW[-1]
        X_0_b2_post = mbcs[1].bodyPosW[-1]
      else:
        X_0_b1_post = list(mbcs[0].bodyPosW)[-1]
        X_0_b2_post = list(mbcs[1].bodyPosW)[-1]
      X_b1_b2_post = X_0_b2*X_0_b1.inv()
      self.assertAlmostEqual((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 0, delta = 1e-5)

      # Force in the world frame must be the same
      f1 = contVec[0].force(solver.lambdaVec(0), contVec[0].r1Cone)
      f2 = contVec[0].force(solver.lambdaVec(0), contVec[0].r2Cone)
      self.assertAlmostEqual((f1+f2).norm(), 0, delta = 1e-5)

    plCstr.removeFromSolver(solver)
    motion2.removeFromSolver(solver)
    motion1.removeFromSolver(solver)
    contCstrSpeed.removeFromSolver(solver)

    solver.removeTask(posture1Task)
    solver.removeTask(posture2Task)

# Test the MultiCoMTask
# We try to move the CoM of two arm at a specific position
class TestTwoArmMultiCoM(unittest.TestCase):
  def test(self):
    mb1, mbc1Init = arms.makeZXZArm(True, sva.PTransformd(eigen.Vector3d(-0.5, 0, 0)))
    rbdyn.forwardKinematics(mb1, mbc1Init)
    rbdyn.forwardVelocity(mb1, mbc1Init)

    mb2, mbc2Init = arms.makeZXZArm(True, sva.PTransformd(eigen.Vector3d(0.5, 0, 0)))
    rbdyn.forwardKinematics(mb2, mbc2Init)
    rbdyn.forwardVelocity(mb2, mbc2Init)

    if not LEGACY:
      X_0_b1 = sva.PTransformd(mbc1Init.bodyPosW[-1])
      X_0_b2 = sva.PTransformd(mbc2Init.bodyPosW[-1])
    else:
      X_0_b1 = sva.PTransformd(list(mbc1Init.bodyPosW)[-1])
      X_0_b2 = sva.PTransformd(list(mbc2Init.bodyPosW)[-1])
    X_b1_b2 = X_0_b2*X_0_b1.inv()

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb1, mb2])
      mbcs = rbdyn.MultiBodyConfigVector([mbc1Init, mbc2Init])
    else:
      mbs = [mb1, mb2]
      mbcs = [rbdyn.MultiBodyConfig(mbc1Init), rbdyn.MultiBodyConfig(mbc2Init)]

    nrGen = 3
    solver = tasks.qp.QPSolver()

    contVec = [ tasks.qp.UnilateralContact(0, 1, "b3", "b3", [eigen.Vector3d.Zero()], sva.RotX(math.pi/2), X_b1_b2, nrGen, 0.7) ]

    if not LEGACY:
      posture1Task = tasks.qp.PostureTask(mbs, 0, mbc1Init.q, 2, 1)
      posture2Task = tasks.qp.PostureTask(mbs, 1, mbc2Init.q, 2, 1)
    else:
      posture1Task = tasks.qp.PostureTask(mbs, 0, rbdList(mbc1Init.q), 2, 1)
      posture2Task = tasks.qp.PostureTask(mbs, 1, rbdList(mbc2Init.q), 2, 1)
    comD = (rbdyn.computeCoM(mb1, mbc1Init) + rbdyn.computeCoM(mb2, mbc2Init))/2 + eigen.Vector3d(0, 0, 0.5)
    multiCoM = tasks.qp.MultiCoMTask(mbs, [0,1], comD, 10, 500)
    multiCoM.updateInertialParameters(mbs)

    contCstrSpeed = tasks.qp.ContactSpeedConstr(0.001)

    solver.addTask(posture1Task)
    solver.addTask(posture2Task)

    solver.nrVars(mbs, contVec, [])

    solver.addTask(mbs, multiCoM)
    contCstrSpeed.addToSolver(mbs, solver)

    solver.updateConstrSize()

    self.assertEqual(solver.nrVars(), 3 + 3 + 1*nrGen)

    for i in range(2000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
        solver.updateMbc(mbcs[1], 1)
      for i in range(2):
        rbdyn.eulerIntegration(mbs[i], mbcs[i], 0.001)
        rbdyn.forwardKinematics(mbs[i], mbcs[i])
        rbdyn.forwardVelocity(mbs[i], mbcs[i])
      # Check that the link hold
      if not LEGACY:
        X_0_b1_post = mbcs[0].bodyPosW[-1]
        X_0_b2_post = mbcs[1].bodyPosW[-1]
      else:
        X_0_b1_post = list(mbcs[0].bodyPosW)[-1]
        X_0_b2_post = list(mbcs[1].bodyPosW)[-1]
      X_b1_b2_post = X_0_b2*X_0_b1.inv()
      self.assertAlmostEqual((X_b1_b2.matrix() - X_b1_b2_post.matrix()).norm(), 0, delta = 1e-5)

    self.assertAlmostEqual(multiCoM.speed().norm(), 0, delta = 1e-3)

    contCstrSpeed.removeFromSolver(solver)
    solver.removeTask(posture1Task)
    solver.removeTask(posture2Task)
    solver.removeTask(multiCoM)

# Test the MultiRobotTransformTask
# We try to set the two end effector at the same frame
class TestMultiRobotTransform(unittest.TestCase):
  def test(self):
    mb1, mbc1Init = arms.makeZXZArm(True, sva.PTransformd(sva.RotZ(-math.pi/4), eigen.Vector3d(-0.5, 0, 0)))
    rbdyn.forwardKinematics(mb1, mbc1Init)
    rbdyn.forwardVelocity(mb1, mbc1Init)

    mb2, mbc2Init = arms.makeZXZArm(False, sva.PTransformd(sva.RotZ(math.pi/2), eigen.Vector3d(0.5, 0, 0)))
    rbdyn.forwardKinematics(mb2, mbc2Init)
    rbdyn.forwardVelocity(mb2, mbc2Init)

    if not LEGACY:
      mbs = rbdyn.MultiBodyVector([mb1, mb2])
      mbcs = rbdyn.MultiBodyConfigVector([mbc1Init, mbc2Init])
    else:
      mbs = [mb1, mb2]
      mbcs = [rbdyn.MultiBodyConfig(mbc1Init), rbdyn.MultiBodyConfig(mbc2Init)]

    solver = tasks.qp.QPSolver()

    if not LEGACY:
      posture1Task = tasks.qp.PostureTask(mbs, 0, mbc1Init.q, 0.1, 10)
      posture2Task = tasks.qp.PostureTask(mbs, 1, mbc2Init.q, 0.1, 10)
    else:
      posture1Task = tasks.qp.PostureTask(mbs, 0, rbdList(mbc1Init.q), 2, 1)
      posture2Task = tasks.qp.PostureTask(mbs, 1, rbdList(mbc2Init.q), 2, 1)
    mrtt = tasks.qp.MultiRobotTransformTask(mbs, 0, 1, "b3", "b3", sva.PTransformd(sva.RotZ(-math.pi/8)), sva.PTransformd.Identity(), 100, 1000)
    if not LEGACY:
      mrtt.dimWeight(eigen.VectorXd(0, 0, 1, 1, 1, 0))
    else:
      mrtt.dimWeight(eigen.Vector6d(0, 0, 1, 1, 1, 0))

    solver.addTask(posture1Task)
    solver.addTask(posture2Task)
    solver.addTask(mrtt)

    solver.nrVars(mbs, [], [])
    solver.updateConstrSize()
    # 3 dof + 9 dof
    self.assertEqual(solver.nrVars(), 3 + 9)
    for i in range(2000):
      if not LEGACY:
        self.assertTrue(solver.solve(mbs, mbcs))
      else:
        self.assertTrue(solver.solveNoMbcUpdate(mbs, mbcs))
        solver.updateMbc(mbcs[0], 0)
        solver.updateMbc(mbcs[1], 1)
      for i in range(2):
        rbdyn.eulerIntegration(mbs[i], mbcs[i], 0.001)
        rbdyn.forwardKinematics(mbs[i], mbcs[i])
        rbdyn.forwardVelocity(mbs[i], mbcs[i])
    self.assertAlmostEqual(mrtt.eval().norm(), 0, delta = 1e-3)

    solver.removeTask(posture1Task)
    solver.removeTask(posture2Task)
    solver.removeTask(mrtt)

if __name__ == "__main__":
  if not LEGACY:
    print("Running tests with new bindings")
  else:
    print("Running tests with legacy bindings")
  suite = unittest.TestSuite()
  suite.addTest(TestTwoArmContact('test'))
  suite.addTest(TestTwoArmDDynamicContact('test'))
  suite.addTest(TestTwoArmMultiCoM('test'))
  suite.addTest(TestMultiRobotTransform('test'))
  unittest.TextTestRunner(verbosity=2).run(suite)
