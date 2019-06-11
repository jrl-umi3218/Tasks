#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

try:
  import eigen
except ImportError:
  import eigen3 as eigen
try:
  import sva
except ImportError:
  import spacevecalg as sva
import rbdyn

def makeZXZArm(isFixed = True, X_base = sva.PTransformd.Identity()):
  mbg = rbdyn.MultiBodyGraph()

  mass = 1.0
  I = eigen.Matrix3d.Identity()
  h = eigen.Vector3d.Zero()

  rbi = sva.RBInertiad(mass, h, I)

  b0 = rbdyn.Body(rbi, "b0")
  b1 = rbdyn.Body(rbi, "b1")
  b2 = rbdyn.Body(rbi, "b2")
  b3 = rbdyn.Body(rbi, "b3")

  mbg.addBody(b0)
  mbg.addBody(b1)
  mbg.addBody(b2)
  mbg.addBody(b3)

  j0 = rbdyn.Joint(rbdyn.Joint.Rev, eigen.Vector3d.UnitZ(), True, "j0")
  j1 = rbdyn.Joint(rbdyn.Joint.Rev, eigen.Vector3d.UnitX(), True, "j1")
  j2 = rbdyn.Joint(rbdyn.Joint.Rev, eigen.Vector3d.UnitZ(), True, "j2")

  mbg.addJoint(j0)
  mbg.addJoint(j1)
  mbg.addJoint(j2)

  to = sva.PTransformd(eigen.Vector3d(0, 0.5, 0))
  _from = sva.PTransformd(eigen.Vector3d(0, 0, 0))

  mbg.linkBodies("b0", sva.PTransformd.Identity(), "b1", _from, "j0")
  mbg.linkBodies("b1", to, "b2", _from, "j1")
  mbg.linkBodies("b2", to, "b3", _from, "j2")

  mb = mbg.makeMultiBody("b0", isFixed, X_base)

  mbc = rbdyn.MultiBodyConfig(mb)
  mbc.zero(mb)

  return mb,mbc

def makeEnv():
  mbg = rbdyn.MultiBodyGraph()

  mass = 1.0
  I = eigen.Matrix3d.Identity()
  h = eigen.Vector3d.Zero()

  rbi = sva.RBInertiad(mass, h, I)

  b0 = rbdyn.Body(rbi, "b0")

  mbg.addBody(b0)

  mb = mbg.makeMultiBody("b0", True)

  mbc = rbdyn.MultiBodyConfig(mb)
  mbc.zero(mb)

  return mb,mbc
