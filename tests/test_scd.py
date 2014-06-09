import numpy as np

from mayavi import mlab
from tvtk.api import tvtk

from eigen3 import *
import spacevecalg as sva

import sch

def setTransform(actor, transform):
  R = transform.rotation()
  T = transform.translation()

  actor.user_transform.set_matrix((R[0], R[1], R[2], T[0],
                                   R[3], R[4], R[5], T[1],
                                   R[6], R[7], R[8], T[2],
                                     0.,   0.,   0.,   1.))


@mlab.animate(delay=33)
def anim():
  x = 0.2
  y = 0.3
  z = 0.4

  cb1 = sch.Box(x, y, z)
  cb2 = sch.Box(x, y, z)

  pair = sch.CD_Pair(cb1, cb2)

  b1s = tvtk.CubeSource(x_length=x, y_length=y, z_length=z)
  b2s = tvtk.CubeSource(x_length=x, y_length=y, z_length=z)

  b1m = tvtk.PolyDataMapper(input=b1s.output)
  b2m = tvtk.PolyDataMapper(input=b2s.output)

  b1a = tvtk.Actor(mapper=b1m)
  b2a = tvtk.Actor(mapper=b2m)

  line = tvtk.LineSource()
  lineM = tvtk.PolyDataMapper(input=line.output)
  lineA = tvtk.Actor(mapper=lineM)

  b1a.user_transform = tvtk.Transform()
  b2a.user_transform = tvtk.Transform()

  b1T = sva.PTransformd.Identity()
  b2T = sva.PTransformd(Vector3d.UnitZ()*2.)

  setTransform(b1a, b1T)
  setTransform(b2a, b2T)

  viewer = mlab.gcf()
  viewer.scene.add_actors([b1a, b2a, lineA])

  rx = ry = rz = 0.
  t = 0
  while True:
    b1T = sva.PTransformd(sva.RotZ(rz)*sva.RotY(ry)*sva.RotX(rx))
    b2T = sva.PTransformd(Vector3d(0., 0., 2. + np.sin(t)))

    setTransform(b1a, b1T)
    setTransform(b2a, b2T)

    cb1.transform(b1T)
    cb2.transform(b2T)

    p1 = Vector3d()
    p2 = Vector3d()
    pair.distance(p1, p2)

    line.point1 = list(p1)
    line.point2 = list(p2)


    rx += 0.01
    ry += 0.005
    rz += 0.002
    t += 0.05
    viewer.scene.render()
    yield



if __name__ == '__main__':
  a = anim()


