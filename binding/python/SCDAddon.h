

// This file is part of Tasks.
//
// Tasks is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Tasks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Tasks.  If not, see <http://www.gnu.org/licenses/>.

namespace SCD
{

void transform(S_Object& obj, const sva::PTransform& t)
{
  SCD::Matrix4x4 m;
  const Eigen::Matrix3d& rot = t.rotation();
  const Eigen::Vector3d& tran = t.translation();

  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      m(i,j) = rot(i,j);
    }
  }

  m(0,3) = tran(0);
  m(1,3) = tran(1);
  m(2,3) = tran(2);

  obj.setTransformation(m);
}



S_Object* Sphere(double radius)
{
  return new S_Sphere(radius);
}



S_Object* Box(double x, double y, double z)
{
  return new S_Box(x, y, z);
}



S_Object* STPBV(const std::string& filename)
{
  STP_BV* s = new STP_BV;
  s->constructFromFile(filename);
  return s;
}

}

