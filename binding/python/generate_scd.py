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

from pybindgen import *
import sys


def import_eigen3_types(mod):
  mod.add_class('Vector3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Matrix3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('MatrixXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('VectorXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')



def import_sva_types(mod):
  mod.add_class('PTransform', foreign_cpp_namespace='sva', import_from_module='spacevecalg')



def build_SCD(scd):
  obj = scd.add_class('S_Object')
  pair = scd.add_class('CD_Pair')

  obj.add_function_as_method('transform', None, [param('SCD::S_Object&', 'obj'),
                                                 param('const sva::PTransform&', 'trans')],
                             custom_name='transform')


  scd.add_function('Sphere', retval('SCD::S_Object*', caller_owns_return=True),
                   [param('double', 'radius')])
  scd.add_function('Box', retval('SCD::S_Object*', caller_owns_return=True),
                   [param('double', 'x'),
                    param('double', 'y'),
                    param('double', 'z')])
  scd.add_function('STPBV', retval('SCD::S_Object*', caller_owns_return=True),
                   [param('const std::string&', 'filename')])
  scd.add_function('Polyhedron', retval('SCD::S_Object*', caller_owns_return=True),
                   [param('const std::string&', 'filename')])


  pair.add_constructor([param('SCD::S_Object*', 'obj1', transfer_ownership=False),
                        param('SCD::S_Object*', 'obj2', transfer_ownership=False)])

  pair.add_method('getDistance', retval('double'), [], custom_name='distance')

  pair.add_function_as_method('distance', retval('double'),
                              [param('SCD::CD_Pair&', 'pair'),
                               param('Eigen::Vector3d&', 'p1'),
                               param('Eigen::Vector3d&', 'p2')],
                              custom_name='distance')




if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  scd = Module('_scd', cpp_namespace='::SCD')
  scd.add_include('<SpaceVecAlg>')
  scd.add_include('<SCD/S_Object/S_Object.h>')
  scd.add_include('<SCD/S_Object/S_Sphere.h>')
  scd.add_include('<SCD/S_Object/S_Box.h>')
  scd.add_include('<SCD/S_Polyhedron/S_Polyhedron.h>')
  scd.add_include('<SCD/STP-BV/STP_BV.h>')
  scd.add_include('<SCD/CD/CD_Pair.h>')

  scd.add_include('"SCDAddon.h"')

  dom_ex = scd.add_exception('std::domain_error', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')
  out_ex = scd.add_exception('std::out_of_range', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')

  # import Eigen3, and sva
  import_eigen3_types(scd)
  import_sva_types(scd)

  # SCD
  build_SCD(scd)

  with open(sys.argv[1], 'w') as f:
    scd.generate(f)

