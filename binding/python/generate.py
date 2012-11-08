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

def import_SCD_types(mod):
  mod.add_class('S_Object', foreign_cpp_namespace='SCD', import_from_module='scd')
  mod.add_class('CD_Pair', foreign_cpp_namespace='SCD', import_from_module='scd')



def import_rbd_types(mod):
  mod.add_class('Body', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('Joint', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('MultiBody', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('MultiBodyConfig', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('Jacobian', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('CoMJacobianDummy', foreign_cpp_namespace='rbd', import_from_module='rbdyn')



def import_sva_types(mod):
  mod.add_class('MotionVec', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('ForceVec', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('RBInertia', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('ABInertia', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('PTransform', foreign_cpp_namespace='sva', import_from_module='spacevecalg')



def import_eigen3_types(mod):
  mod.add_class('Vector3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Vector6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Matrix3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Matrix6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('MatrixXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('VectorXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Quaterniond', foreign_cpp_namespace='Eigen', import_from_module='eigen3')



def build_tasks(posTask, oriTask, positionTask, comTask):
  def add_std_func(cls):
    cls.add_method('update', None,
                   [param('const rbd::MultiBody&', 'mb'),
                    param('const rbd::MultiBodyConfig&', 'mbc')])

    cls.add_method('updateDot', None,
                   [param('const rbd::MultiBody&', 'mb'),
                    param('const rbd::MultiBodyConfig&', 'mbc')])

    cls.add_method('eval', retval('Eigen::VectorXd'), [], is_const=True)
    cls.add_method('jac', retval('Eigen::MatrixXd'), [], is_const=True)
    cls.add_method('jacDot', retval('Eigen::MatrixXd'), [], is_const=True)


  # PositionTask
  posTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('int', 'bodyId'),
                           param('const Eigen::Vector3d&', 'pos'),
                           param('const Eigen::Vector3d&', 'bodyPoint',
                                 default_value='Eigen::Vector3d::Zero()')])

  posTask.add_method('position', None, [param('const Eigen::Vector3d&', 'pos')])
  posTask.add_method('position', retval('Eigen::Vector3d'), [], is_const=True)
  posTask.add_method('bodyPoint', None, [param('const Eigen::Vector3d&', 'point')])
  posTask.add_method('bodyPoint', retval('Eigen::Vector3d'), [], is_const=True)
  add_std_func(posTask)

  # OrientationTask
  oriTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('int', 'bodyId'),
                           param('const Eigen::Quaterniond&', 'ori')])
  oriTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('int', 'bodyId'),
                           param('const Eigen::Matrix3d&', 'ori')])

  oriTask.add_method('orientation', None, [param('const Eigen::Matrix3d&', 'ori')])
  oriTask.add_method('orientation', None, [param('const Eigen::Quaterniond&', 'ori')])
  oriTask.add_method('orientation', retval('Eigen::Matrix3d'), [], is_const=True)
  add_std_func(oriTask)

  # PostureTask
  postureTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                               param('std::vector<std::vector<double> >', 'q')])

  postureTask.add_method('posture', None,
                         [param('std::vector<std::vector<double> >', 'q')])
  postureTask.add_method('posture',
                         retval('std::vector<std::vector<double> >','q'), [],
                         is_const=True)
  add_std_func(postureTask)

  # CoMTask
  comTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('const Eigen::Vector3d&', 'com')])

  comTask.add_method('com', None, [param('const Eigen::Vector3d&', 'com')])
  comTask.add_method('com', retval('const Eigen::Vector3d&', 'com'), [],
                     is_const=True)
  add_std_func(comTask)



def build_qp(tasks):
  qp = tasks.add_cpp_namespace('qp')

  sol = qp.add_class('QPSolver')

  frictionCone = qp.add_struct('FrictionCone')
  contact = qp.add_struct('Contact')

  constr = qp.add_class('Constraint')
  eqConstr = qp.add_class('EqualityConstraint')
  ineqConstr = qp.add_class('InequalityConstraint')
  boundConstr = qp.add_class('BoundConstraint')

  task = qp.add_class('Task')
  hlTask = qp.add_class('HighLevelTask')

  spTask = qp.add_class('SetPointTask', parent=task)

  posTask = qp.add_class('PositionTask', parent=hlTask)
  oriTask = qp.add_class('OrientationTask', parent=hlTask)
  postureTask = qp.add_class('PostureTask', parent=task)
  comTask = qp.add_class('CoMTask', parent=hlTask)

  motionConstr = qp.add_class('MotionConstr', parent=[eqConstr, boundConstr,
                                                     constr])
  contactAccConstr = qp.add_class('ContactAccConstr', parent=[eqConstr, constr])

  selfCollisionConstr = qp.add_class('SelfCollisionConstr', parent=[ineqConstr, constr])
  seCollisionConstr = qp.add_class('StaticEnvCollisionConstr', parent=[ineqConstr, constr])

  jointLimitsConstr = qp.add_class('JointLimitsConstr', parent=[boundConstr, constr])
  torqueLimitsConstr = qp.add_class('TorqueLimitsConstr', parent=[boundConstr, constr])


  # build list type
  tasks.add_container('std::vector<tasks::qp::Contact>', 'tasks::qp::Contact', 'vector')
  tasks.add_container('std::vector<Eigen::Vector3d>', 'Eigen::Vector3d', 'vector')


  # QPSolver
  def add_std_solver_add_rm_nr(name, type):
    for t in type:
      ptr = '%s*' % t
      sol.add_method('add%s' % name, None,
                     [param(ptr, 'ptr', transfer_ownership=False)])
      sol.add_method('remove%s' % name, None,
                     [param(ptr, 'ptr', transfer_ownership=False)])

    sol.add_method('nr%ss' % name, retval('int'), [], is_const=True)

  sol.add_constructor([])
  sol.add_method('update', retval('bool'),
                 [param('const rbd::MultiBody&', 'mb'),
                  param('const rbd::MultiBodyConfig&', 'mbc')])

  sol.add_method('updateEqConstrSize', None, [])
  sol.add_method('updateInEqConstrSize', None, [])

  sol.add_method('nrVars', None,
                 [param('const rbd::MultiBody&', 'mb'),
                  param('std::vector<tasks::qp::Contact>&', 'cont')])
  sol.add_method('nrVars', retval('int'), [], is_const=True)

  constrName = ['MotionConstr', 'ContactAccConstr', 'SelfCollisionConstr', 'JointLimitsConstr',
                'StaticEnvCollisionConstr', 'TorqueLimitsConstr']
  eqConstrName = ['MotionConstr', 'ContactAccConstr']
  ineqConstrName = ['SelfCollisionConstr', 'StaticEnvCollisionConstr']
  boundConstrName = ['MotionConstr', 'JointLimitsConstr', 'TorqueLimitsConstr']
  taskName = ['SetPointTask', 'tasks::qp::PostureTask']

  add_std_solver_add_rm_nr('EqualityConstraint', eqConstrName)
  add_std_solver_add_rm_nr('InequalityConstraint', ineqConstrName)
  add_std_solver_add_rm_nr('BoundConstraint', boundConstrName)
  add_std_solver_add_rm_nr('Constraint', constrName)
  add_std_solver_add_rm_nr('Task', taskName)
  sol.add_method('resetTasks', None, [])


  # FrictionCone
  frictionCone.add_constructor([])
  frictionCone.add_constructor([param('Eigen::Matrix3d', 'frame'), param('int', 'nrGen'),
                                param('double', 'angle')])

  frictionCone.add_instance_attribute('generators', 'std::vector<Eigen::Vector3d>')

  # Contact
  contact.add_constructor([])
  contact.add_constructor([param('int', 'bodyId'),
                           param('const std::vector<Eigen::Vector3d>&', 'points'),
                           param('Eigen::Matrix3d', 'frame'),
                           param('int', 'nrGen'), param('double', 'angle')])

  contact.add_instance_attribute('bodyId', 'int')
  contact.add_instance_attribute('point', 'Eigen::Vector3d')
  contact.add_instance_attribute('cone', 'tasks::qp::FrictionCone')


  # Constraint
  constr.add_method('updateNrVars', None,
                    [param('const rbd::MultiBody&', 'mb'),
                     param('int', 'alphaD'),
                     param('int', 'lambda'),
                     param('int', 'torque'),
                     param('const std::vector<tasks::qp::Contact>', 'cont')])

  constr.add_method('update', None,
                    [param('const rbd::MultiBody&', 'mb'),
                     param('const rbd::MultiBodyConfig&', 'mbc')])

  # EqualityConstraint
  eqConstr.add_method('nrEqLine', retval('int'), [])
  # eqConstr.add_method('AEq', retval('Eigen::MatrixXd'), [],
  #                     is_virtual=True, is_pure_virtual=True)
  # eqConstr.add_method('BEq', retval('Eigen::MatrixXd'), [],
  #                     is_virtual=True, is_pure_virtual=True)

  # InequalityConstraint
  ineqConstr.add_method('nrInEqLine', retval('int'), [])
  # ineqConstr.add_method('AInEq', retval('Eigen::MatrixXd'), [],
  #                       is_virtual=True, is_pure_virtual=True)
  # ineqConstr.add_method('BInEq', retval('Eigen::MatrixXd'), [],
  #                       is_virtual=True, is_pure_virtual=True)

  # BoundConstraint
  boundConstr.add_method('beginVar', retval('int'), [])
  # boundConstr.add_method('Lower', retval('Eigen::MatrixXd'), [],
  #                        is_virtual=True, is_pure_virtual=True)
  # boundConstr.add_method('Upper', retval('Eigen::VectorXd'), [],
  #                        is_virtual=True, is_pure_virtual=True)

  # Task
  # task.add_constructor([param('double', 'weight')])

  task.add_method('weight', retval('double'), [], is_const=True)
  task.add_method('weight', None, [param('double', 'weight')])

  # task.add_method('update', None,
  #                 [param('const rbd::MultiBody&', 'mb'),
  #                  param('const rbd::MultiBodyConfig&', 'mbc')],
  #                 is_virtual=True, is_pure_virtual=True)

  # task.add_method('Q', retval('Eigen::MatrixXd'), [],
  #                 is_virtual=True, is_pure_virtual=True, is_const=True)
  # task.add_method('C', retval('Eigen::VectorXd'), [],
  #                 is_virtual=True, is_pure_virtual=True, is_const=True)

  # HighLevelTask
  hlTask.add_method('dim', retval('int'), [])

  hlTask.add_method('update', None,
                    [param('const rbd::MultiBody&', 'mb'),
                     param('const rbd::MultiBodyConfig&', 'mbc')])

  hlTask.add_method('jac', retval('Eigen::MatrixXd'), [])
  hlTask.add_method('jacDot', retval('Eigen::MatrixXd'), [])
  hlTask.add_method('eval', retval('Eigen::VectorXd'), [])

  # SetPointTask
  def spConstructor(hlTaskName):
    for t in hlTaskName:
      name = 'tasks::qp::%s *' % t
      spTask.add_constructor([param('const MultiBody&', 'mb'),
                              param(name, 'hlTask',
                                    transfer_ownership=False),
                              param('double', 'stiffness'),
                              param('double', 'weight')])

  spConstructor(['PositionTask', 'OrientationTask', 'CoMTask'])

  spTask.add_method('stiffness', retval('double'), [], is_const=True)
  spTask.add_method('stiffness', None, [param('double', 'weight')])

  spTask.add_method('update', None,
                    [param('const rbd::MultiBody&', 'mb'),
                     param('const rbd::MultiBodyConfig&', 'mbc')])

  spTask.add_method('Q', retval('Eigen::MatrixXd'), [], is_const=True)
  spTask.add_method('C', retval('Eigen::VectorXd'), [], is_const=True)

  # PositionTask
  posTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('int', 'bodyId'),
                           param('const Eigen::Vector3d&', 'pos'),
                           param('const Eigen::Vector3d&', 'bodyPoint',
                                 default_value='Eigen::Vector3d::Zero()')])

  posTask.add_method('position', None, [param('const Eigen::Vector3d&', 'pos')])
  posTask.add_method('position', retval('Eigen::Vector3d'), [], is_const=True)

  posTask.add_method('bodyPoint', None, [param('const Eigen::Vector3d&', 'point')])
  posTask.add_method('bodyPoint', retval('Eigen::Vector3d'), [], is_const=True)

  # OrientationTask
  oriTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('int', 'bodyId'),
                           param('const Eigen::Quaterniond&', 'ori')])
  oriTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('int', 'bodyId'),
                           param('const Eigen::Matrix3d&', 'ori')])

  oriTask.add_method('orientation', None, [param('const Eigen::Matrix3d&', 'ori')])
  oriTask.add_method('orientation', None, [param('const Eigen::Quaterniond&', 'ori')])
  oriTask.add_method('orientation', retval('Eigen::Matrix3d'), [], is_const=True)

  # PostureTask
  postureTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                               param('std::vector<std::vector<double> >', 'q'),
                               param('double', 'stiffness'),
                               param('double', 'weight')])


  postureTask.add_method('stiffness', retval('double'), [], is_const=True)
  postureTask.add_method('stiffness', None, [param('double', 'weight')])

  postureTask.add_method('posture', None,
                         [param('std::vector<std::vector<double> >', 'q')])
  postureTask.add_method('posture',
                         retval('std::vector<std::vector<double> >','q'), [],
                         is_const=True)

  # CoMTask
  comTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('const Eigen::Vector3d&', 'com')])

  comTask.add_method('com', None, [param('const Eigen::Vector3d&', 'com')])
  comTask.add_method('com', retval('const Eigen::Vector3d&', 'com'), [],
                     is_const=True)

  # MotionConstr
  motionConstr.add_constructor([param('const rbd::MultiBody', 'mb')])

  # ContactAccConstr
  contactAccConstr.add_constructor([param('const rbd::MultiBody', 'mb')])

  # SelfCollisionConstr
  selfCollisionConstr.add_constructor([param('const rbd::MultiBody', 'mb'),
                                       param('double', 'step')])
  selfCollisionConstr.add_method('addCollision', None,
                                 [param('const rbd::MultiBody&', 'mb'),
                                  param('int', 'body1Id'),
                                  param('SCD::S_Object*', 'body1', transfer_ownership=False),
                                  param('const sva::PTransform&', 'body1T'),
                                  param('int', 'body2Id'),
                                  param('SCD::S_Object*', 'body2', transfer_ownership=False),
                                  param('const sva::PTransform&', 'body2T'),
                                  param('double', 'di'),
                                  param('double', 'ds'),
                                  param('double', 'damping')])

  selfCollisionConstr.add_method('rmCollision', None, [param('int', 'body1Id'),
                                                       param('int', 'body2Id')])

  selfCollisionConstr.add_method('reset', None, []),

  # StaticEnvCollisionConstr
  seCollisionConstr.add_constructor([param('const rbd::MultiBody', 'mb'),
                                     param('double', 'step')])
  seCollisionConstr.add_method('addCollision', None,
                               [param('const rbd::MultiBody&', 'mb'),
                                param('int', 'bodyId'),
                                param('SCD::S_Object*', 'body', transfer_ownership=False),
                                param('const sva::PTransform&', 'bodyT'),
                                param('int', 'envId'),
                                param('SCD::S_Object*', 'env', transfer_ownership=False),
                                param('double', 'di'),
                                param('double', 'ds'),
                                param('double', 'damping')])

  seCollisionConstr.add_method('rmCollision', None, [param('int', 'bodyId'),
                                                     param('int', 'envId')])

  seCollisionConstr.add_method('reset', None, []),

  # JointLimitsConstr
  jointLimitsConstr.add_constructor([param('const rbd::MultiBody&', 'mb'),
                                    param('const std::vector<std::vector<double> >&', 'lbound'),
                                    param('const std::vector<std::vector<double> >&', 'ubound'),
                                    param('double', 'step')])

  # TorqueLimitsConstr
  torqueLimitsConstr.add_constructor([param('const rbd::MultiBody&', 'mb'),
                                      param('const std::vector<std::vector<double> >&', 'lbound'),
                                      param('const std::vector<std::vector<double> >&', 'ubound')])



if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  tasks = Module('_tasks', cpp_namespace='::tasks')
  tasks.add_include('<Tasks.h>')
  tasks.add_include('<QPSolver.h>')
  tasks.add_include('<QPTasks.h>')
  tasks.add_include('<QPConstr.h>')

  tasks.add_include('<MultiBodyConfig.h>')

  tasks.add_include('<SCD/S_Object/S_Object.h>')
  tasks.add_include('<SCD/S_Object/S_Sphere.h>')
  tasks.add_include('<SCD/CD/CD_Pair.h>')

  dom_ex = tasks.add_exception('std::domain_error', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')
  out_ex = tasks.add_exception('std::out_of_range', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')

  # import Eigen3, sva and rbd types
  import_eigen3_types(tasks)
  import_sva_types(tasks)
  import_rbd_types(tasks)
  import_SCD_types(tasks)

  posTask = tasks.add_class('PositionTask')
  oriTask = tasks.add_class('OrientationTask')
  postureTask = tasks.add_class('PostureTask')
  comTask = tasks.add_class('CoMTask')

  # build list type
  tasks.add_container('std::vector<double>', 'double', 'vector')
  tasks.add_container('std::vector<std::vector<double> >', 'std::vector<double>', 'vector')

  build_tasks(posTask, oriTask, postureTask, comTask)

  # qp
  build_qp(tasks)




  with open(sys.argv[1], 'w') as f:
    tasks.generate(f)

