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

def import_sch_types(mod):
  mod.add_class('S_Object', foreign_cpp_namespace='sch', import_from_module='sch')
  mod.add_class('CD_Pair', foreign_cpp_namespace='sch', import_from_module='sch')



def import_rbd_types(mod):
  mod.add_class('Body', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('Joint', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('MultiBody', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('MultiBodyConfig', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('Jacobian', foreign_cpp_namespace='rbd', import_from_module='rbdyn')
  mod.add_class('CoMJacobianDummy', foreign_cpp_namespace='rbd', import_from_module='rbdyn')



def import_sva_types(mod):
  mod.add_class('MotionVecd', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('ForceVecd', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('RBInertiad', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('ABInertiad', foreign_cpp_namespace='sva', import_from_module='spacevecalg')
  mod.add_class('PTransformd', foreign_cpp_namespace='sva', import_from_module='spacevecalg')



def import_eigen3_types(mod):
  mod.add_class('Vector2d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Vector3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Vector6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Matrix3d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('Matrix6d', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('MatrixXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')
  mod.add_class('VectorXd', foreign_cpp_namespace='Eigen', import_from_module='eigen3')

  mod.add_class('Quaterniond', foreign_cpp_namespace='Eigen', import_from_module='eigen3')


def build_boost_timer(tasks):
  cpuTime = tasks.add_struct('cpu_times', foreign_cpp_namespace='boost::timer')
  cpuTime.add_instance_attribute('wall', 'long long int')
  cpuTime.add_instance_attribute('user', 'long long int')
  cpuTime.add_instance_attribute('system', 'long long int')


def build_tasks(posTask, oriTask, surfOriTask, gazeTask, positionTask, comTask,
                multiCoMTask, momTask, linVelTask, oriTrackTask,
                multiRobotTransformTask):
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

  # SurfaceOrientationTask
  surfOriTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                               param('int', 'bodyId'),
                               param('const Eigen::Quaterniond&', 'ori'),
                               param('const sva::PTransformd&', 'X_b_s')])
  surfOriTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                               param('int', 'bodyId'),
                               param('const Eigen::Matrix3d&', 'ori'),
                               param('const sva::PTransformd&', 'X_b_s')])

  surfOriTask.add_method('orientation', None, [param('const Eigen::Matrix3d&', 'ori')])
  surfOriTask.add_method('orientation', None, [param('const Eigen::Quaterniond&', 'ori')])
  surfOriTask.add_method('orientation', retval('Eigen::Matrix3d'), [], is_const=True)
  add_std_func(surfOriTask)

  # GazeTask
  gazeTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                            param('int', 'bodyId'),
                            param('const Eigen::Vector2d&', 'point2d'),
                            param('double', 'depthEstimate'),
                            param('const sva::PTransformd&', 'X_b_gaze'),
                            param('const Eigen::Vector2d&', 'point2d_ref',
                                  default_value='Eigen::Vector2d::Zero()')])
  gazeTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                            param('int', 'bodyId'),
                            param('const Eigen::Vector3d&', 'point3d'),
                            param('const sva::PTransformd&', 'X_b_gaze'),
                            param('const Eigen::Vector2d&', 'point2d_ref',
                                  default_value='Eigen::Vector2d::Zero()')])

  gazeTask.add_method('error', None, [param('const Eigen::Vector2d&', 'point2d'),
                                      param('const Eigen::Vector2d&', 'point2d_ref',
                                            default_value='Eigen::Vector2d::Zero()')])
  gazeTask.add_method('error', None, [param('const Eigen::Vector3d&', 'point3d'),
                                      param('const Eigen::Vector2d&', 'point2d_ref',
                                            default_value='Eigen::Vector2d::Zero()')])
  gazeTask.add_method('update', None, [param('const rbd::MultiBody&', 'mb'),
                                       param('const rbd::MultiBodyConfig&', 'mbc'),
                                       param('const std::vector<sva::MotionVecd>&', 'normalAccB')])
  gazeTask.add_method('eval', retval('Eigen::VectorXd'), [], is_const=True)
  gazeTask.add_method('speed', retval('Eigen::VectorXd'), [], is_const=True)
  gazeTask.add_method('normalAcc', retval('Eigen::VectorXd'), [], is_const=True)
  gazeTask.add_method('jac', retval('Eigen::MatrixXd'), [], is_const=True)
  gazeTask.add_method('jacDot', retval('Eigen::MatrixXd'), [], is_const=True)


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
  comTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('const Eigen::Vector3d&', 'com'),
                           param('std::vector<double>', 'weight')],
                           throw=[dom_ex])

  comTask.add_method('com', None, [param('const Eigen::Vector3d&', 'com')])
  comTask.add_method('com', retval('const Eigen::Vector3d&', 'com'), [],
                     is_const=True)
  comTask.add_method('updateInertialParameters', None,
                     [param('const rbd::MultiBody&', 'mb')])
  add_std_func(comTask)

  # MultiCoMTask
  multiCoMTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                param('std::vector<int>', 'robotIndexes'),
                                param('const Eigen::Vector3d&', 'com')])

  multiCoMTask.add_method('robotIndexes', retval('std::vector<int>'), [],
                          is_const=True)
  multiCoMTask.add_method('com', None, [param('const Eigen::Vector3d&', 'com')])
  multiCoMTask.add_method('com', retval('const Eigen::Vector3d&'), [],
                          is_const=True)
  multiCoMTask.add_method('updateInertialParameters', None,
                     [param('const std::vector<rbd::MultiBody>&', 'mbs')])
  multiCoMTask.add_method('update', None,
                          [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                           param('const std::vector<rbd::MultiBodyConfig>&', 'mbcs')])

  multiCoMTask.add_method('eval', retval('const Eigen::VectorXd&'), [],
                          is_const=True)
  multiCoMTask.add_method('speed', retval('const Eigen::VectorXd&'), [],
                          is_const=True)
  multiCoMTask.add_method('normalAcc', retval('const Eigen::VectorXd&'), [],
                          is_const=True)
  multiCoMTask.add_method('jac', retval('const Eigen::MatrixXd&'),
                          [param('int', 'index')], is_const=True)

  # MultiRobotTransformTask
  multiRobotTransformTask.add_constructor(
    [param('const std::vector<rbd::MultiBody>&', 'mbs'),
     param('int', 'r1Index'), param('int', 'r2Index'),
     param('int', 'r1BodyIndex'), param('int', 'r2BodyIndex'),
     param('const sva::PTransformd&', 'X_r1b_r1s'),
     param('const sva::PTransformd&', 'X_r2b_r2s')])

  multiRobotTransformTask.add_method('r1Index', retval('int'), [], is_const=True)
  multiRobotTransformTask.add_method('r2Index', retval('int'), [], is_const=True)
  multiRobotTransformTask.add_method('X_r1b_r1s', None,
    [param('const sva::PTransformd&', 'X_r1b_r1s')])
  multiRobotTransformTask.add_method('X_r1b_r1s', retval('sva::PTransformd'),
                                     [], is_const=True)
  multiRobotTransformTask.add_method('X_r2b_r2s', None,
                                     [param('const sva::PTransformd&', 'X_r2b_r2s')])
  multiRobotTransformTask.add_method('X_r2b_r2s', retval('sva::PTransformd'),
                                     [], is_const=True)
  multiRobotTransformTask.add_method('update', None,
    [param('const std::vector<rbd::MultiBody>&', 'mbs'),
      param('const std::vector<rbd::MultiBodyConfig>&', 'mbcs'),
      param('const std::vector<std::vector<sva::MotionVecd> >&', 'normalAccB')])

  multiRobotTransformTask.add_method('eval', retval('const Eigen::VectorXd&'), [],
                                     is_const=True)
  multiRobotTransformTask.add_method('speed', retval('const Eigen::VectorXd&'), [],
                                     is_const=True)
  multiRobotTransformTask.add_method('normalAcc', retval('const Eigen::VectorXd&'), [],
                                     is_const=True)
  multiRobotTransformTask.add_method('jac', retval('const Eigen::MatrixXd&'),
                                     [param('int', 'index')], is_const=True)

  # MomentumTask
  momTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('const sva::ForceVecd&', 'mom')])

  momTask.add_method('momentum', None, [param('const sva::ForceVecd&', 'mom')])
  momTask.add_method('momentum', retval('const sva::ForceVecd&', 'momentum'), [],
                     is_const=True)
  add_std_func(momTask)

  # LinVelocityTask
  linVelTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                           param('int', 'bodyId'),
                           param('const Eigen::Vector3d&', 'pos'),
                           param('const Eigen::Vector3d&', 'bodyPoint',
                                 default_value='Eigen::Vector3d::Zero()')])

  linVelTask.add_method('velocity', None, [param('const Eigen::Vector3d&', 'pos')])
  linVelTask.add_method('velocity', retval('Eigen::Vector3d'), [], is_const=True)
  linVelTask.add_method('bodyPoint', None, [param('const Eigen::Vector3d&', 'point')])
  linVelTask.add_method('bodyPoint', retval('Eigen::Vector3d'), [], is_const=True)
  add_std_func(linVelTask)

  # OrientationTrackingTask
  oriTrackTask.add_constructor([param('const rbd::MultiBody&', 'mb'),
                                param('int', 'bodyId'),
                                param('const Eigen::Vector3d&', 'bodyPoint'),
                                param('const Eigen::Vector3d&', 'bodyAxis'),
                                param('const std::vector<int>&', 'trackingJointId'),
                                param('const Eigen::Vector3d&', 'trackedPoint')])

  oriTrackTask.add_method('trackedPoint', None, [param('const Eigen::Vector3d&', 'bPoint')])
  oriTrackTask.add_method('trackedPoint', retval('Eigen::Vector3d'), [], is_const=True)
  oriTrackTask.add_method('bodyPoint', None, [param('const Eigen::Vector3d&', 'tPoint')])
  oriTrackTask.add_method('bodyPoint', retval('Eigen::Vector3d'), [], is_const=True)
  oriTrackTask.add_method('bodyAxis', None, [param('const Eigen::Vector3d&', 'axis')])
  oriTrackTask.add_method('bodyAxis', retval('Eigen::Vector3d'), [], is_const=True)
  add_std_func(oriTrackTask)



def build_qp(tasks):
  qp = tasks.add_cpp_namespace('qp')

  sol = qp.add_class('QPSolver')
  solData = qp.add_class('SolverData')

  frictionCone = qp.add_struct('FrictionCone')
  contactId = qp.add_struct('ContactId')
  unilateralContact = qp.add_struct('UnilateralContact')
  bilateralContact = qp.add_struct('BilateralContact')
  jointStiffness = qp.add_struct('JointStiffness')
  springJoint = qp.add_struct('SpringJoint')
  qBound = tasks.add_struct('QBound')
  alphaBound = tasks.add_struct('AlphaBound')
  torqueBound = tasks.add_struct('TorqueBound')
  polyTorqueBound = tasks.add_struct('PolyTorqueBound')

  constr = qp.add_class('Constraint')
  eqConstr = qp.add_class('Equality')
  ineqConstr = qp.add_class('Inequality')
  genineqConstr = qp.add_class('GenInequality')
  boundConstr = qp.add_class('Bound')

  task = qp.add_class('Task')
  hlTask = qp.add_class('HighLevelTask')

  spTask = qp.add_class('SetPointTask', parent=task)
  pidTask = qp.add_class('PIDTask', parent=task)
  toTask = qp.add_class('TargetObjectiveTask', parent=task)

  posTask = qp.add_class('PositionTask', parent=hlTask)
  oriTask = qp.add_class('OrientationTask', parent=hlTask)
  surfOriTask = qp.add_class('SurfaceOrientationTask', parent=hlTask)
  gazeTask = qp.add_class('GazeTask', parent=hlTask)
  postureTask = qp.add_class('PostureTask', parent=task)
  comTask = qp.add_class('CoMTask', parent=hlTask)
  multiCoMTask = qp.add_class('MultiCoMTask', parent=task)
  multiRobotTransformTask = qp.add_class('MultiRobotTransformTask', parent=task)
  momTask = qp.add_class('MomentumTask', parent=hlTask)
  contactTask = qp.add_class('ContactTask', parent=task)
  gripperTorqueTask = qp.add_class('GripperTorqueTask', parent=task)
  linVelTask = qp.add_class('LinVelocityTask', parent=hlTask)
  oriTrackTask = qp.add_class('OrientationTrackingTask', parent=hlTask)
  jointsSelector = qp.add_class('JointsSelector', parent=hlTask)

  motionConstr = qp.add_class('MotionConstr', parent=[genineqConstr, constr])
  motionPolyConstr = qp.add_class('MotionPolyConstr', parent=[genineqConstr,
                                                              constr])
  motionSpringConstr = qp.add_class('MotionSpringConstr', parent=[genineqConstr,
                                                                  constr])
  positiveLambdaConstr = qp.add_class('PositiveLambda', parent=[boundConstr,
                                                                constr])
  contactConstrCommon = qp.add_class('ContactConstrCommon')
  contactAccConstr = qp.add_class('ContactAccConstr', parent=[eqConstr, constr, contactConstrCommon])
  contactSpeedConstr = qp.add_class('ContactSpeedConstr', parent=[eqConstr, contactConstrCommon])
  contactPosConstr = qp.add_class('ContactPosConstr', parent=[eqConstr, contactConstrCommon])

  collisionConstr = qp.add_class('CollisionConstr', parent=[ineqConstr, constr])
  comIncPlaneConstr = qp.add_class('CoMIncPlaneConstr', parent=[ineqConstr, constr])

  jointLimitsConstr = qp.add_class('JointLimitsConstr', parent=[boundConstr, constr])
  damperJointLimitsConstr = qp.add_class('DamperJointLimitsConstr', parent=[boundConstr, constr])

  gripperTorqueConstr = qp.add_class('GripperTorqueConstr', parent=[ineqConstr, constr])
  boundedSpeedConstr = qp.add_class('BoundedSpeedConstr', parent=[genineqConstr, constr])


  constrName = ['MotionConstr', 'MotionPolyConstr', 'ContactAccConstr', 'ContactSpeedConstr',
                'CollisionConstr', 'JointLimitsConstr', 'DamperJointLimitsConstr',
                'MotionSpringConstr', 'GripperTorqueConstr', 'BoundedSpeedConstr',
                'CoMIncPlaneConstr', 'PositiveLambda', 'ContactPosConstr']
  eqConstrName = ['ContactAccConstr', 'ContactSpeedConstr', 'ContactPosConstr']
  ineqConstrName = ['CollisionConstr', 'GripperTorqueConstr', 'CoMIncPlaneConstr']
  genineqConstrName = ['MotionConstr', 'MotionPolyConstr', 'MotionSpringConstr']
  boundConstrName = ['PositiveLambda', 'JointLimitsConstr', 'DamperJointLimitsConstr']
  taskName = ['SetPointTask', 'PIDTask', 'TargetObjectiveTask',
              'tasks::qp::PostureTask', 'tasks::qp::ContactTask',
              'tasks::qp::GripperTorqueTask', 'tasks::qp::MultiCoMTask',
              'tasks::qp::MultiRobotTransformTask']
  hlTaskName = ['PositionTask', 'OrientationTask', 'SurfaceOrientationTask', 'GazeTask', 'CoMTask', 'LinVelocityTask',
                'OrientationTrackingTask', 'MomentumTask', 'JointsSelector']
  constrList = [motionConstr, motionPolyConstr, contactAccConstr, contactSpeedConstr,
                collisionConstr, jointLimitsConstr, damperJointLimitsConstr,
                motionSpringConstr, gripperTorqueConstr, boundedSpeedConstr,
                comIncPlaneConstr, positiveLambdaConstr, contactPosConstr]

  # build list type
  tasks.add_container('std::vector<tasks::qp::FrictionCone>',
                      'tasks::qp::FrictionCone', 'vector')
  tasks.add_container('std::vector<tasks::qp::UnilateralContact>',
                      'tasks::qp::UnilateralContact', 'vector')
  tasks.add_container('std::vector<tasks::qp::BilateralContact>',
                      'tasks::qp::BilateralContact', 'vector')
  tasks.add_container('std::vector<tasks::qp::JointStiffness>',
                      'tasks::qp::JointStiffness', 'vector')
  tasks.add_container('std::vector<tasks::qp::SpringJoint>',
                      'tasks::qp::SpringJoint', 'vector')
  tasks.add_container('std::vector<Eigen::Vector3d>', 'Eigen::Vector3d', 'vector')
  tasks.add_container('std::vector<Eigen::Matrix3d>', 'Eigen::Matrix3d', 'vector')
  tasks.add_container('std::vector<Eigen::VectorXd>', 'Eigen::VectorXd', 'vector')
  tasks.add_container('std::vector<std::vector<Eigen::VectorXd> >',
                      'std::vector<Eigen::VectorXd>', 'vector')


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
  sol.add_method('solve', retval('bool'),
                 [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                  param('std::vector<rbd::MultiBodyConfig>&', 'mbcs')])
  sol.add_method('solveNoMbcUpdate', retval('bool'),
                 [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                  param('const std::vector<rbd::MultiBodyConfig>&', 'mbcs')])
  sol.add_method('updateMbc', None,
                 [param('rbd::MultiBodyConfig&', 'mbc'),
                  param('int', 'robotIndex')], is_const=True)

  sol.add_method('updateConstrSize', None, [])

  sol.add_method('nrVars', None,
                 [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                  param('std::vector<tasks::qp::UnilateralContact>&', 'uni'),
                  param('std::vector<tasks::qp::BilateralContact>&', 'bi')])
  sol.add_method('nrVars', retval('int'), [], is_const=True)

  sol.add_method('updateTasksNrVars', None, [param('const std::vector<rbd::MultiBody>&', 'mb')],
                 is_const=True)
  sol.add_method('updateConstrsNrVars', None, [param('const std::vector<rbd::MultiBody>&', 'mb')],
                 is_const=True)
  sol.add_method('updateNrVars', None, [param('const std::vector<rbd::MultiBody>&', 'mb')],
                 is_const=True)

  add_std_solver_add_rm_nr('EqualityConstraint', eqConstrName)
  add_std_solver_add_rm_nr('InequalityConstraint', ineqConstrName)
  add_std_solver_add_rm_nr('GenInequalityConstraint', genineqConstrName)
  add_std_solver_add_rm_nr('BoundConstraint', boundConstrName)
  add_std_solver_add_rm_nr('Constraint', constrName)
  add_std_solver_add_rm_nr('Task', taskName)
  sol.add_method('resetTasks', None, [])

  sol.add_method('solver', None, [param('const std::string&', 'name')])

  sol.add_method('result', retval('const Eigen::VectorXd&'), [], is_const=True)
  sol.add_method('alphaDVec', retval('Eigen::VectorXd'), [], is_const=True)
  sol.add_method('alphaDVec', retval('Eigen::VectorXd'),
                 [param('int', 'robotIndex')], is_const=True)
  sol.add_method('lambdaVec', retval('Eigen::VectorXd'), [], is_const=True)
  sol.add_method('lambdaVec', retval('Eigen::VectorXd'),
                 [param('int', 'contactIndex')], is_const=True)
  sol.add_method('lambdaVec', retval('Eigen::VectorXd'), [], is_const=True)

  sol.add_method('contactLambdaPosition', retval('int'),
                 [param('const tasks::qp::ContactId&', 'contactId')], is_const=True)
  sol.add_method('data', retval('tasks::qp::SolverData'), [], is_const=True)

  sol.add_method('solveTime', retval('boost::timer::cpu_times'),
                 [], is_const=True)
  sol.add_method('solveAndBuildTime', retval('boost::timer::cpu_times'),
                 [], is_const=True)

  # SolverData
  solData.add_method('nrVars', retval('int'), [], is_const=True)
  solData.add_method('totalAlphaD', retval('int'), [], is_const=True)
  solData.add_method('totalLambda', retval('int'), [], is_const=True)
  solData.add_method('alphaD', retval('int'), [param('int','robotIndex')],
                     is_const=True)
  solData.add_method('lambda', retval('int'), [param('int','contactIndex')],
                     is_const=True)
  solData.add_method('alphaDBegin', retval('int'), [], is_const=True)
  solData.add_method('alphaDBegin', retval('int'), [param('int','robotIndex')],
                     is_const=True)
  solData.add_method('lambdaBegin', retval('int'), [], is_const=True)
  solData.add_method('lambdaBegin', retval('int'), [param('int','contactIndex')],
                     is_const=True)
  solData.add_method('nrUniLambda', retval('int'), [], is_const=True)
  solData.add_method('nrBiLambda', retval('int'), [], is_const=True)
  solData.add_method('unilateralBegin', retval('int'), [], is_const=True)
  solData.add_method('bilateralBegin', retval('int'), [], is_const=True)
  solData.add_method('nrContacts', retval('int'), [], is_const=True)
  solData.add_method('unilateralContacts',
                     retval('std::vector<tasks::qp::UnilateralContact>'),
                     [], is_const=True)
  solData.add_method('bilateralContacts',
                     retval('std::vector<tasks::qp::BilateralContact>'),
                     [], is_const=True)
  solData.add_method('allContacts',
                     retval('std::vector<tasks::qp::BilateralContact>'),
                     [], is_const=True)
  solData.add_method('computeNormalAccB', None,
                     [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                      param('const std::vector<rbd::MultiBodyConfig>&', 'mbcs')])
  solData.add_method('normalAccB',
                     retval('std::vector<sva::MotionVecd>'),
                     [param('int', 'robotIndex')], is_const=True)

  # FrictionCone
  frictionCone.add_constructor([])
  frictionCone.add_constructor([param('Eigen::Matrix3d', 'frame'), param('int', 'nrGen'),
                                param('double', 'mu'),
                                param('double', 'direction', default_value='1.')])

  frictionCone.add_instance_attribute('generators', 'std::vector<Eigen::Vector3d>')

  # ContactId
  contactId.add_constructor([])
  contactId.add_constructor([param('int', 'r1Index'), param('int', 'r2Index'),
                             param('int', 'r1BodyId'), param('int', 'r2BodyId'),
                             param('int', 'nSurf', default_value='-1')])
  contactId.add_instance_attribute('r1Index', 'int')
  contactId.add_instance_attribute('r2Index', 'int')
  contactId.add_instance_attribute('r1BodyId', 'int')
  contactId.add_instance_attribute('r2BodyId', 'int')
  contactId.add_instance_attribute('ambiguityId', 'int')
  contactId.add_binary_comparison_operator('==')
  contactId.add_binary_comparison_operator('!=')
  contactId.add_binary_comparison_operator('<')

  # UnilateralContact
  unilateralContact.add_constructor([])
  unilateralContact.add_constructor([
    param('int', 'r1Index'), param('int', 'r2Index'),
    param('int', 'r1BodyId'), param('int', 'r2BodyId'),
    param('const std::vector<Eigen::Vector3d>&', 'r1Points'),
    param('const Eigen::Matrix3d&', 'r1Frame'),
    param('const sva::PTransformd&', 'X_b1_b2'),
    param('int', 'nrGen'), param('double', 'mu'),
    param('const sva::PTransformd&','X_b1_cf',
           default_value='sva::PTransformd::Identity()')])
  unilateralContact.add_constructor([
    param('int', 'r1Index'), param('int', 'r2Index'),
    param('int', 'r1BodyId'), param('int', 'r2BodyId'),
    param('int', 'ambiguityId'),
    param('const std::vector<Eigen::Vector3d>&', 'r1Points'),
    param('const Eigen::Matrix3d&', 'r1Frame'),
    param('const sva::PTransformd&', 'X_b1_b2'),
    param('int', 'nrGen'), param('double', 'mu'),
    param('const sva::PTransformd&','X_b1_cf',
           default_value='sva::PTransformd::Identity()')])
  unilateralContact.add_constructor([param('const tasks::qp::ContactId&', 'contactId'),
    param('const std::vector<Eigen::Vector3d>&', 'r1Points'),
    param('const Eigen::Matrix3d&', 'r1Frame'),
    param('const sva::PTransformd&', 'X_b1_b2'),
    param('int', 'nrGen'), param('double', 'mu'),
    param('const sva::PTransformd&', 'X_b1_cf',
           default_value='sva::PTransformd::Identity()')])

  unilateralContact.add_instance_attribute('contactId', 'tasks::qp::ContactId')
  unilateralContact.add_instance_attribute('r1Points', 'std::vector<Eigen::Vector3d>')
  unilateralContact.add_instance_attribute('r2Points', 'std::vector<Eigen::Vector3d>')
  unilateralContact.add_instance_attribute('r1Cone', 'tasks::qp::FrictionCone')
  unilateralContact.add_instance_attribute('r2Cone', 'tasks::qp::FrictionCone')
  unilateralContact.add_instance_attribute('X_b1_b2', 'sva::PTransformd')
  unilateralContact.add_instance_attribute('X_b1_cf', 'sva::PTransformd')
  unilateralContact.add_method('sForce', retval('Eigen::Vector3d'),
                               [param('const Eigen::VectorXd&', 'lambda'),
                                param('int', 'point'),
                                param('const tasks::qp::FrictionCone&', 'c')],
                               throw=[dom_ex], custom_name='force', is_const=True)
  unilateralContact.add_method('sForce', retval('Eigen::Vector3d'),
                               [param('const Eigen::VectorXd&', 'lambda'),
                                param('const tasks::qp::FrictionCone&', 'c')],
                               throw=[dom_ex], custom_name='force', is_const=True)
  unilateralContact.add_method('sNrLambda', retval('int'), [param('int', 'point')],
                               is_const=True, throw=[dom_ex], custom_name='nrLambda')
  unilateralContact.add_method('nrLambda', retval('int'), [],
                               is_const=True, throw=[dom_ex])

  # BilateralContact
  bilateralContact.add_constructor([])
  bilateralContact.add_constructor([
    param('int', 'r1Index'), param('int', 'r2Index'),
    param('int', 'r1BodyId'), param('int', 'r2BodyId'),
    param('const std::vector<Eigen::Vector3d>&', 'r1Points'),
    param('const std::vector<Eigen::Matrix3d>&', 'r1Frames'),
    param('const sva::PTransformd&', 'X_b1_b2'),
    param('int', 'nrGen'), param('double', 'mu'),
    param('const sva::PTransformd&','X_b1_cf',
           default_value='sva::PTransformd::Identity()')])
  bilateralContact.add_constructor([
    param('int', 'r1Index'), param('int', 'r2Index'),
    param('int', 'r1BodyId'), param('int', 'r2BodyId'),
    param('int', 'ambiguityId'),
    param('const std::vector<Eigen::Vector3d>&', 'r1Points'),
    param('const std::vector<Eigen::Matrix3d>&', 'r1Frames'),
    param('const sva::PTransformd&', 'X_b1_b2'),
    param('int', 'nrGen'), param('double', 'mu'),
    param('const sva::PTransformd&','X_b1_cf',
           default_value='sva::PTransformd::Identity()')])
  bilateralContact.add_constructor([
    param('const tasks::qp::ContactId&', 'contactId'),
    param('const std::vector<Eigen::Vector3d>&', 'r1Points'),
    param('const std::vector<Eigen::Matrix3d>&', 'r1Frames'),
    param('const sva::PTransformd&', 'X_b1_b2'),
    param('int', 'nrGen'), param('double', 'mu'),
    param('const sva::PTransformd&','X_b1_cf',
          default_value='sva::PTransformd::Identity()')])
  bilateralContact.add_constructor([
    param('const tasks::qp::UnilateralContact&', 'uc')])

  bilateralContact.add_instance_attribute('contactId', 'tasks::qp::ContactId')
  bilateralContact.add_instance_attribute('r1Points', 'std::vector<Eigen::Vector3d>')
  bilateralContact.add_instance_attribute('r2Points', 'std::vector<Eigen::Vector3d>')
  bilateralContact.add_instance_attribute('r1Cones', 'std::vector<tasks::qp::FrictionCone>')
  bilateralContact.add_instance_attribute('r2Cones', 'std::vector<tasks::qp::FrictionCone>')
  bilateralContact.add_instance_attribute('X_b1_b2', 'sva::PTransformd')
  bilateralContact.add_instance_attribute('X_b1_cf', 'sva::PTransformd')
  bilateralContact.add_method('sForce', retval('Eigen::Vector3d'),
                              [param('const Eigen::VectorXd&', 'lambda'),
                               param('int', 'point'),
                               param('const std::vector<tasks::qp::FrictionCone>&', 'c')],
                              throw=[dom_ex], custom_name='force', is_const=True)
  bilateralContact.add_method('sForce', retval('Eigen::Vector3d'),
                              [param('const Eigen::VectorXd&', 'lambda'),
                               param('const std::vector<tasks::qp::FrictionCone>&', 'c')],
                              throw=[dom_ex], custom_name='force', is_const=True)
  bilateralContact.add_method('sNrLambda', retval('int'), [param('int', 'point')],
                              is_const=True, throw=[dom_ex], custom_name='nrLambda')
  bilateralContact.add_method('nrLambda', retval('int'), [],
                              is_const=True, throw=[dom_ex])

  # JointStiffness
  jointStiffness.add_constructor([])
  jointStiffness.add_constructor([param('int', 'jointId'),
                                  param('double', 'stiffness')])
  jointStiffness.add_instance_attribute('jointId', 'int')
  jointStiffness.add_instance_attribute('stiffness', 'std::vector<Eigen::Vector3d>')

  # SpringJoint
  springJoint.add_constructor([])
  springJoint.add_constructor([param('int', 'jointId'),
                               param('double', 'K'), param('double', 'C'),
                               param('double', 'O')])
  springJoint.add_instance_attribute('jointId', 'int')
  springJoint.add_instance_attribute('K', 'double')
  springJoint.add_instance_attribute('C', 'double')
  springJoint.add_instance_attribute('O', 'double')

  # QBound
  qBound.add_constructor([])
  qBound.add_constructor([param('std::vector<std::vector<double> >', 'lQB'),
                          param('std::vector<std::vector<double> >', 'uQB')])
  qBound.add_instance_attribute('lQBound', 'std::vector<std::vector<double> >')
  qBound.add_instance_attribute('uQBound', 'std::vector<std::vector<double> >')

  # AlphaBound
  alphaBound.add_constructor([])
  alphaBound.add_constructor([param('std::vector<std::vector<double> >', 'lAB'),
                              param('std::vector<std::vector<double> >', 'uAB')])
  alphaBound.add_instance_attribute('lAlphaBound', 'std::vector<std::vector<double> >')
  alphaBound.add_instance_attribute('uAlphaBound', 'std::vector<std::vector<double> >')

  # TorqueBound
  torqueBound.add_constructor([])
  torqueBound.add_constructor([param('std::vector<std::vector<double> >', 'lTB'),
                               param('std::vector<std::vector<double> >', 'uTB')])
  torqueBound.add_instance_attribute('lTorqueBound', 'std::vector<std::vector<double> >')
  torqueBound.add_instance_attribute('uTorqueBound', 'std::vector<std::vector<double> >')

  # PolyTorqueBound
  polyTorqueBound.add_constructor([])
  polyTorqueBound.add_constructor([param('std::vector<std::vector<Eigen::VectorXd> >', 'lPTB'),
                                   param('std::vector<std::vector<Eigen::VectorXd> >', 'uPTB')])
  polyTorqueBound.add_instance_attribute('lPolyTorqueBound',
                                         'std::vector<std::vector<Eigen::VectorXd> >')
  polyTorqueBound.add_instance_attribute('uPolyTorqueBound',
                                         'std::vector<std::vector<Eigen::VectorXd> >')

  # Constraint
  constr.add_method('updateNrVars', None,
                    [param('const std::vector<rbd::MultiBody>&', 'mb'),
                     param('tasks::qp::SolverData', 'data')])

  constr.add_method('update', None,
                    [param('const std::vector<rbd::MultiBody>&', 'mb'),
                     param('const std::vector<rbd::MultiBodyConfig>&', 'mbc'),
                     param('const tasks::qp::SolverData&', 'data')])

  # InequalityConstraint
  eqConstr.add_method('maxEq', retval('int'), [])
  eqConstr.add_method('nrEq', retval('int'), [])
  eqConstr.add_method('AEq', retval('Eigen::MatrixXd'), [])
  eqConstr.add_method('bEq', retval('Eigen::VectorXd'), [])

  # InequalityConstraint
  ineqConstr.add_method('maxInEq', retval('int'), [])
  ineqConstr.add_method('nrInEq', retval('int'), [])
  ineqConstr.add_method('AInEq', retval('Eigen::MatrixXd'), [])
  ineqConstr.add_method('bInEq', retval('Eigen::VectorXd'), [])

  # GenInequalityConstraint
  genineqConstr.add_method('maxGenInEq', retval('int'), [])
  genineqConstr.add_method('nrGenInEq', retval('int'), [])
  genineqConstr.add_method('AGenInEq', retval('Eigen::MatrixXd'), [])
  genineqConstr.add_method('LowerGenInEq', retval('Eigen::VectorXd'), [])
  genineqConstr.add_method('UpperGenInEq', retval('Eigen::VectorXd'), [])

  # BoundConstraint
  boundConstr.add_method('beginVar', retval('int'), [])
  boundConstr.add_method('Lower', retval('Eigen::MatrixXd'), [])
  boundConstr.add_method('Upper', retval('Eigen::VectorXd'), [])

  # Task
  # task.add_constructor([param('double', 'weight')])

  task.add_method('weight', retval('double'), [], is_const=True)
  task.add_method('weight', None, [param('double', 'weight')])

  # task.add_method('update', None,
  #                 [param('const std::vector<rbd::MultiBody>&', 'mb'),
  #                  param('const std::vector<rbd::MultiBodyConfig>&', 'mbc')],
  #                 is_virtual=True, is_pure_virtual=True)

  # task.add_method('Q', retval('Eigen::MatrixXd'), [],
  #                 is_virtual=True, is_pure_virtual=True, is_const=True)
  # task.add_method('C', retval('Eigen::VectorXd'), [],
  #                 is_virtual=True, is_pure_virtual=True, is_const=True)

  # HighLevelTask
  hlTask.add_method('dim', retval('int'), [])

  hlTask.add_method('update', None,
                    [param('const std::vector<rbd::MultiBody>&', 'mb'),
                     param('const std::vector<rbd::MultiBodyConfig>&', 'mbc'),
                     param('const tasks::qp::SolverData&', 'data')])

  hlTask.add_method('jac', retval('Eigen::MatrixXd'), [])
  hlTask.add_method('eval', retval('Eigen::VectorXd'), [])
  hlTask.add_method('speed', retval('Eigen::VectorXd'), [])
  hlTask.add_method('normalAcc', retval('Eigen::VectorXd'), [])

  # SetPointTask
  def spConstructor(hlTaskName):
    for t in hlTaskName:
      name = 'tasks::qp::%s *' % t
      spTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                              param('int', 'robotIndex'),
                              param(name, 'hlTask',
                                    transfer_ownership=False),
                              param('double', 'stiffness'),
                              param('double', 'weight')])

      spTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                              param('int', 'robotIndex'),
                              param(name, 'hlTask',
                                    transfer_ownership=False),
                              param('double', 'stiffness'),
                              param('const Eigen::VectorXd&', 'dimWeight'),
                              param('double', 'weight')])

  spConstructor(hlTaskName)

  spTask.add_method('stiffness', retval('double'), [], is_const=True)
  spTask.add_method('stiffness', None, [param('double', 'weight')])

  spTask.add_method('dimWeight', retval('const Eigen::VectorXd&'), [], is_const=True)
  spTask.add_method('dimWeight', None, [param('const Eigen::VectorXd&', 'dim')])

  spTask.add_method('update', None,
                    [param('const std::vector<rbd::MultiBody>&', 'mb'),
                     param('const std::vector<rbd::MultiBodyConfig>&', 'mbc'),
                     param('const tasks::qp::SolverData&', 'data')])

  spTask.add_method('Q', retval('Eigen::MatrixXd'), [], is_const=True)
  spTask.add_method('C', retval('Eigen::VectorXd'), [], is_const=True)

  # TargetObjectiveTask
  def toConstructor(hlTaskName):
    for t in hlTaskName:
      name = 'tasks::qp::%s *' % t
      toTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                              param('int', 'robotIndex'),
                              param(name, 'hlTask',
                                    transfer_ownership=False),
                              param('double', 'timeStep'),
                              param('double', 'duration'),
                              param('const Eigen::VectorXd&', 'objDot'),
                              param('double', 'weight')])

      toTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                              param('int', 'robotIndex'),
                              param(name, 'hlTask',
                                    transfer_ownership=False),
                              param('double', 'timeStep'),
                              param('double', 'duration'),
                              param('const Eigen::VectorXd&', 'objDot'),
                              param('const Eigen::VectorXd&', 'dimWeight'),
                              param('double', 'weight')])

  toConstructor(hlTaskName)

  toTask.add_method('duration', retval('double'), [], is_const=True)
  toTask.add_method('duration', None, [param('double', 'd')])

  toTask.add_method('iter', retval('int'), [], is_const=True)
  toTask.add_method('iter', None, [param('int', 'i')])

  toTask.add_method('nrIter', retval('int'), [], is_const=True)
  toTask.add_method('nrIter', None, [param('int', 'i')])

  toTask.add_method('objDot', retval('const Eigen::VectorXd&'), [], is_const=True)
  toTask.add_method('objDot', None, [param('const Eigen::VectorXd&', 'obj')])

  toTask.add_method('dimWeight', retval('const Eigen::VectorXd&'), [], is_const=True)
  toTask.add_method('dimWeight', None, [param('const Eigen::VectorXd&', 'w')])

  toTask.add_method('phi', retval('const Eigen::VectorXd&'), [], is_const=True)
  toTask.add_method('psi', retval('const Eigen::VectorXd&'), [], is_const=True)

  toTask.add_method('update', None,
                    [param('const std::vector<rbd::MultiBody>&', 'mb'),
                     param('const std::vector<rbd::MultiBodyConfig>&', 'mbc'),
                     param('const tasks::qp::SolverData&', 'data')])

  toTask.add_method('Q', retval('Eigen::MatrixXd'), [], is_const=True)
  toTask.add_method('C', retval('Eigen::VectorXd'), [], is_const=True)

  # PIDTask
  def pidConstructor(hlTaskName):
    for t in hlTaskName:
      name = 'tasks::qp::%s *' % t
      pidTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                              param('int', 'robotIndex'),
                              param(name, 'hlTask',
                                    transfer_ownership=False),
                              param('double', 'P'),
                              param('double', 'I'),
                              param('double', 'D'),
                              param('double', 'weight')])

      pidTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                              param('int', 'robotIndex'),
                              param(name, 'hlTask',
                                    transfer_ownership=False),
                              param('double', 'P'),
                              param('double', 'I'),
                              param('double', 'D'),
                              param('const Eigen::VectorXd&', 'dimWeight'),
                              param('double', 'weight')])

  pidConstructor(hlTaskName)

  pidTask.add_method('P', retval('double'), [], is_const=True)
  pidTask.add_method('P', None, [param('double', 'weight')])
  pidTask.add_method('I', retval('double'), [], is_const=True)
  pidTask.add_method('I', None, [param('double', 'weight')])
  pidTask.add_method('D', retval('double'), [], is_const=True)
  pidTask.add_method('D', None, [param('double', 'weight')])

  pidTask.add_method('error', None, [param('const Eigen::VectorXd&', 'error')])
  pidTask.add_method('errorD', None, [param('const Eigen::VectorXd&', 'errorD')])
  pidTask.add_method('errorI', None, [param('const Eigen::VectorXd&', 'errorI')])

  pidTask.add_method('dimWeight', retval('const Eigen::VectorXd&'), [], is_const=True)
  pidTask.add_method('dimWeight', None, [param('const Eigen::VectorXd&', 'dim')])

  pidTask.add_method('update', None,
                    [param('const std::vector<rbd::MultiBody>&', 'mb'),
                     param('const std::vector<rbd::MultiBodyConfig>&', 'mbc'),
                     param('const tasks::qp::SolverData&', 'data')])

  pidTask.add_method('Q', retval('Eigen::MatrixXd'), [], is_const=True)
  pidTask.add_method('C', retval('Eigen::VectorXd'), [], is_const=True)

  # PositionTask
  posTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                           param('int', 'robotIndex'),
                           param('int', 'bodyId'),
                           param('const Eigen::Vector3d&', 'pos'),
                           param('const Eigen::Vector3d&', 'bodyPoint',
                                 default_value='Eigen::Vector3d::Zero()')])

  posTask.add_method('position', None, [param('const Eigen::Vector3d&', 'pos')])
  posTask.add_method('position', retval('Eigen::Vector3d'), [], is_const=True)

  posTask.add_method('bodyPoint', None, [param('const Eigen::Vector3d&', 'point')])
  posTask.add_method('bodyPoint', retval('Eigen::Vector3d'), [], is_const=True)

  # OrientationTask
  oriTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                           param('int', 'robotIndex'),
                           param('int', 'bodyId'),
                           param('const Eigen::Quaterniond&', 'ori')])
  oriTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                           param('int', 'robotIndex'),
                           param('int', 'bodyId'),
                           param('const Eigen::Matrix3d&', 'ori')])

  oriTask.add_method('orientation', None, [param('const Eigen::Matrix3d&', 'ori')])
  oriTask.add_method('orientation', None, [param('const Eigen::Quaterniond&', 'ori')])
  oriTask.add_method('orientation', retval('Eigen::Matrix3d'), [], is_const=True)

  # SurfaceOrientationTask
  surfOriTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                               param('int', 'robotIndex'),
                               param('int', 'bodyId'),
                               param('const Eigen::Quaterniond&', 'ori'),
                               param('const sva::PTransformd&', 'X_b_s')])
  surfOriTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                               param('int', 'robotIndex'),
                               param('int', 'bodyId'),
                               param('const Eigen::Matrix3d&', 'ori'),
                               param('const sva::PTransformd&', 'X_b_s')])

  surfOriTask.add_method('orientation', None, [param('const Eigen::Matrix3d&', 'ori')])
  surfOriTask.add_method('orientation', None, [param('const Eigen::Quaterniond&', 'ori')])
  surfOriTask.add_method('orientation', retval('Eigen::Matrix3d'), [], is_const=True)

  # GazeTask
  gazeTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                            param('int', 'robotIndex'),
                            param('int', 'bodyId'),
                            param('const Eigen::Vector2d&', 'point2d'),
                            param('double', 'depthEstimate'),
                            param('const sva::PTransformd&', 'X_b_gaze'),
                            param('const Eigen::Vector2d&', 'point2d_ref',
                                  default_value='Eigen::Vector2d::Zero()')])
  gazeTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                            param('int', 'robotIndex'),
                            param('int', 'bodyId'),
                            param('const Eigen::Vector3d&', 'point3d'),
                            param('const sva::PTransformd&', 'X_b_gaze'),
                            param('const Eigen::Vector2d&', 'point2d_ref',
                                  default_value='Eigen::Vector2d::Zero()')])

  gazeTask.add_method('error', None, [param('const Eigen::Vector2d&', 'point2d'),
                                      param('const Eigen::Vector2d&', 'point2d_ref',
                                            default_value='Eigen::Vector2d::Zero()')])
  gazeTask.add_method('error', None, [param('const Eigen::Vector3d&', 'point3d'),
                                      param('const Eigen::Vector2d&', 'point2d_ref',
                                            default_value='Eigen::Vector2d::Zero()')])

  # PostureTask
  postureTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                               param('int', 'robotIndex'),
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

  postureTask.add_method('jointsStiffness', None,
                         [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                          param('std::vector<tasks::qp::JointStiffness>', 'js')])

  postureTask.add_method('update', None,
                    [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                     param('const std::vector<rbd::MultiBodyConfig>&', 'mbcs'),
                     param('const tasks::qp::SolverData&', 'data')])
  postureTask.add_method('eval', retval('Eigen::VectorXd'), [])

  # CoMTask
  comTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                           param('int', 'robotIndex'),
                           param('const Eigen::Vector3d&', 'com')])
  comTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                           param('int', 'robotIndex'),
                           param('const Eigen::Vector3d&', 'com'),
                           param('std::vector<double>', 'weight')],
                           throw=[dom_ex])

  comTask.add_method('com', None, [param('const Eigen::Vector3d&', 'com')])
  comTask.add_method('com', retval('const Eigen::Vector3d&', 'com'), [],
                     is_const=True)
  comTask.add_method('updateInertialParameters', None,
                     [param('const std::vector<rbd::MultiBody>&', 'mbs')])

  # MultiCoMTask
  multiCoMTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                param('std::vector<int>', 'robotIndexex'),
                                param('const Eigen::Vector3d&', 'com'),
                                param('double', 'stiffness'),
                                param('double', 'weight')])
  multiCoMTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                param('std::vector<int>', 'robotIndexex'),
                                param('const Eigen::Vector3d&', 'com'),
                                param('double', 'stiffness'),
                                param('const Eigen::Vector3d&', 'dimWeight'),
                                param('double', 'weight')])

  multiCoMTask.add_method('com', None, [param('const Eigen::Vector3d&', 'com')])
  multiCoMTask.add_method('com', retval('const Eigen::Vector3d&'), [],
                          is_const=True)
  multiCoMTask.add_method('updateInertialParameters', None,
                          [param('const std::vector<rbd::MultiBody>&', 'mbs')])
  multiCoMTask.add_method('stiffness', None, [param('double', 'stiffness')])
  multiCoMTask.add_method('stiffness', retval('double'), [],
                          is_const=True)
  multiCoMTask.add_method('dimWeight', None,
                          [param('const Eigen::Vector3d&', 'dimWeight')])
  multiCoMTask.add_method('dimWeight', retval('Eigen::Vector3d'), [],
                          is_const=True)
  multiCoMTask.add_method('eval', retval('Eigen::VectorXd'), [],
                          is_const=True)
  multiCoMTask.add_method('speed', retval('Eigen::VectorXd'), [],
                          is_const=True)
  multiCoMTask.add_method('update', None,
                          [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                           param('const std::vector<rbd::MultiBodyConfig>&', 'mbcs'),
                           param('const tasks::qp::SolverData&', 'data')])

  # MultiRobotTransformTask
  multiRobotTransformTask.add_constructor(
      [param('const std::vector<rbd::MultiBody>&', 'mbs'),
      param('int', 'r1Index'), param('int', 'r2Index'),
      param('int', 'r1BodyIndex'), param('int', 'r2BodyIndex'),
      param('const sva::PTransformd&', 'X_r1b_r1s'),
      param('const sva::PTransformd&', 'X_r2b_r2s'),
      param('double', 'stiffness'), param('double', 'weight')])

  multiRobotTransformTask.add_method('X_r1b_r1s', None,
    [param('const sva::PTransformd&', 'X_r1b_r1s')])
  multiRobotTransformTask.add_method('X_r1b_r1s', retval('sva::PTransformd'),
                                     [], is_const=True)
  multiRobotTransformTask.add_method('X_r2b_r2s', None,
                                     [param('const sva::PTransformd&', 'X_r2b_r2s')])
  multiRobotTransformTask.add_method('X_r2b_r2s', retval('sva::PTransformd'),
                                     [], is_const=True)

  multiRobotTransformTask.add_method('eval', retval('const Eigen::VectorXd&'), [],
                                     is_const=True)
  multiRobotTransformTask.add_method('speed', retval('const Eigen::VectorXd&'), [],
                                     is_const=True)

  multiRobotTransformTask.add_method('stiffness', None, [param('double', 'stiffness')])
  multiRobotTransformTask.add_method('stiffness', retval('double'), [],
                                     is_const=True)
  multiRobotTransformTask.add_method('dimWeight', None,
                                     [param('const Eigen::Vector6d&', 'dimWeight')])
  multiRobotTransformTask.add_method('dimWeight', retval('Eigen::Vector6d'), [],
                                     is_const=True)
  multiRobotTransformTask.add_method('update', None,
    [param('const std::vector<rbd::MultiBody>&', 'mbs'),
      param('const std::vector<rbd::MultiBodyConfig>&', 'mbcs'),
      param('const tasks::qp::SolverData&', 'data')])

  # MomentumTask
  momTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                           param('int', 'robotIndex'),
                           param('const sva::ForceVecd&', 'mom')])

  momTask.add_method('momentum', None, [param('const sva::ForceVecd&', 'mom')])
  momTask.add_method('momentum', retval('const sva::ForceVecd&', 'momentum'), [],
                     is_const=True)

  # ContactTask
  contactTask.add_constructor([param('const tasks::qp::ContactId&', 'contactId'),
                               param('double', 'stiffness'),
                               param('double', 'weight')])

  contactTask.add_method('error', None, [param('const Eigen::Vector3d&', 'error')])
  contactTask.add_method('errorD', None, [param('const Eigen::Vector3d&', 'errorD')])

  # GripperTorqueTask
  gripperTorqueTask.add_constructor([param('const tasks::qp::ContactId&', 'contactId'),
                                     param('const Eigen::Vector3d&', 'origin'),
                                     param('const Eigen::Vector3d&', 'axis'),
                                     param('double', 'weight')])

  # LinVelocityTask
  linVelTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                              param('int', 'robotIndex'),
                              param('int', 'bodyId'),
                              param('const Eigen::Vector3d&', 'pos'),
                              param('const Eigen::Vector3d&', 'bodyPoint',
                                    default_value='Eigen::Vector3d::Zero()')])

  linVelTask.add_method('velocity', None, [param('const Eigen::Vector3d&', 'pos')])
  linVelTask.add_method('velocity', retval('Eigen::Vector3d'), [], is_const=True)

  linVelTask.add_method('bodyPoint', None, [param('const Eigen::Vector3d&', 'point')])
  linVelTask.add_method('bodyPoint', retval('Eigen::Vector3d'), [], is_const=True)

  # OrientationTrackingTask
  oriTrackTask.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                param('int', 'robotIndex'),
                                param('int', 'bodyId'),
                                param('const Eigen::Vector3d&', 'bodyPoint'),
                                param('const Eigen::Vector3d&', 'bodyAxis'),
                                param('const std::vector<int>&', 'trackingJointId'),
                                param('const Eigen::Vector3d&', 'trackedPoint')])

  oriTrackTask.add_method('trackedPoint', None, [param('const Eigen::Vector3d&', 'bPoint')])
  oriTrackTask.add_method('trackedPoint', retval('Eigen::Vector3d'), [], is_const=True)
  oriTrackTask.add_method('bodyPoint', None, [param('const Eigen::Vector3d&', 'tPoint')])
  oriTrackTask.add_method('bodyPoint', retval('Eigen::Vector3d'), [], is_const=True)
  oriTrackTask.add_method('bodyAxis', None, [param('const Eigen::Vector3d&', 'axis')])
  oriTrackTask.add_method('bodyAxis', retval('Eigen::Vector3d'), [], is_const=True)

  # JointsSelector
  jointsSelector.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                  param('int', 'robotIndex'),
                                  param('tasks::qp::HighLevelTask*', 'hl',
                                        transfer_ownership=False),
                                  param('const std::vector<int>&', 'selectedJointsId')])
  jointsSelector.add_method('ActiveJoints', retval('tasks::qp::JointsSelector'),
                            [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                             param('int', 'robotIndex'),
                             param('tasks::qp::HighLevelTask*', 'hl',
                                   transfer_ownership=False),
                             param('const std::vector<int>&', 'activeJointsId')],
                            is_static=True)
  jointsSelector.add_method('UnactiveJoints', retval('tasks::qp::JointsSelector'),
                            [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                             param('int', 'robotIndex'),
                             param('tasks::qp::HighLevelTask*', 'hl',
                                   transfer_ownership=False),
                             param('const std::vector<int>&', 'unactiveJointsId')],
                            is_static=True)


  # MotionConstr
  def addMotionDefault(motion):
    motion.add_method('computeTorque', None, [param('const Eigen::VectorXd&', 'alphaD'),
                                              param('const Eigen::VectorXd&', 'lambda')])
    motion.add_method('torque', retval('Eigen::VectorXd'), [], is_const=True)
    motion.add_method('torque', None, [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                       param('std::vector<rbd::MultiBodyConfig>&', 'mbcs')], is_const=True)

  motionConstr.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                param('int', 'robotIndex'),
                                param('const tasks::TorqueBound&', 'tb')])
  addMotionDefault(motionConstr)

  # MotionPolyConstr
  motionPolyConstr.add_constructor([param('const std::vector<rbd::MultiBody>', 'mbs'),
                                    param('int', 'robotIndex'),
                                    param('const tasks::PolyTorqueBound&', 'ptb')])
  addMotionDefault(motionPolyConstr)

  # MotionSpringConstr
  motionSpringConstr.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                      param('int', 'robotIndex'),
                                      param('const tasks::TorqueBound&', 'tb'),
                                      param('const std::vector<tasks::qp::SpringJoint>&', 'springs')])
  addMotionDefault(motionSpringConstr)

  # PositiveLambda
  positiveLambdaConstr.add_constructor([])

  # ContactConstrCommon
  contactConstrCommon.add_method('addVirtualContact', retval('bool'),
                                 [param('const tasks::qp::ContactId&', 'contactId')])
  contactConstrCommon.add_method('removeVirtualContact', retval('bool'),
                                 [param('const tasks::qp::ContactId&', 'contactId')])
  contactConstrCommon.add_method('resetVirtualContacts', None, [])

  contactConstrCommon.add_method('addDofContact', retval('bool'),
                                 [param('const tasks::qp::ContactId&', 'contactId'),
                                  param('const Eigen::MatrixXd&', 'dof')])
  contactConstrCommon.add_method('removeDofContact', retval('bool'),
                                 [param('const tasks::qp::ContactId&', 'contactId')])
  contactConstrCommon.add_method('resetDofContacts', None, [])

  # ContactAccConstr
  contactAccConstr.add_constructor([])
  contactAccConstr.add_method('updateDofContacts', None, [])

  # ContactSpeedConstr
  contactSpeedConstr.add_constructor([param('double', 'timeStep')])
  contactSpeedConstr.add_method('updateDofContacts', None, [])

  # ContactPosConstr
  contactPosConstr.add_constructor([param('double', 'timeStep')])
  contactPosConstr.add_method('updateDofContacts', None, [])

  # CollisionConstr
  collisionConstr.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                   param('double', 'step')])
  collisionConstr.add_method('addCollision', None,
                             [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                             param('int', 'collId'),
                             param('int', 'r1Index'), param('int', 'r1BodyId'),
                             param('sch::S_Object*', 'body1', transfer_ownership=False),
                             param('const sva::PTransformd&', 'X_op1_o1'),
                             param('int', 'r2Index'), param('int', 'r2BodyId'),
                             param('sch::S_Object*', 'body2', transfer_ownership=False),
                             param('const sva::PTransformd&', 'X_op2_o2'),
                             param('double', 'di'),
                             param('double', 'ds'),
                             param('double', 'damping'),
                             param('double', 'dampingOff', default_value='0.')])

  collisionConstr.add_method('rmCollision', retval('bool'),
                             [param('int', 'collId')])
  collisionConstr.add_method('nrCollisions', retval('int'),
                             [], is_const=True)
  collisionConstr.add_method('reset', None, []),

  collisionConstr.add_method('updateNrCollisions', None, []),


  # CoMIncPlaneConstr
  comIncPlaneConstr.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                     param('int', 'robotIndex'),
                                     param('double', 'step')])
  comIncPlaneConstr.add_method('addPlane', None,
                               [param('int', 'planeId'),
                                param('const Eigen::Vector3d&', 'normal'),
                                param('double', 'offset'),
                                param('double', 'di'),
                                param('double', 'ds'),
                                param('double', 'damping'),
                                param('double', 'dampingOff', default_value='0.')])

  comIncPlaneConstr.add_method('rmPlane', retval('bool'), [param('int', 'planeId')])
  comIncPlaneConstr.add_method('nrPlanes', retval('int'),
                               [], is_const=True)
  comIncPlaneConstr.add_method('reset', None, []),

  comIncPlaneConstr.add_method('updateNrPlanes', None, []),


  # JointLimitsConstr
  jointLimitsConstr.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                     param('int', 'robotIndex'),
                                     param('const tasks::QBound&', 'qb'),
                                     param('double', 'step')])

  # DamperJointLimitsConstr
  damperJointLimitsConstr.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                           param('int', 'robotIndex'),
                                           param('const tasks::QBound&', 'qb'),
                                           param('const tasks::AlphaBound&', 'ab'),
                                           param('double', 'interPercent'),
                                           param('double', 'securityPercent'),
                                           param('double', 'damperOffset'),
                                           param('double', 'step')])

  # GripperTorqueConstr
  gripperTorqueConstr.add_constructor([])
  gripperTorqueConstr.add_method('addGripper', None,
                                 [param('const tasks::qp::ContactId&', 'contactId'),
                                  param('double', 'torqueLimit'),
                                  param('const Eigen::Vector3d&', 'origin'),
                                  param('const Eigen::Vector3d&', 'axis')]);
  gripperTorqueConstr.add_method('rmGripper', retval('bool'),
                                 [param('const tasks::qp::ContactId&', 'contactId')])
  gripperTorqueConstr.add_method('reset', None, [])

  # BoundedSpeedConstr
  boundedSpeedConstr.add_constructor([param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                      param('int', 'robotIndex'),
                                      param('double', 'timeStep')])
  boundedSpeedConstr.add_method('addBoundedSpeed', None,
                                [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                 param('int', 'bodyId'),
                                 param('const Eigen::Vector3d&', 'bodyPoint'),
                                 param('const Eigen::MatrixXd&', 'dof'),
                                 param('const Eigen::VectorXd&', 'speed')])
  boundedSpeedConstr.add_method('addBoundedSpeed', None,
                                [param('const std::vector<rbd::MultiBody>&', 'mbs'),
                                 param('int', 'bodyId'),
                                 param('const Eigen::Vector3d&', 'bodyPoint'),
                                 param('const Eigen::MatrixXd&', 'dof'),
                                 param('const Eigen::VectorXd&', 'lowerSpeed'),
                                 param('const Eigen::VectorXd&', 'upperSpeed')])
  boundedSpeedConstr.add_method('removeBoundedSpeed', retval('bool'),
                                [param('int', 'bodyId')])
  boundedSpeedConstr.add_method('resetBoundedSpeeds', None, [])
  boundedSpeedConstr.add_method('nrBoundedSpeeds', retval('int'), [])

  boundedSpeedConstr.add_method('updateBoundedSpeeds', None, [])


  def add_add_remove_solver(constr):
    for c in constr:
      c.add_method('addToSolver', None, [param('tasks::qp::QPSolver&', 'solver')])
      c.add_method('removeFromSolver', None, [param('tasks::qp::QPSolver&',  'solver')])
  add_add_remove_solver(constrList)




if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.exit(1)

  tasks = Module('_tasks', cpp_namespace='::tasks')
  tasks.add_include('<Tasks.h>')
  tasks.add_include('<QPSolver.h>')
  tasks.add_include('<QPTasks.h>')
  tasks.add_include('<QPConstr.h>')
  tasks.add_include('<QPContactConstr.h>')
  tasks.add_include('<QPMotionConstr.h>')
  tasks.add_include('<Bounds.h>')

  tasks.add_include('<RBDyn/MultiBodyConfig.h>')

  tasks.add_include('<sch/S_Object/S_Object.h>')
  tasks.add_include('<sch/S_Object/S_Sphere.h>')
  tasks.add_include('<sch/CD/CD_Pair.h>')

  dom_ex = tasks.add_exception('std::domain_error', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')
  out_ex = tasks.add_exception('std::out_of_range', foreign_cpp_namespace=' ',
                               message_rvalue='%(EXC)s.what()')

  build_boost_timer(tasks)

  # import Eigen3, sva and rbd types
  import_eigen3_types(tasks)
  import_sva_types(tasks)
  import_rbd_types(tasks)
  import_sch_types(tasks)

  posTask = tasks.add_class('PositionTask')
  oriTask = tasks.add_class('OrientationTask')
  surfOriTask = tasks.add_class('SurfaceOrientationTask')
  gazeTask = tasks.add_class('GazeTask')
  postureTask = tasks.add_class('PostureTask')
  comTask = tasks.add_class('CoMTask')
  multiCoMTask = tasks.add_class('MultiCoMTask')
  momTask = tasks.add_class('MomentumTask')
  linVelTask = tasks.add_class('LinVelocityTask')
  oriTrackTask = tasks.add_class('OrientationTrackingTask')
  multiRobotTransformTask = tasks.add_class('MultiRobotTransformTask')

  # build list type
  tasks.add_container('std::vector<int>', 'int', 'vector')
  tasks.add_container('std::vector<double>', 'double', 'vector')
  tasks.add_container('std::vector<std::vector<double> >', 'std::vector<double>', 'vector')
  tasks.add_container('std::vector<rbd::MultiBody>',
                      'rbd::MultiBody', 'vector')
  tasks.add_container('std::vector<rbd::MultiBodyConfig>',
                      'rbd::MultiBodyConfig', 'vector')
  tasks.add_container('std::vector<sva::MotionVecd>', 'sva::MotionVecd', 'vector')
  tasks.add_container('std::vector<std::vector<sva::MotionVecd> >',
                      'std::vector<sva::MotionVecd>', 'vector')

  build_tasks(posTask, oriTask, surfOriTask, gazeTask, postureTask, comTask, multiCoMTask,
              momTask, linVelTask, oriTrackTask, multiRobotTransformTask)

  # qp
  build_qp(tasks)




  with open(sys.argv[1], 'w') as f:
    tasks.generate(f)

