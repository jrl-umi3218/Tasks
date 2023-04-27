#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

make_cython_bindings(
  tasks
  TARGETS
  Tasks
  MODULES
  tasks.tasks
  tasks.qp.qp
  EXPORT_SOURCES
  tasks/c_tasks.pxd
  tasks/tasks.pxd
  tasks/__init__.py
  tasks/qp/c_qp.pxd
  tasks/qp/qp.pxd
  tasks/qp/__init__.py
  PRIVATE_SOURCES
  tasks/tasks.pyx
  tasks/qp/c_qp_private.pxd
  tasks/qp/qp.pyx
  include/qp_wrapper.hpp
  tests/arms.py
  tests/test_qp_multi_robot.py
  tests/test_qp_solver.py
)
