# Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
#
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

MACRO(GET_TASKS_PROPERTY PROPERTY)
  GET_TARGET_PROPERTY(TASKS_${PROPERTY} Tasks ${PROPERTY})
  IF(NOT TASKS_${PROPERTY})
    SET(TASKS_${PROPERTY} "")
  ENDIF()
ENDMACRO()
GET_TASKS_PROPERTY(COMPILE_FLAGS)
GET_TASKS_PROPERTY(INCLUDE_DIRECTORIES)
GET_TASKS_PROPERTY(LINK_FLAGS)
GET_TASKS_PROPERTY(LINK_LIBRARIES)
GET_TASKS_PROPERTY(LOCATION)
CONFIGURE_FILE(setup.in.py ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)

# Build the bindings locally at build time for test purposes
ADD_CUSTOM_TARGET(tasks-python-bindings ALL
  COMMAND python setup.py build_ext --inplace
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating local Tasks Python bindings"
  DEPENDS include/qp_wrapper.hpp tasks/c_tasks.pxd tasks/tasks.pxd tasks/tasks.pyx
          tasks/qp/c_qp.pxd tasks/qp/qp.pxd tasks/qp/qp.pyx
)
ADD_DEPENDENCIES(tasks-python-bindings Tasks)

if(NOT ${DISABLE_TESTS})
  if(NOT DEFINED PYTHON_DEB_ROOT)
    ADD_TEST(NAME PythonBindingsTest
      COMMAND nosetests
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
  endif()
endif()

SET(PIP_EXTRA_OPTIONS "")
IF(${PYTHON_BINDING_USER_INSTALL})
  SET(PIP_EXTRA_OPTIONS "--user")
ENDIF()
# Install the bindings
IF(DEFINED PYTHON_DEB_ROOT)
  ADD_CUSTOM_TARGET(install-tasks-python-bindings
    COMMAND python setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb && python3 setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Install Tasks Python bindings (Debian package)"
  )
ELSE()
  ADD_CUSTOM_TARGET(install-tasks-python-bindings
    COMMAND pip install . ${PIP_EXTRA_OPTIONS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Install Tasks Python bindings"
  )
ENDIF()
INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_MAKE_PROGRAM} install-tasks-python-bindings)")
