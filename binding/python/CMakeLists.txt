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

SET(SETUP_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION>/")
IF(DEFINED CMAKE_BUILD_TYPE)
  FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE}")
ELSE()
  FOREACH(CFG ${CMAKE_CONFIGURATION_TYPES})
    FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CFG}")
  ENDFOREACH()
ENDIF()
CONFIGURE_FILE(setup.in.py "${CMAKE_CURRENT_BINARY_DIR}/setup.in.py")
FILE(GENERATE
      OUTPUT "${SETUP_LOCATION}/setup.py"
      INPUT  "${CMAKE_CURRENT_BINARY_DIR}/setup.in.py")

# Build the bindings locally at build time for test purposes
ADD_CUSTOM_TARGET(tasks-python-bindings ALL
  COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" python setup.py build_ext --inplace
  COMMENT "Generating local Tasks Python bindings"
  DEPENDS include/qp_wrapper.hpp tasks/c_tasks.pxd tasks/tasks.pxd tasks/tasks.pyx
          tasks/qp/c_qp.pxd tasks/qp/qp.pxd tasks/qp/qp.pyx
  SOURCES include/qp_wrapper.hpp tasks/c_tasks.pxd tasks/tasks.pxd tasks/tasks.pyx
          tasks/qp/c_qp.pxd tasks/qp/qp.pxd tasks/qp/qp.pyx
)
ADD_DEPENDENCIES(tasks-python-bindings Tasks)

if(NOT ${DISABLE_TESTS})
  if(NOT DEFINED PYTHON_DEB_ROOT)
    ADD_TEST(NAME PythonBindingsTest
      COMMAND ${CMAKE_COMMAND} -E chdir "$<TARGET_FILE_DIR:Tasks>" nosetests ${SETUP_LOCATION}
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
    COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" python setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb
    COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" python3 setup.py install --root=${PYTHON_DEB_ROOT} --install-layout=deb
    COMMENT "Install Tasks Python bindings (Debian package)"
  )
ELSE()
  ADD_CUSTOM_TARGET(install-tasks-python-bindings
    COMMAND ${CMAKE_COMMAND} -E chdir "${SETUP_LOCATION}" pip install . ${PIP_EXTRA_OPTIONS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Install Tasks Python bindings"
  )
ENDIF()
INSTALL(CODE "EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" --build \"${CMAKE_BINARY_DIR}\" --config \${CMAKE_INSTALL_CONFIG_NAME} --target install-tasks-python-bindings)")
