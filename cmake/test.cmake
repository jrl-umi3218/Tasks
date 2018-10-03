# Copyright (C) 2008-2014,2018 LAAS-CNRS, JRL AIST-CNRS.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#.rst:
# .. variable:: DISABLE_TESTS
#
#   Boolean variable to configure unit test compilation declared with
#   :command:`ADD_UNIT_TEST`.

#   A target *build_tests* is added to compile the tests and
#   a test that run target build_tests is run before all other tests.
#   So command ``make test`` compiles and runs the unit-tests.
#
#   * if ``OFF`` (default), nothing special is done.
#   * if ``ON``, the unit-test is not compiled with target *all*.
IF(NOT DEFINED DISABLE_TESTS)
  SET(DISABLE_TESTS OFF)
ENDIF(NOT DEFINED DISABLE_TESTS)

ADD_CUSTOM_TARGET(build_tests)
IF(DISABLE_TESTS)
   ADD_TEST(ctest_build_tests "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target build_tests)
ENDIF()

#.rst:
# .. command:: ADD_UNIT_TEST (NAME SOURCE)
#
#   The behaviour of this function depends on :variable:`DISABLE_TESTS` option.
#
MACRO(ADD_UNIT_TEST NAME SOURCE)
  IF(DISABLE_TESTS)
    ADD_EXECUTABLE(${NAME} EXCLUDE_FROM_ALL ${SOURCE})
  ELSE(DISABLE_TESTS)
    ADD_EXECUTABLE(${NAME} ${SOURCE})
  ENDIF(DISABLE_TESTS)

  ADD_DEPENDENCIES(build_tests ${NAME})

  ADD_TEST(${NAME} ${RUNTIME_OUTPUT_DIRECTORY}/${NAME})
  # Support definition of DYLD_LIBRARY_PATH for OSX systems
  IF(APPLE)
    SET_TESTS_PROPERTIES(${NAME} PROPERTIES ENVIRONMENT "DYLD_LIBRARY_PATH=$ENV{DYLD_LIBRARY_PATH}")
  ENDIF(APPLE)

  IF(DISABLE_TESTS)
    SET_TESTS_PROPERTIES(${NAME} PROPERTIES DEPENDS ctest_build_tests)
  ENDIF(DISABLE_TESTS)
ENDMACRO(ADD_UNIT_TEST NAME SOURCE)

#.rst:
# .. command:: ADD_PYTHON_UNIT_TEST (NAME SOURCE [MODULE_PATH])
#
#   Add a test called `NAME` that runs an equivalent of ``python ${SOURCE}``,
#   optionnaly with a `PYTHONPATH` set to `CMAKE_BINARY_DIR/MODULE_PATH` if `MODULE_PATH` is set.
#   `SOURCE` is relative to `PROJECT_SOURCE_DIR`
#
#   .. note:: :command:`FINDPYTHON` should have been called first.
#
MACRO(ADD_PYTHON_UNIT_TEST NAME SOURCE)
  SET(MODULE_PATH "${ARGN}")  # ARGN is not a variable
  ADD_TEST(${NAME} ${PYTHON_EXECUTABLE} "${PROJECT_SOURCE_DIR}/${SOURCE}")
  IF(MODULE_PATH)
    SET(PYTHONPATH "${CMAKE_BINARY_DIR}/${MODULE_PATH}")
    IF(DEFINED ENV{PYTHONPATH})
      SET(PYTHONPATH "${PYTHONPATH}:$ENV{PYTHONPATH}")
    ENDIF(DEFINED ENV{PYTHONPATH})
    SET_PROPERTY(TEST ${NAME} PROPERTY ENVIRONMENT "PYTHONPATH=${PYTHONPATH}")
  ENDIF(MODULE_PATH)
ENDMACRO(ADD_PYTHON_UNIT_TEST NAME SOURCE)

# DEFINE_UNIT_TEST(NAME LIB)
# ----------------------
#
# Compile a program and add it as a test
#
MACRO(DEFINE_UNIT_TEST NAME LIB)
  ADD_UNIT_TEST(${NAME} ${NAME}.cc)
  TARGET_LINK_LIBRARIES(${NAME} ${LIB})
ENDMACRO(DEFINE_UNIT_TEST)
