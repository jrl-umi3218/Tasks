# Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
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

ENABLE_TESTING()

set(BOOST_COMPONENTS unit_test_framework timer system)
search_for_boost()
add_definitions(-DBOOST_TEST_DYN_LINK)
IF(WIN32)
  #This is one of the way to avoid link errors related to static variables in program_options
  ADD_DEFINITIONS( -DBOOST_ALL_DYN_LINK )
ENDIF(WIN32)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

set(HEADERS arms.h)

macro(addUnitTest name)
  add_executable(${name} ${name}.cpp ${HEADERS})
  target_link_libraries(${name} ${Boost_LIBRARIES} Tasks)
  PKG_CONFIG_USE_DEPENDENCY(${name} sch-core)
  PKG_CONFIG_USE_DEPENDENCY(${name} SpaceVecAlg)
  PKG_CONFIG_USE_DEPENDENCY(${name} RBDyn)
  PKG_CONFIG_USE_DEPENDENCY(${name} eigen-qld)
  if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(${name} eigen-lssol)
  endif()
  add_test(${name}Unit ${name})
  # Adding a project configuration file (for MSVC only)
  GENERATE_MSVC_DOT_USER_FILE(${name})
endmacro(addUnitTest)

addUnitTest(QPSolverTest)
addUnitTest(QPMultiRobotTest)
addUnitTest(TasksTest)
addUnitTest(AllocationTest)
