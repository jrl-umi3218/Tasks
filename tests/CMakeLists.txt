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

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_definitions(-DBOOST_TEST_DYN_LINK)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

set(HEADERS arms.h)

macro(addUnitTest name)
  if(${UNIT_TESTS})
    add_executable(${name} ${name}.cpp ${HEADERS})
    target_link_libraries(${name} ${Boost_LIBRARIES} RBDyn Tasks)
    PKG_CONFIG_USE_DEPENDENCY(${name} sch-core)
    PKG_CONFIG_USE_DEPENDENCY(${name} SpaceVecAlg)
    PKG_CONFIG_USE_DEPENDENCY(${name} RBDyn)
    PKG_CONFIG_USE_DEPENDENCY(${name} eigen-qld)
    if(${EIGEN_LSSOL_FOUND})
      PKG_CONFIG_USE_DEPENDENCY(${name} eigen-lssol)
    endif()
    add_test(${name}Unit ${name})
  endif()
endmacro(addUnitTest)

addUnitTest(QPSolverTest)
addUnitTest(QPMultiRobotTest)
addUnitTest(TasksTest)
