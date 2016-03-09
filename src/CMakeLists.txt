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
# along with Tasks.  If not, see <http://www.gnu.org/licenses/>.

set(SOURCES Tasks.cpp QPSolver.cpp QPTasks.cpp QPConstr.cpp
            QPContacts.cpp QPSolverData.cpp QPMotionConstr.cpp
            GenQPSolver.cpp QPContactConstr.cpp QLDQPSolver.cpp)
set(HEADERS Tasks.h QPSolver.h QPTasks.h QPConstr.h
            QPContacts.h QPSolverData.h QPMotionConstr.h
            GenQPSolver.h Bounds.h QPContactConstr.h QLDQPSolver.h)
set(PRIVATE_HEADERS utils.h GenQPUtils.h)

if(${EIGEN_LSSOL_FOUND})
  list(APPEND SOURCES LSSOLQPSolver.cpp)
  list(APPEND HEADERS LSSOLQPSolver.h)
endif()

set(BOOST_COMPONENTS timer system)
SEARCH_FOR_BOOST()

include_directories(${Boost_INCLUDE_DIRS})

add_library(Tasks SHARED ${SOURCES} ${HEADERS} ${PRIVATE_HEADERS})
set_target_properties(Tasks PROPERTIES COMPILE_FLAGS "-Dtasks_EXPORTS")

PKG_CONFIG_USE_DEPENDENCY(Tasks sch-core)
PKG_CONFIG_USE_DEPENDENCY(Tasks SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(Tasks RBDyn)
PKG_CONFIG_USE_DEPENDENCY(Tasks eigen-qld)

if(${EIGEN_LSSOL_FOUND})
  PKG_CONFIG_USE_DEPENDENCY(Tasks eigen-lssol)
  ADD_DEFINITIONS(-DLSSOL_SOLVER_FOUND)
endif()

target_link_libraries(Tasks ${Boost_LIBRARIES})


set(INSTALL_PATH include/Tasks)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS Tasks 
                      RUNTIME DESTINATION bin
                      LIBRARY DESTINATION lib
                      ARCHIVE DESTINATION lib
)
