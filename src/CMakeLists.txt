#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

set(SOURCES Tasks.cpp QPSolver.cpp QPTasks.cpp QPConstr.cpp
            QPContacts.cpp QPSolverData.cpp QPMotionConstr.cpp
            GenQPSolver.cpp QPContactConstr.cpp QLDQPSolver.cpp)
set(HEADERS Tasks/Tasks.h Tasks/QPSolver.h Tasks/QPTasks.h Tasks/QPConstr.h
            Tasks/QPContacts.h Tasks/QPSolverData.h Tasks/QPMotionConstr.h
            Tasks/GenQPSolver.h Tasks/Bounds.h Tasks/QPContactConstr.h)
set(PRIVATE_HEADERS utils.h GenQPUtils.h QLDQPSolver.h)

if(${EIGEN_LSSOL_FOUND})
  list(APPEND SOURCES LSSOLQPSolver.cpp)
  list(APPEND PRIVATE_HEADERS LSSOLQPSolver.h)
endif()

set(BOOST_COMPONENTS timer)
SEARCH_FOR_BOOST()
IF(WIN32)
  #This is one of the way to avoid link errors related to static variables in program_options
  ADD_DEFINITIONS( -DBOOST_ALL_DYN_LINK )
ENDIF(WIN32)


add_library(Tasks SHARED ${SOURCES} ${HEADERS} ${PRIVATE_HEADERS})
target_include_directories(Tasks BEFORE PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:${include_install_destination}>)
target_include_directories(Tasks SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
set_target_properties(Tasks PROPERTIES COMPILE_FLAGS "-Dtasks_EXPORTS")
set_target_properties(Tasks PROPERTIES SOVERSION 0 VERSION 0.9.0)

PKG_CONFIG_USE_DEPENDENCY(Tasks sch-core)
PKG_CONFIG_USE_DEPENDENCY(Tasks SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(Tasks RBDyn)
PKG_CONFIG_USE_DEPENDENCY(Tasks eigen-qld)

if(${EIGEN_LSSOL_FOUND})
  PKG_CONFIG_USE_DEPENDENCY(Tasks eigen-lssol)
  ADD_DEFINITIONS(-DLSSOL_SOLVER_FOUND)
endif()

target_link_libraries(Tasks PUBLIC ${Boost_LIBRARIES})

# Targets:
#   * <prefix>/lib/libbar.a
#   * <prefix>/lib/libbaz.a
#   * header location after install: <prefix>/include/foo/Bar.hpp
#   * headers can be included by C++ code `#include <foo/Bar.hpp>`
install(
    TARGETS Tasks
    EXPORT "${TARGETS_EXPORT_NAME}"
    LIBRARY DESTINATION "lib"
    ARCHIVE DESTINATION "lib"
    RUNTIME DESTINATION "bin"
    INCLUDES DESTINATION "${INCLUDE_INSTALL_DIR}"
)
install(FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DESTINATION})

