#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

set(SOURCES
    Tasks.cpp
    QPSolver.cpp
    QPTasks.cpp
    QPConstr.cpp
    QPContacts.cpp
    QPSolverData.cpp
    QPMotionConstr.cpp
    GenQPSolver.cpp
    QPContactConstr.cpp
    QLDQPSolver.cpp
)
set(HEADERS
    Tasks/Tasks.h
    Tasks/QPSolver.h
    Tasks/QPTasks.h
    Tasks/QPConstr.h
    Tasks/QPContacts.h
    Tasks/QPSolverData.h
    Tasks/QPMotionConstr.h
    Tasks/GenQPSolver.h
    Tasks/Bounds.h
    Tasks/QPContactConstr.h
)
set(PRIVATE_HEADERS utils.h GenQPUtils.h QLDQPSolver.h)

if(${eigen-lssol_FOUND})
  list(APPEND SOURCES LSSOLQPSolver.cpp)
  list(APPEND PRIVATE_HEADERS LSSOLQPSolver.h)
endif()

add_project_dependency(Boost REQUIRED COMPONENTS timer)

add_library(Tasks SHARED ${SOURCES} ${HEADERS} ${PRIVATE_HEADERS})
target_link_libraries(
  Tasks PUBLIC RBDyn::RBDyn sch-core::sch-core eigen-qld::eigen-qld Boost::timer
               Boost::disable_autolinking Boost::dynamic_linking
)
if(${eigen-lssol_FOUND})
  target_link_libraries(Tasks PUBLIC eigen-lssol::eigen-lssol)
  target_compile_definitions(Tasks PRIVATE -DLSSOL_SOLVER_FOUND)
endif()
target_include_directories(
  Tasks
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../include>
         $<INSTALL_INTERFACE:include>
)
set_target_properties(Tasks PROPERTIES COMPILE_FLAGS "-Dtasks_EXPORTS")
set_target_properties(
  Tasks PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR} VERSION ${PROJECT_VERSION}
)
set_target_properties(Tasks PROPERTIES CXX_STANDARD 17)

pkg_config_append_libs(Tasks)

install(
  TARGETS Tasks
  EXPORT "${TARGETS_EXPORT_NAME}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  INCLUDES
  DESTINATION "${INCLUDE_INSTALL_DIR}"
)
install(FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DESTINATION})
