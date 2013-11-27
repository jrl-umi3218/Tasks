set(SOURCES Tasks.cpp QPSolver.cpp QPTasks.cpp QPConstr.cpp
            QPContacts.cpp QPSolverData.cpp)
set(HEADERS Tasks.h QPSolver.h QPTasks.h QPConstr.h
            QPContacts.h QPSolverData.h)

if(${POSTURE_GENERATOR})
  list(APPEND SOURCES PGJacobian.cpp PostureGenerator.cpp PGConstr.cpp)
  list(APPEND HEADERS PGJacobian.h PostureGenerator.h PGConstr.h)
endif()


add_library(Tasks SHARED ${SOURCES} ${HEADERS})


PKG_CONFIG_USE_DEPENDENCY(Tasks scd)
PKG_CONFIG_USE_DEPENDENCY(Tasks SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(Tasks RBDyn)
PKG_CONFIG_USE_DEPENDENCY(Tasks EigenQP)

if(${POSTURE_GENERATOR})
  PKG_CONFIG_USE_DEPENDENCY(Tasks ipopt)
endif()


set(INSTALL_PATH include/Tasks)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS Tasks DESTINATION "lib")

