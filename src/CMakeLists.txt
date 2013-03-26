enable_language(Fortran)

set(SOURCES Tasks.cpp QPSolver.cpp QPTasks.cpp QPConstr.cpp QL.cpp QL.f
            QPContacts.cpp QPSolverData.cpp
            PGJacobian.cpp PostureGenerator.cpp PGConstr.cpp)
set(HEADERS Tasks.h QPSolver.h QPTasks.h QPConstr.h QL.h
            QPContacts.h QPSolverData.h
            PGJacobian.h PostureGenerator.h PGConstr.h)

add_definitions("-DWITH_OPENGL")

add_library(Tasks SHARED ${SOURCES} ${HEADERS})


PKG_CONFIG_USE_DEPENDENCY(Tasks ipopt)
PKG_CONFIG_USE_DEPENDENCY(Tasks scd)
PKG_CONFIG_USE_DEPENDENCY(Tasks SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(Tasks RBDyn)


set(INSTALL_PATH include/Tasks)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS Tasks DESTINATION "lib")

