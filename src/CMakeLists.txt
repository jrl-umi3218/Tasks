enable_language(Fortran)

set(SOURCES Tasks.cpp QPSolver.cpp QPTasks.cpp QPConstr.cpp QL.cpp QL.f
            PGJacobian.cpp)
set(HEADERS Tasks.h QPSolver.h QPTasks.h QPConstr.h QL.h
            PGJacobian.h)

pkg_check_modules(EIGEN3 REQUIRED eigen3)
find_path(SVA_INCLUDE_DIRS SpaceVecAlg/SpaceVecAlg
                           /usr/include
                           /usr/local/include)

find_path(RBD_INCLUDE_DIRS RBDyn/
                           /usr/include
                           /usr/local/include)

find_path(SCD_INCLUDE_DIRS SCD/
                           /usr/include
                           /usr/local/include)

include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${SVA_INCLUDE_DIRS}/SpaceVecAlg)
include_directories(${RBD_INCLUDE_DIRS}/RBDyn)
include_directories(${SCD_INCLUDE_DIRS})
add_definitions("-DWITH_OPENGL")

add_library(Tasks SHARED ${SOURCES} ${HEADERS})
target_link_libraries(Tasks RBDyn scd)

set(INSTALL_PATH include/Tasks)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS Tasks DESTINATION "lib")
