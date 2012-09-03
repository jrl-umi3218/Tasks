set(SOURCES Tasks.cpp)
set(HEADERS Tasks.h)

pkg_check_modules(EIGEN3 REQUIRED eigen3)
find_path(SVA_INCLUDE_DIRS SpaceVecAlg/SpaceVecAlg
                           /usr/include
                           /usr/local/include)

find_path(RBD_INCLUDE_DIRS RBDyn/
                           /usr/include
                           /usr/local/include)

include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${SVA_INCLUDE_DIRS}/SpaceVecAlg)
include_directories(${RBD_INCLUDE_DIRS}/RBDyn)

add_library(Tasks SHARED ${SOURCES} ${HEADERS})

set(INSTALL_PATH include/Tasks)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS Tasks DESTINATION "lib")
