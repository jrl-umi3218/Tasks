ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})
add_definitions("-DWITH_OPENGL")

add_executable(QPSolverTest QPSolverTest.cpp)
target_link_libraries(QPSolverTest ${Boost_LIBRARIES} RBDyn Tasks)
PKG_CONFIG_USE_DEPENDENCY(QPSolverTest scd)
PKG_CONFIG_USE_DEPENDENCY(QPSolverTest SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(QPSolverTest RBDyn)
add_test(QPSolverTestUnit QPSolverTest)

add_executable(PGJacobianTest PGJacobianTest.cpp)
target_link_libraries(PGJacobianTest ${Boost_LIBRARIES} RBDyn Tasks)
PKG_CONFIG_USE_DEPENDENCY(PGJacobianTest SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(PGJacobianTest RBDyn)
add_test(PGJacobiantUnit PGJacobianTest)

add_executable(PGTest PGTest.cpp)
target_link_libraries(PGTest ${Boost_LIBRARIES} RBDyn Tasks)
PKG_CONFIG_USE_DEPENDENCY(PGTest SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(PGTest RBDyn)
add_test(PGUnit PGTest)

