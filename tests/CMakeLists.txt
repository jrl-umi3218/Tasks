ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

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
    PKG_CONFIG_USE_DEPENDENCY(${name} EigenQP)
    add_test(${name}Unit ${name})
  endif()
endmacro(addUnitTest)

addUnitTest(QPSolverTest)
addUnitTest(QPMultiRobotTest)
addUnitTest(TasksTest)


if(${POSTURE_GENERATOR})
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
endif()
