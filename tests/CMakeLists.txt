#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

enable_testing()

find_package(Boost REQUIRED COMPONENTS unit_test_framework timer system)

set(HEADERS arms.h)

macro(addUnitTest name)
  add_executable(${name} ${name}.cpp ${HEADERS})
  target_link_libraries(${name} PRIVATE Tasks Boost::unit_test_framework)
  add_test(${name}Unit ${name})
  # Adding a project configuration file (for MSVC only)
  generate_msvc_dot_user_file(${name})
endmacro(addUnitTest)

addUnitTest(QPSolverTest)
addUnitTest(QPMultiRobotTest)
addUnitTest(TasksTest)
addUnitTest(AllocationTest)
