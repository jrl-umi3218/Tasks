#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 2.8)

include(cmake/base.cmake)
include(cmake/cython/cython.cmake)
include(cmake/msvc-specific.cmake)

set(PROJECT_NAME Tasks)
set(PROJECT_DESCRIPTION "...")
set(PROJECT_URL "https://github.com/jrl-umi3218/Tasks")
set(PROJECT_DEBUG_POSTFIX "_d")

# Disable -Werror on Unix for now.
set(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

setup_project(FALSE)
add_definitions(-mfpmath=sse -msse2)

option(DISABLE_TESTS "Disable unit tests." OFF)

if(NOT WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

add_project_dependency(RBDyn REQUIRED)
add_project_dependency(sch-core REQUIRED)
add_project_dependency(eigen-qld REQUIRED)
add_project_dependency(eigen-lssol QUIET)

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif(MSVC)

add_subdirectory(src)

if(NOT ${DISABLE_TESTS})
  add_subdirectory(tests)
endif()

if(${PYTHON_BINDING})
 add_subdirectory(binding/python)
endif()

# Add dependency towards the library in the pkg-config file.
pkg_config_append_libs(Tasks)
setup_project_finalize()
setup_project_package_finalize()
