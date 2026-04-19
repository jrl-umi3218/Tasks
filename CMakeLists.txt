#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 3.14)

set(PROJECT_NAME Tasks)
set(PROJECT_DESCRIPTION "...")
set(PROJECT_URL "https://github.com/jrl-umi3218/Tasks")
set(PROJECT_DEBUG_POSTFIX "_d")
set(PROJECT_VERSION 1.8.3)
set(PROJECT_USE_CMAKE_EXPORT TRUE)

set(CMAKE_CXX_STANDARD 17)
# Disable -Werror on Unix for now.
set(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

# Check if the submodule cmake have been initialized
set(JRL_CMAKE_MODULES "${CMAKE_CURRENT_LIST_DIR}/cmake")
if(EXISTS "${JRL_CMAKE_MODULES}/base.cmake")
  message(STATUS "JRL cmakemodules found in 'cmake/' git submodule")
else()
  find_package(jrl-cmakemodules QUIET CONFIG)
  if(jrl-cmakemodules_FOUND)
    get_property(
      JRL_CMAKE_MODULES
      TARGET jrl-cmakemodules::jrl-cmakemodules
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    )
    message(STATUS "JRL cmakemodules found on system at ${JRL_CMAKE_MODULES}")
  else()
    message(STATUS "JRL cmakemodules not found. Let's fetch it.")
    include(FetchContent)
    FetchContent_Declare(
      "jrl-cmakemodules"
      GIT_REPOSITORY "https://github.com/jrl-umi3218/jrl-cmakemodules.git"
    )
    FetchContent_MakeAvailable("jrl-cmakemodules")
    FetchContent_GetProperties("jrl-cmakemodules" SOURCE_DIR JRL_CMAKE_MODULES)
  endif()
endif()

include(${JRL_CMAKE_MODULES}/base.cmake)
include(${JRL_CMAKE_MODULES}/cython/cython.cmake)
include(${JRL_CMAKE_MODULES}/msvc-specific.cmake)

project(Tasks CXX)

if(UNIX
   AND NOT EMSCRIPTEN
   AND NOT (APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm")
)
  include(CheckCXXSourceCompiles)
  check_cxx_compiler_flag("-msse2" SSE2_SUPPORTED)
  if(SSE2_SUPPORTED)
    add_definitions(-mfpmath=sse -msse2)
  endif()
endif()

if(NOT WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

add_project_dependency(RBDyn REQUIRED NO_MODULE)
add_project_dependency(sch-core REQUIRED NO_MODULE)
add_project_dependency(eigen-qld REQUIRED NO_MODULE)
find_package(eigen-lssol QUIET NO_MODULE)
if(${eigen-lssol_FOUND})
  add_project_dependency(eigen-lssol REQUIRED NO_MODULE)
endif()

# For MSVC, set local environment variable to enable finding the built dll of the main
# library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif(MSVC)

add_subdirectory(src)

if(${BUILD_TESTING})
  add_subdirectory(tests)
endif()

if(${PYTHON_BINDING})
  add_subdirectory(binding/python)
endif()
