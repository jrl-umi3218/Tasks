cmake_minimum_required(VERSION 2.8)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME Tasks)
SET(PROJECT_DESCRIPTION "...")
SET(PROJECT_URL "https://github.com/jorisv//Tasks")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)

project(${PROJECT_NAME} CXX)
SETUP_PROJECT(FALSE)

option(PYTHON_BINDING "Generate python binding." ON)
option(POSTURE_GENERATOR "Build the non linear posture generator." OFF)
option(UNIT_TESTS "Generate unit tests." ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic")

#########################
# External dependencies
SEARCH_FOR_EIGEN()

ADD_REQUIRED_DEPENDENCY("sch-core")
ADD_REQUIRED_DEPENDENCY("SpaceVecAlg")
ADD_REQUIRED_DEPENDENCY("RBDyn")
ADD_REQUIRED_DEPENDENCY("eigen-qld")
ADD_OPTIONAL_DEPENDENCY("eigen-lssol")

if(${POSTURE_GENERATOR})
  ADD_REQUIRED_DEPENDENCY("ipopt")
endif()

#########################
# 

add_subdirectory(src)
add_subdirectory(tests)

if(${PYTHON_BINDING})
 add_subdirectory(binding/python)
endif()

# Add dependency towards the library in the pkg-config file.
PKG_CONFIG_APPEND_LIBS(Tasks)
SETUP_PROJECT_FINALIZE()
