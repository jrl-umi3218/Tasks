cmake_minimum_required(VERSION 2.8)

message(STATUS "Conan FlatBuffers Wrapper")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/conan/CMakeModules)
include(${CMAKE_SOURCE_DIR}/CMakeListsOriginal.txt)
