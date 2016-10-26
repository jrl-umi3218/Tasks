# Copyright (C) 2008-2014 LAAS-CNRS, JRL AIST-CNRS.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

 # ------ #
 # README #
 # ------ #

# This file factorizes all rules to define a project for JRL or LAAS.
# It supposes that some variables have already been defined:
# - PROJECT_NAME	Project name.
#			Please keep respect our coding style and choose a name
#			which respects the following regexp: [a-z][a-z0-9-]*
#			I.e. a lower-case letter then one or more lower-case
#			letter, number or hyphen ``-''.
# - PROJECT_VERSION     Project version (X.Y.Z where X, Y, Z are unsigned
#                       integers). If not defined, it will automatically
#                       be computed through `git describe`.
#                       See BASE_COMPUTE_VERSION for more information.
# - PROJECT_DESCRIPTION One line summary of the package goal.
# - PROJECT_URL		Project's website.
#
# Please note that functions starting with an underscore are internal
# functions and should not be used directly.

 # ---- #
 # TODO #
 # ---- #

# - make install should trigger make doc
# - unit tests should be tagged as
#   EXCLUDE_FROM_ALL and make test should trigger their compilation.

# Include base features.
INCLUDE(cmake/logging.cmake)
INCLUDE(cmake/portability.cmake)
INCLUDE(cmake/compiler.cmake)
INCLUDE(cmake/debian.cmake)
INCLUDE(cmake/dist.cmake)
INCLUDE(cmake/distcheck.cmake)
INCLUDE(cmake/doxygen.cmake)
INCLUDE(cmake/header.cmake)
INCLUDE(cmake/pkg-config.cmake)
INCLUDE(cmake/uninstall.cmake)
INCLUDE(cmake/install-data.cmake)
INCLUDE(cmake/release.cmake)
INCLUDE(cmake/version.cmake)

 # --------- #
 # Constants #
 # --------- #

# Variables requires by SETUP_PROJECT.
SET(REQUIRED_VARIABLES PROJECT_NAME PROJECT_DESCRIPTION PROJECT_URL)

 # --------------------- #
 # Project configuration #
 # --------------------- #

# _ADD_TO_LIST LIST VALUE
# -------------
#
# Add a value to a comma-separated list.
#
# LIST		: the list.
# VALUE		: the value to be appended.
# SEPARATOR	: the separation symol.
#
MACRO(_ADD_TO_LIST LIST VALUE SEPARATOR)
  IF("${${LIST}}" STREQUAL "")
    SET(${LIST} "${VALUE}")
  ELSE("${${LIST}}" STREQUAL "")
    IF(NOT "${VALUE}" STREQUAL "")
      SET(${LIST} "${${LIST}}${SEPARATOR} ${VALUE}")
    ENDIF(NOT "${VALUE}" STREQUAL "")
  ENDIF("${${LIST}}" STREQUAL "")
ENDMACRO(_ADD_TO_LIST LIST VALUE)

# _CONCATE_ARGUMENTS
# -------------
#
# Concatenate all arguments into the output variable.
#
# OUTPUT	: the output variable.
# SEPARTOR	: the list separation symbol.
# ARG1...ARGN	: the values to be concatenated.
#
MACRO(_CONCATENATE_ARGUMENTS OUTPUT SEPARATOR)
  FOREACH(I RANGE 2 ${ARGC})
    _ADD_TO_LIST("${OUTPUT}" "${ARGV${I}}" "${SEPARATOR}")
  ENDFOREACH(I RANGE 2 ${ARGC})
  MESSAGE(${${OUTPUT}})
ENDMACRO(_CONCATENATE_ARGUMENTS OUTPUT)

# _SETUP_PROJECT_PACKAGE_INIT
# -------------
#
# Initialize the PackageConfig installation configuration.
# https://github.com/forexample/package-example
#
# This function does not take any argument but check that some
# variables are defined (see documentation at the beginning of this
# file).
#
MACRO(_SETUP_PROJECT_PACKAGE_INIT)
####
# Installation (https://github.com/forexample/package-example)

# Layout. This works for all platforms:
#   * <prefix>/lib/cmake/<PROJECT-NAME>
#   * <prefix>/lib/
#   * <prefix>/include/
set(config_install_dir "lib/cmake/${PROJECT_NAME}")
set(include_install_dir "include")
set(include_install_destination "${include_install_dir}/${PROJECT_NAME}")

set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")

# Configuration
set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
set(targets_export_name "${PROJECT_NAME}Targets")
set(namespace "${PROJECT_NAME}::")
ENDMACRO(_SETUP_PROJECT_PACKAGE_INIT)


# SETUP_PROJECT
# -------------
#
# Initialize the project. Should be called first in the root
# CMakeList.txt.
#
# This function does not take any argument but check that some
# variables are defined (see documentation at the beginning of this
# file).
#
MACRO(SETUP_PROJECT)
  INCLUDE(cmake/GNUInstallDirs.cmake)
  SET(CMAKE_INSTALL_FULL_PKGLIBDIR ${CMAKE_INSTALL_FULL_LIBDIR}/${PROJECT_NAME})
  SET(CMAKE_INSTALL_PKGLIBDIR ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME})

  # Check that required variables are defined.
  FOREACH(VARIABLE ${REQUIRED_VARIABLES})
    IF (NOT DEFINED ${VARIABLE})
      MESSAGE(FATAL_ERROR
	"Required variable ``${VARIABLE}'' has not been defined.")
    ENDIF(NOT DEFINED ${VARIABLE})
  ENDFOREACH(VARIABLE)

  # Define project name.
  PROJECT(${PROJECT_NAME} CXX)

  # If the project version number is not set, compute it automatically.
  IF(NOT DEFINED PROJECT_VERSION)
    VERSION_COMPUTE()
  ENDIF()

  IF(DEFINED PROJECT_DEBUG_POSTFIX)
    SET(CMAKE_DEBUG_POSTFIX ${PROJECT_DEBUG_POSTFIX})
    STRING(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type)
    IF(${cmake_build_type} MATCHES debug)
      SET(PKGCONFIG_POSTFIX ${PROJECT_DEBUG_POSTFIX})
    ELSE()
      SET(PKGCONFIG_POSTFIX "")
    ENDIF()
    IF(DEFINED CMAKE_CONFIGURATION_TYPES)
      SET(PKGCONFIG_POSTFIX ${PROJECT_DEBUG_POSTFIX})
    ENDIF()
  ENDIF()

  IF(${ARGC})
    SET(CMAKE_VERBOSE_MAKEFILE ${ARGV0})
  ELSE(${ARGC})
    # Be verbose by default.
    SET(CMAKE_VERBOSE_MAKEFILE TRUE)
  ENDIF(${ARGC})
  
  # If the variable INSTALL_DOCUMENTATION is not set, then the project will install the doc. Otherwise, 
  # the doc will be installed if INSTALL_DOCUMENTATION is set to TRUE.
  IF(DEFINED INSTALL_DOCUMENTATION)
    SET(_INSTALL_DOC ${INSTALL_DOCUMENTATION} CACHE INTERNAL "")
  ELSE(DEFINED INSTALL_DOCUMENTATION)
    SET(_INSTALL_DOC TRUE CACHE INTERNAL "")
  ENDIF(DEFINED INSTALL_DOCUMENTATION)

  ENABLE_TESTING()

  LOGGING_INITIALIZE()

  #FIXME: normalize naming to <MODULE>_SETUP()
  _SETUP_PROJECT_WARNINGS()
  _SETUP_PROJECT_HEADER()
  _SETUP_PROJECT_DIST()
  DISTCHECK_SETUP()
  RELEASE_SETUP()
  _SETUP_PROJECT_DEB()
  _SETUP_PROJECT_UNINSTALL()
  _SETUP_PROJECT_PKG_CONFIG()
  _SETUP_PROJECT_DOCUMENTATION()
  _SETUP_PROJECT_PACKAGE_INIT()
ENDMACRO(SETUP_PROJECT)

# _SETUP_PROJECT_PACKAGE_FINALIZE
# -------------
#
# Initialize the PackageConfig installation configuration.
# https://github.com/forexample/package-example
#
# This function does not take any argument but check that some
# variables are defined (see documentation at the beginning of this
# file).
#
# assumes SETUP_PROJECT() was called
# internally the real requirement is that _SETUP_PROJECT_PACKAGE_INIT() was called
MACRO(_SETUP_PROJECT_PACKAGE_FINALIZE)

####
# Installation (https://github.com/forexample/package-example)

# Layout. This works for all platforms:
#   * <prefix>/lib/cmake/<PROJECT-NAME>
#   * <prefix>/lib/
#   * <prefix>/include/
set(config_install_dir "lib/cmake/${PROJECT_NAME}")
set(include_install_dir "include")
set(include_install_destination "${include_install_dir}/${PROJECT_NAME}")

set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")

# Configuration
set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
set(targets_export_name "${PROJECT_NAME}Targets")
set(namespace "${PROJECT_NAME}::")

# Include module with fuction 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
# Note: PROJECT_VERSION is used as a VERSION
write_basic_package_version_file(
    "${version_config}" COMPATIBILITY SameMajorVersion
)

# Configure '<PROJECT-NAME>Config.cmake'
# Use variables:
#   * targets_export_name
#   * PROJECT_NAME
configure_package_config_file(
    "cmake/Config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${config_install_dir}"
)

# # Targets:
# #   * <prefix>/lib/libbar.a
# #   * <prefix>/lib/libbaz.a
# #   * header location after install: <prefix>/include/foo/Bar.hpp
# #   * headers can be included by C++ code `#include <foo/Bar.hpp>`
# install(
#     TARGETS ${SETUP_PROJECT_PACKAGE_TARGETS}
#     EXPORT "${targets_export_name}"
#     LIBRARY DESTINATION "lib"
#     ARCHIVE DESTINATION "lib"
#     RUNTIME DESTINATION "bin"
#     INCLUDES DESTINATION "${include_install_dir}"
# )

# # Headers:
# #   * Source/foo/Bar.hpp -> <prefix>/include/foo/Bar.hpp
# #   * Source/foo/Baz.hpp -> <prefix>/include/foo/Baz.hpp
# install(
#     DIRECTORY "Source/foo"
#     DESTINATION "${include_install_dir}"
#     FILES_MATCHING PATTERN "*.h*"
# )

# # Export headers:
# #   * ${CMAKE_CURRENT_BINARY_DIR}/bar_export.h -> <prefix>/include/bar_export.h
# #   * ${CMAKE_CURRENT_BINARY_DIR}/baz_export.h -> <prefix>/include/baz_export.h
# install(
#     FILES
#         "${CMAKE_CURRENT_BINARY_DIR}/bar_export.h"
#         "${CMAKE_CURRENT_BINARY_DIR}/baz_export.h"
#     DESTINATION "${include_install_dir}"
# )

# Config
#   * <prefix>/lib/cmake/Foo/FooConfig.cmake
#   * <prefix>/lib/cmake/Foo/FooConfigVersion.cmake
install(
    FILES "${project_config}" "${version_config}"
    DESTINATION "${config_install_dir}"
)

# Config
#   * <prefix>/lib/cmake/Foo/FooTargets.cmake
install(
    EXPORT "${targets_export_name}"
    NAMESPACE "${namespace}"
    DESTINATION "${config_install_dir}"
)
ENDMACRO(_SETUP_PROJECT_PACKAGE_FINALIZE)


# SETUP_PROJECT_FINALIZE
# ----------------------
#
# To be called at the end of the CMakeLists.txt to
# finalize the project setup.
#
MACRO(SETUP_PROJECT_FINALIZE)
  _SETUP_PROJECT_PKG_CONFIG_FINALIZE()
  _SETUP_PROJECT_DOCUMENTATION_FINALIZE()
  _SETUP_PROJECT_HEADER_FINAlIZE()
  _SETUP_PROJECT_PACKAGE_FINALIZE()
  _SETUP_DEBIAN()
  # Install data if needed
  _INSTALL_PROJECT_DATA()

  LOGGING_FINALIZE()
ENDMACRO(SETUP_PROJECT_FINALIZE)
