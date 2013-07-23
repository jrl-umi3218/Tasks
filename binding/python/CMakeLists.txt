INCLUDE(../../cmake/python.cmake)
FINDPYTHON()

## Define PYTHON_DISTLIB
EXECUTE_PROCESS(
  COMMAND "${PYTHON_EXECUTABLE}" "-c"
  "import sys, os; print os.sep.join(['lib', 'python' + sys.version[:3], 'dist-packages'])"
  OUTPUT_VARIABLE PYTHON_DISTLIB
  ERROR_QUIET)
# Remove final \n of the variable PYTHON_DISTLIB
STRING(REPLACE "\n" "" PYTHON_DISTLIB "${PYTHON_DISTLIB}")


# Tasks

# create the package in build dir for testing purpose
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/tasks)
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/tasks/__init__.py COPYONLY)

set(OUTPUT_BINDING ${CMAKE_CURRENT_BINARY_DIR}/tasks.cpp)

# generate python binding code
add_custom_command (
  OUTPUT ${OUTPUT_BINDING}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate.py ${OUTPUT_BINDING}
  DEPENDS generate.py
)

# build the library
set(SOURCES ${OUTPUT_BINDING})
include_directories(.)
include_directories(../../src)
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(_tasks SHARED ${SOURCES})
PKG_CONFIG_USE_DEPENDENCY(_tasks scd) # usefull because of WITH_OPENGL flags
PKG_CONFIG_USE_DEPENDENCY(_tasks SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(_tasks RBDyn)
PKG_CONFIG_USE_DEPENDENCY(_tasks EigenQP)

target_link_libraries(_tasks Tasks)
set_target_properties(_tasks PROPERTIES PREFIX "")



# SCD
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/scd)
configure_file(__scd_init__.py ${CMAKE_CURRENT_BINARY_DIR}/scd/__init__.py COPYONLY)

set(OUTPUT_BINDING ${CMAKE_CURRENT_BINARY_DIR}/scd.cpp)

# generate python binding code
add_custom_command (
  OUTPUT ${OUTPUT_BINDING}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate_scd.py ${OUTPUT_BINDING}
  DEPENDS generate_scd.py
)

# build the library
set(SOURCES ${OUTPUT_BINDING})
include_directories(.)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${SVA_INCLUDE_DIRS}/SpaceVecAlg)
include_directories(${SCD_INCLUDE_DIRS})

add_library(_scd SHARED ${SOURCES})
PKG_CONFIG_USE_DEPENDENCY(_scd scd)
set_target_properties(_scd PROPERTIES PREFIX "")


# install rules
set(INSTALL_PATH "${PYTHON_DISTLIB}/tasks/")
install(TARGETS _tasks DESTINATION ${INSTALL_PATH})
install(FILES __init__.py DESTINATION ${INSTALL_PATH})

set(INSTALL_PATH "${PYTHON_DISTLIB}/scd/")
install(TARGETS _scd DESTINATION ${INSTALL_PATH})
install(FILES __scd_init__.py DESTINATION ${INSTALL_PATH} RENAME __init__.py)

