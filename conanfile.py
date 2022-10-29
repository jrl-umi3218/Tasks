# -*- coding: utf-8 -*-
#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from conans import python_requires, tools
import os

base = python_requires("Eigen3ToPython/latest@multi-contact/dev")

class TasksConan(base.Eigen3ToPythonConan):
    name = "Tasks"
    version = "1.7.1"
    description = "Real time control of robots using constrained optimization"
    topics = ("robotics", "control", "optimization", "python")
    url = "https://github.com/jrl-umi3218/Tasks"
    homepage = "https://github.com/jrl-umi3218/Tasks"
    author = "Pierre Gergondet <pierre.gergondet@gmail.com>"
    license = "BSD-2-Clause"
    exports = ["LICENSE"]
    exports_sources = ["CMakeLists.txt", "conan/CMakeLists.txt", "conan/FindBoost.cmake", "binding/*", "cmake/*", "src/*"]
    generators = ["cmake_paths"]
    settings = "os", "arch", "compiler", "build_type"

    requires = (
        "RBDyn/latest@multi-contact/dev",
        "sch-core-python/latest@multi-contact/dev",
        "eigen-qld/latest@multi-contact/dev"
    )

    def package_info(self):
        self.cpp_info.libs = tools.collect_libs(self)
        self.env_info.LD_LIBRARY_PATH.append(os.path.join(self.package_folder, 'lib'))
        self.env_info.PYTHONPATH.append(self._extra_python_path())

    def source(self):
        base.Eigen3ToPythonConan.source(self)
        # Make sure we find conan's Boost not system Boost
        pattern = 'include(CMakeFindDependencyMacro)'
        replacement = '''set(BOOST_ROOT "${{PACKAGE_PREFIX_DIR}}")
set(Boost_NO_SYSTEM_PATHS ON)
list(APPEND CMAKE_MODULE_PATH "${{CMAKE_CURRENT_LIST_DIR}}")
{}'''.format(pattern)
        tools.replace_in_file('cmake/Config.cmake.in', pattern, replacement)
        # Install the up-to-date FindBoost.cmake
        pattern = 'add_subdirectory(src)'
        replacement = '''list(APPEND CMAKE_MODULE_PATH "${{CMAKE_CURRENT_LIST_DIR}}/conan")
{}
install(FILES conan/FindBoost.cmake DESTINATION lib/cmake/RBDyn)'''.format(pattern)
        tools.replace_in_file('CMakeListsOriginal.txt', pattern, replacement)

    def package(self):
        cmake = self._configure_cmake()
        cmake.definitions['INSTALL_DOCUMENTATION'] = False
        cmake.configure()
        cmake.install()
