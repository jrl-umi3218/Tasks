# -*- coding: utf-8 -*-
#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from conans import python_requires, tools
import os

base = python_requires("Eigen3ToPython/latest@multi-contact/dev")

class TasksConan(base.Eigen3ToPythonConan):
    name = "Tasks"
    version = "1.2.0"
    description = "Real time control of robots using constrained optimization"
    topics = ("robotics", "control", "optimization", "python")
    url = "https://github.com/jrl-umi3218/Tasks"
    homepage = "https://github.com/jrl-umi3218/Tasks"
    author = "Pierre Gergondet <pierre.gergondet@gmail.com>"
    license = "BSD-2-Clause"
    exports = ["LICENSE"]
    exports_sources = ["CMakeLists.txt", "conan/CMakeLists.txt", "conan/FindBoost.cmake", "binding/*", "cmake/*", "src/*"]
    generators = ["cmake_find_package", "cmake_paths"]
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
        replacement = '''{}
install(FILES conan/FindBoost.cmake DESTINATION lib/cmake/Tasks)'''.format(pattern)
        tools.replace_in_file('CMakeListsOriginal.txt', pattern, replacement)
        # Link with Boost::Boost if consumed by conan
        pattern = "Boost::timer Boost::disable_autolinking Boost::dynamic_linking"
        replacement = "$<BUILD_INTERFACE:Boost::Boost>$<INSTALL_INTERFACE:$<IF:$<BOOL:CONAN_BOOST_ROOT>,Boost::Boost,{}>>".format(pattern.replace(' ', '$<SEMICOLON>'))
        tools.replace_in_file('src/CMakeLists.txt', pattern, replacement)

    def package(self):
        cmake = self._configure_cmake()
        pattern = "BOOL:CONAN_BOOST_ROOT"
        replacement = "BOOL:${CONAN_BOOST_ROOT}"
        generated = os.path.join(self.build_folder, 'CMakeFiles', 'Export', 'lib', 'cmake', self.name, '{}Targets.cmake'.format(self.name))
        tools.replace_in_file(generated, pattern, replacement)
        cmake.install()
