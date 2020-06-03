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
    exports_sources = ["CMakeLists.txt", "conan/CMakeLists.txt", "binding/*", "cmake/*", "src/*"]
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
