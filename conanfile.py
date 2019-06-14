# -*- coding: utf-8 -*-
#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from conans import python_requires, tools
import os

base = python_requires("RBDyn/1.1.0@gergondet/stable")


class TasksConan(base.RBDynConan):
    name = "Tasks"
    version = "0.9.0"
    description = "Real time control of robots using constrained optimization"
    topics = ("robotics", "control", "optimization", "python")
    url = "https://github.com/jrl-umi3218/Tasks"
    homepage = "https://github.com/jrl-umi3218/Tasks"
    author = "Pierre Gergondet <pierre.gergondet@gmail.com>"
    license = "BSD-2-Clause"
    exports = ["LICENSE"]
    exports_sources = ["CMakeLists.txt", "conan/CMakeLists.txt", "binding/*", "cmake/*", "conan/CMakeModules/*", "doc/*", "src/*"]
    generators = "cmake"
    settings = "os", "arch", "compiler", "build_type"
    options = { "python_version": ["2.7", "3.3", "3.4", "3.5", "3.6", "3.7"] }
    default_options = { "python_version": base.get_python_version() }

    requires = (
        "RBDyn/1.1.0@gergondet/stable",
        "sch-core-python/1.0.0@gergondet/stable",
        "eigen-qld/1.0.0@gergondet/stable"
    )
