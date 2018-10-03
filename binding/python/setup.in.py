# Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
#
# This file is part of Tasks.
#
# Tasks is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Tasks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Tasks.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
try:
  from setuptools import setup
  from setuptools import Extension
except ImportError:
  from distutils.core import setup
  from distutils.extension import Extension

from Cython.Build import cythonize

import filecmp
import hashlib
import os
import shutil
import subprocess

from numpy import get_include as numpy_get_include

win32_build = os.name == 'nt'
debug_build = "$<CONFIGURATION>".lower() == "debug"

def exists_or_create(path):
  if not os.path.exists(path):
    os.makedirs(path)

def copy_if_different(f):
  in_ = '@CMAKE_CURRENT_SOURCE_DIR@/{}'.format(f)
  out_ = '{}/{}'.format(this_path, f)
  if not os.path.exists(out_) or not filecmp.cmp(in_, out_):
    shutil.copyfile(in_, out_)

this_path  = os.path.dirname(os.path.realpath(__file__))
exists_or_create(this_path + '/tasks/qp/')
exists_or_create(this_path + '/tests/')
with open(this_path + '/tasks/__init__.py', 'w') as fd:
    fd.write('from . tasks import *\n')
    fd.write('from . import qp\n')
with open(this_path + '/tasks/qp/__init__.py', 'w') as fd:
    fd.write('from .qp import *\n')
[copy_if_different('tasks/{}'.format(f)) for f in ['c_tasks.pxd', 'tasks.pyx', 'tasks.pxd']]
[copy_if_different('tasks/qp/{}'.format(f)) for f in ['c_qp_private.pxd', 'c_qp.pxd', 'qp.pyx', 'qp.pxd']]
[copy_if_different('tests/{}'.format(f)) for f in ['arms.py', 'TestQPMultiRobot.py', 'TestQPSolver.py', 'utils.py']]

sha512 = hashlib.sha512()
src_files = ['tasks/tasks.pyx', 'tasks/c_tasks.pxd', 'tasks/tasks.pxd']
src_files.extend( ['tasks/qp/qp.pyx', 'tasks/qp/c_qp.pxd', 'tasks/qp/qp.pxd', 'include/qp_wrapper.hpp'] )
src_files = [ '{}/{}'.format('@CMAKE_CURRENT_SOURCE_DIR@', f) for f in src_files ]
for f in src_files:
  chunk = 2**16
  with open(f, 'r') as fd:
    while True:
      data = fd.read(chunk)
      if data:
        sha512.update(data.encode('ascii'))
      else:
        break
version_hash = sha512.hexdigest()[:7]

class pkg_config(object):
  def __init__(self):
    self.compile_args = []
    self.include_dirs = [ x for x in '$<TARGET_PROPERTY:Tasks,INCLUDE_DIRECTORIES>'.split(';') if len(x) ]
    self.include_dirs.append('@CMAKE_CURRENT_SOURCE_DIR@/include')
    self.library_dirs = [ x for x in '$<TARGET_PROPERTY:Tasks,LINK_FLAGS>'.split(';') if len(x) ]
    if debug_build:
      self.libraries = ['Tasks_d']
    else:
      self.libraries = ['Tasks']
    tasks_location = '$<TARGET_FILE:Tasks>'
    self.library_dirs.append(os.path.dirname(tasks_location))
    self.found = True

config = pkg_config()

config.compile_args.append('-std=c++11')
config.include_dirs.append(os.getcwd() + "/include")
if win32_build:
  config.compile_args.append("-DWIN32")

def GenExtension(name, pkg, ):
  pyx_src = name.replace('.', '/')
  pyx_src = pyx_src + '.pyx'
  ext_src = pyx_src
  if pkg.found:
    return Extension(name, [ext_src], extra_compile_args = pkg.compile_args, include_dirs = pkg.include_dirs + [numpy_get_include()], library_dirs = pkg.library_dirs, libraries = pkg.libraries)
  else:
    print("Failed to find {}".format(pkg.name))
    return None

extensions = [
  GenExtension('tasks.tasks', config),
  GenExtension('tasks.qp.qp', config)
]

extensions = [ x for x in extensions if x is not None ]
packages = ['tasks']
data = ['__init__.py', 'c_tasks.pxd', 'tasks.pxd', 'qp/__init__.py', 'qp/c_qp.pxd', 'qp/qp.pxd']

extensions = cythonize(extensions)

setup(
    name = 'tasks',
    version='1.0.0-{}'.format(version_hash),
    ext_modules = extensions,
    packages = packages,
    package_data = { 'tasks': data }
)
