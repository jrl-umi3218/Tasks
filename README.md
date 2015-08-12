# Tasks

[![Build Status](https://travis-ci.org/jorisv/Tasks.svg?branch=master)](https://travis-ci.org/jorisv/Tasks)

Tasks is a library tha allow to control a robot or to compute the inverse kinematics.

Pulling git subtree
------

To update sync cmake or .travis directory with their upstream git repository:

	git fetch git://github.com/jrl-umi3218/jrl-cmakemodules.git master
	git subtree pull --prefix cmake git://github.com/jrl-umi3218/jrl-cmakemodules.git master --squash


## Install

```sh
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX="/your/install/prefix" -DCMAKE_BUILD_TYPE=RelWithDebInfo
make
make install
```

Note: if you are on a Debian-based distribution (e.g. Ubuntu), you may want to
add the `-DPYTHON_DEB_LAYOUT` flag to the `cmake` command in order to install
this package with the specific Debian layout.
