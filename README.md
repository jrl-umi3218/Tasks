# Tasks

[![License LGPL 3](https://img.shields.io/badge/license-LGPLv3-green.svg)](http://www.gnu.org/licenses/lgpl-3.0.txt)
[![Build Status](https://travis-ci.org/jrl-umi3218/Tasks.svg?branch=master)](https://travis-ci.org/jrl-umi3218/Tasks)
[![AppVeyor status](https://ci.appveyor.com/api/projects/status/kteqpch13y0ac3wq/branch/master?svg=true)](https://ci.appveyor.com/project/gergondet/tasks/branch/master)

Tasks is library for real time control of robots and kinematic trees using constrained optimization.
It has been used extensively to control humanoid robots such as HOAP-3, HRP-2, HRP-4 and Atlas.

## Documentation

Features:
 * Support Kinematics Tree with Revolute/Prismatic/Spherical/Free/Planar/Cylindrical joints
 * Dynamic motion (motion must fulfill the equation of motion)
 * Contact forces in friction cones
 * Static contacts
 * Articular position, speed and torque limits
 * Collision avoidance
 * Multi-robot contact (can solve problems involving multiple robots in contact)
 * Tasks:
  * Posture target (articular position target, mandatory if you want avoid singularity issues)
  * Link Position/Orientation target
  * Link Velocity target
  * Center of Mass (CoM) target
  * Momentum target
  * Contact force target

To make sure that Tasks works as intended, unit tests are available for each algorithm.

The [SpaceVecAlg and RBDyn tutorial](https://github.com/jorisv/sva_rbdyn_tutorials) is also a big resources to understand how to use Tasks by providing a lot of IPython Notebook that will present real use case.

An online documentation can be found [online](https://jrl-umi3218.github.io/Tasks).

## Installing

### Ubuntu 14.04 and 16.04 binary ppa install

Use the [multi-contact-unstable](https://launchpad.net/~pierre-gergondet+ppa/+archive/ubuntu/multi-contact-unstable) ppa:
```bash
sudo add-apt-repository ppa:pierre-gergondet+ppa/multi-contact-unstable
sudo apt-get update
sudo apt-get install libtasks-dev libtasks-qld-doc
```

### Homebrew OS X install

Install from the command line using [Homebrew](brew.sh):

```bash
# install homebrew package manager
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
# install caskroom application manager
brew install caskroom/cask/brew-cask
# tap homebrew-science package repository
brew tap homebrew/science
# tap ahundt-robotics repository
brew tap ahundt/robotics
# install tasks and all its dependencies
brew install tasks
```

### Manually build from source

#### Dependencies

To compile you need the following tools:

 * [Git]()
 * [CMake]() >= 2.8
 * [pkg-config]()
 * [doxygen]()
 * [g++]() >= 4.7 (for C++11 support)
 * [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.49
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
 * [SpaceVecAlg](https://github.com/jrl-umi3218/SpaceVecAlg)
 * [RBDyn](https://github.com/jrl-umi3218/RBDyn)
 * [eigen-qld](https://github.com/jrl-umi3218/eigen-qld)
 * [eigen-lssol]() (Optional, if you have the LSSOL licence ask us this library)
 * [sch-core](https://github.com/jrl-umi3218/sch-core)

For Python bindings:
 * [Cython](http://cython.org/) = 0.25
 * [Eigen3ToPython](https://github.com/jrl-umi3218/Eigen3ToPython)
 * [sch-core-python](https://github.com/jrl-umi3218/sch-core-python)

#### Building

```sh
git clone --recursive https://github.com/jrl-umi3218/Tasks
cd Tasks
mkdir _build
cd _build
cmake [options] ..
make && make intall
```

Where the main options are:

 * `-DCMAKE_BUIlD_TYPE=Release` Build in Release mode
 * `-DCMAKE_INSTALL_PREFIX=some/path/to/install` default is `/usr/local`
 * `-DPYTHON_BINDING=ON` Build the python binding
 * `-DUNIT_TESTS=ON` Build unit tests.
 * `-DPYTHON_DEB_LAYOUT=OFF` install python library in `site-packages` (ON will install in `dist-packages`)


## Pulling git subtree

To update cmake or .travis directory with their upstream git repository:

```
git fetch git://github.com/jrl-umi3218/jrl-cmakemodules.git master
git subtree pull --prefix cmake git://github.com/jrl-umi3218/jrl-cmakemodules.git master --squash

git fetch git://github.com/jrl-umi3218/jrl-travis.git master
git subtree pull --prefix .travis git://github.com/jrl-umi3218/jrl-travis.git master --squash
```
