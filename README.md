# Tasks

[![Build Status](https://travis-ci.org/jorisv/Tasks.svg?branch=master)](https://travis-ci.org/jorisv/Tasks)

The Tasks library is designed to make real-time control for kinematics tree and list of kinematics tree.

## Documentation

Features:
 * Support Kinematics Tree with Revolute/Prismatic/Spherical/Free/Planar/Cylindrical joints
 * Dynamic motion (motion must fulfill the equation of motion)
 * Contacts forces in friction cone
 * Static contacts
 * Articular position, speed and torque limits
 * Collision avoidance
 * Multi-robot (can solve a problem involving multiple robots in contact)
 * Tasks:
  * Posture target (articular position target, mandatory if you want avoid singularity issues)
  * Link Position/Orientation target
  * Link Velocity target
  * CoM target
  * Momentum target
  * Contact force target

To make sure that Tasks works as intended, unit tests are available for each algorithm.
Besides, the library has been used extensively to control humanoid robots such as HOAP-3, HRP-2, HRP-4 and Atlas.

The [SpaceVecAlg and RBDyn tutorial](https://github.com/jorisv/sva_rbdyn_tutorials) is also a big resources to understand how to use Tasks by providing a lot of IPython Notebook that will present real use case.

## Installing

### Manual

#### Dependencies

To compile you need the following tools:

 * [Git]()
 * [CMake]() >= 2.8
 * [pkg-config]()
 * [doxygen]()
 * [g++]() >= 4.7 (for C++11 support)
 * [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.49
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
 * [PyBindGen](https://launchpad.net/pybindgen) = 0.16
 * [Eigen3ToPython](https://github.com/jorisv/Eigen3ToPython) (to use the python binding)
 * [SpaceVecAlg](https://github.com/jorisv/SpaceVecAlg)
 * [RBDyn](https://github.com/jorisv/RBDyn)
 * [eigen-qld](https://github.com/jorisv/eigen-qld)
 * [eigen-lssol]() (Optional, if you have the LSSOL licence ask us this library)
 * [sch-core](https://github.com/jrl-umi3218/sch-core)

#### Building

```sh
git clone --recursive https://github.com/jorisv/RBDyn
cd RBDyn
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
