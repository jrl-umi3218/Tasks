#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from eigen.c_eigen cimport *
from sva.c_sva cimport *
from rbdyn.c_rbdyn cimport *
from sch.c_sch import *
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "<boost/timer/timer.hpp>" namespace "boost::timer":
  cdef cppclass cpu_times:
    long long int wall
    long long int user
    long long int system

cdef extern from "<Tasks/Tasks.h>" namespace "tasks":
  cdef cppclass PositionTask:
    PositionTask(const MultiBody&, const string&, const Vector3d&, const Vector3d&)
    void position(const Vector3d&)
    Vector3d position() const
    void bodyPoint(const Vector3d&)
    Vector3d bodyPoint() const
    # Common to *Task
    void update(const MultiBody&, const MultiBodyConfig&)
    void updateDot(const MultiBody&, const MultiBodyConfig&)
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass OrientationTask:
    OrientationTask(const MultiBody&, const string&, const Quaterniond&)
    OrientationTask(const MultiBody&, const string&, const Matrix3d&)
    void orientation(const Quaterniond&)
    void orientation(const Matrix3d&)
    Matrix3d orientation() const
    # Common to *Task
    void update(const MultiBody&, const MultiBodyConfig&)
    void updateDot(const MultiBody&, const MultiBodyConfig&)
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass SurfaceOrientationTask:
    SurfaceOrientationTask(const MultiBody&, const string&, const Quaterniond&, const PTransformd&)
    SurfaceOrientationTask(const MultiBody&, const string&, const Matrix3d&, const PTransformd&)
    void orientation(const Quaterniond&)
    void orientation(const Matrix3d&)
    Matrix3d orientation() const
    # Common to *Task
    void update(const MultiBody&, const MultiBodyConfig&)
    void updateDot(const MultiBody&, const MultiBodyConfig&)
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass GazeTask:
    GazeTask(const MultiBody&, const string&, const Vector2d&, double, const PTransformd&, const Vector2d&)
    GazeTask(const MultiBody&, const string&, const Vector3d&, const PTransformd&, const Vector2d&)
    void error(const Vector2d&, const Vector2d&)
    void error(const Vector3d&, const Vector2d&)
    VectorXd speed() const
    VectorXd normalAcc() const
    void update(const MultiBody&, const MultiBodyConfig&, const vector[MotionVecd]&)
    # Common to *Task
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass  PositionBasedVisServoTask:
    PositionBasedVisServoTask(const MultiBody&, const string&, const PTransformd&, const PTransformd&)
    void error(const PTransformd&)
    void update(const MultiBody&, const MultiBodyConfig&, const vector[MotionVecd]&)
    VectorXd speed() const
    VectorXd normalAcc() const
    # Common to *Task
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass PostureTask:
    PostureTask(const MultiBody&, vector[vector[double]])
    void posture(vector[vector[double]])
    vector[vector[double]] posture()
    # Common to *Task
    void update(const MultiBody&, const MultiBodyConfig&)
    void updateDot(const MultiBody&, const MultiBodyConfig&)
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass CoMTask:
    CoMTask(const MultiBody&, const Vector3d&)
    CoMTask(const MultiBody&, const Vector3d&, vector[double]) except +
    void com(const Vector3d&)
    Vector3d com() const
    void updateInertialParameters(const MultiBody&)
    # Common to *Task
    void update(const MultiBody&, const MultiBodyConfig&)
    void updateDot(const MultiBody&, const MultiBodyConfig&)
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass MultiCoMTask:
    MultiCoMTask(const vector[MultiBody]&, vector[int], const Vector3d& com)
    vector[int] robotIndexes() const
    void com(const Vector3d&)
    Vector3d com() const
    void updateInertialParameters(const vector[MultiBody]&)
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&)
    VectorXd eval() const
    VectorXd speed() const
    VectorXd normalAcc() const
    MatrixXd jac(int) const

  cdef cppclass MultiRobotTransformTask:
    MultiRobotTransformTask(const vector[MultiBody]&, int, int, const string&, const string&, const PTransformd&, const PTransformd&)
    int r1Index() const
    int r2Index() const
    void X_r1b_r1s(const PTransformd&)
    PTransformd X_r1b_r1s()
    void X_r2b_r2s(const PTransformd&)
    PTransformd X_r2b_r2s()
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const vector[vector[MotionVecd]]&)
    VectorXd eval() const
    VectorXd speed() const
    VectorXd normalAcc() const
    MatrixXd jac(int) const

  cdef cppclass MomentumTask:
    MomentumTask(const MultiBody&, const ForceVecd&)
    void momentum(const ForceVecd&)
    ForceVecd momentum()
    # Common to *Task
    void update(const MultiBody&, const MultiBodyConfig&)
    void updateDot(const MultiBody&, const MultiBodyConfig&)
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass LinVelocityTask:
    LinVelocityTask(const MultiBody&, const string&, const Vector3d&, const Vector3d&)
    void velocity(const Vector3d&)
    Vector3d velocity()
    void bodyPoint(const Vector3d&)
    Vector3d bodyPoint()
    # Common to *Task
    void update(const MultiBody&, const MultiBodyConfig&)
    void updateDot(const MultiBody&, const MultiBodyConfig&)
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass OrientationTrackingTask:
    OrientationTrackingTask(const MultiBody&, const string&, const Vector3d&, const Vector3d&, const vector[string]&, const Vector3d&)
    void trackedPoint(const Vector3d&)
    Vector3d trackedPoint()
    void bodyPoint(const Vector3d&)
    Vector3d bodyPoint()
    void bodyAxis(const Vector3d&)
    Vector3d bodyAxis()
    # Common to *Task
    void update(const MultiBody&, const MultiBodyConfig&)
    void updateDot(const MultiBody&, const MultiBodyConfig&)
    VectorXd eval() const
    MatrixXd jac() const
    MatrixXd jacDot() const

  cdef cppclass TransformTask:
    TransformTask(const MultiBody&, const string&, const PTransformd&, const PTransformd&, const Matrix3d&)
    void E_0_c(const Matrix3d&)
    Matrix3d E_0_c() const
    # Common *TransformTask
    void target(const PTransformd&)
    PTransformd target() const
    void X_b_p(const PTransformd&)
    PTransformd X_b_p() const
    void update(const MultiBody&, const MultiBodyConfig&, const vector[MotionVecd]&)
    VectorXd eval() const
    MatrixXd jac() const
    VectorXd speed() const
    VectorXd normalAcc() const

  cdef cppclass SurfaceTransformTask:
    SurfaceTransformTask(const MultiBody&, const string&, const PTransformd&, const PTransformd&)
    # Common *TransformTask
    void target(const PTransformd&)
    PTransformd target() const
    void X_b_p(const PTransformd&)
    PTransformd X_b_p() const
    void update(const MultiBody&, const MultiBodyConfig&, const vector[MotionVecd]&)
    VectorXd eval() const
    MatrixXd jac() const
    VectorXd speed() const
    VectorXd normalAcc() const

cdef extern from "<Tasks/Bounds.h>" namespace "tasks":
  cdef cppclass QBound:
    QBound()
    QBound(vector[vector[double]], vector[vector[double]])
    vector[vector[double]] lQBound
    vector[vector[double]] uQBound

  cdef cppclass AlphaBound:
    AlphaBound()
    AlphaBound(vector[vector[double]], vector[vector[double]])
    vector[vector[double]] lAlphaBound
    vector[vector[double]] uAlphaBound

  cdef cppclass TorqueBound:
    TorqueBound()
    TorqueBound(vector[vector[double]], vector[vector[double]])
    vector[vector[double]] lTorqueBound
    vector[vector[double]] uTorqueBound

  cdef cppclass PolyTorqueBound:
    PolyTorqueBound()
    PolyTorqueBound(vector[vector[VectorXd]], vector[vector[VectorXd]])
    vector[vector[VectorXd]] lPolyTorqueBound
    vector[vector[VectorXd]] uPolyTorqueBound

