#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from eigen.c_eigen cimport *
from sva.c_sva cimport *
from rbdyn.c_rbdyn cimport *
from sch.c_sch cimport *
cimport tasks.c_tasks as c_tasks
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "<Tasks/QPContacts.h>" namespace "tasks::qp":
  cdef cppclass FrictionCone:
    FrictionCone()
    FrictionCone(Matrix3d, int, double, double)
    vector[Vector3d] generators

  cdef cppclass ContactId:
    ContactId()
    ContactId(int, int, const string&, const string&, int)
    int r1Index
    int r2Index
    string r1BodyName
    string r2BodyName
    int ambiguityId
    bool operator==(const ContactId&)
    bool operator!=(const ContactId&)
    bool operator<(const ContactId&)

  cdef cppclass UnilateralContact:
    UnilateralContact()
    UnilateralContact(int, int, const string&, const string&, const vector[Vector3d]&, const Matrix3d&, const PTransformd&, int, double, const PTransformd&)
    UnilateralContact(int, int, const string&, const string&, int, const vector[Vector3d]&, const Matrix3d&, const PTransformd&, int, double, const PTransformd&)
    UnilateralContact(const ContactId&, const vector[Vector3d]&, const Matrix3d&, const PTransformd&, int, double, const PTransformd&)
    ContactId contactId
    vector[Vector3d] r1Points
    vector[Vector3d] r2Points
    FrictionCone r1Cone
    FrictionCone r2Cone
    PTransformd X_b1_b2
    PTransformd X_b1_cf
    Vector3d sForce(const VectorXd&, int, const FrictionCone&) except +
    Vector3d sForce(const VectorXd&, const FrictionCone&) except +
    ForceVecd sForce(const VectorXd&, const vector[Vector3d]&, const FrictionCone&) except +
    int sNrLambda(int) except +
    int nrLambda() except +

  cdef cppclass BilateralContact:
    BilateralContact()
    BilateralContact(int, int, const string&, const string&, const vector[Vector3d]&, const vector[Matrix3d]&, const PTransformd&, int, double, const PTransformd&)
    BilateralContact(int, int, const string&, const string&, int, const vector[Vector3d]&, const vector[Matrix3d]&, const PTransformd&, int, double, const PTransformd&)
    BilateralContact(const ContactId&, const vector[Vector3d]&, const vector[Matrix3d]&, const PTransformd&, int, double, const PTransformd&)
    ContactId contactId
    vector[Vector3d] r1Points
    vector[Vector3d] r2Points
    vector[FrictionCone] r1Cones
    vector[FrictionCone] r2Cones
    PTransformd X_b1_b2
    PTransformd X_b1_cf
    Vector3d sForce(const VectorXd&, int, const vector[FrictionCone]&) except +
    Vector3d sForce(const VectorXd&, const vector[FrictionCone]&) except +
    ForceVecd sForce(const VectorXd&, const vector[Vector3d]&, const vector[FrictionCone]&) except +
    int sNrLambda(int) except +
    int nrLambda() except +

  cdef cppclass SolverData:
    int nrVars() const
    int totalAlphaD() const
    int totalLambda() const
    int alphaD(int) const
    int _lambda "lambda"(int) const
    int alphaDBegin() const
    int alphaDBegin(int) const
    int lambdaBegin() const
    int lambdaBegin(int) const
    int nrUniLambda() const
    int nrBiLambda() const
    int unilateralBegin() const
    int bilateralBegin() const
    int nrContacts() const
    vector[UnilateralContact] unilateralContacts() const
    vector[BilateralContact] bilateralContacts() const
    vector[BilateralContact] allContacts() const
    void computeNormalAccB(const vector[MultiBody]&, const vector[MultiBodyConfig]&)
    vector[MotionVecd] normalAccB(int) const

cdef extern from "<Tasks/QPTasks.h>" namespace "tasks::qp":
  cdef cppclass JointStiffness:
    JointStiffness()
    JointStiffness(const string&, double)
    string jointName
    double stiffness

  cdef cppclass JointGains:
    JointGains()
    JointGains(const string&, double)
    JointGains(const string&, double, double)
    string jointName
    double stiffness
    double damping

cdef extern from "<Tasks/QPMotionConstr.h>" namespace "tasks::qp":
  cdef cppclass SpringJoint:
    SpringJoint()
    SpringJoint(const string&, double, double, double)
    string jointName
    double K
    double C
    double O

cdef extern from "<Tasks/QPSolver.h>" namespace "tasks::qp":
  cdef cppclass Constraint:
    void updateNrVars(const vector[MultiBody]&, SolverData)
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)

  cdef cppclass Equality:
    int maxEq()
    int nrEq()
    MatrixXd AEq()
    VectorXd bEq()

  cdef cppclass Inequality:
    int maxInEq()
    int nrInEq()
    MatrixXd AInEq()
    VectorXd bInEq()

  cdef cppclass GenInequality:
    int maxGenInEq()
    int nrGenInEq()
    MatrixXd AGenInEq()
    VectorXd LowerGenInEq()
    VectorXd UpperGenInEq()

  cdef cppclass Bound:
    int beginVar()
    MatrixXd Lower()
    VectorXd Upper()

  cdef cppclass ConstraintFunction[T](Constraint):
    pass

  cdef cppclass Task:
    double weight() const
    void weight(double)

  cdef cppclass HighLevelTask:
    int dim()
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)
    MatrixXd jac()
    VectorXd eval()
    VectorXd speed()
    VectorXd normalAcc()

cdef extern from "<Tasks/QPTasks.h>" namespace "tasks::qp":
  cdef cppclass PositionTask(HighLevelTask):
    PositionTask(const vector[MultiBody]&, int, const string&, const Vector3d&, const Vector3d&)
    void position(const Vector3d&)
    Vector3d position() const
    void bodyPoint(const Vector3d&)
    Vector3d bodyPoint() const

  cdef cppclass VectorOrientationTask(HighLevelTask):
    pass

  cdef cppclass OrientationTask(HighLevelTask):
    OrientationTask(const vector[MultiBody]&, int, const string&, const Quaterniond&)
    OrientationTask(const vector[MultiBody]&, int, const string&, const Matrix3d&)
    void orientation(const Matrix3d&)
    void orientation(const Quaterniond&)
    Matrix3d orientation() const

  cdef cppclass TransformTask(HighLevelTask):
    TransformTask(const vector[MultiBody]&, int, const string&, const PTransformd&, const PTransformd&, const Matrix3d&)
    void E_0_c(const Matrix3d&)
    Matrix3d E_0_c() const
    # Common to *TransformTask
    void target(const PTransformd&)
    PTransformd target() const
    void X_b_p(const PTransformd&)
    PTransformd X_b_p() const

  cdef cppclass SurfaceTransformTask(HighLevelTask):
    SurfaceTransformTask(const vector[MultiBody]&, int, const string&, const PTransformd&, const PTransformd&)
    # Common to *TransformTask
    void target(const PTransformd&)
    PTransformd target() const
    void X_b_p(const PTransformd&)
    PTransformd X_b_p() const

  cdef cppclass SurfaceOrientationTask(HighLevelTask):
    SurfaceOrientationTask(const vector[MultiBody]&, int, const string&, const Quaterniond&, const PTransformd&)
    SurfaceOrientationTask(const vector[MultiBody]&, int, const string&, const Matrix3d&, const PTransformd&)
    void orientation(const Quaterniond&)
    void orientation(const Matrix3d&)
    Matrix3d orientation() const

  cdef cppclass GazeTask(HighLevelTask):
    GazeTask(const vector[MultiBody]&, int, const string&, const Vector2d&, double, const PTransformd&, const Vector2d&)
    GazeTask(const vector[MultiBody]&, int, const string&, const Vector3d&, const PTransformd&, const Vector2d&)
    void error(const Vector2d&, const Vector2d&)
    void error(const Vector3d&, const Vector2d&)

  cdef cppclass PositionBasedVisServoTask(HighLevelTask):
    PositionBasedVisServoTask(const vector[MultiBody]&, int, const string&, const PTransformd&, const PTransformd&)
    void error(const PTransformd&)

  cdef cppclass PostureTask(Task):
    PostureTask(const vector[MultiBody]&, int, vector[vector[double]], double, double)
    double stiffness() const
    void stiffness(double)
    double damping() const
    void gains(double)
    void gains(double, double)
    void posture(vector[vector[double]])
    vector[vector[double]] posture()
    void jointsStiffness(const vector[MultiBody]&, vector[JointStiffness])
    void jointsGains(const vector[MultiBody]&, vector[JointGains])
    void ptUpdate "update"(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)
    VectorXd ptEval "eval"()

  cdef cppclass CoMTask(HighLevelTask):
    CoMTask(const vector[MultiBody]&, int, const Vector3d&)
    CoMTask(const vector[MultiBody]&, int, const Vector3d&, vector[double]) except +
    void com(const Vector3d&)
    Vector3d com() const
    void updateInertialParameters(const vector[MultiBody]&)

  cdef cppclass MomentumTask(HighLevelTask):
    MomentumTask(const vector[MultiBody]&, int, const ForceVecd&)
    void momentum(const ForceVecd&)
    ForceVecd momentum() const

  cdef cppclass LinVelocityTask(HighLevelTask):
    LinVelocityTask(const vector[MultiBody]&, int, const string&, const Vector3d&, const Vector3d&)
    void velocity(const Vector3d&)
    Vector3d velocity() const
    void bodyPoint(const Vector3d&)
    Vector3d bodyPoint() const

  cdef cppclass OrientationTrackingTask(HighLevelTask):
    OrientationTrackingTask(const vector[MultiBody]&, int, string&, const Vector3d&, const Vector3d&, const vector[string]&, const Vector3d&)
    void trackedPoint(const Vector3d&)
    Vector3d trackedPoint() const
    void bodyPoint(const Vector3d&)
    Vector3d bodyPoint() const
    void bodyAxis(const Vector3d&)
    Vector3d bodyAxis() const

  cdef cppclass JointsSelector(HighLevelTask):
    JointsSelector(const vector[MultiBody]&, int, HighLevelTask*, const vector[string]&)
    vector[JointsSelector_SelectedData] selectedJoints() const

  cdef cppclass JointsSelector_SelectedData "tasks::qp::JointsSelector::SelectedData":
    int posInDof
    int dof


  cdef cppclass SetPointTask(Task):
    SetPointTask(const vector[MultiBody]&, int, PositionTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, PositionTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, OrientationTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, OrientationTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, GazeTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, GazeTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, CoMTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, CoMTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, LinVelocityTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, MomentumTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, MomentumTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, JointsSelector*, double, double)
    SetPointTask(const vector[MultiBody]&, int, JointsSelector*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, TransformTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, TransformTask*, double, VectorXd, double)
    SetPointTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double)
    SetPointTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, VectorXd, double)

    double stiffness() const
    void stiffness(double)

    # SetPointTaskCommon
    VectorXd dimWeight() const
    void dimWeight(const VectorXd&)
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)
    MatrixXd Q() const
    VectorXd C() const

  cdef cppclass TrackingTask(Task):
    TrackingTask(const vector[MultiBody]&, int, PositionTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, PositionTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, OrientationTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, OrientationTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, GazeTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, GazeTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, CoMTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, CoMTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, MomentumTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, MomentumTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, JointsSelector*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, JointsSelector*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, TransformTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, TransformTask*, double, double, VectorXd, double)
    TrackingTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double, double)
    TrackingTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double, VectorXd, double)

    void setGains(double, double)
    void errorPos(const VectorXd&)
    void errorVel(const VectorXd&)
    void refAccel(const VectorXd&)

    # SetPointTaskCommon
    VectorXd dimWeight() const
    void dimWeight(const VectorXd&)
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)
    MatrixXd Q() const
    VectorXd C() const

  cdef cppclass TrajectoryTask(Task):
    TrajectoryTask(const vector[MultiBody]&, int, PositionTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, PositionTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, OrientationTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, OrientationTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, GazeTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, GazeTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, CoMTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, CoMTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, MomentumTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, MomentumTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, JointsSelector*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, JointsSelector*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, TransformTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, TransformTask*, double, double, VectorXd, double)
    TrajectoryTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double, double)
    TrajectoryTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double, VectorXd, double)

    void setGains(double, double)
    void refVel(const VectorXd&)
    void refAccel(const VectorXd&)

    # SetPointTaskCommon
    VectorXd dimWeight() const
    void dimWeight(const VectorXd&)
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)
    MatrixXd Q() const
    VectorXd C() const

  cdef cppclass TargetObjectiveTask(Task):
    TargetObjectiveTask(const vector[MultiBody]&, int, PositionTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, PositionTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, OrientationTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, OrientationTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, GazeTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, GazeTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, CoMTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, CoMTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, MomentumTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, MomentumTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, JointsSelector*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, JointsSelector*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, TransformTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, TransformTask*, double, double, VectorXd, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double, VectorXd, double)
    TargetObjectiveTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double, VectorXd, VectorXd, double)

    double duration() const
    void duration(double)
    int iter() const
    void iter(int)
    int nrIter() const
    void nrIter(int)
    VectorXd objDot() const
    void objDot(const VectorXd&)
    VectorXd phi() const
    VectorXd psi() const

    # SetPointTaskCommon
    VectorXd dimWeight() const
    void dimWeight(const VectorXd&)
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)
    MatrixXd Q() const
    VectorXd C() const

  cdef cppclass PIDTask(Task):
    PIDTask(const vector[MultiBody]&, int, PositionTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, PositionTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, OrientationTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, OrientationTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, SurfaceOrientationTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, GazeTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, GazeTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, PositionBasedVisServoTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, CoMTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, CoMTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, LinVelocityTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, OrientationTrackingTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, MomentumTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, MomentumTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, JointsSelector*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, JointsSelector*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, TransformTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, TransformTask*, double, double, double, VectorXd, double)
    PIDTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double, double, double)
    PIDTask(const vector[MultiBody]&, int, SurfaceTransformTask*, double, double, double, VectorXd, double)

    double P() const
    void P(double)
    double I() const
    void I(double)
    double D() const
    void D(double)
    void error(const VectorXd&)
    void errorI(const VectorXd&)
    void errorD(const VectorXd&)

    # SetPointTaskCommon
    VectorXd dimWeight() const
    void dimWeight(const VectorXd&)
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)
    MatrixXd Q() const
    VectorXd C() const

  cdef cppclass MultiCoMTask(Task):
    MultiCoMTask(const vector[MultiBody]&, vector[int], const Vector3d&, double, double)
    MultiCoMTask(const vector[MultiBody]&, vector[int], const Vector3d&, double, const Vector3d&, double)
    void com(const Vector3d&)
    Vector3d com() const
    void updateInertialParameters(const vector[MultiBody]&)
    void stiffness(double)
    double stiffness() const
    VectorXd dimWeight() const
    void dimWeight(const VectorXd&)
    VectorXd eval() const
    VectorXd speed() const
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)

  cdef cppclass MultiRobotTransformTask(Task):
    MultiRobotTransformTask(const vector[MultiBody]&, int, int, const string&, const string&, const PTransformd&, const PTransformd&, double, double)
    void X_r1b_r1s(const PTransformd&)
    PTransformd X_r1b_r1s() const
    void X_r2b_r2s(const PTransformd&)
    PTransformd X_r2b_r2s() const
    VectorXd eval() const
    VectorXd speed() const
    void stiffness(double)
    double stiffness() const
    VectorXd dimWeight() const
    void dimWeight(const VectorXd&)
    void update(const vector[MultiBody]&, const vector[MultiBodyConfig]&, const SolverData&)

  cdef cppclass ContactTask(Task):
    ContactTask(const ContactId&, double, double)
    void error(const Vector3d&)
    void errorD(const Vector3d&)

  cdef cppclass GripperTorqueTask(Task):
    GripperTorqueTask(const ContactId&, const Vector3d&, const Vector3d&, double)

cdef extern from "<Tasks/QPMotionConstr.h>" namespace "tasks::qp":
  cdef cppclass MotionPolyConstr(ConstraintFunction[GenInequality], GenInequality, Constraint):
    MotionPolyConstr(const vector[MultiBody]&, int, const c_tasks.PolyTorqueBound&)
    # Motion default
    void computeTorque(const VectorXd&, const VectorXd&)
    VectorXd torque() const
    void torque(const vector[MultiBody]&, const vector[MultiBodyConfig]&) const
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass MotionConstr(ConstraintFunction[GenInequality], GenInequality, Constraint):
    MotionConstr(const vector[MultiBody]&, int, const c_tasks.TorqueBound&)
    # Motion default
    void computeTorque(const VectorXd&, const VectorXd&)
    VectorXd torque() const
    void torque(const vector[MultiBody]&, const vector[MultiBodyConfig]&) const
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass MotionSpringConstr(ConstraintFunction[GenInequality], GenInequality, Constraint):
    MotionSpringConstr(const vector[MultiBody]&, int, const c_tasks.TorqueBound&, const vector[SpringJoint]&)
    # Motion default
    void computeTorque(const VectorXd&, const VectorXd&)
    VectorXd torque() const
    void torque(const vector[MultiBody]&, const vector[MultiBodyConfig]&) const
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass PositiveLambda(ConstraintFunction[Bound], Bound, Constraint):
    PositiveLambda()
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

cdef extern from "<Tasks/QPContactConstr.h>" namespace "tasks::qp":
  cdef cppclass ContactConstrCommon:
    bool addVirtualContact(const ContactId&)
    bool removeVirtualContact(const ContactId&)
    void resetVirtualContacts()
    bool addDofContact(const ContactId&, const MatrixXd&)
    bool removeDofContact(const ContactId&)
    void resetDofContacts()

  cdef cppclass ContactAccConstr(ContactConstrCommon, ConstraintFunction[Equality], Equality, Constraint):
    ContactAccConstr()
    void updateDofContacts()
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass ContactSpeedConstr(ContactConstrCommon, ConstraintFunction[Equality], Equality, Constraint):
    ContactSpeedConstr(double)
    void updateDofContacts()
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass ContactPosConstr(ContactConstrCommon, ConstraintFunction[Equality], Equality, Constraint):
    ContactPosConstr(double)
    void updateDofContacts()
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

cdef extern from "<Tasks/QPConstr.h>" namespace "tasks::qp":
  cdef cppclass CollisionConstr(ConstraintFunction[Inequality], Inequality, Constraint):
    CollisionConstr(const vector[MultiBody]&, double)
    void addCollision(const vector[MultiBody]&, int, int, const string&, S_Object*, const PTransformd&, int, const string&, S_Object*, const PTransformd&, double, double, double, double)
    bool rmCollision(int)
    int nrCollisions() const
    void reset()
    void updateNrCollisions()
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass CoMIncPlaneConstr(ConstraintFunction[Inequality], Inequality, Constraint):
    CoMIncPlaneConstr(const vector[MultiBody]&, int, double)
    void addPlane(int, const Vector3d&, double, double, double, double, double)
    bool rmPlane(int)
    int nrPlanes()
    void reset()
    void updateNrPlanes()
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass JointLimitsConstr(ConstraintFunction[Bound], Bound, Constraint):
    JointLimitsConstr(const vector[MultiBody]&, int, const c_tasks.QBound&, double)
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass DamperJointLimitsConstr(ConstraintFunction[Bound], Bound, Constraint):
    DamperJointLimitsConstr(const vector[MultiBody]&, int, const c_tasks.QBound&, const c_tasks.AlphaBound&, double, double, double, double)
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass GripperTorqueConstr(ConstraintFunction[Inequality], Inequality, Constraint):
    GripperTorqueConstr()
    void addGripper(const ContactId&, double, const Vector3d&, const Vector3d&)
    bool rmGripper(const ContactId&)
    void reset()
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass BoundedSpeedConstr(ConstraintFunction[GenInequality], GenInequality, Constraint):
    BoundedSpeedConstr(const vector[MultiBody]&, int, double)
    void addBoundedSpeed(const vector[MultiBody]&, const string&, const Vector3d&, const MatrixXd&, const VectorXd&)
    void addBoundedSpeed(const vector[MultiBody]&, const string&, const Vector3d&, const MatrixXd&, const VectorXd&, const VectorXd&)
    bool removeBoundedSpeed(const string&)
    void resetBoundedSpeeds()
    int nrBoundedSpeeds()
    void updateBoundedSpeeds()
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

  cdef cppclass ImageConstr(ConstraintFunction[Inequality], Inequality, Constraint):
    ImageConstr(const vector[MultiBody]&, int, const string&, const PTransformd&, double, double)
    void setLimits(const Vector2d&, const Vector2d&, const double, const double, const double, const double)
    int addPoint(const Vector2d&, const double)
    int addPoint(const Vector3d&)
    void addPoint(const vector[MultiBody]&, const string&, const PTransformd&)
    void reset()
    void updatePoint(const int, const Vector2d&)
    void updatePoint(const int, const Vector2d&, const double)
    void updatePoint(const int, const Vector3d&)
    void computeComponents(const MultiBody&, const MultiBodyConfig&, const SolverData&, const Vector2d&, const double, Jacobian&, const int, const PTransformd&, MatrixXd&, Vector2d&)
    void addToSolver(QPSolver &)
    void addToSolver(const vector[MultiBody]&, QPSolver &)
    void removeFromSolver(QPSolver &)

cdef extern from "<Tasks/QPSolver.h>" namespace "tasks::qp":
  cdef cppclass QPSolver:
    QPSolver()
    bool solve(const vector[MultiBody]&, vector[MultiBodyConfig]&)
    bool solveNoMbcUpdate(const vector[MultiBody]&, const vector[MultiBodyConfig]&)
    void updateMbc(MultiBodyConfig&, int) const
    void updateConstrSize()
    void nrVars(const vector[MultiBody]&, vector[UnilateralContact]&, vector[BilateralContact]&)
    int nrVars() const
    void updateTasksNrVars(const vector[MultiBody]&) const
    void updateConstrsNrVars(const vector[MultiBody]&) const
    void updateNrVars(const vector[MultiBody]&) const

    # EqualityConstraint
    void addEqualityConstraint(Equality*)
    void removeEqualityConstraint(Equality*)
    int nrEqualityConstraints() const
    # InequalityConstraint
    void addInequalityConstraint(Inequality*)
    void removeInequalityConstraint(Inequality*)
    int nrInequalityConstraints() const
    # GenInequalityConstraint
    void addGenInequalityConstraint(GenInequality*)
    void removeGenInequalityConstraint(GenInequality*)
    int nrGenInequalityConstraints() const
    # BoundConstraint
    void addBoundConstraint(Bound*)
    void removeBoundConstraint(Bound*)
    int nrBoundConstraints() const
    # Constraint
    void addConstraint(Constraint*)
    void addConstraint(const vector[MultiBody]&, Constraint*)
    void removeConstraint(Constraint*)
    int nrConstraints()
    # Task
    void addTask(Task*)
    void removeTask(Task*)
    int nrTasks() const
    void addTask(const vector[MultiBody]&, Task*)
    void resetTasks()

    void solver(const string&)
    string solver()
    VectorXd result() const
    VectorXd alphaDVec() const
    VectorXd alphaDVec(int) const
    VectorXd lambdaVec() const
    VectorXd lambdaVec(int) const
    int contactLambdaPosition(const ContactId&) const
    SolverData data() const
    c_tasks.cpu_times solveTime() const
    c_tasks.cpu_times solveAndBuildTime() const

