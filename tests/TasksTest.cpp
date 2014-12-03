// This file is part of Tasks.
//
// Tasks is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Tasks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Tasks.  If not, see <http://www.gnu.org/licenses/>.

// includes
// std
#include <fstream>
#include <iostream>
#include <tuple>

// boost
#define BOOST_TEST_MODULE TasksTest
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// Eigen
#include <unsupported/Eigen/Polynomials>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include <RBDyn/EulerIntegration.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>

// Tasks
#include "Tasks.h"

// Arms
#include "arms.h"


/// run the task update(mb, mbc) method
template<typename Task>
struct ClassicUpdater
{
	void operator()(Task& task, const rbd::MultiBody& mb,
		const rbd::MultiBodyConfig& mbc)
	{
		task.update(mb, mbc);
	}
};


/// Compute normal acceleration (like QPSolverData::computeNormalAccB)
void computeNormalAccB(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, std::vector<sva::MotionVecd>& normalAccB)
{
	const std::vector<int>& pred = mb.predecessors();
	const std::vector<int>& succ = mb.successors();

	for(int i = 0; i < mb.nrJoints(); ++i)
	{
		const sva::PTransformd& X_p_i = mbc.parentToSon[i];
		const sva::MotionVecd& vj_i = mbc.jointVelocity[i];
		const sva::MotionVecd& vb_i = mbc.bodyVelB[i];

		if(pred[i] != -1)
			normalAccB[succ[i]] = X_p_i*normalAccB[pred[i]] + vb_i.cross(vj_i);
		else
			normalAccB[succ[i]] = vb_i.cross(vj_i);
	}
}


/// run the task update(mb, mbc, bodyNormalAcc) method
template<typename Task>
struct NormalAccUpdater
{
	NormalAccUpdater(const rbd::MultiBody& mb):
		normalAccB(mb.nrBodies())
	{}

	void operator()(Task& task, const rbd::MultiBody& mb,
		const rbd::MultiBodyConfig& mbc)
	{
		computeNormalAccB(mb, mbc, normalAccB);
		task.update(mb, mbc, normalAccB);
	}

	std::vector<sva::MotionVecd> normalAccB;
};


/// run the task update(mb, mbc, com, bodyNormalAcc) method
template<typename Task>
struct NormalAccCoMUpdater
{
	NormalAccCoMUpdater(const rbd::MultiBody& mb):
		normalAccB(mb.nrBodies())
	{}

	void operator()(Task& task, const rbd::MultiBody& mb,
		const rbd::MultiBodyConfig& mbc)
	{
		computeNormalAccB(mb, mbc, normalAccB);
		Eigen::Vector3d com = rbd::computeCoM(mb, mbc);
		task.update(mb, mbc, com, normalAccB);
	}

	std::vector<sva::MotionVecd> normalAccB;
};


/// Test position task (eval, speed and acc are defined)
struct PosTester
{
	void operator()(const Eigen::VectorXd& speedCur,
		const Eigen::VectorXd& accCur, const Eigen::VectorXd& speedDiff,
		const Eigen::VectorXd& accDiff, double tol)
	{
		BOOST_CHECK_SMALL((speedCur - speedDiff).norm(), tol);
		BOOST_CHECK_SMALL((accCur - accDiff).norm(), tol);
	}
};


/// Test position task with orientation error (speedDiff is not realiable)
struct OriTaskTester
{
	void operator()(const Eigen::VectorXd& /* speedCur */,
		const Eigen::VectorXd& accCur, const Eigen::VectorXd& /* speedDiff */,
		const Eigen::VectorXd& accDiff, double tol)
	{
		BOOST_CHECK_SMALL((accCur - accDiff).norm(), tol);
	}
};


/// Test velocity task (eval, and acc are defined)
struct VelTester
{
	void operator()(const Eigen::VectorXd& /* speedCur */,
		const Eigen::VectorXd& accCur, const Eigen::VectorXd& speedDiff,
		const Eigen::VectorXd& /* accDiff */, double tol)
	{
		BOOST_CHECK_SMALL((accCur - speedDiff).norm(), tol);
	}
};


/**
	* Use finite difference to test a task.
	* Task is the task to test, Updater is a function to call the Task::update
	* method, Tester test the task speed and acceleration against the finite
	* difference speed and accelartion.
	*/
template<typename Task, typename Updater, typename Tester>
void testTaskNumDiff(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	Task& task, Updater updater, Tester tester,
	int nrIter=100, double diffStep=1e-6, double tol=1e-4)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;

	MultiBodyConfig mbcPost, mbcCur;

	Eigen::VectorXd q(mb.nrParams());
	Eigen::VectorXd alpha(mb.nrDof());
	Eigen::VectorXd alphaD(mb.nrDof());

	for(int i = 0; i < nrIter; ++i)
	{
		q.setRandom();
		alpha.setRandom();
		alphaD.setRandom();

		mbcCur = mbc;
		vectorToParam(q, mbcCur.q);
		vectorToParam(alpha, mbcCur.alpha);
		vectorToParam(alphaD, mbcCur.alphaD);

		mbcPost = mbcCur;

		eulerIntegration(mb, mbcPost, diffStep);

		forwardKinematics(mb, mbcCur);
		forwardKinematics(mb, mbcPost);
		forwardVelocity(mb, mbcCur);
		forwardVelocity(mb, mbcPost);

		updater(task, mb, mbcCur);
		VectorXd evalCur = task.eval();
		VectorXd speedCur = -task.speed();
		VectorXd accCur = -task.normalAcc() - task.jac()*alphaD;

		updater(task, mb, mbcPost);
		VectorXd evalPost = task.eval();
		VectorXd speedPost = -task.speed();

		VectorXd speedDiff = (evalPost - evalCur)/diffStep;
		VectorXd accDiff = (speedPost - speedCur)/diffStep;

		tester(speedCur, accCur, speedDiff, accDiff, tol);
	}
}


BOOST_AUTO_TEST_CASE(PositionTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::PositionTask pt(mb, 3, Vector3d::Random(), Vector3d::Random());

	testTaskNumDiff(mb, mbc, pt, ClassicUpdater<tasks::PositionTask>(),
		PosTester());
	testTaskNumDiff(mb, mbc, pt, NormalAccUpdater<tasks::PositionTask>(mb),
		PosTester());
}


BOOST_AUTO_TEST_CASE(OrientationTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::OrientationTask ot(mb, 3, Quaterniond(Vector4d::Random().normalized()));

	testTaskNumDiff(mb, mbc, ot, ClassicUpdater<tasks::OrientationTask>(),
		OriTaskTester());
	testTaskNumDiff(mb, mbc, ot, NormalAccUpdater<tasks::OrientationTask>(mb),
		OriTaskTester());
}


BOOST_AUTO_TEST_CASE(SurfaceOrientationTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::SurfaceOrientationTask sot(mb, 3,
		Quaterniond(Vector4d::Random().normalized()),
		sva::PTransformd(Quaterniond(Vector4d::Random().normalized()), Vector3d::Random()));

	testTaskNumDiff(mb, mbc, sot, ClassicUpdater<tasks::SurfaceOrientationTask>(),
		OriTaskTester());
	testTaskNumDiff(mb, mbc, sot, NormalAccUpdater<tasks::SurfaceOrientationTask>(mb),
		OriTaskTester());
}


BOOST_AUTO_TEST_CASE(CoMTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::CoMTask ct(mb, Vector3d::Random());

	testTaskNumDiff(mb, mbc, ct, ClassicUpdater<tasks::CoMTask>(),
		PosTester());
	testTaskNumDiff(mb, mbc, ct, NormalAccCoMUpdater<tasks::CoMTask>(mb),
		PosTester());
}


BOOST_AUTO_TEST_CASE(MomentumTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::MomentumTask mt(mb, sva::ForceVecd(Vector6d::Random()));

	testTaskNumDiff(mb, mbc, mt, ClassicUpdater<tasks::MomentumTask>(),
		VelTester());
	testTaskNumDiff(mb, mbc, mt, NormalAccUpdater<tasks::MomentumTask>(mb),
		VelTester());
}


BOOST_AUTO_TEST_CASE(LinVelocityTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::LinVelocityTask lvt(mb, 3, Vector3d::Random());

	testTaskNumDiff(mb, mbc, lvt, ClassicUpdater<tasks::LinVelocityTask>(),
		VelTester());
	testTaskNumDiff(mb, mbc, lvt, NormalAccUpdater<tasks::LinVelocityTask>(mb),
		VelTester());
}
