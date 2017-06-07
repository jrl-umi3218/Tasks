// Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
//
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

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include <RBDyn/EulerIntegration.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>

// Tasks
#include "Tasks/Tasks.h"

// Arms
#include "arms.h"


template<typename Task>
struct TanAccel
{
	Eigen::VectorXd tanAcc(Task& task, const std::vector<Eigen::VectorXd>& alphaD)
	{
		return Eigen::VectorXd(task.jac()*alphaD[0]);
	}
};


template<typename Task>
struct MRTanAccel
{
	MRTanAccel(int taskDim):
		tanAccV(taskDim)
	{}

	Eigen::VectorXd tanAcc(Task& task, const std::vector<Eigen::VectorXd>& alphaD)
	{
		tanAccV.setZero();
		for(std::size_t i = 0; i < alphaD.size(); ++i)
		{
			tanAccV += task.jac(int(i))*alphaD[i];
		}
		return tanAccV;
	}

	Eigen::VectorXd tanAccV;
};


/// run the task update(mb, mbc) method
template<typename Task>
struct ClassicUpdater : public TanAccel<Task>
{
	void operator()(Task& task, const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs)
	{
		task.update(mbs[0], mbcs[0]);
	}
};


/// run the task update(mbs, mbcs) method
template<typename Task>
struct MRClassicUpdater : public MRTanAccel<Task>
{
	MRClassicUpdater(int taskDim):
		MRTanAccel<Task>(taskDim)
	{}

	void operator()(Task& task, const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs)
	{
		task.update(mbs, mbcs);
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


/// Multi-robot version of computeNormalAccB
void computeNormalAccB(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	std::vector<std::vector<sva::MotionVecd>>& normalAccB)
{
	for(std::size_t i = 0; i < mbs.size(); ++i)
	{
		computeNormalAccB(mbs[i], mbcs[i], normalAccB[i]);
	}
}


/// run the task update(mb, mbc, bodyNormalAcc) method
template<typename Task>
struct NormalAccUpdater : public TanAccel<Task>
{
	NormalAccUpdater(const rbd::MultiBody& mb):
		normalAccB(mb.nrBodies())
	{}

	void operator()(Task& task, const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs)
	{
		computeNormalAccB(mbs[0], mbcs[0], normalAccB);
		task.update(mbs[0], mbcs[0], normalAccB);
	}

	std::vector<sva::MotionVecd> normalAccB;
};


/// run the task update(mbs, mbcs, bodyNormalAccs) method
template<typename Task>
struct MRNormalAccUpdater : public MRTanAccel<Task>
{
	MRNormalAccUpdater(const std::vector<rbd::MultiBody>& mbs, int taskDim):
		MRTanAccel<Task>(taskDim),
		normalAccBs(mbs.size())
	{
		for(std::size_t i = 0; i < mbs.size(); ++i)
		{
			normalAccBs[i].resize(mbs[i].nrBodies());
		}
	}

	void operator()(Task& task, const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs)
	{
		computeNormalAccB(mbs, mbcs, normalAccBs);
		task.update(mbs, mbcs, normalAccBs);
	}

	std::vector<std::vector<sva::MotionVecd>> normalAccBs;
};


/// run the task update(mb, mbc, com, bodyNormalAcc) method
template<typename Task>
struct NormalAccCoMUpdater : public TanAccel<Task>
{
	NormalAccCoMUpdater(const rbd::MultiBody& mb):
		normalAccB(mb.nrBodies())
	{}

	void operator()(Task& task, const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs)
	{
		computeNormalAccB(mbs[0], mbcs[0], normalAccB);
		Eigen::Vector3d com = rbd::computeCoM(mbs[0], mbcs[0]);
		task.update(mbs[0], mbcs[0], com, normalAccB);
	}

	std::vector<sva::MotionVecd> normalAccB;
};


/// run the task update(mbs, mbcs, coms, bodyNormalAccs) method
template<typename Task>
struct MRNormalAccCoMUpdater : public MRTanAccel<Task>
{
	MRNormalAccCoMUpdater(const std::vector<rbd::MultiBody>& mbs, int taskDim):
		MRTanAccel<Task>(taskDim),
		normalAccBs(mbs.size()),
		coms(mbs.size())
	{
		for(std::size_t i = 0; i < mbs.size(); ++i)
		{
			normalAccBs[i].resize(mbs[i].nrBodies());
		}
	}

	void operator()(Task& task, const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs)
	{
		computeNormalAccB(mbs, mbcs, normalAccBs);
		for(std::size_t i = 0; i < mbs.size(); ++i)
		{
			coms[i] = rbd::computeCoM(mbs[i], mbcs[i]);
		}
		task.update(mbs, mbcs, coms, normalAccBs);
	}

	std::vector<std::vector<sva::MotionVecd>> normalAccBs;
	std::vector<Eigen::Vector3d> coms;
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


/// Test TransformTask, position is reliable but only rotational acceleration
/// is reliable
struct PosTTTester
{
	void operator()(const Eigen::VectorXd& speedCur,
		const Eigen::VectorXd& accCur, const Eigen::VectorXd& speedDiff,
		const Eigen::VectorXd& accDiff, double tol)
	{
		BOOST_CHECK_SMALL((speedCur.tail<3>() - speedDiff.tail<3>()).norm(), tol);
		BOOST_CHECK_SMALL((accCur - accDiff).norm(), tol);
	}
};


/// Test MultiRobotTransformTask, the rotation part is not reliable since
/// he use rotationVelocity
struct PosMRTTTester
{
	void operator()(const Eigen::VectorXd& speedCur,
		const Eigen::VectorXd& accCur, const Eigen::VectorXd& speedDiff,
		const Eigen::VectorXd& accDiff, double tol)
	{
		BOOST_CHECK_SMALL((speedCur.tail<3>() - speedDiff.tail<3>()).norm(), tol);
		BOOST_CHECK_SMALL((accCur.tail<3>() - accDiff.tail<3>()).norm(), tol);
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

/// Test position task (eval, speed and acc are defined)
struct VectOriTester
{
	void operator()(const Eigen::VectorXd& speedCur,
		const Eigen::VectorXd& accCur, const Eigen::VectorXd& speedDiff,
		const Eigen::VectorXd& accDiff, double tol)
	{
		BOOST_CHECK_SMALL((speedCur - speedDiff).norm(), tol);
		BOOST_CHECK_SMALL((accCur - accDiff).norm(), tol);
	}
};

/**
	* Use finite difference to test a task.
	* Task is the task to test, Updater is a function to call the Task::update
	* method, Tester test the task speed and acceleration against the finite
	* difference speed and accelartion.
	*/
template<typename Task, typename Updater, typename Tester>
void testTaskNumDiff(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	Task& task, Updater updater, Tester tester,
	int nrIter=100, double diffStep=1e-6, double tol=1e-4)
{
	using namespace Eigen;
	using namespace sva;
	using namespace rbd;

	std::vector<MultiBodyConfig> mbcsPost(mbcs), mbcsCur(mbcs);

	std::vector<Eigen::VectorXd> q(mbs.size());
	std::vector<Eigen::VectorXd> alpha(mbs.size());
	std::vector<Eigen::VectorXd> alphaD(mbs.size());

	for(int i = 0; i < nrIter; ++i)
	{
		for(std::size_t r = 0; r < mbs.size(); ++r)
		{
			const rbd::MultiBody& mb = mbs[r];
			const rbd::MultiBodyConfig& mbc = mbcs[r];
			rbd::MultiBodyConfig& mbcPost = mbcsPost[r];
			rbd::MultiBodyConfig& mbcCur = mbcsCur[r];

			q[r].setRandom(mb.nrParams());
			alpha[r].setRandom(mb.nrDof());
			alphaD[r].setRandom(mb.nrDof());

			mbcCur = mbc;
			vectorToParam(q[r], mbcCur.q);
			vectorToParam(alpha[r], mbcCur.alpha);
			vectorToParam(alphaD[r], mbcCur.alphaD);

			mbcPost = mbcCur;

			eulerIntegration(mb, mbcPost, diffStep);

			forwardKinematics(mb, mbcCur);
			forwardKinematics(mb, mbcPost);
			forwardVelocity(mb, mbcCur);
			forwardVelocity(mb, mbcPost);
		}

		updater(task, mbs, mbcsCur);
		VectorXd evalCur = task.eval();
		VectorXd speedCur = -task.speed();
		VectorXd accCur = -task.normalAcc();
		accCur -= updater.tanAcc(task, alphaD);

		updater(task, mbs, mbcsPost);
		VectorXd evalPost = task.eval();
		VectorXd speedPost = -task.speed();

		VectorXd speedDiff = (evalPost - evalCur)/diffStep;
		VectorXd accDiff = (speedPost - speedCur)/diffStep;

		tester(speedCur, accCur, speedDiff, accDiff, tol);
	}
}


template<typename Task, typename Updater, typename Tester>
void testTaskNumDiff(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	Task& task, Updater updater, Tester tester,
	int nrIter=100, double diffStep=1e-6, double tol=1e-4)
{
	testTaskNumDiff(std::vector<rbd::MultiBody>{mb},
		std::vector<rbd::MultiBodyConfig>{mbc},
		task, updater, tester, nrIter, diffStep, tol);
}


BOOST_AUTO_TEST_CASE(PositionTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::PositionTask pt(mb, "b3", Vector3d::Random(), Vector3d::Random());

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

	tasks::OrientationTask ot(mb, "b3", Quaterniond(Vector4d::Random().normalized()));

	testTaskNumDiff(mb, mbc, ot, ClassicUpdater<tasks::OrientationTask>(),
		OriTaskTester());
	testTaskNumDiff(mb, mbc, ot, NormalAccUpdater<tasks::OrientationTask>(mb),
		OriTaskTester());
}


BOOST_AUTO_TEST_CASE(TransformTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	sva::PTransformd X_b_s(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());
	sva::PTransformd X_0_t(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());
	Eigen::Quaterniond E_0_c(Vector4d::Random().normalized());

	tasks::TransformTask tt(mb, "b3", X_0_t, X_b_s, E_0_c.matrix());

	testTaskNumDiff(mb, mbc, tt,
		NormalAccUpdater<tasks::TransformTask>(mb), PosTTTester(), 100);
}


BOOST_AUTO_TEST_CASE(SurfaceTransformTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	sva::PTransformd X_b_s(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());
	sva::PTransformd X_0_t(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());

	tasks::SurfaceTransformTask tt(mb, "b3", X_0_t, X_b_s);

	testTaskNumDiff(mb, mbc, tt,
		NormalAccUpdater<tasks::SurfaceTransformTask>(mb), PosMRTTTester(), 100);
}


BOOST_AUTO_TEST_CASE(MultiRobotTransformTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb1, mb2;
	MultiBodyConfig mbc1, mbc2;

	std::tie(mb1, mbc1) = makeZXZArm();
	std::tie(mb2, mbc2) = makeZXZArm();

	std::vector<MultiBody> mbs{mb1, mb2};
	std::vector<MultiBodyConfig> mbcs{mbc1, mbc2};

	sva::PTransformd X_r1b_r1s(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());
	sva::PTransformd X_r2b_r2s(Quaterniond(Vector4d::Random().normalized()),
		Vector3d::Random());

	tasks::MultiRobotTransformTask mrtt(mbs, 0, 1, "b3", "b3", X_r1b_r1s, X_r2b_r2s);

	testTaskNumDiff(mbs, mbcs, mrtt,
		MRNormalAccUpdater<tasks::MultiRobotTransformTask>(mbs, 6), PosMRTTTester());
}


BOOST_AUTO_TEST_CASE(SurfaceOrientationTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::SurfaceOrientationTask sot(mb, "b3",
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


BOOST_AUTO_TEST_CASE(MultiCoMTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb1, mb2;
	MultiBodyConfig mbc1, mbc2;

	std::tie(mb1, mbc1) = makeZXZArm();
	std::tie(mb2, mbc2) = makeZXZArm();

	std::vector<MultiBody> mbs{mb1, mb2};
	std::vector<MultiBodyConfig> mbcs{mbc1, mbc2};

	tasks::MultiCoMTask mct(mbs, {0, 1}, Vector3d::Random());

	testTaskNumDiff(mbs, mbcs, mct, MRClassicUpdater<tasks::MultiCoMTask>(3),
		PosTester());
	testTaskNumDiff(mbs, mbcs, mct, MRNormalAccUpdater<tasks::MultiCoMTask>(mbs, 3),
		PosTester());
	testTaskNumDiff(mbs, mbcs, mct, MRNormalAccCoMUpdater<tasks::MultiCoMTask>(mbs, 3),
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

	tasks::LinVelocityTask lvt(mb, "b3", Vector3d::Random());

	testTaskNumDiff(mb, mbc, lvt, ClassicUpdater<tasks::LinVelocityTask>(),
		VelTester());
	testTaskNumDiff(mb, mbc, lvt, NormalAccUpdater<tasks::LinVelocityTask>(mb),
		VelTester());
}

BOOST_AUTO_TEST_CASE(VectorOrientationTaskTest)
{
	using namespace Eigen;
	using namespace rbd;

	MultiBody mb;
	MultiBodyConfig mbc;

	std::tie(mb, mbc) = makeZXZArm();

	tasks::VectorOrientationTask vot(mb, "b3", Vector3d::Random(), Vector3d::Random());

	testTaskNumDiff(mb, mbc, vot, NormalAccUpdater<tasks::VectorOrientationTask>(mb),
		VectOriTester());
}
