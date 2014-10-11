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

#pragma once

// includes
// Eigen
#include <Eigen/Core>

// RBDyn
#include <RBDyn/CoM.h>
#include <RBDyn/Jacobian.h>
#include <RBDyn/Momentum.h>

// forward declarations
// RBDyn
namespace rbd
{
	class MultiBody;
	class MultiBodyConfig;
}

namespace tasks
{


class PositionTask
{
public:
	PositionTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Vector3d& pos,
		const Eigen::Vector3d& bodyPoint=Eigen::Vector3d::Zero());

	void position(const Eigen::Vector3d& pos);
	const Eigen::Vector3d& position() const;

	void bodyPoint(const Eigen::Vector3d& point);
	const Eigen::Vector3d& bodyPoint() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const std::vector<sva::MotionVecd>& normalAccB);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::VectorXd& speed() const;
	const Eigen::VectorXd& normalAcc() const;

	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	Eigen::Vector3d pos_;
	sva::PTransformd point_;
	int bodyIndex_;
	rbd::Jacobian jac_;

	Eigen::VectorXd eval_;
	Eigen::VectorXd speed_;
	Eigen::VectorXd normalAcc_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};



class OrientationTask
{
public:
	OrientationTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Quaterniond& ori);
	OrientationTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Matrix3d& ori);

	void orientation(const Eigen::Quaterniond& ori);
	void orientation(const Eigen::Matrix3d& ori);
	const Eigen::Matrix3d& orientation() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const std::vector<sva::MotionVecd>& normalAccB);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::VectorXd& speed() const;
	const Eigen::VectorXd& normalAcc() const;

	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	Eigen::Matrix3d ori_;
	int bodyIndex_;
	rbd::Jacobian jac_;

	Eigen::VectorXd eval_;
	Eigen::VectorXd speed_;
	Eigen::VectorXd normalAcc_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};



class PostureTask
{
public:
	PostureTask(const rbd::MultiBody& mb, std::vector<std::vector<double> > q);

	void posture(std::vector<std::vector<double> > q);
	const std::vector<std::vector<double> > posture() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	std::vector<std::vector<double> > q_;

	Eigen::VectorXd eval_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};



class CoMTask
{
public:
	CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com);
	CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com,
				 std::vector<double> weight);

	void com(const Eigen::Vector3d& com);
	const Eigen::Vector3d com() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const Eigen::Vector3d& com, const std::vector<sva::MotionVecd>& normalAccB);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::VectorXd& speed() const;
	const Eigen::VectorXd& normalAcc() const;

	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	Eigen::Vector3d com_;
	rbd::CoMJacobian jac_;

	Eigen::VectorXd eval_;
	Eigen::VectorXd speed_;
	Eigen::VectorXd normalAcc_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};



class MultiCoMTask
{
public:
	MultiCoMTask(const std::vector<rbd::MultiBody>& mbs,
		std::vector<int> robotIndexes, const Eigen::Vector3d& com);

	const std::vector<int>& robotIndexes() const;

	void com(const Eigen::Vector3d& com);
	const Eigen::Vector3d com() const;

	void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs);
	void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const std::vector<std::vector<sva::MotionVecd>>& normalAccB);
	void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const std::vector<Eigen::Vector3d>& coms,
		const std::vector<std::vector<sva::MotionVecd>>& normalAccB);

	const Eigen::VectorXd& eval() const;
	const Eigen::VectorXd& speed() const;
	const Eigen::VectorXd& normalAcc() const;

	const Eigen::MatrixXd& jac(int index) const;

private:
	void update(const std::vector<rbd::MultiBody>& mbs,
		const std::vector<rbd::MultiBodyConfig>& mbcs,
		const std::vector<std::vector<sva::MotionVecd>>& normalAccB,
		const Eigen::Vector3d& multiCoM);

private:
	Eigen::Vector3d com_;
	std::vector<int> robotIndexes_;
	std::vector<double> robotsWeight_;
	std::vector<rbd::CoMJacobian> jac_;

	Eigen::VectorXd eval_;
	Eigen::VectorXd speed_;
	Eigen::VectorXd normalAcc_;
	std::vector<Eigen::MatrixXd> jacMat_;
};



class MomentumTask
{
public:
	MomentumTask(const rbd::MultiBody& mb, const sva::ForceVecd mom);

	void momentum(const sva::ForceVecd& mom);
	const sva::ForceVecd momentum() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const std::vector<sva::MotionVecd>& normalAccB);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::VectorXd& speed() const;
	const Eigen::VectorXd& normalAcc() const;

	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:

	sva::ForceVecd momentum_;
	rbd::CentroidalMomentumMatrix momentumMatrix_;
	Eigen::VectorXd eval_;
	Eigen::VectorXd speed_;
	Eigen::VectorXd normalAcc_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};



class LinVelocityTask
{
public:
	LinVelocityTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Vector3d& vel,
		const Eigen::Vector3d& bodyPoint=Eigen::Vector3d::Zero());

	void velocity(const Eigen::Vector3d& s);
	const Eigen::Vector3d& velocity() const;

	void bodyPoint(const Eigen::Vector3d& point);
	const Eigen::Vector3d& bodyPoint() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
		const std::vector<sva::MotionVecd>& normalAccB);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::VectorXd& speed() const;
	const Eigen::VectorXd& normalAcc() const;

	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	Eigen::Vector3d vel_;
	sva::PTransformd point_;
	int bodyIndex_;
	rbd::Jacobian jac_;

	Eigen::VectorXd eval_;
	Eigen::VectorXd speed_;
	Eigen::VectorXd normalAcc_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};



class OrientationTrackingTask
{
public:
	OrientationTrackingTask(const rbd::MultiBody& mb, int bodyId,
		const Eigen::Vector3d& bodyPoint, const Eigen::Vector3d& bodyAxis,
		const std::vector<int>& trackingJointsId,
		const Eigen::Vector3d& trackedPoint);

	void trackedPoint(const Eigen::Vector3d& tp);
	const Eigen::Vector3d& trackedPoint() const;

	void bodyPoint(const Eigen::Vector3d& bp);
	const Eigen::Vector3d& bodyPoint() const;

	void bodyAxis(const Eigen::Vector3d& ba);
	const Eigen::Vector3d& bodyAxis() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	virtual const Eigen::MatrixXd& jac();
	virtual const Eigen::MatrixXd& jacDot();
	virtual const Eigen::VectorXd& eval();

private:
	void zeroJacobian(Eigen::MatrixXd& jac) const;

private:
	int bodyIndex_;
	sva::PTransformd bodyPoint_;
	Eigen::Vector3d bodyAxis_;
	std::vector<int> zeroJacIndex_;
	Eigen::Vector3d trackedPoint_;
	rbd::Jacobian jac_;

	Eigen::VectorXd eval_;
	Eigen::MatrixXd shortJacMat_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};

} // namespace tasks
