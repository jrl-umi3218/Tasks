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
#include <CoM.h>
#include <Jacobian.h>

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
	PositionTask(const rbd::MultiBody& mb, int bodyId, const Eigen::Vector3d& pos);

	void position(const Eigen::Vector3d& pos);
	const Eigen::Vector3d& position() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	Eigen::Vector3d pos_;
	int bodyIndex_;
	rbd::Jacobian jac_;

	Eigen::VectorXd eval_;
	Eigen::MatrixXd shortJacMat_;
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
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	Eigen::Matrix3d ori_;
	int bodyIndex_;
	rbd::Jacobian jac_;

	Eigen::VectorXd eval_;
	Eigen::MatrixXd shortJacMat_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};



class PostureTask
{
public:
	PostureTask(const rbd::MultiBody& mb, std::vector<std::vector<double>> q);

	void posture(std::vector<std::vector<double>> q);
	const std::vector<std::vector<double>> posture() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	std::vector<std::vector<double>> q_;

	Eigen::VectorXd eval_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};



class CoMTask
{
public:
	CoMTask(const rbd::MultiBody& mb, Eigen::Vector3d& com);

	void com(const Eigen::Vector3d& com);
	const Eigen::Vector3d com() const;

	void update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);
	void updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc);

	const Eigen::VectorXd& eval() const;
	const Eigen::MatrixXd& jac() const;
	const Eigen::MatrixXd& jacDot() const;

private:
	Eigen::Vector3d com_;
	rbd::CoMJacobianDummy jac_;

	Eigen::VectorXd eval_;
	Eigen::MatrixXd jacMat_;
	Eigen::MatrixXd jacDotMat_;
};

} // namespace tasks
