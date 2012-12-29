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
// std
#include <vector>

// Eigen
#include <Eigen/Core>

// forward declaration
// rbd
namespace rbd
{
class Jacobian;
class MultiBody;
class MultiBodyConfig;
}


namespace tasks
{

namespace pg
{

Eigen::Matrix<double, 4, 3> angularVelToQuatVel(const std::vector<double>& q);
Eigen::Matrix<double, 3, 4> quatVelToAngularVel(const std::vector<double>& q);

class PGJacobian
{
public:
	PGJacobian(const rbd::MultiBody& mb, const rbd::Jacobian& jac);

	const Eigen::MatrixXd& jacobian(const rbd::MultiBody& mb,
		const rbd::MultiBodyConfig& mbc,
		const Eigen::MatrixXd& jacobian);

	void fullJacobian(const rbd::MultiBody& mb, const rbd::Jacobian& jac,
		const Eigen::MatrixXd& jacobian, Eigen::MatrixXd& res) const;

	int params() const
	{
		return jac_.cols();
	}

private:
	struct Joint
	{
		Joint(int index, int pos, int posPG):
			index(index),
			pos(pos),
			posPG(posPG)
		{}

		int index;
		int pos, posPG;
	};
private:
	std::vector<Joint> free_;
	std::vector<Joint> spherical_;
	std::vector<Joint> other_;

	Eigen::MatrixXd jac_;
};

} // namespace pg

} // namespace tasks
