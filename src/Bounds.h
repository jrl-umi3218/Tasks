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
#include <utility>
#include <vector>

// Eigen
#include <Eigen/Core>


namespace tasks
{

/// General position vector bound.
struct QBound
{
	QBound() {}
	QBound(std::vector<std::vector<double>> lQB,
		std::vector<std::vector<double>> uQB):
		lQBound(std::move(lQB)),
		uQBound(std::move(uQB))
	{ }

	std::vector<std::vector<double>> lQBound;
	std::vector<std::vector<double>> uQBound;
};


/// General velocity vector bound.
struct AlphaBound
{
	AlphaBound() {}
	AlphaBound(std::vector<std::vector<double>> lAB,
		std::vector<std::vector<double>> uAB):
		lAlphaBound(std::move(lAB)),
		uAlphaBound(std::move(uAB))
	{ }

	std::vector<std::vector<double>> lAlphaBound;
	std::vector<std::vector<double>> uAlphaBound;
};


/// General force vector bound.
struct TorqueBound
{
	TorqueBound() {}
	TorqueBound(std::vector<std::vector<double>> lTB,
		std::vector<std::vector<double>> uTB):
		lTorqueBound(std::move(lTB)),
		uTorqueBound(std::move(uTB))
	{ }

	std::vector<std::vector<double>> lTorqueBound;
	std::vector<std::vector<double>> uTorqueBound;
};


/// General force vector bound.
struct PolyTorqueBound
{
	PolyTorqueBound() {}
	PolyTorqueBound(std::vector<std::vector<Eigen::VectorXd>> lPTB,
		std::vector<std::vector<Eigen::VectorXd>> uPTB):
		lPolyTorqueBound(std::move(lPTB)),
		uPolyTorqueBound(std::move(uPTB))
	{ }


	std::vector<std::vector<Eigen::VectorXd>> lPolyTorqueBound;
	std::vector<std::vector<Eigen::VectorXd>> uPolyTorqueBound;
};

} // namespace tasks
