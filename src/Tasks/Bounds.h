/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
// std
#include <utility>
#include <vector>

// Eigen
#include <Eigen/Core>


namespace tasks
{


/**
	* General position vector bounds
	* \f$ \underline{q} \f$ and \f$ \overline{q} \f$.
	*/
struct QBound
{
	QBound() {}
	QBound(std::vector<std::vector<double>> lQB,
		std::vector<std::vector<double>> uQB):
		lQBound(std::move(lQB)),
		uQBound(std::move(uQB))
	{ }

	/// \f$ \underline{q} \f$
	std::vector<std::vector<double>> lQBound;
	/// \f$ \overline{q} \f$
	std::vector<std::vector<double>> uQBound;
};


/**
	* General velocity vector bounds
	* \f$ \underline{\alpha} \f$ and \f$ \overline{\alpha} \f$.
	*/
struct AlphaBound
{
	AlphaBound() {}
	AlphaBound(std::vector<std::vector<double>> lAB,
		std::vector<std::vector<double>> uAB):
		lAlphaBound(std::move(lAB)),
		uAlphaBound(std::move(uAB))
	{ }

	/// \f$ \underline{\alpha} \f$
	std::vector<std::vector<double>> lAlphaBound;
	/// \f$ \overline{\alpha} \f$
	std::vector<std::vector<double>> uAlphaBound;
};


/**
	* General force vector bounds
	* \f$ \underline{\tau} \f$ and \f$ \overline{\tau} \f$.
	*/
struct TorqueBound
{
	TorqueBound() {}
	TorqueBound(std::vector<std::vector<double>> lTB,
		std::vector<std::vector<double>> uTB):
		lTorqueBound(std::move(lTB)),
		uTorqueBound(std::move(uTB))
	{ }

	/// \f$ \underline{\tau} \f$
	std::vector<std::vector<double>> lTorqueBound;
	/// \f$ \overline{\tau} \f$
	std::vector<std::vector<double>> uTorqueBound;
};


/**
	* General force vector bounds in function of articular position
	* \f$ \underline{\tau}(q) \f$ and \f$ \overline{\tau}(q) \f$.
	* The upper and lower bound function is represented
	* by a polynomial of any degree (coefficients are ordered by degree).
	*/
struct PolyTorqueBound
{
	PolyTorqueBound() {}
	PolyTorqueBound(std::vector<std::vector<Eigen::VectorXd>> lPTB,
		std::vector<std::vector<Eigen::VectorXd>> uPTB):
		lPolyTorqueBound(std::move(lPTB)),
		uPolyTorqueBound(std::move(uPTB))
	{ }

	/// \f$ \underline{\tau}(q) \f$
	std::vector<std::vector<Eigen::VectorXd>> lPolyTorqueBound;
	/// \f$ \overline{\tau}(q) \f$
	std::vector<std::vector<Eigen::VectorXd>> uPolyTorqueBound;
};


} // namespace tasks
