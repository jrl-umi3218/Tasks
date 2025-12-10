/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// associated header
#include "Tasks/GenQPSolver.h"

// includes
#include <cassert>
// std
#include <map>

// Tasks
#include "QLDQPSolver.h"

#ifdef LSSOL_SOLVER_FOUND
#  include "LSSOLQPSolver.h"
#endif

namespace tasks
{

namespace qp
{

#ifdef LSSOL_SOLVER_FOUND
const std::string GenQPSolver::default_qp_solver("LSSOL");
#else
const std::string GenQPSolver::default_qp_solver("QLD");
#endif

template<typename T>
T * allocateQP()
{
  return new T;
}

static const std::map<std::string, std::function<GenQPSolver *(void)>> qpFactory = {
#ifdef LSSOL_SOLVER_FOUND
    {"LSSOL", allocateQP<LSSOLQPSolver>},
#endif
    {"QLD", allocateQP<QLDQPSolver>}};

GenQPSolver * createQPSolver(const std::string & name)
{
  return qpFactory.at(name)();
}

void GenQPSolver::setDependencies(int nrVars, std::vector<std::tuple<int, int, double>> dependencies)
{
  assert(static_cast<size_t>(nrVars) > dependencies.size());
  dependencies_ = dependencies;
  using dependency_t = std::tuple<int, int, double>;
  /* Sort dependencies by the index of removed variables */
  std::sort(dependencies_.begin(), dependencies_.end(),
            [](const dependency_t & lhs, const dependency_t & rhs) { return std::get<1>(lhs) < std::get<1>(rhs); });
  fullToReduced_.resize(static_cast<size_t>(nrVars), -1);
  reducedToFull_.resize(static_cast<size_t>(nrVars) - dependencies_.size(), -1);
  /* Initialize the multipliers and offset variable */
  multipliers_ = Eigen::SparseMatrix<double>(nrVars, nrVars - static_cast<int>(dependencies.size()));
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<size_t>(nrVars));
  /* Number of removed variables encountered so far */
  size_t shift = 0;
  for(size_t i = 0; i < fullToReduced_.size(); ++i)
  {
    if(shift < dependencies_.size() && std::get<1>(dependencies_[shift]) == static_cast<int>(i))
    {
      shift++;
      continue;
    }
    fullToReduced_[i] = static_cast<int>(i - shift);
    triplets.push_back({static_cast<Eigen::SparseMatrix<double>::StorageIndex>(i), fullToReduced_[i], 1.0});
    reducedToFull_[i - shift] = static_cast<int>(i);
  }
  for(const auto & d : dependencies_)
  {
    auto leader_idx = std::get<0>(d);
    auto mimic_idx = std::get<1>(d);
    auto mult = std::get<2>(d);
    triplets.push_back({mimic_idx, fullToReduced_[static_cast<size_t>(leader_idx)], mult});
  }
  multipliers_.setFromTriplets(triplets.begin(), triplets.end());
}

} // namespace qp

} // namespace tasks
