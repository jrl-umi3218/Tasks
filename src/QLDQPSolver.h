/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
// eigen-qld
#include <eigen-qld/QLD.h>

// Tasks
#include "Tasks/GenQPSolver.h"

namespace tasks
{

namespace qp
{

/**
 * GenQPSolver interface implementation with the QLD QP solver.
 */
class TASKS_DLLAPI QLDQPSolver : public GenQPSolver
{
public:
  QLDQPSolver();

  virtual void updateSize(int nrVars, int nrEq, int nrInEq, int nrGenInEq) override;
  virtual void updateMatrix(const std::vector<Task *> & tasks,
                            const std::vector<Equality *> & eqConstr,
                            const std::vector<Inequality *> & inEqConstr,
                            const std::vector<GenInequality *> & genInEqConstr,
                            const std::vector<Bound *> & boundConstr) override;
  virtual bool solve() override;
  virtual const Eigen::VectorXd & result() const override;
  virtual std::ostream & errorMsg(const std::vector<rbd::MultiBody> & mbs,
                                  const std::vector<Task *> & tasks,
                                  const std::vector<Equality *> & eqConstr,
                                  const std::vector<Inequality *> & inEqConstr,
                                  const std::vector<GenInequality *> & genInEqConstr,
                                  const std::vector<Bound *> & boundConstr,
                                  std::ostream & out) const override;

private:
  Eigen::QLD qld_;

  Eigen::MatrixXd Aeq_, Aineq_;
  Eigen::VectorXd beq_, bineq_;

  Eigen::MatrixXd AeqFull_, AineqFull_;

  Eigen::VectorXd XL_;
  Eigen::VectorXd XU_;

  Eigen::VectorXd XLFull_;
  Eigen::VectorXd XUFull_;

  Eigen::MatrixXd Q_;
  Eigen::VectorXd C_;

  Eigen::MatrixXd QFull_;
  Eigen::VectorXd CFull_;

  Eigen::VectorXd XFull_;

  int nrAeqLines_;
  int nrAineqLines_;
};

} // namespace qp

} // namespace tasks
