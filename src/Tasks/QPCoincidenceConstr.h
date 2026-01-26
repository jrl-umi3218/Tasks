#pragma once

// include
// std
#include <vector>

// Eigen
#include <Eigen/Core>

// Tasks
#include "QPSolver.h"

namespace tasks
{

namespace qp
{

class TASKS_DLLAPI CoincidenceConstr : public ConstraintFunction<Equality>
{
public:
  CoincidenceConstr(int robotIndex,
                    const std::string & body1Name,
                    const std::string & body2Name,
                    const Eigen::VectorXd & jointSelector)
  : robotIndex_(robotIndex), body1Name_(body1Name), body2Name_(body2Name), body1Index_(-1), body2Index_(-1),
    jointSelector_(jointSelector)
  {
  }
  virtual int maxEq() const override;
  virtual const Eigen::MatrixXd & AEq() const override;
  virtual const Eigen::VectorXd & bEq() const override;
  void setJointSelector(const Eigen::VectorXd & selector);

protected:
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;
  int robotIndex_;
  std::string body1Name_, body2Name_;
  int body1Index_, body2Index_;
  Eigen::VectorXd jointSelector_;
};

class FixedCoincidenceConstr : public CoincidenceConstr
{
public:
  /**
   * Constructor
   * @param robotIndex
   * @param body1Name
   * @param body2Name
   * @param point1
   * @param point2
   */
  FixedCoincidenceConstr(int robotIndex = 0,
                         const std::string & body1Name = "",
                         const std::string & body2Name = "",
                         const Eigen::Vector3d & point1 = Eigen::Vector3d::Zero(),
                         const Eigen::Vector3d & point2 = Eigen::Vector3d::Zero(),
                         const Eigen::VectorXd & jointSelector = Eigen::VectorXd());

  virtual std::string nameEq() const override;
  virtual std::string descEq(const std::vector<rbd::MultiBody> & mbs, int i) override;
  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;

private:
  Eigen::Vector3d point1_, point2_;
};

class RotationalCoincidenceConstr : public CoincidenceConstr
{
public:
  /**
   * Constructor
   * @param robotIndex
   * @param body1Name
   * @param body2Name
   */
  RotationalCoincidenceConstr(int robotIndex = 0,
                              const std::string & body1Name = "",
                              const std::string & body2Name = "",
                              const Eigen::VectorXd & jointSelector = Eigen::VectorXd());

  virtual std::string nameEq() const override;
  virtual std::string descEq(const std::vector<rbd::MultiBody> & mbs, int i) override;
  virtual void updateNrVars(const std::vector<rbd::MultiBody> & mbs, const SolverData & data) override;
  virtual void update(const std::vector<rbd::MultiBody> & mbs,
                      const std::vector<rbd::MultiBodyConfig> & mbcs,
                      const SolverData & data) override;
};

} // namespace qp

} // namespace tasks
