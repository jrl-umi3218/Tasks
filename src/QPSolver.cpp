/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// associated header
#include "Tasks/QPSolver.h"

// includes
// std
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

// Tasks
#include "Tasks/GenQPSolver.h"
#include "Tasks/QPMotionConstr.h"

// included by Rafa as a test
#include "LSSOLQPSolver.h"

namespace tasks
{

namespace qp
{

/**
 *													QPSolver
 */

QPSolver::QPSolver()
: constr_(), eqConstr_(), inEqConstr_(), genInEqConstr_(), boundConstr_(), tasks_(), maxEqLines_(0), maxInEqLines_(0),
  maxGenInEqLines_(0), solver_(createQPSolver(GenQPSolver::default_qp_solver))
{
}

// must declare it in cpp because of GenQPSolver fwd declarition
QPSolver::~QPSolver() {}

bool QPSolver::solve(const std::vector<rbd::MultiBody> & mbs, std::vector<rbd::MultiBodyConfig> & mbcs)
{
  bool success = solveNoMbcUpdate(mbs, mbcs);

  postUpdate(mbs, mbcs, success);

  return success;
}

bool QPSolver::solveNoMbcUpdate(const std::vector<rbd::MultiBody> & mbs, const std::vector<rbd::MultiBodyConfig> & mbcs)
{
  solverAndBuildTimer_.start();
  preUpdate(mbs, mbcs);

  solverTimer_.start();
  bool success = solver_->solve();
  solverTimer_.stop();

  data_.lambdaVecPrev_ = lambdaVec();

  std::cout << "Rafa, in QPSolver::solveNoMbcUpdate, lambdaVec() = "
            << lambdaVec().transpose() << std::endl;

  //volatile bool check = data_.lambdaVecPrev_.isZero(0);
  //volatile int dummy = 0;

  //if (check) {
  //  dummy = 1;
  //  std::cout << "Rafa, in QPSolver::solveNoMbcUpdate, lambdaVec() = 0" << std::endl;
  //}
  
  // if (data_.lambdaVecPrev_.isZero(0))
  //   *(int*)0 = 0;
  
  if(!success)
  {
    solver_->errorMsg(mbs, tasks_, eqConstr_, inEqConstr_, genInEqConstr_, boundConstr_, std::cerr) << std::endl;
  }
  solverAndBuildTimer_.stop();

  return success;
}

void QPSolver::updateMbc(rbd::MultiBodyConfig & mbc, int rI) const
{
  rbd::vectorToParam(solver_->result().segment(data_.alphaDBegin_[rI], data_.alphaD_[rI]), mbc.alphaD);
}

template<typename T>
int accumMaxLines(int acc, T * constr)
{
  return acc + constr_traits<T>::maxLines(constr);
}

void QPSolver::updateConstrSize()
{
  maxEqLines_ = std::accumulate(eqConstr_.begin(), eqConstr_.end(), 0, accumMaxLines<Equality>);
  maxInEqLines_ = std::accumulate(inEqConstr_.begin(), inEqConstr_.end(), 0, accumMaxLines<Inequality>);
  maxGenInEqLines_ = std::accumulate(genInEqConstr_.begin(), genInEqConstr_.end(), 0, accumMaxLines<GenInequality>);

  std::cout << "Rafa, in tasks::qp::QPSolver::updateConstrSize, maxEqLines_ = " << maxEqLines_ << ", maxInEqLines_ = " << maxInEqLines_ << ", maxGenInEqLines_ = " << maxGenInEqLines_ << std::endl;

  solver_->updateSize(data_.nrVars_, maxEqLines_, maxInEqLines_, maxGenInEqLines_);
}

void QPSolver::nrVars(const std::vector<rbd::MultiBody> & mbs,
                      std::vector<UnilateralContact> uni,
                      std::vector<BilateralContact> bi)
{
  std::vector<std::tuple<int, int, double>> dependencies;
  data_.alphaD_.resize(mbs.size());
  data_.alphaDBegin_.resize(mbs.size());

  data_.uniCont_ = std::move(uni);
  data_.biCont_ = std::move(bi);

  int nrContacts = data_.nrContacts();

  data_.lambda_.resize(nrContacts);
  data_.lambdaBegin_.resize(nrContacts);

  data_.mobileRobotIndex_.clear();
  data_.normalAccB_.resize(mbs.size());

  int cumAlphaD = 0;
  for(std::size_t r = 0; r < mbs.size(); ++r)
  {
    const rbd::MultiBody & mb = mbs[r];
    data_.alphaD_[r] = mb.nrDof();
    data_.alphaDBegin_[r] = cumAlphaD;
    data_.normalAccB_[r].resize(mb.nrBodies(), sva::MotionVecd(Eigen::Vector6d::Zero()));
    cumAlphaD += mb.nrDof();
    if(mb.nrDof() > 0)
    {
      data_.mobileRobotIndex_.push_back(int(r));
      for(const auto & j : mb.joints())
      {
        if(j.isMimic())
        {
          dependencies.emplace_back(data_.alphaDBegin_[r] + mb.jointPosInDof(mb.jointIndexByName(j.mimicName())),
                                    data_.alphaDBegin_[r] + mb.jointPosInDof(mb.jointIndexByName(j.name())),
                                    j.mimicMultiplier());
        }
      }
    }
  }
  data_.totalAlphaD_ = cumAlphaD;

  int cumLambda = cumAlphaD;
  int cIndex = 0;
  data_.allCont_.clear();
  // counting unilateral contact
  for(const UnilateralContact & c : data_.uniCont_)
  {
    data_.lambdaBegin_[cIndex] = cumLambda;
    int lambda = 0;
    for(std::size_t p = 0; p < c.r1Points.size(); ++p)
    {
      lambda += c.nrLambda(int(p));
    }
    data_.lambda_[cIndex] = lambda;
    cumLambda += lambda;
    ++cIndex;

    data_.allCont_.emplace_back(c);
  }
  data_.nrUniLambda_ = cumLambda - cumAlphaD;

  // counting bilateral contact
  for(const BilateralContact & c : data_.biCont_)
  {
    data_.lambdaBegin_[cIndex] = cumLambda;
    int lambda = 0;
    for(std::size_t p = 0; p < c.r1Points.size(); ++p)
    {
      lambda += c.nrLambda(int(p));
    }
    data_.lambda_[cIndex] = lambda;
    cumLambda += lambda;
    ++cIndex;

    data_.allCont_.emplace_back(c);
  }
  data_.nrBiLambda_ = cumLambda - data_.nrUniLambda_ - cumAlphaD;

  data_.totalLambda_ = data_.nrUniLambda_ + data_.nrBiLambda_;
  data_.nrVars_ = data_.totalAlphaD_ + data_.totalLambda_;

  for(Task * t : tasks_)
  {
    t->updateNrVars(mbs, data_);
  }

  for(Constraint * c : constr_)
  {
    c->updateNrVars(mbs, data_);
  }

  solver_->setDependencies(data_.nrVars_, dependencies);
  solver_->updateSize(data_.nrVars_, maxEqLines_, maxInEqLines_, maxGenInEqLines_);
}

int QPSolver::nrVars() const
{
  return data_.nrVars_;
}

void QPSolver::updateTasksNrVars(const std::vector<rbd::MultiBody> & mbs) const
{
  for(Task * t : tasks_)
  {
    t->updateNrVars(mbs, data_);
  }
}

void QPSolver::updateConstrsNrVars(const std::vector<rbd::MultiBody> & mbs) const
{
  for(Constraint * c : constr_)
  {
    c->updateNrVars(mbs, data_);
  }
}

void QPSolver::updateNrVars(const std::vector<rbd::MultiBody> & mbs) const
{
  updateTasksNrVars(mbs);
  updateConstrsNrVars(mbs);
}

bool QPSolver::hasConstraint(const Constraint* co)
{
  return std::find(constr_.begin(), constr_.end(), co) != constr_.end(); 
}

/*
std::shared_ptr<tasks::qp::MotionConstr> QPSolver::getMotionConstr()
{
  std::shared_ptr<MotionConstr> motionConstr = NULL;
  
  for(Constraint* c: constr_)
  {
    std::shared_ptr<Constraint> c_prime = std::make_shared<Constraint>().reset(c);
    motionConstr = std::dynamic_pointer_cast<MotionConstr>(c_prime);
    
    if (motionConstr)
    {
      break;
    }
  }
}
*/

void QPSolver::addEqualityConstraint(Equality * co)
{
  eqConstr_.push_back(co);
}

void QPSolver::removeEqualityConstraint(Equality * co)
{
  eqConstr_.erase(std::find(eqConstr_.begin(), eqConstr_.end(), co));
}

int QPSolver::nrEqualityConstraints() const
{
  return static_cast<int>(eqConstr_.size());
}

void QPSolver::addInequalityConstraint(Inequality * co)
{
  inEqConstr_.push_back(co);
}

void QPSolver::removeInequalityConstraint(Inequality * co)
{
  inEqConstr_.erase(std::find(inEqConstr_.begin(), inEqConstr_.end(), co));
}

int QPSolver::nrInequalityConstraints() const
{
  return static_cast<int>(inEqConstr_.size());
}

void QPSolver::addGenInequalityConstraint(GenInequality * co)
{
  genInEqConstr_.push_back(co);
}

void QPSolver::removeGenInequalityConstraint(GenInequality * co)
{
  genInEqConstr_.erase(std::find(genInEqConstr_.begin(), genInEqConstr_.end(), co));
}

int QPSolver::nrGenInequalityConstraints() const
{
  return static_cast<int>(genInEqConstr_.size());
}

void QPSolver::addBoundConstraint(Bound * co)
{
  boundConstr_.push_back(co);
}

void QPSolver::removeBoundConstraint(Bound * co)
{
  boundConstr_.erase(std::find(boundConstr_.begin(), boundConstr_.end(), co));
}

int QPSolver::nrBoundConstraints() const
{
  return static_cast<int>(boundConstr_.size());
}

void QPSolver::addConstraint(Constraint * co)
{
  if(std::find(constr_.begin(), constr_.end(), co) == constr_.end())
  {
    constr_.push_back(co);
  }
}

void QPSolver::addConstraint(const std::vector<rbd::MultiBody> & mbs, Constraint * co)
{
  if(std::find(constr_.begin(), constr_.end(), co) == constr_.end())
  {
    constr_.push_back(co);
    // check if nrVars has been call at least one
    if(data_.nrVars_ > 0)
    {
      co->updateNrVars(mbs, data_);
    }
  }
}

void QPSolver::removeConstraint(Constraint * co)
{
  auto it = std::find(constr_.begin(), constr_.end(), co);
  if(it != constr_.end())
  {
    constr_.erase(it);
  }
}

int QPSolver::nrConstraints() const
{
  return static_cast<int>(constr_.size());
}

void QPSolver::addTask(Task * task)
{
  if(std::find(tasks_.begin(), tasks_.end(), task) == tasks_.end())
  {
    tasks_.push_back(task);
  }
}

void QPSolver::addTask(const std::vector<rbd::MultiBody> & mbs, Task * task)
{
  if(std::find(tasks_.begin(), tasks_.end(), task) == tasks_.end())
  {
    tasks_.push_back(task);
    // check if nrVars has been call at least one
    if(data_.nrVars_ > 0)
    {
      task->updateNrVars(mbs, data_);
    }
  }
}

void QPSolver::removeTask(Task * task)
{
  auto it = std::find(tasks_.begin(), tasks_.end(), task);
  if(it != tasks_.end())
  {
    tasks_.erase(it);
  }
}

int QPSolver::nrTasks() const
{
  return static_cast<int>(tasks_.size());
}

void QPSolver::solver(const std::string & name)
{
  solver_ = std::unique_ptr<GenQPSolver>(createQPSolver(name));
  solver_->updateSize(data_.nrVars_, maxEqLines_, maxInEqLines_, maxGenInEqLines_);
}

std::string QPSolver::solver() const
{
  return solver_->name();
}

void QPSolver::resetTasks()
{
  tasks_.clear();
}

const SolverData & QPSolver::data() const
{
  return data_;
}

SolverData & QPSolver::data()
{
  return data_;
}

const Eigen::VectorXd & QPSolver::result() const
{
  return solver_->result();
}

Eigen::VectorXd QPSolver::alphaDVec() const
{
  return solver_->result().head(data_.totalAlphaD_);
}

Eigen::VectorXd QPSolver::alphaDVec(int rIndex) const
{
  return solver_->result().segment(data_.alphaDBegin_[rIndex], data_.alphaD_[rIndex]);
}

Eigen::VectorXd QPSolver::lambdaVec() const
{
  // Rafa added this
  bool nantest = false;
  bool zerotest = false;
  //for (int i = 0; i < solver_->result().size(); i++) {
  //  nantest |= std::isnan(solver_->result()[i]);
  //}
  
  nantest = (isnan(solver_->result().array())).count() > 0;
  zerotest = solver_->result().segment(data_.lambdaBegin(), data_.totalLambda_).isZero(0);

  // std::cout << "Rafa, in QPSolver::lambdaVec, solver_->result() = " << solver_->result().transpose() << std::endl;
  
  if (nantest || zerotest)
  {
    if (zerotest)
      std::cout << "Rafa, in QPSolver::lambdaVec, zerotest" << std::endl;
    
    std::cout << "Rafa, in QPSolver::lambdaVec, solver_->result() = " << solver_->result().transpose() << std::endl;

    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->AFull_.size() = " << static_cast<LSSOLQPSolver*>(solver_.get())->AFull_.size() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->AL_.size() = " << static_cast<LSSOLQPSolver*>(solver_.get())->AL_.size() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->AU_.size() = " << static_cast<LSSOLQPSolver*>(solver_.get())->AU_.size() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->XLFull_.size() = " << static_cast<LSSOLQPSolver*>(solver_.get())->XLFull_.size() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->XUFull_.size() = " << static_cast<LSSOLQPSolver*>(solver_.get())->XUFull_.size() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->QFull_.size() = " << static_cast<LSSOLQPSolver*>(solver_.get())->QFull_.size() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->CFull_.size() = " << static_cast<LSSOLQPSolver*>(solver_.get())->CFull_.size() << std::endl;

    /*
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->AFull_ = " << static_cast<LSSOLQPSolver*>(solver_.get())->AFull_ << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->AL_ = " << static_cast<LSSOLQPSolver*>(solver_.get())->AL_.transpose() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->AU_ = " << static_cast<LSSOLQPSolver*>(solver_.get())->AU_.transpose() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->XLFull_ = " << static_cast<LSSOLQPSolver*>(solver_.get())->XLFull_.transpose() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->XUFull_ = " << static_cast<LSSOLQPSolver*>(solver_.get())->XUFull_.transpose() << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->QFull_ = " << static_cast<LSSOLQPSolver*>(solver_.get())->QFull_ << std::endl;
    std::cout << "Rafa, in QPSolver::lambdaVec, static_cast<LSSOLQPSolver*>(solver_.get())->CFull_ = " << static_cast<LSSOLQPSolver*>(solver_.get())->CFull_.transpose() << std::endl;
    */
    
    if (nantest)
      *(int*)0 = 0;
  }
  
  return solver_->result().segment(data_.lambdaBegin(), data_.totalLambda_);
}

Eigen::VectorXd QPSolver::lambdaVec(int cIndex) const
{
  return solver_->result().segment(data_.lambdaBegin_[cIndex], data_.lambda_[cIndex]);
}

int QPSolver::contactLambdaPosition(const ContactId & cId) const
{
  int pos = 0;

  for(const BilateralContact & bc : data_.allContacts())
  {
    if(bc.contactId == cId)
    {
      return pos;
    }

    for(std::size_t i = 0; i < bc.r1Points.size(); ++i)
    {
      pos += bc.nrLambda(int(i));
    }
  }

  return -1;
}

boost::timer::cpu_times QPSolver::solveTime() const
{
  return solverTimer_.elapsed();
}

boost::timer::cpu_times QPSolver::solveAndBuildTime() const
{
  return solverAndBuildTimer_.elapsed();
}

void QPSolver::preUpdate(const std::vector<rbd::MultiBody> & mbs, const std::vector<rbd::MultiBodyConfig> & mbcs)
{
  data_.computeNormalAccB(mbs, mbcs);
  for(std::size_t i = 0; i < constr_.size(); ++i)
  {
    constr_[i]->update(mbs, mbcs, data_);
  }

  for(std::size_t i = 0; i < tasks_.size(); ++i)
  {
    tasks_[i]->update(mbs, mbcs, data_);
  }

  solver_->updateMatrix(tasks_, eqConstr_, inEqConstr_, genInEqConstr_, boundConstr_);
}

void QPSolver::postUpdate(const std::vector<rbd::MultiBody> & /* mbs */,
                          std::vector<rbd::MultiBodyConfig> & mbcs,
                          bool success)
{
  if(success)
  {
    for(std::size_t r = 0; r < mbcs.size(); ++r)
    {
      updateMbc(mbcs[r], int(r));
    }

    // don't write contact force to the structure since contact force are used
    // to compute C vector.
  }
}

/**
 *													PassivityPIDTerm_QPSolver
 */

PassivityPIDTerm_QPSolver::PassivityPIDTerm_QPSolver()
: QPSolver()
{
}

bool PassivityPIDTerm_QPSolver::solve(const std::vector<rbd::MultiBody> & mbs,
                                      std::vector<rbd::MultiBodyConfig> & mbcs_real,
                                      std::vector<rbd::MultiBodyConfig> & mbcs_calc)
{
  bool success = solveNoMbcUpdate(mbs, mbcs_real, mbcs_calc);
  
  postUpdate(mbs, mbcs_real, success);
  postUpdate(mbs, mbcs_calc, success);
  
  return success;
}

bool PassivityPIDTerm_QPSolver::solveNoMbcUpdate(const std::vector<rbd::MultiBody> & mbs,
                                                 const std::vector<rbd::MultiBodyConfig> & mbcs_real,
                                                 const std::vector<rbd::MultiBodyConfig> & mbcs_calc)
{
  solverAndBuildTimer_.start();
  preUpdate(mbs, mbcs_real, mbcs_calc);
  
  solverTimer_.start();
  bool success = solver_->solve();
  solverTimer_.stop();
  
  data_.lambdaVecPrev_ = lambdaVec();
  
  if(!success)
  {
    solver_->errorMsg(mbs, tasks_, eqConstr_, inEqConstr_, genInEqConstr_, boundConstr_, std::cerr) << std::endl;
  }
  solverAndBuildTimer_.stop();
  
  return success;
}

void PassivityPIDTerm_QPSolver::preUpdate(const std::vector<rbd::MultiBody>& mbs,
                                          const std::vector<rbd::MultiBodyConfig>& mbcs_real,
                                          const std::vector<rbd::MultiBodyConfig>& mbcs_calc)
{
  data_.computeNormalAccB(mbs, mbcs_real);
  for(std::size_t i = 0; i < constr_.size(); ++i)
  {
    constr_[i]->update(mbs, mbcs_real, data_);
  }
  
  data_.computeNormalAccB(mbs, mbcs_calc);
  for(std::size_t i = 0; i < tasks_.size(); ++i)
  {
    tasks_[i]->update(mbs, mbcs_calc, data_);
  }
  
  solver_->updateMatrix(tasks_, eqConstr_, inEqConstr_, genInEqConstr_, boundConstr_);
}


} // namespace qp

} // namespace tasks
