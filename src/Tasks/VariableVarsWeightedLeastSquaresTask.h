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

#pragma once

// includes
// Eigen
#include <Eigen/Core>
#include <Eigen/StdVector>

// RBDyn
//#include <RBDyn/Jacobian.h>
//#include <RBDyn/CoM.h>

// sch
//#include <sch/Matrix/SCH_Types.h>

// Tasks
#include "QPSolver.h"

// unique_ptr
#include <memory>


namespace tasks
{
namespace qp
{

class TASKS_DLLAPI VariableVarsWeightedLeastSquaresTask : public Task
{
private:
    //const std::vector<rbd::MultiBody>& mbs_;
    Eigen::VectorXd dimWeight_;
    int robotIndex_;
    int variableVarsPos, variableVarsBegin_;
    
    Eigen::MatrixXd Q_; // < diagonal dimWeight_
    Eigen::VectorXd C_; // < zero
    
public:
    VariableVarsWeightedLeastSquaresTask(const std::vector<rbd::MultiBody>& mbs, int robotIndex, double weight)
    : Task(weight),
    robotIndex_(robotIndex)
    {
    }
    
    virtual ~VariableVarsWeightedLeastSquaresTask(){}
    
    virtual void addNrVariableVars( SolverData& data ) override
    {
        variableVarsPos = data.variableVars()[robotIndex_];
        data.variableVars()[robotIndex_] += dimWeight_.rows();
    }

    virtual std::pair<int, int> begin() const override
    {
        return std::make_pair(variableVarsBegin_, variableVarsBegin_);
    }
    
    void dimWeight(const Eigen::VectorXd& dim)
    {
        dimWeight_ = dim;
    }
    
    const Eigen::VectorXd& dimWeight() const
    {
        return dimWeight_;
    }
    
    virtual void update(const std::vector<rbd::MultiBody>& mbs,
                        const std::vector<rbd::MultiBodyConfig>& mbcs,
                        const SolverData& data) override
    {
        // don't need to do this every solver update.
        Q_ = dimWeight_.asDiagonal();
        C_ = Eigen::VectorXd::Zero( dimWeight_.rows() );
    }
    
    virtual void updateNrVars(const std::vector<rbd::MultiBody>& mbs,
                              const SolverData& data) override
    {
        variableVarsBegin_ = data.variableVarsBegin(robotIndex_) + variableVarsPos;
    }
    
    virtual const Eigen::MatrixXd& Q() const override { return Q_; }
    virtual const Eigen::VectorXd& C() const override { return C_; }
};
   

} // namespace qp

} // namespace tasks

