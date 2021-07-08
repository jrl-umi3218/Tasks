/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// associated header
#include "Tasks/QPSolverData.h"

// includes
// RBDyn
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

namespace tasks
{

namespace qp
{

SolverData::SolverData()
: alphaD_(), alphaDBegin_(), lambda_(), totalAlphaD_(0), totalLambda_(0), nrUniLambda_(0), nrBiLambda_(0), nrVars_(0),
  uniCont_(), biCont_(), allCont_(), mobileRobotIndex_(), normalAccB_(), lambdaVecPrev_()
{
}

void SolverData::computeNormalAccB(const std::vector<rbd::MultiBody> & mbs,
                                   const std::vector<rbd::MultiBodyConfig> & mbcs)
{
  // we just need to update mobile robot normal acceleration
  for(int r : mobileRobotIndex_)
  {
    const rbd::MultiBody & mb = mbs[r];
    const rbd::MultiBodyConfig & mbc = mbcs[r];
    std::vector<sva::MotionVecd> & normalAccBr = normalAccB_[r];

    const std::vector<int> & pred = mb.predecessors();
    const std::vector<int> & succ = mb.successors();

    for(int i = 0; i < mb.nrJoints(); ++i)
    {
      const sva::PTransformd & X_p_i = mbc.parentToSon[i];
      const sva::MotionVecd & vj_i = mbc.jointVelocity[i];
      const sva::MotionVecd & vb_i = mbc.bodyVelB[i];

      if(pred[i] != -1)
        normalAccBr[succ[i]] = X_p_i * normalAccBr[pred[i]] + vb_i.cross(vj_i);
      else
        normalAccBr[succ[i]] = vb_i.cross(vj_i);
    }
  }
}

  
} // namespace qp

} // namespace tasks
