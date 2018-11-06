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


// associated header
#include "Tasks/Tasks.h"

// includes
// std
#include <numeric>
#include <set>

// rbd
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

namespace tasks
{


/**
	*													PositionTask
	*/


PositionTask::PositionTask(const rbd::MultiBody& mb, const std::string& bodyName,
	const Eigen::Vector3d& pos, const Eigen::Vector3d& bodyPoint):
	pos_(pos),
	point_(bodyPoint),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName, bodyPoint),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


void PositionTask::position(const Eigen::Vector3d& pos)
{
	pos_ = pos;
}


const Eigen::Vector3d& PositionTask::position() const
{
	return pos_;
}


void PositionTask::bodyPoint(const Eigen::Vector3d& point)
{
	jac_.point(point);
}


const Eigen::Vector3d& PositionTask::bodyPoint() const
{
	return jac_.point();
}


void PositionTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = pos_ - (point_*mbc.bodyPosW[bodyIndex_]).translation();
	speed_ = jac_.velocity(mb, mbc).linear();
	normalAcc_ = jac_.normalAcceleration(mb, mbc).linear();

	const auto& shortJacMat =
		jac_.jacobian(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void PositionTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = pos_ - (point_*mbc.bodyPosW[bodyIndex_]).translation();
	speed_ = jac_.velocity(mb, mbc).linear();
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB).linear();

	const auto& shortJacMat =
		jac_.jacobian(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void PositionTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	const auto& shortJacMat =
		jac_.jacobianDot(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacDotMat_);
}


const Eigen::VectorXd& PositionTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& PositionTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& PositionTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& PositionTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& PositionTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													OrientationTask
	*/


OrientationTask::OrientationTask(const rbd::MultiBody& mb, const std::string& bodyName, const Eigen::Quaterniond& ori):
	ori_(ori.matrix()),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


OrientationTask::OrientationTask(const rbd::MultiBody& mb, const std::string& bodyName, const Eigen::Matrix3d& ori):
	ori_(ori),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


void OrientationTask::orientation(const Eigen::Quaterniond& ori)
{
	ori_ = ori.matrix();
}


void OrientationTask::orientation(const Eigen::Matrix3d& ori)
{
	ori_ = ori;
}


const Eigen::Matrix3d& OrientationTask::orientation() const
{
	return ori_;
}


void OrientationTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = sva::rotationError(mbc.bodyPosW[bodyIndex_].rotation(), ori_);
	speed_ = jac_.velocity(mb, mbc).angular();
	normalAcc_ = jac_.normalAcceleration(mb, mbc).angular();

	const auto& shortJacMat = jac_.jacobian(mb, mbc).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void OrientationTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = sva::rotationError(mbc.bodyPosW[bodyIndex_].rotation(), ori_);
	speed_ = jac_.velocity(mb, mbc).angular();
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB).angular();

	const auto& shortJacMat = jac_.jacobian(mb, mbc).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void OrientationTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	const auto& shortJacMat = jac_.jacobianDot(mb, mbc).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacDotMat_);
}


const Eigen::VectorXd& OrientationTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& OrientationTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& OrientationTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& OrientationTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& OrientationTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													TransformTaskCommon
	*/


TransformTaskCommon::TransformTaskCommon(const rbd::MultiBody& mb, const std::string& bodyName,
		const sva::PTransformd& X_0_t, const sva::PTransformd& X_b_p):
	X_0_t_(X_0_t),
	X_b_p_(X_b_p),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	eval_(6),
	speed_(6),
	normalAcc_(6),
	jacMat_(6, mb.nrDof())
{
}


void TransformTaskCommon::target(const sva::PTransformd& X_0_t)
{
	X_0_t_ = X_0_t;
}


const sva::PTransformd& TransformTaskCommon::target() const
{
	return X_0_t_;
}


void TransformTaskCommon::X_b_p(const sva::PTransformd& X_b_p)
{
	X_b_p_ = X_b_p;
}


const sva::PTransformd& TransformTaskCommon::X_b_p() const
{
	return X_b_p_;
}


const Eigen::VectorXd& TransformTaskCommon::eval() const
{
	return eval_;
}


const Eigen::VectorXd& TransformTaskCommon::speed() const
{
	return speed_;
}


const Eigen::VectorXd& TransformTaskCommon::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& TransformTaskCommon::jac() const
{
	return jacMat_;
}


/**
	*													SurfaceTransformTask
	*/


SurfaceTransformTask::SurfaceTransformTask(const rbd::MultiBody& mb, const std::string& bodyName,
		const sva::PTransformd& X_0_t, const sva::PTransformd& X_b_p):
	TransformTaskCommon(mb, bodyName, X_0_t, X_b_p),
	jacMatTmp_(6, jac_.dof())
{
}


void SurfaceTransformTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const std::vector<sva::MotionVecd>& normalAccB)
{
	sva::PTransformd X_0_p = X_b_p_*mbc.bodyPosW[bodyIndex_];
	sva::PTransformd X_p_t = X_0_t_*X_0_p.inv();

	// see Section 4.2.6 of Joris Vaillant's PhD thesis (French) for details
	sva::MotionVecd err_p = sva::transformVelocity(X_p_t);
	sva::MotionVecd V_p_p = jac_.velocity(mb, mbc, X_b_p_);
	sva::MotionVecd w_p_p = sva::MotionVecd(V_p_p.angular(), Eigen::Vector3d::Zero());
	sva::MotionVecd AN_p_p = jac_.normalAcceleration(mb, mbc, normalAccB, X_b_p_,
		sva::MotionVecd(Eigen::Vector6d::Zero()));
	sva::MotionVecd wAN_p_p = sva::MotionVecd(AN_p_p.angular(), Eigen::Vector3d::Zero());
	sva::MotionVecd V_err_p = err_p.cross(w_p_p) - V_p_p;

	eval_ = err_p.vector();
	speed_ = -V_err_p.vector();
	normalAcc_ = -(V_err_p.cross(w_p_p) + err_p.cross(wAN_p_p) - AN_p_p).vector();

	jacMatTmp_ = jac_.jacobian(mb, mbc, X_0_p);

	for(int i = 0; i < jac_.dof(); ++i)
	{
		jacMatTmp_.col(i).head<6>() -= err_p.cross(
			sva::MotionVecd(jacMatTmp_.col(i).head<3>(), Eigen::Vector3d::Zero())).vector();
	}

	jac_.fullJacobian(mb, jacMatTmp_, jacMat_);
}


/**
	*													TransformTask
	*/


TransformTask::TransformTask(const rbd::MultiBody& mb, const std::string& bodyName,
		const sva::PTransformd& X_0_t, const sva::PTransformd& X_b_p,
		const Eigen::Matrix3d& E_0_c):
	TransformTaskCommon(mb, bodyName, X_0_t, X_b_p),
	E_0_c_(E_0_c)
{
}


void TransformTask::E_0_c(const Eigen::Matrix3d& E_0_c)
{
	E_0_c_ = E_0_c;
}


const Eigen::Matrix3d& TransformTask::E_0_c() const
{
	return E_0_c_;
}


void TransformTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const std::vector<sva::MotionVecd>& normalAccB)
{
	sva::PTransformd X_0_p(X_b_p_*mbc.bodyPosW[bodyIndex_]);
	sva::PTransformd E_p_c(Eigen::Matrix3d(E_0_c_*X_0_p.rotation().transpose()));
	sva::PTransformd X_b_p_c(E_p_c*X_b_p_);
	sva::MotionVecd V_p_c(jac_.velocity(mb, mbc, X_b_p_c));
	sva::MotionVecd w_p_c(V_p_c.angular(), Eigen::Vector3d::Zero());

	eval_ = (sva::PTransformd(E_0_c_)*sva::transformError(X_0_p, X_0_t_)).vector();
	speed_ = V_p_c.vector();
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB, X_b_p_c, w_p_c).vector();
	const auto& shortJacMat = jac_.jacobian(mb, mbc, E_p_c*X_0_p);

	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


/**
	*													MultiRobotTransformTask
	*/


MultiRobotTransformTask::MultiRobotTransformTask(
	const std::vector<rbd::MultiBody>& mbs,
	int r1Index, int r2Index,
	const std::string& r1BodyName, const std::string& r2BodyName,
	const sva::PTransformd& X_r1b_r1s, const sva::PTransformd& X_r2b_r2s):
	r1Index_(r1Index),
	r2Index_(r2Index),
	r1BodyIndex_(mbs[r1Index].bodyIndexByName(r1BodyName)),
	r2BodyIndex_(mbs[r2Index].bodyIndexByName(r2BodyName)),
	X_r1b_r1s_(X_r1b_r1s),
	X_r2b_r2s_(X_r2b_r2s),
	jacR1B_(mbs[r1Index], r1BodyName),
	jacR2B_(mbs[r2Index], r2BodyName),
	eval_(6),
	speed_(6),
	normalAcc_(6),
	jacMat1_(6, jacR1B_.dof()),
	jacMat2_(6, jacR2B_.dof()),
	fullJacMat_({Eigen::MatrixXd::Zero(6,mbs[r1Index].nrDof()),
							Eigen::MatrixXd::Zero(6,mbs[r2Index].nrDof())})
{}


int MultiRobotTransformTask::r1Index() const
{
	return r1Index_;
}


int MultiRobotTransformTask::r2Index() const
{
	return r2Index_;
}


void MultiRobotTransformTask::X_r1b_r1s(const sva::PTransformd& X_r1b_r1s)
{
	X_r1b_r1s_ = X_r1b_r1s;
}


const sva::PTransformd& MultiRobotTransformTask::X_r1b_r1s() const
{
	return X_r1b_r1s_;
}


void MultiRobotTransformTask::X_r2b_r2s(const sva::PTransformd& X_r2b_r2s)
{
	X_r2b_r2s_ = X_r2b_r2s;
}


const sva::PTransformd& MultiRobotTransformTask::X_r2b_r2s() const
{
	return X_r2b_r2s_;
}


void MultiRobotTransformTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const std::vector<std::vector<sva::MotionVecd>>& normalAccB)
{
	using namespace Eigen;

	const rbd::MultiBody& mb1 = mbs[r1Index_];
	const rbd::MultiBody& mb2 = mbs[r2Index_];
	const rbd::MultiBodyConfig& mbc1 = mbcs[r1Index_];
	const rbd::MultiBodyConfig& mbc2 = mbcs[r2Index_];
	const sva::PTransformd& X_0_r1b = mbc1.bodyPosW[r1BodyIndex_];
	const sva::PTransformd& X_0_r2b = mbc2.bodyPosW[r2BodyIndex_];
	const std::vector<sva::MotionVecd>& normalAccBR1 = normalAccB[r1Index_];
	const std::vector<sva::MotionVecd>& normalAccBR2 = normalAccB[r2Index_];

	sva::PTransformd X_0_r1s = X_r1b_r1s_*X_0_r1b;
	sva::PTransformd X_0_r2s = X_r2b_r2s_*X_0_r2b;
	sva::PTransformd X_r1s_r2s = X_0_r2s*X_0_r1s.inv();

	sva::PTransformd E_r2s_r1s(Matrix3d(X_r1s_r2s.rotation().transpose()));
	sva::PTransformd X_r2b_r2s_r1s(E_r2s_r1s*X_r2b_r2s_);

	// see Section 4.2.6 of Joris Vaillant's PhD thesis (French) for details
	sva::MotionVecd err_r1s(sva::transformVelocity(X_r1s_r2s));

	sva::MotionVecd V_r1s_r1s = jacR1B_.velocity(mb1, mbc1, X_r1b_r1s_);
	sva::MotionVecd V_r2s_r1s = jacR2B_.velocity(mb2, mbc2, X_r2b_r2s_r1s);

	sva::MotionVecd V_err_s = V_r2s_r1s - V_r1s_r1s;
	sva::MotionVecd w_r1s(V_r1s_r1s.angular(), Vector3d::Zero());

	sva::MotionVecd V_err_r1s = err_r1s.cross(w_r1s) + V_err_s;

	sva::MotionVecd AN_r1s_r1s = jacR1B_.normalAcceleration(mb1, mbc1,
		normalAccBR1, X_r1b_r1s_, sva::MotionVecd(Vector6d::Zero()));
	sva::MotionVecd wAN_r1s_r1s(AN_r1s_r1s.angular(), Vector3d::Zero());
	sva::MotionVecd AN_r2s_r1s = jacR2B_.normalAcceleration(mb2, mbc2, normalAccBR2,
		X_r2b_r2s_r1s, sva::MotionVecd(V_err_s.angular(), Vector3d::Zero()));
	sva::MotionVecd AN_err_s = AN_r2s_r1s - AN_r1s_r1s;

	sva::MotionVecd AN_err_r1s = V_err_r1s.cross(w_r1s) +
		err_r1s.cross(wAN_r1s_r1s) + AN_err_s;

	eval_ = err_r1s.vector();
	speed_ = -V_err_r1s.vector();
	normalAcc_ = -AN_err_r1s.vector();

	jacMat1_.noalias() = jacR1B_.jacobian(mb1, mbc1, X_0_r1s);
	for(int i = 0; i < jacR1B_.dof(); ++i)
	{
		jacMat1_.col(i).head<6>() -= err_r1s.cross(
			sva::MotionVecd(jacMat1_.col(i).head<3>(), Vector3d::Zero())).vector();
	}
	jacMat2_.noalias() = -jacR2B_.jacobian(mb2, mbc2, E_r2s_r1s*X_0_r2s);

	jacR1B_.fullJacobian(mb1, jacMat1_, fullJacMat_[0]);
	jacR2B_.fullJacobian(mb2, jacMat2_, fullJacMat_[1]);
}


const Eigen::VectorXd& MultiRobotTransformTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& MultiRobotTransformTask::speed() const
{
	return speed_;
}

const Eigen::VectorXd& MultiRobotTransformTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& MultiRobotTransformTask::jac(int index) const
{
	return fullJacMat_[index];
}


/**
	*													SurfaceOrientationTask
	*/


SurfaceOrientationTask::SurfaceOrientationTask(const rbd::MultiBody& mb,
	const std::string& bodyName, const Eigen::Quaterniond& ori,
	const sva::PTransformd& X_b_s):
	ori_(ori.matrix()),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	X_b_s_(X_b_s),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


SurfaceOrientationTask::SurfaceOrientationTask(const rbd::MultiBody& mb,
	const std::string& bodyName, const Eigen::Matrix3d& ori, const sva::PTransformd& X_b_s):
	ori_(ori),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	X_b_s_(X_b_s),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
}


void SurfaceOrientationTask::orientation(const Eigen::Quaterniond& ori)
{
	ori_ = ori.matrix();
}


void SurfaceOrientationTask::orientation(const Eigen::Matrix3d& ori)
{
	ori_ = ori;
}


const Eigen::Matrix3d& SurfaceOrientationTask::orientation() const
{
	return ori_;
}


void SurfaceOrientationTask::update(const rbd::MultiBody& mb,
																	const rbd::MultiBodyConfig& mbc)
{
	eval_ = sva::rotationVelocity<double>
		(ori_*mbc.bodyPosW[bodyIndex_].rotation().transpose()*X_b_s_.rotation().transpose());
	speed_ = jac_.velocity(mb, mbc, X_b_s_).angular();
	// since X_b_s is constant, the X_b_s velocity
	// (last argument of normalAccelation) is a 0 velocity vector
	normalAcc_ = jac_.normalAcceleration(mb, mbc, X_b_s_,
		sva::MotionVecd(Eigen::Vector6d::Zero())).angular();

	const auto& shortJacMat =
		jac_.jacobian(mb, mbc, X_b_s_*mbc.bodyPosW[bodyIndex_]).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void SurfaceOrientationTask::update(const rbd::MultiBody& mb,
																	const rbd::MultiBodyConfig& mbc,
																	const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = sva::rotationVelocity<double>
		(ori_*mbc.bodyPosW[bodyIndex_].rotation().transpose()*X_b_s_.rotation().transpose());
	speed_ = jac_.velocity(mb, mbc, X_b_s_).angular();
	// since X_b_s is constant, the X_b_s velocity
	// (third argument of normalAccelation) is a 0 velocity vector
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB, X_b_s_,
		sva::MotionVecd(Eigen::Vector6d::Zero())).angular();

	const auto& shortJacMat =
		jac_.jacobian(mb, mbc, X_b_s_*mbc.bodyPosW[bodyIndex_]).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void SurfaceOrientationTask::updateDot(const rbd::MultiBody& mb,
																		 const rbd::MultiBodyConfig& mbc)
{
	const auto& shortJacMat = jac_.bodyJacobianDot(mb, mbc).block(0, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacDotMat_);
}


const Eigen::VectorXd& SurfaceOrientationTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& SurfaceOrientationTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& SurfaceOrientationTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& SurfaceOrientationTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& SurfaceOrientationTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													GazeTask
	*/


GazeTask::GazeTask(const rbd::MultiBody &mb, const std::string& bodyName,
	const Eigen::Vector2d& point2d, double depthEstimate, const sva::PTransformd& X_b_gaze,
	const Eigen::Vector2d& point2d_ref):
	point2d_(new Eigen::Vector2d(point2d)),
	point2d_ref_(new Eigen::Vector2d(point2d_ref)),
	depthEstimate_(depthEstimate),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	X_b_gaze_(X_b_gaze),
	L_img_(new Eigen::Matrix<double, 2, 6>(Eigen::Matrix<double, 2, 6>::Zero())),
	surfaceVelocity_(new Eigen::Matrix<double,6,1>(Eigen::Matrix<double,6,1>::Zero())),
	L_Z_dot_(new Eigen::Matrix<double,1,6>(Eigen::Matrix<double,1,6>::Zero())),
	L_img_dot_(new Eigen::Matrix<double,2,6>(Eigen::Matrix<double,2,6>::Zero())),
	eval_(2),
	speed_(2),
	normalAcc_(2),
	jacMat_(2, mb.nrDof()),
	jacDotMat_(2, mb.nrDof())
{
}


GazeTask::GazeTask(const rbd::MultiBody &mb, const std::string& bodyName,
	const Eigen::Vector3d &point3d, const sva::PTransformd& X_b_gaze,
	const Eigen::Vector2d& point2d_ref):
	point2d_(new Eigen::Vector2d()),
	point2d_ref_(new Eigen::Vector2d(point2d_ref)),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	X_b_gaze_(X_b_gaze),
	L_img_(new Eigen::Matrix<double, 2, 6>(Eigen::Matrix<double, 2, 6>::Zero())),
	surfaceVelocity_(new Eigen::Matrix<double,6,1>(Eigen::Matrix<double,6,1>::Zero())),
	L_Z_dot_(new Eigen::Matrix<double,1,6>(Eigen::Matrix<double,1,6>::Zero())),
	L_img_dot_(new Eigen::Matrix<double,2,6>(Eigen::Matrix<double,2,6>::Zero())),
	eval_(2),
	speed_(2),
	normalAcc_(2),
	jacMat_(2, mb.nrDof()),
	jacDotMat_(2, mb.nrDof())
{
	*point2d_ << point3d[0]/point3d[2], point3d[1]/point3d[2];
	depthEstimate_ = point3d[2];
}

GazeTask::GazeTask(const GazeTask& rhs):
	point2d_(new Eigen::Vector2d(*rhs.point2d_)),
	point2d_ref_(new Eigen::Vector2d(*rhs.point2d_ref_)),
	depthEstimate_(rhs.depthEstimate_),
	bodyIndex_(rhs.bodyIndex_),
	jac_(rhs.jac_),
	X_b_gaze_(rhs.X_b_gaze_),
	L_img_(new Eigen::Matrix<double, 2, 6>(*rhs.L_img_)),
	surfaceVelocity_(new Eigen::Matrix<double, 6, 1>(*rhs.surfaceVelocity_)),
	L_Z_dot_(new Eigen::Matrix<double, 1, 6>(*rhs.L_Z_dot_)),
	L_img_dot_(new Eigen::Matrix<double, 2, 6>(*rhs.L_img_dot_)),
	eval_(rhs.eval_), speed_(rhs.speed_),
	normalAcc_(rhs.normalAcc_), jacMat_(rhs.jacMat_),
	jacDotMat_(rhs.jacDotMat_)
{
}

GazeTask& GazeTask::operator=(const GazeTask& rhs)
{
	if(&rhs != this)
	{
		*point2d_ = *rhs.point2d_;
		*point2d_ref_ = *rhs.point2d_ref_;
		depthEstimate_ = rhs.depthEstimate_;
		bodyIndex_ = rhs.bodyIndex_;
		jac_ = rhs.jac_;
		X_b_gaze_ = rhs.X_b_gaze_;
		*L_img_ = *rhs.L_img_;
		*surfaceVelocity_ = *rhs.surfaceVelocity_;
		*L_Z_dot_ = *rhs.L_Z_dot_;
		*L_img_dot_ = *rhs.L_img_dot_;
		eval_ = rhs.eval_;
		speed_ = rhs.speed_;
		normalAcc_ = rhs.normalAcc_;
		jacMat_ = rhs.jacMat_;
		jacDotMat_ = rhs.jacDotMat_;
	}
	return *this;
}

void GazeTask::error(const Eigen::Vector2d& point2d, const Eigen::Vector2d& point2d_ref)
{
	*point2d_ = point2d;
	*point2d_ref_ = point2d_ref;
}


void GazeTask::error(const Eigen::Vector3d& point3d, const Eigen::Vector2d& point2d_ref)
{
	*point2d_ << point3d[0]/point3d[2], point3d[1]/point3d[2];
	depthEstimate_ = point3d[2];
	*point2d_ref_ = point2d_ref;
}


void GazeTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	const std::vector<sva::MotionVecd>& normalAccB)
{
	// compute eval term
	eval_ = *point2d_ref_ - *point2d_;

	// compute speed term
	rbd::imagePointJacobian(*point2d_, depthEstimate_, *L_img_);
	*surfaceVelocity_ = (jac_.velocity(mb, mbc, X_b_gaze_)).vector();
	speed_ = (*L_img_)*(*surfaceVelocity_);

	// compute norm accel term
	rbd::depthDotJacobian(speed_, depthEstimate_, *L_Z_dot_);
	rbd::imagePointJacobianDot(*point2d_, speed_, depthEstimate_, (*L_Z_dot_)*(*surfaceVelocity_), *L_img_dot_);
	normalAcc_ = (*L_img_)*(jac_.normalAcceleration(mb, mbc, normalAccB, X_b_gaze_,
		sva::MotionVecd(Eigen::Vector6d::Zero()))).vector() + (*L_img_dot_)*(*surfaceVelocity_);

	// compute the task Jacobian
	Eigen::MatrixXd shortJacMat = (*L_img_)*jac_.jacobian(mb, mbc, X_b_gaze_*mbc.bodyPosW[bodyIndex_]).block(0, 0, 6, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


const Eigen::VectorXd& GazeTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& GazeTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& GazeTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& GazeTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& GazeTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													PositionBasedVisServoTask
	*/


PositionBasedVisServoTask::PositionBasedVisServoTask(const rbd::MultiBody &mb,
	const std::string& bodyName, const sva::PTransformd& X_t_s,
	const sva::PTransformd& X_b_s):
	X_t_s_(X_t_s),
	X_b_s_(X_b_s),
	angle_(0.),
	axis_(Eigen::Vector3d::Zero()),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	L_pbvs_(new Eigen::Matrix<double, 6, 6>(Eigen::Matrix<double, 6, 6>::Zero())),
	surfaceVelocity_(new Eigen::Matrix<double, 6, 1>(Eigen::Matrix<double, 6, 1>::Zero())),
	omegaSkew_(Eigen::Matrix3d::Zero()),
	L_pbvs_dot_(new Eigen::Matrix<double, 6, 6>(Eigen::Matrix<double, 6, 6>::Zero())),
	eval_(6),
	speed_(6),
	normalAcc_(6),
	jacMat_(6, mb.nrDof()),
	jacDotMat_(6, mb.nrDof())
{
}

PositionBasedVisServoTask::PositionBasedVisServoTask(const PositionBasedVisServoTask& rhs):
	X_t_s_(rhs.X_t_s_),
	X_b_s_(rhs.X_b_s_),
	angle_(rhs.angle_),
	axis_(rhs.axis_),
	bodyIndex_(rhs.bodyIndex_),
	jac_(rhs.jac_),
	L_pbvs_(new Eigen::Matrix<double, 6, 6>(*rhs.L_pbvs_)),
	surfaceVelocity_(new Eigen::Matrix<double, 6, 1>(*rhs.surfaceVelocity_)),
	omegaSkew_(rhs.omegaSkew_),
	L_pbvs_dot_(new Eigen::Matrix<double, 6, 6>(*rhs.L_pbvs_dot_)),
	eval_(rhs.eval_), speed_(rhs.speed_), normalAcc_(rhs.normalAcc_),
	jacMat_(rhs.jacMat_), jacDotMat_(rhs.jacDotMat_)
{
}

PositionBasedVisServoTask& PositionBasedVisServoTask::operator=(const PositionBasedVisServoTask& rhs)
{
	if(&rhs != this)
	{
		X_t_s_ = rhs.X_t_s_;
		X_b_s_ = rhs.X_b_s_;
		angle_ = rhs.angle_;
		axis_ = rhs.axis_;
		bodyIndex_ = rhs.bodyIndex_;
		jac_ = rhs.jac_;
		*L_pbvs_ = *rhs.L_pbvs_;
		*surfaceVelocity_ = *rhs.surfaceVelocity_;
		omegaSkew_ = rhs.omegaSkew_;
		*L_pbvs_dot_ = *rhs.L_pbvs_dot_;
		eval_ = rhs.eval_;
		speed_ = rhs.speed_;
		normalAcc_ = rhs.normalAcc_;
		jacMat_ = rhs.jacMat_;
		jacDotMat_ = rhs.jacDotMat_;
	}
	return *this;
}

void PositionBasedVisServoTask::error(const sva::PTransformd& X_t_s)
{
	X_t_s_ = X_t_s;
}

void PositionBasedVisServoTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	const std::vector<sva::MotionVecd>& normalAccB)
{
	// compute eval term
	rbd::getAngleAxis(X_t_s_.rotation().transpose(), angle_, axis_);
	eval_.tail(3) = -X_t_s_.translation();
	eval_.head(3) = angle_*axis_;

	// compute speed term
	rbd::poseJacobian(X_t_s_.rotation(), *L_pbvs_);
	*surfaceVelocity_ = (jac_.velocity(mb, mbc, X_b_s_)).vector();
	speed_ = (*L_pbvs_)*(*surfaceVelocity_);

	// compute norm accel term
	rbd::getSkewSym(surfaceVelocity_->head(3), omegaSkew_);
	*L_pbvs_dot_ << Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero(),
								 Eigen::Matrix3d::Zero(), -X_t_s_.rotation().transpose()*omegaSkew_;
	normalAcc_ = (*L_pbvs_)*(jac_.normalAcceleration(mb, mbc, normalAccB, X_b_s_,
		sva::MotionVecd(Eigen::Vector6d::Zero()))).vector() + (*L_pbvs_dot_)*(*surfaceVelocity_);

	// compute the task Jacobian
	Eigen::MatrixXd shortJacMat = (*L_pbvs_)*jac_.jacobian(mb, mbc, X_b_s_*mbc.bodyPosW[bodyIndex_]).block(0, 0, 6, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


const Eigen::VectorXd& PositionBasedVisServoTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& PositionBasedVisServoTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& PositionBasedVisServoTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& PositionBasedVisServoTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& PositionBasedVisServoTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													PostureTask
	*/


PostureTask::PostureTask(const rbd::MultiBody& mb, std::vector<std::vector<double> > q):
	q_(q),
	eval_(mb.nrDof()),
	jacMat_(mb.nrDof(), mb.nrDof()),
	jacDotMat_(mb.nrDof(), mb.nrDof())
{
	eval_.setZero();
	jacMat_.setIdentity();
	jacDotMat_.setZero();

	if(mb.nrDof() > 0 && mb.joint(0).type() == rbd::Joint::Free)
	{
		for(int i = 0; i < 6; ++i)
		{
			jacMat_(i, i) = 0;
		}
	}
}


void PostureTask::posture(std::vector<std::vector<double> > q)
{
	q_ = q;
}


const std::vector<std::vector<double> > PostureTask::posture() const
{
	return q_;
}


void PostureTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;

	int pos = mb.jointPosInDof(1);

	// we drop the first joint (fixed or free flyier).
	for(int i = 1; i < mb.nrJoints(); ++i)
	{
		// if dof == 1 is a prismatic/revolute joint
		// else if dof == 4 is a spherical one
		// else is a fixed one
		if(mb.joint(i).dof() == 1)
		{
			eval_(pos) = q_[i][0] - mbc.q[i][0];
			++pos;
		}
		else if(mb.joint(i).dof() == 4)
		{
			Matrix3d orid(
				Quaterniond(q_[i][0], q_[i][1], q_[i][2], q_[i][3]).matrix());

			Vector3d err = sva::rotationError(mbc.jointConfig[i].rotation(), orid);

			eval_.segment(pos, 3) = err;
			pos += 3;
		}
	}
}


void PostureTask::updateDot(const rbd::MultiBody& /* mb */, const rbd::MultiBodyConfig& /* mbc */)
{}


const Eigen::VectorXd& PostureTask::eval() const
{
	return eval_;
}


const Eigen::MatrixXd& PostureTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& PostureTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													CoMTask
	*/


CoMTask::CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com):
	com_(com),
	jac_(mb),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{}


CoMTask::CoMTask(const rbd::MultiBody& mb, const Eigen::Vector3d& com,
								std::vector<double> weight):
	com_(com),
	jac_(mb, std::move(weight)),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{}


void CoMTask::com(const Eigen::Vector3d& com)
{
	com_ = com;
}


const Eigen::Vector3d CoMTask::com() const
{
	return com_;
}


void CoMTask::updateInertialParameters(const rbd::MultiBody& mb)
{
	jac_.updateInertialParameters(mb);
}


void CoMTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = com_ - rbd::computeCoM(mb, mbc);

	speed_ = jac_.velocity(mb, mbc);
	normalAcc_ = jac_.normalAcceleration(mb, mbc);
	jacMat_ = jac_.jacobian(mb, mbc);
}


void CoMTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const Eigen::Vector3d& com,
	const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = com_ - com;

	speed_ = jac_.velocity(mb, mbc);
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB);
	jacMat_ = jac_.jacobian(mb, mbc);
}


void CoMTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	jacDotMat_ = jac_.jacobianDot(mb, mbc);
}


const Eigen::VectorXd& CoMTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& CoMTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& CoMTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& CoMTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& CoMTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													MultiCoMTask
	*/


MultiCoMTask::MultiCoMTask(const std::vector<rbd::MultiBody>& mbs,
	std::vector<int> robotIndexes, const Eigen::Vector3d& com):
	com_(com),
	robotIndexes_(std::move(robotIndexes)),
	robotsWeight_(),
	jac_(robotIndexes_.size()),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(robotIndexes_.size())
{
	computeRobotsWeight(mbs);

	// create CoMJacobian and jacobian matrix
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];

		jac_[i] = rbd::CoMJacobian(mb,
			std::vector<double>(mb.nrBodies(), robotsWeight_[i]));
		jacMat_[i].resize(3, mb.nrDof());
	}
}


void MultiCoMTask::com(const Eigen::Vector3d& com)
{
	com_ = com;
}


const Eigen::Vector3d MultiCoMTask::com() const
{
	return com_;
}


const std::vector<int>& MultiCoMTask::robotIndexes() const
{
	return robotIndexes_;
}


void MultiCoMTask::updateInertialParameters(const std::vector<rbd::MultiBody>& mbs)
{
	computeRobotsWeight(mbs);

	// upadte CoMJacobian per body weight
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];

		jac_[i].weight(mb, std::vector<double>(mb.nrBodies(), robotsWeight_[i]));
	}
}


void MultiCoMTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs)
{
	eval_ = com_;
	speed_.setZero();
	normalAcc_.setZero();
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];
		const rbd::MultiBodyConfig& mbc = mbcs[r];

		eval_ -= rbd::computeCoM(mbs[r], mbcs[r])*robotsWeight_[i];
		speed_ += jac_[i].velocity(mb, mbc);
		normalAcc_ += jac_[i].normalAcceleration(mb, mbc);
		jacMat_[i] = jac_[i].jacobian(mb, mbc);
	}
}


void MultiCoMTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const std::vector<std::vector<sva::MotionVecd>>& normalAccB)
{
	eval_ = com_;
	speed_.setZero();
	normalAcc_.setZero();
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];
		const rbd::MultiBodyConfig& mbc = mbcs[r];

		eval_ -= rbd::computeCoM(mbs[r], mbcs[r])*robotsWeight_[i];
		speed_ += jac_[i].velocity(mb, mbc);
		normalAcc_ += jac_[i].normalAcceleration(mb, mbc, normalAccB[r]);
		jacMat_[i] = jac_[i].jacobian(mb, mbc);
	}
}


void MultiCoMTask::update(const std::vector<rbd::MultiBody>& mbs,
	const std::vector<rbd::MultiBodyConfig>& mbcs,
	const std::vector<Eigen::Vector3d>& coms,
	const std::vector<std::vector<sva::MotionVecd>>& normalAccB)
{
	eval_ = com_;
	speed_.setZero();
	normalAcc_.setZero();
	for(std::size_t i = 0; i < robotIndexes_.size(); ++i)
	{
		int r = robotIndexes_[i];
		const rbd::MultiBody& mb = mbs[r];
		const rbd::MultiBodyConfig& mbc = mbcs[r];

		eval_ -= coms[r]*robotsWeight_[i];
		speed_ += jac_[i].velocity(mb, mbc);
		normalAcc_ += jac_[i].normalAcceleration(mb, mbc, normalAccB[r]);
		jacMat_[i] = jac_[i].jacobian(mb, mbc);
	}
}


void MultiCoMTask::computeRobotsWeight(const std::vector<rbd::MultiBody>& mbs)
{
	double totalMass = 0.;

	robotsWeight_.clear();
	robotsWeight_.reserve(robotIndexes_.size());
	// compute the total mass and the weight of each robot
	for(int r: robotIndexes_)
	{
		double rm = std::accumulate(mbs[r].bodies().begin(),
			mbs[r].bodies().end(), 0., [](double ac, const rbd::Body& b)
				{return ac + b.inertia().mass();}
		);
		robotsWeight_.push_back(rm);
		totalMass += rm;
	}

	// divide all robotsWeight values by the total mass
	std::transform(robotsWeight_.begin(), robotsWeight_.end(),
		robotsWeight_.begin(),
		[totalMass](double w){return w/totalMass;});
}


const Eigen::VectorXd& MultiCoMTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& MultiCoMTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& MultiCoMTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& MultiCoMTask::jac(int index) const
{
	return jacMat_[index];
}


/**
	*													MomentumTask
	*/


MomentumTask::MomentumTask(const rbd::MultiBody& mb, const sva::ForceVecd mom):
	momentum_(mom),
	momentumMatrix_(mb),
	eval_(6),
	speed_(6),
	normalAcc_(6),
	jacMat_(6,mb.nrDof()),
	jacDotMat_(6,mb.nrDof())
{
	speed_.setZero();
}


void MomentumTask::momentum(const sva::ForceVecd& mom)
{
	momentum_ = mom;
}


const sva::ForceVecd MomentumTask::momentum() const
{
	return momentum_;
}


void MomentumTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	Eigen::Vector3d com = rbd::computeCoM(mb, mbc);

	eval_ = momentum_.vector() - rbd::computeCentroidalMomentum(mb, mbc, com).vector();
	normalAcc_ = momentumMatrix_.normalMomentumDot(mb, mbc,
		com,  rbd::computeCoMVelocity(mb, mbc)).vector();

	momentumMatrix_.computeMatrix(mb, mbc, com);
	jacMat_ = momentumMatrix_.matrix();
}


void MomentumTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	const std::vector<sva::MotionVecd>& normalAccB)
{
	Eigen::Vector3d com = rbd::computeCoM(mb, mbc);

	eval_ = momentum_.vector() - rbd::computeCentroidalMomentum(mb, mbc, com).vector();
	normalAcc_ = momentumMatrix_.normalMomentumDot(mb, mbc,
		com,  rbd::computeCoMVelocity(mb, mbc), normalAccB).vector();

	momentumMatrix_.computeMatrix(mb, mbc, com);
	jacMat_ = momentumMatrix_.matrix();
}


void MomentumTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	momentumMatrix_.computeMatrixDot(mb, mbc, rbd::computeCoM(mb, mbc),
				rbd::computeCoMVelocity(mb, mbc));
	jacDotMat_ = momentumMatrix_.matrixDot();
}


const Eigen::VectorXd& MomentumTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& MomentumTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& MomentumTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& MomentumTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& MomentumTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													LinVelocityTask
	*/


LinVelocityTask::LinVelocityTask(const rbd::MultiBody& mb, const std::string& bodyName,
	const Eigen::Vector3d& v, const Eigen::Vector3d& bodyPoint):
	vel_(v),
	point_(bodyPoint),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName, bodyPoint),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
	// this task don't have any derivative
	speed_.setZero();
}


void LinVelocityTask::velocity(const Eigen::Vector3d& v)
{
	vel_ = v;
}


const Eigen::Vector3d& LinVelocityTask::velocity() const
{
	return vel_;
}


void LinVelocityTask::bodyPoint(const Eigen::Vector3d& point)
{
	jac_.point(point);
}


const Eigen::Vector3d& LinVelocityTask::bodyPoint() const
{
	return jac_.point();
}


void LinVelocityTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	eval_ = vel_ - jac_.velocity(mb, mbc).linear();
	normalAcc_ = jac_.normalAcceleration(mb, mbc).linear();

	const auto& shortJacMat = jac_.jacobian(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void LinVelocityTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc, const std::vector<sva::MotionVecd>& normalAccB)
{
	eval_ = vel_ - jac_.velocity(mb, mbc).linear();
	normalAcc_ = jac_.normalAcceleration(mb, mbc, normalAccB).linear();

	const auto& shortJacMat = jac_.jacobian(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacMat_);
}


void LinVelocityTask::updateDot(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc)
{
	const auto& shortJacMat = jac_.jacobianDot(mb, mbc).block(3, 0, 3, jac_.dof());
	jac_.fullJacobian(mb, shortJacMat, jacDotMat_);
}


const Eigen::VectorXd& LinVelocityTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& LinVelocityTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& LinVelocityTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& LinVelocityTask::jac() const
{
	return jacMat_;
}


const Eigen::MatrixXd& LinVelocityTask::jacDot() const
{
	return jacDotMat_;
}


/**
	*													OrientationTrackingTask
	*/


OrientationTrackingTask::OrientationTrackingTask(const rbd::MultiBody& mb,
	const std::string& bodyName, const Eigen::Vector3d& bodyPoint,
	const Eigen::Vector3d& bodyAxis,
	const std::vector<std::string>& trackingJointsName,
	const Eigen::Vector3d& trackedPoint):
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	bodyPoint_(bodyPoint),
	bodyAxis_(bodyAxis),
	zeroJacIndex_(),
	trackedPoint_(trackedPoint),
	jac_(mb, bodyName),
	eval_(3),
	shortJacMat_(3, jac_.dof()),
	jacMat_(3, mb.nrDof()),
	jacDotMat_(3, mb.nrDof())
{
	std::set<int> trackingJointsIndex;
	for(const std::string& name: trackingJointsName)
	{
		trackingJointsIndex.insert(mb.jointIndexByName(name));
	}

	int jacPos = 0;
	for(int i: jac_.jointsPath())
	{
		const rbd::Joint& curJoint = mb.joint(i);
		if(trackingJointsIndex.find(i) == std::end(trackingJointsIndex))
		{
			for(int j = 0; j < curJoint.dof(); ++j)
			{
				zeroJacIndex_.push_back(jacPos + j);
			}
		}

		jacPos += curJoint.dof();
	}
}


void OrientationTrackingTask::trackedPoint(const Eigen::Vector3d& tp)
{
	trackedPoint_ = tp;
}


const Eigen::Vector3d& OrientationTrackingTask::trackedPoint() const
{
	return trackedPoint_;
}


void OrientationTrackingTask::bodyPoint(const Eigen::Vector3d& bp)
{
	bodyPoint_ = sva::PTransformd(bp);
}


const Eigen::Vector3d& OrientationTrackingTask::bodyPoint() const
{
	return bodyPoint_.translation();
}


void OrientationTrackingTask::bodyAxis(const Eigen::Vector3d& ba)
{
	bodyAxis_ = ba;
}


const Eigen::Vector3d& OrientationTrackingTask::bodyAxis() const
{
	return bodyAxis_;
}


void OrientationTrackingTask::update(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc)
{
	using namespace Eigen;
	const sva::PTransformd& bodyTf = mbc.bodyPosW[bodyIndex_];
	Vector3d desDir(trackedPoint_ - (bodyPoint_*bodyTf).translation());
	Vector3d curDir(bodyTf.rotation().transpose()*bodyAxis_);
	desDir.normalize();
	curDir.normalize();

	Matrix3d targetOri(
		Quaterniond::FromTwoVectors(curDir, desDir).inverse().matrix());

	eval_ = sva::rotationError<double>(mbc.bodyPosW[bodyIndex_].rotation(),
														targetOri*mbc.bodyPosW[bodyIndex_].rotation());

	shortJacMat_ = jac_.jacobian(mb, mbc).block(0, 0, 3, shortJacMat_.cols());
	zeroJacobian(shortJacMat_);
	jac_.fullJacobian(mb, shortJacMat_, jacMat_);
}


void OrientationTrackingTask::updateDot(const rbd::MultiBody& mb,
	const rbd::MultiBodyConfig& mbc)
{
	shortJacMat_ = jac_.jacobianDot(mb, mbc).block(0, 0, 3, shortJacMat_.cols());
	zeroJacobian(shortJacMat_);
	jac_.fullJacobian(mb, shortJacMat_, jacDotMat_);
}


const Eigen::MatrixXd& OrientationTrackingTask::jac()
{
	return jacMat_;
}


const Eigen::MatrixXd& OrientationTrackingTask::jacDot()
{
	return jacDotMat_;
}


const Eigen::VectorXd& OrientationTrackingTask::eval()
{
	return eval_;
}


void OrientationTrackingTask::zeroJacobian(Eigen::MatrixXd& jac) const
{
	for(int i: zeroJacIndex_)
	{
		jac.col(i).setZero();
	}
}


/**
	*													RelativeDistTask
	*/

RelativeDistTask::RelativeDistTask(const rbd::MultiBody& mb,
	const double timestep, const rbInfo& rbi1, const rbInfo& rbi2,
	const Eigen::Vector3d& u1, const Eigen::Vector3d& u2) :
	timestep_(timestep),
	isVectorFixed_(false),
	eval_(1),
	speed_(1),
	normalAcc_(1),
	jacMat_(1, mb.nrDof())
{
	rbInfo_[0] = RelativeDistTask::RelativeBodiesInfo(mb, rbi1, u1);
	rbInfo_[1] = RelativeDistTask::RelativeBodiesInfo(mb, rbi2, u2);
	if(u1!=Eigen::Vector3d::Zero() && u2!=Eigen::Vector3d::Zero())
	{
		isVectorFixed_ = true;
	}
}

RelativeDistTask::RelativeBodiesInfo::RelativeBodiesInfo(const rbd::MultiBody& mb,
	const rbInfo& rbi, const Eigen::Vector3d& fixedVector) :
	b1Index(mb.bodyIndexByName(std::get<0>(rbi))),
	r_b1_p(std::get<1>(rbi)),
	r_0_b2p(std::get<2>(rbi)),
	jac(mb, std::get<0>(rbi), r_b1_p),
	u(fixedVector),
	offn(Eigen::Vector3d::Zero())
{
}

void RelativeDistTask::robotPoint(const int bIndex, const Eigen::Vector3d& point)
{
	for(RelativeBodiesInfo& rbi : rbInfo_)
	{
		if(rbi.b1Index==bIndex)
		{
			rbi.r_b1_p = point;
			rbi.jac.point(point);
		}
	}
}

void RelativeDistTask::envPoint(const int bIndex, const Eigen::Vector3d& point)
{
	for(RelativeBodiesInfo& rbi : rbInfo_)
	{
		if(rbi.b1Index==bIndex)
		{
			rbi.r_0_b2p = point;
		}
	}
}

void RelativeDistTask::vector(const int bIndex, const Eigen::Vector3d& u)
{
	for(RelativeBodiesInfo& rbi : rbInfo_)
	{
		if(rbi.b1Index==bIndex)
		{
			rbi.u = u;
		}
	}
}

void RelativeDistTask::update(const rbd::MultiBody& mb, const rbd::MultiBodyConfig& mbc,
	const std::vector<sva::MotionVecd>& normalAccB)
{
	double sign = 1;
	jacMat_.setZero();
	eval_.setZero();
	speed_.setZero();
	normalAcc_.setZero();
	for(RelativeBodiesInfo& rbi : rbInfo_)
	{
		//Compute the error
		sva::PTransformd X_0_p1 = sva::PTransformd(rbi.r_b1_p)*mbc.bodyPosW[rbi.b1Index];
		Eigen::Vector3d r_0_p2 = rbi.r_0_b2p;
		Eigen::Vector3d r_p2_p1 = X_0_p1.translation()-r_0_p2;
		double d;
		Eigen::Vector3d n, dn;
		if(isVectorFixed_)
		{
			d = std::abs(r_p2_p1.dot(rbi.u));
			n = -rbi.u;
			dn = Eigen::Vector3d::Zero();
		}
		else
		{
			d = r_p2_p1.norm();
			n = (r_0_p2 - X_0_p1.translation())/d;
			dn = (n - rbi.offn)/timestep_;
		}
		eval_[0] += sign*d;

		//Compute the jacobian matrix
		Eigen::MatrixXd shortMat;
		Eigen::MatrixXd fullJac;
		fullJac.resize(3, mb.nrDof());
		shortMat = rbi.jac.jacobian(mb, mbc);
		rbi.jac.fullJacobian(mb, shortMat.block(3, 0, 3, rbi.jac.dof()), fullJac);
		jacMat_ += sign*n.transpose()*fullJac;

		//Compute the speed
		speed_[0] += sign*rbi.jac.velocity(mb, mbc).linear().dot(n);

		//Compute the normal acceleration (JDot alpha)
		normalAcc_[0] += sign*(fullJac*rbd::dofToVector(mb, mbc.alpha)).dot(dn)
					+ sign*rbi.jac.normalAcceleration(mb, mbc, normalAccB).linear().dot(n);

		//Update offn and sign
		//little hack: sign is +1 for the first pair and -1 for the second
		sign *= -1;
		rbi.offn = n;
	}
}

const Eigen::VectorXd& RelativeDistTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& RelativeDistTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& RelativeDistTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& RelativeDistTask::jac() const
{
	return jacMat_;
}

/**
	*													VectorOrientationTask
	*/

VectorOrientationTask::VectorOrientationTask(const rbd::MultiBody &mb,
	const std::string& bodyName,
	const Eigen::Vector3d &bodyVector, const Eigen::Vector3d &targetVector) :
	actualVector_(Eigen::Vector3d::Zero()),
	bodyVector_(bodyVector),
	targetVector_(targetVector),
	bodyIndex_(mb.bodyIndexByName(bodyName)),
	jac_(mb, bodyName),
	eval_(3),
	speed_(3),
	normalAcc_(3),
	jacMat_(3, mb.nrDof())
{
}

void VectorOrientationTask::update(const rbd::MultiBody &mb,
	const rbd::MultiBodyConfig &mbc,
	const std::vector<sva::MotionVecd> &normalAccB)
{
	//Evaluation of eval
	Eigen::Matrix3d E_0_b = mbc.bodyPosW[bodyIndex_].rotation().transpose();
	actualVector_ = E_0_b*bodyVector_;
	eval_ = targetVector_ - actualVector_;

	//Evaluation of speed and jacMat
	Eigen::MatrixXd shortMat, fullJac(3, mb.nrDof());
	shortMat = jac_.bodyJacobian(mb, mbc);
	jac_.fullJacobian(mb, shortMat.block(0, 0, 3, jac_.dof()), fullJac);
	jacMat_ = -E_0_b*skewMatrix(bodyVector_)*fullJac;
	Eigen::Vector3d w_b_b = jac_.bodyVelocity(mb, mbc).angular();
	speed_ = E_0_b*(w_b_b.cross(bodyVector_));

	//Evaluation of normalAcc
	Eigen::Vector3d bodyNormalAcc = jac_.bodyNormalAcceleration(mb, mbc, normalAccB).angular();
	normalAcc_ = w_b_b.cross(w_b_b.cross(bodyVector_));
	normalAcc_ += bodyNormalAcc.cross(bodyVector_);
	normalAcc_ = E_0_b*normalAcc_;
}

Eigen::Matrix3d VectorOrientationTask::skewMatrix(const Eigen::Vector3d &v)
{
	Eigen::Matrix3d m;
	m << 0., -v[2], v[1],
		 v[2], 0., -v[0],
		 -v[1], v[0], 0.;
	return m;
}

void VectorOrientationTask::bodyVector(const Eigen::Vector3d &vector)
{
	bodyVector_ = vector;
}

const Eigen::Vector3d &VectorOrientationTask::bodyVector() const
{
	return actualVector_;
}

void VectorOrientationTask::target(const Eigen::Vector3d &vector)
{
	targetVector_ = vector;
}

const Eigen::Vector3d& VectorOrientationTask::target() const
{
	return targetVector_;
}

const Eigen::Vector3d& VectorOrientationTask::actual() const
{
	return actualVector_;
}

const Eigen::VectorXd& VectorOrientationTask::eval() const
{
	return eval_;
}


const Eigen::VectorXd& VectorOrientationTask::speed() const
{
	return speed_;
}


const Eigen::VectorXd& VectorOrientationTask::normalAcc() const
{
	return normalAcc_;
}


const Eigen::MatrixXd& VectorOrientationTask::jac() const
{
	return jacMat_;
}

} // namespace tasks
