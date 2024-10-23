#include "update/UpdaterZeroVelocity.h"
#include "feat/FeatureDatabase.h"
#include "feat/FeatureHelper.h"
#include "state/Propagator.h"
#include "state/State.h"
#include "state/StateHelper.h"
#include "utils/Print_Logger.h"
#include "state/Propagator.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "utils/Print_Logger.h"
#include "utils/quat_ops.h"
#include "types/Vec.h"
#include "utils/colors.h"
#include "types/IMU.h"

using namespace  ov_core ; 
using namespace mins ; 
using namespace ov_type ; 

UpdaterZeroVelocity::UpdaterZeroVelocity(OptionsZeroVelocityUpdate& options , OptionsIMU& noises,
                            std::shared_ptr<Propagator> prop ,double gravity_mag)
    : _options(options), _noises(noises), _prop(prop){
      PRINT2("Entered ZVUPT Constructor") ; 
    _gravity << 0.0,0.0,gravity_mag ; 
    for (int i = 1; i < 1000; i++) {
    boost::math::chi_squared chi_squared_dist(i);
    chi_squared_table[i] = boost::math::quantile(chi_squared_dist, 0.95);
  }
}

bool UpdaterZeroVelocity::try_update(std::shared_ptr<State> state, double timestamp , double time_offset){
    // Return if we don't have any imu data yet
  if (imu_data.empty()) {
    last_zupt_state_timestamp = 0.0;
    return false;
  }
  // Return if the state is already at the desired time
  if (state->time == timestamp) {
    last_zupt_state_timestamp = 0.0;
    return false;
  }
  auto pin_ts = imu_data.back().timestamp ; 
  
  // Set the last time offset value if we have just started the system up
  if (!have_last_prop_time_offset) {
    last_prop_time_offset = time_offset ;
    have_last_prop_time_offset = true;
  }
  double t_off_new = time_offset /*state->cam_dt.at((size_t)cam_id)->value()(0)*/ ;
  // First lets construct an IMU vector of measurements we need
  // double time0 = state->_timestamp+t_off_new;
  double time0 = state->time + last_prop_time_offset ;
  double time1 = timestamp + t_off_new ;
  // Select bounding inertial measurements
  std::vector<ov_core::ImuData> imu_recent = imu_data ; //mins::select_imu_readings(imu_data, time0, time1); 
  // Move forward in time
  last_prop_time_offset = t_off_new;
  // Check that we have at least one measurement to propagate with
  if (imu_recent.size() < 2) {
    PRINT2("[ZUPT]: There are no IMU data to check for zero velocity with!!\n");
    last_zupt_state_timestamp = 0.0;
    return false;
  }

  // If we should integrate the acceleration and say the velocity should be zero
  // Also if we should still inflate the bias based on their random walk noises
  bool integrated_accel_constraint = true ; // untested
  bool model_time_varying_bias = true;

  // Order of our Jacobian
  std::vector<std::shared_ptr<ov_type::Type>> Hx_order;
  Hx_order.push_back(state->imu->q());
  Hx_order.push_back(state->imu->bg());
  Hx_order.push_back(state->imu->ba());
  if (integrated_accel_constraint) {
    Hx_order.push_back(state->imu->v());
  }
  // Large final matrices used for update (we will compress these)
  int h_size = (integrated_accel_constraint) ? 12 : 9;
  int m_size = 6 * ((int)imu_recent.size() - 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m_size, h_size);
  Eigen::VectorXd res = Eigen::VectorXd::Zero(m_size);
   // Loop through all our IMU and construct the residual and Jacobian
  // TODO: should add jacobians here in respect to IMU intrinsics!!
  // State order is: [q_GtoI, bg, ba, v_IinG]
  // Measurement order is: [w_true = 0, a_true = 0 or v_k+1 = 0]
  // w_true = w_m - bw - nw
  // a_true = a_m - ba - R*g - na
  // v_true = v_k - g*dt + R^T*(a_m - ba - na)*dt
  // Mins does not support imu intrinsics so the values will be taken raw froom the imu
  // I sincerely hope it won't blow in my face 
  double dt_summed = 0;
  for (size_t i = 0; i < imu_recent.size() - 1; i++) {

    // Precomputed values
    double dt = imu_recent.at(i + 1).timestamp - imu_recent.at(i).timestamp;
    Eigen::Vector3d a_hat = imu_recent.at(i).am - state->imu->bias_a();
    Eigen::Vector3d w_hat =  imu_recent.at(i).wm - state->imu->bias_g() ;

    // Measurement noise (convert from continuous to discrete)
    // NOTE: The dt time might be different if we have "cut" any imu measurements
    // NOTE: We are performing "whittening" thus, we will decompose R_meas^-1 = L*L^t
    // NOTE: This is then multiplied to the residual and Jacobian (equivalent to just updating with R_meas)
    // NOTE: See Maybeck Stochastic Models, Estimation, and Control Vol. 1 Equations (7-21a)-(7-21c)
    double w_omega = std::sqrt(dt) / _noises.sigma_w;
    double w_accel = std::sqrt(dt) / _noises.sigma_a;
    double w_accel_v = 1.0 / (std::sqrt(dt) * _noises.sigma_a);

    // Measurement residual (true value is zero)
    res.block(6 * i + 0, 0, 3, 1) = -w_omega * w_hat;
    if (!integrated_accel_constraint) {
      res.block(6 * i + 3, 0, 3, 1) = -w_accel * (a_hat - state->imu->Rot() * _gravity);
    } else {
      res.block(6 * i + 3, 0, 3, 1) = -w_accel_v * (state->imu->vel() - _gravity * dt + state->imu->Rot().transpose() * a_hat * dt);
    }

    // Measurement Jacobian
    Eigen::Matrix3d R_GtoI_jacob = state->imu->Rot_fej()  ; 
    // Eigen::Matrix3d R_GtoI_jacob = state->imu->Rot(); //Since fej is supposed to be better than raw rot we will take this first
    H.block(6 * i + 0, 3, 3, 3) = -w_omega * Eigen::Matrix3d::Identity();
    if (!integrated_accel_constraint) {
      H.block(6 * i + 3, 0, 3, 3) = -w_accel * skew_x(R_GtoI_jacob * _gravity);
      H.block(6 * i + 3, 6, 3, 3) = -w_accel * Eigen::Matrix3d::Identity();
    } else {
      H.block(6 * i + 3, 0, 3, 3) = -w_accel_v * R_GtoI_jacob.transpose() * skew_x(a_hat) * dt;
      H.block(6 * i + 3, 6, 3, 3) = -w_accel_v * R_GtoI_jacob.transpose() * dt;
      H.block(6 * i + 3, 9, 3, 3) = w_accel_v * Eigen::Matrix3d::Identity();
    }
    dt_summed += dt;
  }
  // Compress the system (we should be over determined)
  StateHelper::measurement_compress_inplace(H, res);
  if (H.rows() < 1) {
    return false;
  }
   // Multiply our noise matrix by a fixed amount
  // We typically need to treat the IMU as being "worst" to detect / not become overconfident
  Eigen::MatrixXd R = _options.zupt_noise_multiplier * Eigen::MatrixXd::Identity(res.rows(), res.rows());

  // Next propagate the biases forward in time
  // NOTE: G*Qd*G^t = dt*Qd*dt = dt*(1/dt*Qc)*dt = dt*Qc
  Eigen::MatrixXd Q_bias = Eigen::MatrixXd::Identity(6, 6);
  Q_bias.block(0, 0, 3, 3) *= dt_summed *  std::pow(_noises.sigma_wb, 2);   //_noises.sigma_wb_2;
  Q_bias.block(3, 3, 3, 3) *= dt_summed *  std::pow(_noises.sigma_ab, 2);   //_noises.sigma_ab_2;

  // Chi2 distance check
  // NOTE: we also append the propagation we "would do before the update" if this was to be accepted (just the bias evolution)
  // NOTE: we don't propagate first since if we fail the chi2 then we just want to return and do normal logic
  Eigen::MatrixXd P_marg = StateHelper::get_marginal_covariance(state, Hx_order);
  if (model_time_varying_bias) {
    P_marg.block(3, 3, 6, 6) += Q_bias;
  }
  Eigen::MatrixXd S = H * P_marg * H.transpose() + R;
  double chi2 = res.dot(S.llt().solve(res));
  // Get our threshold (we precompute up to 1000 but handle the case that it is more)
  double chi2_check;
  if (res.rows() < 1000) {
    chi2_check = chi_squared_table[res.rows()];
  } else {
    boost::math::chi_squared chi_squared_dist(res.rows());
    chi2_check = boost::math::quantile(chi_squared_dist, 0.95);
    PRINT2("[ZUPT]: chi2_check over the residual limit - %d\n" , (int)res.rows());
  }
  // Check if the image disparity
  bool disparity_passed = false;

  // Check if we are currently zero velocity
  // We need to pass the chi2 and not be above our velocity threshold
  if (!disparity_passed && (chi2 > _options.zupt_chi2_multipler  * chi2_check || state->imu->vel().norm() > _options.zupt_max_velocity)) {
    last_zupt_state_timestamp = 0.0;
    last_zupt_count = 0;
    PRINT2(BOLDBLACK "[ZUPT]: rejected |v_IinG| = %.3f (chi2 %.3f > %.3f)\n" RESET , state->imu->vel().norm(), chi2,
                _options.zupt_chi2_multipler * chi2_check);
    return false;
  }
  PRINT2(REDPURPLE"[ZUPT]: accepted |v_IinG| = %.3f (chi2 %.3f < %.3f)\n" RESET, state->imu->vel().norm(), chi2,
             _options.zupt_chi2_multipler * chi2_check);
  // Do our update, only do this update if we have previously detected
  // If we have succeeded, then we should remove the current timestamp feature tracks
  // This is because we will not clone at this timestep and instead do our zero velocity update
  // NOTE: We want to keep the tracks from the second time we have called the zv-upt since this won't have a clone
  // NOTE: All future times after the second call to this function will also *not* have a clone, so we can remove those
  if (last_zupt_count >= 2) {
    // _db->cleanup_measurements_exact(last_zupt_state_timestamp);
  }
  else{
    last_zupt_count++ ; 
    return false ; 
  }

  // Else we are good, update the system
  // 1) update with our IMU measurements directly
  // 2) propagate and then explicitly say that our ori, pos, and vel should be zero
  // Next propagate the biases forward in time
  // NOTE: G*Qd*G^t = dt*Qd*dt = dt*Qc
  if (model_time_varying_bias) {
    Eigen::MatrixXd Phi_bias = Eigen::MatrixXd::Identity(6, 6);
    std::vector<std::shared_ptr<ov_type::Type>> Phi_order;
    Phi_order.push_back(state->imu->bg());
    Phi_order.push_back(state->imu->ba());
    StateHelper::EKFPropagation(state, Phi_order, Phi_order, Phi_bias, Q_bias);
  }
  // Finally move the state time forward
  StateHelper::EKFUpdate(state, Hx_order, H, res, R,"ZUPT");
  state->time = timestamp + t_off_new ;
  // Finally return
  last_zupt_state_timestamp = timestamp;
  last_zupt_count++;
  return true;
}

void UpdaterZeroVelocity::clean_old_imu_measurements(double oldest_time) {
  if (oldest_time < 0)
    return;
  auto it0 = imu_data.begin();
  while (it0 != imu_data.end()) {
    if (it0->timestamp < oldest_time) {
      it0 = imu_data.erase(it0);
    } else {
      it0++;
    }
  }
} 

void UpdaterZeroVelocity::feed_imu(const ov_core::ImuData &message, double oldest_time) {
    // Append it to our vector
    imu_data.emplace_back(message);
    // Loop through and delete IMU messages that are older than the oldest clone time
    // Sort our imu data (handles any out of order measurements)
    std::sort(imu_data.begin(), imu_data.end(), [](const ov_core::ImuData i, const ov_core::ImuData j) {
       return i.timestamp < j.timestamp;
    });
    for (auto data = imu_data.begin(); (*data).timestamp < oldest_time - 1;)
      data = imu_data.erase(data);
} 

std::vector<ov_core::ImuData> mins::select_imu_readings(const std::vector<ov_core::ImuData> &imu_data, double time0, double time1,
                                                           bool warn)
{
    // Ensure we have some measurements in the first place!
    if (imu_data.size() < 2) {
      if(warn) PRINT4("[Prop] Not enough IMU measurements. Requested: %.4f - %.4f\n" , time0, time1);
      return {};
    }
    // Make sure forward request
    if (time1 <= time0) {
      if(warn) (time1 < time0) && warn ? PRINT3(RED "Propagator::select_imu_readings::Backward request. time0: %.4f, time1: %.4f. \n", time0, time1) : void();
      return {};
    }

    // Make sure we have IMU data to process
    if (imu_data.front().timestamp > time0) {
      if(warn) PRINT3( "Propagator::select_imu_readings::Cannot handle request. " );
      if(warn) PRINT3( "time0 %.4f < oldest imu %.4f\n" , time0, imu_data.front().timestamp);
      return {};
    }

    // Make sure we have IMU data to process
    if (imu_data.back().timestamp < time1) {
      if(warn) PRINT3("Propagator::select_imu_readings::Cannot handle request. " );
      if(warn) PRINT3("newest imu %.4f < time1 %.4f\n", imu_data.back().timestamp, time1);
      return {};
    }

    std::vector<ov_core::ImuData> r ; 

    // Add the first data
    size_t i = 0;
    for (; i < imu_data.size() - 1; i++) {
      if (imu_data.at(i).timestamp <= time0 && time0 <= imu_data.at(i + 1).timestamp) {
        r.push_back(Propagator::interpolate_data(imu_data.at(i), imu_data.at(i + 1), time0));
        break;
      }
      assert(imu_data.at(i).timestamp <= time0);
    }

    // Add the middle data
    for (i == 0 ? i = 0 : i--; i < imu_data.size() - 1; i++) {
      if (time0 < imu_data.at(i).timestamp && imu_data.at(i + 1).timestamp < time1) {
        r.push_back(imu_data.at(i));
      }
      if (imu_data.at(i + 1).timestamp > time1)
        break;
    }

    // Add the final data
    for (i == 0 ? i = 0 : i--; i < imu_data.size() - 1; i++) {
      if (imu_data.at(i).timestamp <= time1 && time1 <= imu_data.at(i + 1).timestamp) {
        r.push_back(Propagator::interpolate_data(imu_data.at(i), imu_data.at(i + 1), time1));
        break;
      }
      assert(imu_data.at(i).timestamp <= time1);
    }

    assert(time0 == r.begin()->timestamp);
    assert(time1 == r.rbegin()->timestamp);
    return r ; 
}