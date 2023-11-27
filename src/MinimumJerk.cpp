#include "./MinimumJerk.h"

MinimumJerk::MinimumJerk(double dt_, double mass_)
: dt(dt_), Rs(Matrix36d::Zero()), Rn(Matrix63d::Zero()), u(Vector6d::Zero())
{
  Ri << -6., -3., -0.5, 6., -3., 0.5, 15., 8., 1.5, -15., 7., -1., -10., -6., -1.5, 10., -4., 0.5, 0., 0., 0.5, 0., 0.,
      0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.;

  Ris << 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., -60., -36., -9., 60., -24., 3.;

  Rn(0, 0) = 1;

  mc_rtc::log::info("Initialize Minimum-Jerk trajectory generator with dt {}[s] for body of mass {}[kg]", dt_, mass_);
}

Vector3d MinimumJerk::getTargetAcceleration(void)
{
  return target_acc;
}

Vector3d MinimumJerk::eval(void)
{
  return err;
}

void MinimumJerk::setTargetPos(Vector3d pos_)
{
  X_f = pos_;
  target_acc.setZero();
  acc_mj = 0.0;
  mc_rtc::log::info("[MinJerk Generator] New target position set : \n{}", u);
}

void MinimumJerk::setDuration(double T)
{
  deltaT = T;
}

void MinimumJerk::setJerkLimit(double limit)
{
  J_limit = limit;
}

void MinimumJerk::setAccelerationLimit(double limit)
{
  A_limit = limit;
}

void MinimumJerk::setDampingGain(double gain)
{
  Kd = gain;
}

void MinimumJerk::update(Vector3d pos_, Vector3d vel_, Vector3d acc_, Vector3d acc_dist_)
{
  delta_recompute_source = "None";

  deltaT -= dt;
  if(deltaT < dt)
  {
    delta_recompute_source = "LessThanDT";
    deltaT = dt;
  }

  X_0 = pos_;
  vel = vel_;

  err = X_f - X_0;

  u[3] = err.norm();

  // Project on the error
  V_0 = vel_.transpose() * (err / err.norm());
  A_0 = acc_.transpose() * (err / err.norm());
  A_0_dist = 0.0;

  X_s << 0.0, V_0, A_0;

  Rn(1, 1) = deltaT;
  Rn(2, 2) = pow(deltaT, 2);

  poly_coeffs = Ri * Rn * X_s + Ri * u;

  target_acc = decisionTreeControl();
}

Vector3d MinimumJerk::decisionTreeControl(void)
{
  // Check if current acceleration in any dimension is greater than limit
  if(abs(A_0) > A_limit or abs(A_0_dist) > A_limit)
  {
    delta_recompute_source = (abs(A_0) > A_limit) ? "CrntAccLimit" : "CmdAccLimit";
    return computeTargetAcc(true);
  }

  // Check if constraints are ensured
  double j0_init = computeJerkAtZero(deltaT);
  mc_rtc::log::info("e0 = {}", err.norm());
  mc_rtc::log::info("v0 = {}", V_0);
  mc_rtc::log::info("a0 = {}", A_0);
  mc_rtc::log::info("j0 = {}", j0_init);
  if(abs(j0_init) > J_limit * 1.001)
  {
    delta_recompute_source = "JerkConstraint";
    deltaT = solveNewDelta();

    Rn(1, 1) = deltaT;
    Rn(2, 2) = pow(deltaT, 2);

    poly_coeffs = Ri * Rn * X_s + Ri * u;

    return computeTargetAcc(false);
  }

  return computeTargetAcc(false);
}

Vector3d MinimumJerk::computeTargetAcc(bool onlydDamp)
{
  Rs(0, 4) = 1.0 / deltaT;
  Rs(1, 3) = 2.0 / pow(deltaT, 2);
  Rs(2, 2) = 6.0 / pow(deltaT, 3);

  Vector3d dXs = Rs * poly_coeffs;

  commanded_jerk = (onlydDamp) ? 0.0 : dXs[2];
  acc_mj = acc_mj + dt * commanded_jerk;

  Vector3d proj_back_jerk_acc = acc_mj * (err / err.norm());
  // Vector3d damping_acc = -dt * Kd * (vel - V_0 * (err / err.norm()));
  Vector3d damping_acc = -Kd * (vel - V_0 * (err / err.norm()));
  // mc_rtc::log::info("Amj = {}", proj_back_jerk_acc.transpose());
  // mc_rtc::log::info("Vel = {}", vel.transpose());
  // mc_rtc::log::info("Ortogonal Vel = {}", (vel - V_0 * (err / err.norm())).transpose());
  // mc_rtc::log::info("Admp = {}", damping_acc.transpose());

  return proj_back_jerk_acc + damping_acc;
}

double MinimumJerk::computeJerkAtZero(double deltaT_target)
{
  return (60.0 * err.norm()) / pow(deltaT_target, 3) - (36.0 * V_0) / pow(deltaT_target, 2) - (9 * A_0) / deltaT_target;
}

double MinimumJerk::solveNewDelta(void)
{
  SetD roots;
  SetD roots_j0 = solveJerkCst();
  if(not roots_j0.empty())
  {
    roots.insert(roots_j0.begin(), roots_j0.end());
  }
  // SetD roots_acc = solveAccCst();
  // if(not roots_acc.empty())
  // {
  //   roots.insert(roots_acc.begin(), roots_acc.end());
  // }

  double deltaT_new = -1.0;
  for(auto it = roots.begin(); it != roots.end(); it++)
  {
    auto deltaT_tmp = *it;
    double j0_new = computeJerkAtZero(deltaT_tmp);
    if(abs(j0_new) > J_limit * 1.001)
    {
      continue;
    }
    deltaT_new = deltaT_tmp;
    break;

    // bool allCstEnsured = checkAccCstForDeltaT(deltaT_tmp);
    //
    // if(allCstEnsured)
    // {
    //   deltaT_new = deltaT_tmp;
    //   break;
    // }
  }

  if(deltaT_new >= deltaT)
  {
    assert(not abs(computeJerkAtZero(deltaT_new)) > J_limit * 1.001);
    // assert(checkAccCstForDeltaT(deltaT_new));
    return deltaT_new;
  }
  else
  {
    mc_rtc::log::warning("No solutions found for current state");
    return deltaT;
  }
}

SetD MinimumJerk::solveJerkCst(void)
{
  SetD roots;
  Vector4d coeffs;
  // Compute set of roots for upper limit
  coeffs << -J_limit, -9.0 * A_0, -36.0 * V_0, 60.0 * err.norm();
  auto sub_root_set = RootFinder::solvePolynomial(coeffs, 0.0, INFINITY, 1e-6);
  auto valid_start = sub_root_set.upper_bound(deltaT);
  if(valid_start != sub_root_set.end())
  {
    roots.insert(valid_start, sub_root_set.end());
  }
  // Compute set of roots for lower limit
  coeffs << J_limit, -9.0 * A_0, -36.0 * V_0, 60.0 * err.norm();
  sub_root_set = RootFinder::solvePolynomial(coeffs, 0.0, INFINITY, 1e-6);
  valid_start = sub_root_set.upper_bound(deltaT);
  if(valid_start != sub_root_set.end())
  {
    roots.insert(valid_start, sub_root_set.end());
  }
  return roots;
}

void MinimumJerk::add_to_logger(mc_rtc::Logger & logger)
{
  logger.addLogEntry("MinimumJerk_error_vec", [this]() { return err; });
  logger.addLogEntry("MinimumJerk_error_norm", [this]() { return err.norm(); });
  logger.addLogEntry("MinimumJerk_output_source", [this]() { return output_source; });
  logger.addLogEntry("MinimumJerk_delta_recompute_source", [this]() { return delta_recompute_source; });
  logger.addLogEntry("MinimumJerk_trajectory_time", [this]() { return deltaT; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_pos_vec", [this]() { return X_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_vel_vec", [this]() { return V_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_acc_vec", [this]() { return A_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_dist_acc_vec", [this]() { return A_0_dist; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_upper_limit_acc", [this]() { return A_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_upper_limit_jerk", [this]() { return J_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_lower_limit_acc", [this]() { return -A_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_lower_limit_jerk", [this]() { return -J_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_pos", [this]() { return X_f; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_acc", [this]() { return target_acc; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_mj_acc", [this]() { return acc_mj; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_jerk", [this]() { return commanded_jerk; });
}
