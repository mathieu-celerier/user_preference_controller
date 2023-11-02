#include "./MinimumJerk.h"

MinimumJerk::MinimumJerk(double dt_, double mass_)
: dt(dt_), Rs(Matrix36d::Zero()), Rn(Matrix63d::Zero()), u(Matrix63d::Zero())
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
  u.row(3) = pos_;
  target_acc.setZero();
  mc_rtc::log::info("[MinJerk Generator] New target position set : \n{}", u);
}

void MinimumJerk::setDuration(double T)
{
  deltaT = T;
}

void MinimumJerk::setVelocityLimit(double limit)
{
  V_limit = limit;
}

void MinimumJerk::setAccelerationLimit(double limit)
{
  A_limit = limit;
}

void MinimumJerk::setJerkLimit(double limit)
{
  J_limit = limit;
}

void MinimumJerk::update(Vector3d pos_, Vector3d vel_, Vector3d acc_)
{
  X_0 = pos_;
  V_0 = vel_;
  A_0 = acc_;

  err = X_f - X_0;

  X_s << X_0.transpose(), V_0.transpose(), A_0.transpose();

  Rn(1, 1) = deltaT;
  Rn(2, 2) = pow(deltaT, 2);

  poly_coeffs = Ri * Rn * X_s + Ri * u;

  target_acc = decisionTreeControl();

  deltaT -= dt;
  if(deltaT < dt)
  {
    deltaT = dt;
  }
}

Vector3d MinimumJerk::decisionTreeControl(void)
{
  if(V_0.norm() >= V_limit)
  {
    if(V_0.transpose() * A_0 > 0) // Check if we are trying to slow the motion or not
    {
      // If not, for know for now let's send 0 accleration but may be nice to only set the colinear part to 0
      mc_rtc::log::warning("[MinimumJerk] Velocity too high during acceleration phase. Acceleration set to zero.");
      output_source = "VelSupVelLim";
      return Vector3d::Zero();
    }
    // If we are trying to slow the motion let's keep sending the normal acceleration
  }

  // If the velocity is within limit let's compute all bounds for Delta T
  deltaT = solveNewDeltaT();
  Rn(1, 1) = deltaT;
  Rn(2, 2) = pow(deltaT, 2);
  poly_coeffs = Ri * Rn * X_s + Ri * u;
  return computeTargetAcc();
}

Vector3d MinimumJerk::computeTargetAcc(void)
{
  Rs(0, 4) = 1.0 / deltaT;
  Rs(1, 3) = 2.0 / pow(deltaT, 2);
  Rs(2, 2) = 6.0 / pow(deltaT, 3);

  Matrix3d dX_s = Rs * poly_coeffs;

  commanded_jerk = dX_s.row(2).transpose();

  return target_acc + dt * commanded_jerk;
}

Vector3d MinimumJerk::computeVelAtRoot(double t)
{
  return (1 / deltaT)
         * (5 * poly_coeffs.row(0) * pow(t, 4) + 4 * poly_coeffs.row(1) * pow(t, 3) + 3 * poly_coeffs.row(2) * pow(t, 2)
            + 2 * poly_coeffs.row(3) * t + poly_coeffs.row(4));
}

Vector3d MinimumJerk::computeAccAtRoot(double t)
{
  return (1 / pow(deltaT, 2))
         * (20 * poly_coeffs.row(0) * pow(t, 3) + 12 * poly_coeffs.row(1) * pow(t, 2) + 6 * poly_coeffs.row(2) * t
            + 2 * poly_coeffs.row(3));
}

Vector3d MinimumJerk::computeJerkAtRoot(double t)
{
  return (1 / pow(deltaT, 3))
         * (60 * poly_coeffs.row(0) * pow(t, 2) + 24 * poly_coeffs.row(1) * t + 6 * poly_coeffs.row(2));
}

SetD MinimumJerk::computeRootsOfJerk(double newDeltaT)
{
  SetD jerk_roots;
  for(int i = 0; i < 3; i++)
  {
    Eigen::Vector4d coeffs(60);
    // mc_rtc::log::info("X_0 = {}, V_0 = {}, A_0 = {}, X_f = {}", X_0[i], V_0[i], A_0[i], X_f[0]);
    // coeffs << 60 * (X_f[i] - X_0[i]), -36 * V_0[i], -9 * A_0[i], -J_limit;
    coeffs << -J_limit, -9 * A_0[i], -36 * V_0[i], 60 * (X_f[i] - X_0[i]);
    SetD axis_roots = RootFinder::solvePolynomial(coeffs, 0.0, INFINITY, 0.0);

    if(!axis_roots.empty())
    {
      auto first_root = axis_roots.upper_bound(newDeltaT);
      if(first_root != axis_roots.end())
      {
        jerk_roots.insert(*first_root);
      }
    }
  }

  return jerk_roots;
}

SetD MinimumJerk::velTrajectoryExtremums(void)
{
  SetD acc_roots;
  for(int i = 0; i < 3; i++)
  {
    SetD sub_roots = RootFinder::solvePolynomial(poly_coeffs.col(i), 0.0, 1.0, 0.0);
    acc_roots.insert(sub_roots.begin(), sub_roots.end());
  }

  SetD vel_at_roots;
  for(auto it = acc_roots.begin(); it != acc_roots.end(); it++)
  {
    vel_at_roots.insert(computeVelAtRoot(*it).norm());
  }

  return vel_at_roots;
}

double MinimumJerk::solveNewDeltaT()
{
  auto oldDeltaT = deltaT;
  // auto deltaSet = computeRootsOfJerk(deltaT);
  // if(!deltaSet.empty())
  // {
  //   deltaT = *deltaSet.rbegin();
  //   output_source = "JerkLimit";
  // }

  if(oldDeltaT != deltaT)
  {
    mc_rtc::log::warning("[MinimumJerk] Computed trajectory for trajectory duration generates kinematics out of "
                         "bounds, computing new trajectory for {}s duration",
                         deltaT);
  }
  else
  {
    output_source = "NormalControl";
  }

  return deltaT;
}

void MinimumJerk::add_to_logger(mc_rtc::Logger & logger)
{
  logger.addLogEntry("MinimumJerk_error", [this]() { return err; });
  logger.addLogEntry("MinimumJerk_output_source", [this]() { return output_source; });
  logger.addLogEntry("MinimumJerk_trajectory_time", [this]() { return deltaT; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_pos_vec", [this]() { return X_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_pos_norm", [this]() { return X_0.norm(); });
  logger.addLogEntry("MinimumJerk_trajectory_cur_vel_vec", [this]() { return V_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_vel_norm", [this]() { return V_0.norm(); });
  logger.addLogEntry("MinimumJerk_trajectory_cur_acc_vec", [this]() { return A_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_acc_norm", [this]() { return A_0.norm(); });
  logger.addLogEntry("MinimumJerk_trajectory_cur_norm_limit_vel", [this]() { return V_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_norm_limit_acc", [this]() { return A_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_norm_limit_jerk", [this]() { return J_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_pos", [this]() { return X_f; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_acc", [this]() { return target_acc; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_jerk", [this]() { return commanded_jerk; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_jerk_norm", [this]() { return commanded_jerk.norm(); });
  logger.addLogEntry("MinimumJerk_trajectory_polynomial_coefficient_a",
                     [this]() { return (Eigen::Vector3d)poly_coeffs.row(0); });
  logger.addLogEntry("MinimumJerk_trajectory_polynomial_coefficient_b",
                     [this]() { return (Eigen::Vector3d)poly_coeffs.row(1); });
  logger.addLogEntry("MinimumJerk_trajectory_polynomial_coefficient_c",
                     [this]() { return (Eigen::Vector3d)poly_coeffs.row(2); });
  logger.addLogEntry("MinimumJerk_trajectory_polynomial_coefficient_d",
                     [this]() { return (Eigen::Vector3d)poly_coeffs.row(3); });
  logger.addLogEntry("MinimumJerk_trajectory_polynomial_coefficient_e",
                     [this]() { return (Eigen::Vector3d)poly_coeffs.row(4); });
  logger.addLogEntry("MinimumJerk_trajectory_polynomial_coefficient_f",
                     [this]() { return (Eigen::Vector3d)poly_coeffs.row(5); });
}
