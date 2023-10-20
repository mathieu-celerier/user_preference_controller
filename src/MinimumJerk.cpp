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
  double proj_acc = A_0.transpose() * err;
  // if(proj_acc < 0)
  // {
  //   return computeTargetAcc();
  // }
  // else
  // {
  if(V_0.norm() >= V_limit)
  {
    mc_rtc::log::warning("[MinimumJerk] Velocity too high during acceleration phase. Acceleration set to zero.");
    output_source = "VelSupVelLim";
    return Vector3d::Zero();
  }
  else
  {
    auto V_extremums = velTrajectoryExtremums();
    // std::vector<double> vec_extremum(V_extremums.begin(), V_extremums.end());
    // mc_rtc::log::info("Trajectory velocity's extremums : {}", fmt::join(vec_extremum, ", "));
    if(V_extremums.upper_bound(V_limit) == V_extremums.end())
    {
      output_source = "VelExtremInBounds";
      return computeTargetAcc();
    }
    else
    {
      deltaT = solveNewDeltaT(V_extremums);
      mc_rtc::log::warning("[MinimumJerk] Computed trajectory for trajectory duration generates velocity out of "
                           "bounds, computing new trajectory for {}s duration",
                           deltaT);
      Rn(1, 1) = deltaT;
      Rn(2, 2) = pow(deltaT, 2);
      poly_coeffs = Ri * Rn * X_s + Ri * u;
      output_source = "NewDelta";
      return computeTargetAcc();
    }
  }
  // }
}

Vector3d MinimumJerk::computeTargetAcc(void)
{
  Rs(0, 4) = 1.0 / deltaT;
  Rs(1, 3) = 2.0 / pow(deltaT, 2);
  Rs(2, 2) = 6.0 / pow(deltaT, 3);

  Matrix3d dX_s = Rs * poly_coeffs;

  return A_0 + dt * dX_s.row(2).transpose();
}

Vector3d MinimumJerk::computeVelAtRoot(double t)
{
  return (1 / deltaT)
         * (5 * poly_coeffs.row(0) * pow(t, 4) + 4 * poly_coeffs.row(1) * pow(t, 3) + 3 * poly_coeffs.row(2) * pow(t, 2)
            + 2 * poly_coeffs.row(3) * t + poly_coeffs.row(4));
}

Vector3d MinimumJerk::computeJerkAtRoot(double t)
{
  return (1 / pow(deltaT, 3))
         * (60 * poly_coeffs.row(0) * pow(t, 2) + 24 * poly_coeffs.row(1) * t + 6 * poly_coeffs.row(2));
}

SetD MinimumJerk::velTrajectoryExtremums(void)
{
  SetD acc_roots = RootFinder::solvePolynomial(poly_coeffs.col(0), 0.0, 1.0, 0.0);
  SetD vel_at_roots;
  for(auto it = acc_roots.begin(); it != acc_roots.end(); it++)
  {
    vel_at_roots.insert(computeVelAtRoot(*it).norm());
  }

  return vel_at_roots;
}

double MinimumJerk::solveNewDeltaT(SetD vel_extremums)
{
  double max_delta = -1;
  for(auto it = vel_extremums.begin(); it != vel_extremums.end(); it++)
  {
    auto delta = (deltaT * (*it)) / V_limit;
    if(delta > max_delta) max_delta = delta;
  }
  return max_delta;
}

void MinimumJerk::add_to_logger(mc_rtc::Logger & logger)
{
  logger.addLogEntry("MinimumJerk_error", [this]() { return err; });
  logger.addLogEntry("MinimumJerk_output_source", [this]() { return output_source; });
  logger.addLogEntry("MinimumJerk_trajectory_time", [this]() { return deltaT; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_pos", [this]() { return X_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_vel", [this]() { return V_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_acc", [this]() { return A_0; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_pos", [this]() { return X_f; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_acc", [this]() { return target_acc; });
}
