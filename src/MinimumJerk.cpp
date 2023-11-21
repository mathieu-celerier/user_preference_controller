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

void MinimumJerk::setAccelerationLimit(double limit)
{
  A_limit = limit;
}

void MinimumJerk::setJerkLimit(double limit)
{
  J_limit = limit;
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
  V_0 = vel_;
  A_0 = acc_;
  A_0_dist = target_acc;

  err = X_f - X_0;

  X_s << X_0.transpose(), V_0.transpose(), A_0.transpose();

  Rn(1, 1) = deltaT;
  Rn(2, 2) = pow(deltaT, 2);

  poly_coeffs = Ri * Rn * X_s + Ri * u;

  target_acc = decisionTreeControl();
}

Vector3d MinimumJerk::decisionTreeControl(void)
{
  // Check if current acceleration in any dimension is greater than limit
  if((A_0.array().abs() > A_limit).any() or (target_acc.array().abs() > A_limit).any())
  {
    delta_recompute_source = ((A_0.array().abs() > A_limit).any()) ? "CrntAccLimit" : "CmdAccLimit";
    return Vector3d::Zero();
  }

  // Check if constraints are ensured
  Vector3d j0_init = computeJerkAtZero(deltaT);
  // mc_rtc::log::info("j0 = {}", j0_init.transpose());
  if((j0_init.array().abs() > J_limit).any())
  {
    delta_recompute_source = "JerkConstraint";
    deltaT = solveNewDelta();

    Rn(1, 1) = deltaT;
    Rn(2, 2) = pow(deltaT, 2);

    poly_coeffs = Ri * Rn * X_s + Ri * u;

    return computeTargetAcc();
  }

  if(not checkAccCstForDeltaT(deltaT))
  {
    delta_recompute_source = "AccelerationConstraint";
    deltaT = solveNewDelta();

    Rn(1, 1) = deltaT;
    Rn(2, 2) = pow(deltaT, 2);

    poly_coeffs = Ri * Rn * X_s + Ri * u;

    return computeTargetAcc();
  }

  return computeTargetAcc();
}

Vector3d MinimumJerk::computeTargetAcc(void)
{
  Rs(0, 4) = 1.0 / deltaT;
  Rs(1, 3) = 2.0 / pow(deltaT, 2);
  Rs(2, 2) = 6.0 / pow(deltaT, 3);

  Matrix3d dXs = Rs * poly_coeffs;

  commanded_jerk = dXs.row(2).transpose();

  return target_acc + dt * commanded_jerk;
}

Vector3d MinimumJerk::computeJerkAtZero(double deltaT_target)
{
  return (60.0 * err) / pow(deltaT_target, 3) - (36.0 * V_0) / pow(deltaT_target, 2) - (9 * A_0) / deltaT_target;
}

double MinimumJerk::computeAccAtRoot(double dim, double deltaT_target, double root)
{
  double err_dim = err[dim];
  double V0 = V_0[dim];
  double A0 = A_0_dist[dim];

  double a = 6 * err_dim - 3 * V0 * deltaT_target - 0.5 * A0 * pow(deltaT_target, 2);
  double b = -15 * err_dim + 8 * V0 * deltaT_target + 1.5 * A0 * pow(deltaT_target, 2);
  double c = 10 * err_dim - 6 * V0 * deltaT_target - 1.5 * A0 * pow(deltaT_target, 2);
  double d = 0.5 * A0 * pow(deltaT_target, 2);

  return (1 / pow(deltaT_target, 2)) * (20 * a * pow(root, 3) + 12 * b * pow(root, 2) + 6 * c * (root) + 2 * d);
}

double MinimumJerk::computeRootX1(double dim, double deltaT_target)
{
  double err_dim = err[dim];
  double V0 = V_0[dim];
  double A0 = A_0[dim];

  return (-sqrt(216 * pow(A0, 2) * pow(deltaT_target, 4) + 3024 * A0 * V0 * pow(deltaT_target, 3)
                + (10944 * pow(V0, 2) - 5760 * A0 * err_dim) * pow(deltaT_target, 2)
                - 43200 * V0 * err_dim * deltaT_target + 43200 * pow(err_dim, 2))
          - 24 * (1.5 * A0 * pow(deltaT_target, 2) + 8 * V0 * deltaT_target - 15 * err_dim))
         / (120 * (-0.5 * A0 * pow(deltaT_target, 2) - 3 * V0 * deltaT_target + 6 * err_dim));
}

double MinimumJerk::computeRootX2(double dim, double deltaT_target)
{
  double err_dim = err[dim];
  double V0 = V_0[dim];
  double A0 = A_0[dim];

  return (sqrt(216 * pow(A0, 2) * pow(deltaT_target, 4) + 3024 * A0 * V0 * pow(deltaT_target, 3)
               + (10944 * pow(V0, 2) - 5760 * A0 * err_dim) * pow(deltaT_target, 2)
               - 43200 * V0 * err_dim * deltaT_target + 43200 * pow(err_dim, 2))
          - 24 * (1.5 * A0 * pow(deltaT_target, 2) + 8 * V0 * deltaT_target - 15 * err_dim))
         / (120 * (-0.5 * A0 * pow(deltaT_target, 2) - 3 * V0 * deltaT_target + 6 * err_dim));
}

Vector3d MinimumJerk::computeDiscriminant(double deltaT_target)
{
  Vector3d A_0_sq = A_0_dist.array().square();
  Vector3d V_0_sq = V_0.array().square();
  Vector3d err_sq = err.array().square();
  Vector3d A_0_V_0 = A_0_dist.array() * V_0.array();
  Vector3d A_0_err = A_0_dist.array() * err.array();
  Vector3d V_0_err = V_0.array() * err.array();

  return 216.0 * A_0_sq * pow(deltaT_target, 4) + 3024.0 * A_0_V_0 * pow(deltaT_target, 3)
         + (10944.0 * V_0_sq - 5760.0 * A_0_err) * pow(deltaT_target, 2) - 43200 * V_0_err * deltaT_target
         + 43200 * err_sq;
}

Vector3d MinimumJerk::computeCstCoeffX1(double dim, double limit)
{
  double A0 = A_0_dist[dim];
  double V0 = V_0[dim];
  double err_dim = err[dim];
  Eigen::VectorXd coeffs(13);
  coeffs << -25 * pow(A0, 4) * pow(limit, 2) - 4 * pow(A0, 5) * limit + 2 * pow(A0, 6),
      (-600 * pow(A0, 3) * pow(limit, 2) - 132 * pow(A0, 4) * limit + 84 * pow(A0, 5)) * V0,
      (1200 * pow(A0, 3) * pow(limit, 2) + 156 * pow(A0, 4) * limit - 168 * pow(A0, 5)) * err_dim
          + (-5400 * pow(A0, 2) * pow(limit, 2) - 1800 * pow(A0, 3) * limit + 1476 * pow(A0, 4)) * pow(V0, 2),
      (21600 * pow(A0, 2) * pow(limit, 2) + 4632 * pow(A0, 3) * limit - 5928 * pow(A0, 4)) * V0 * err_dim
          + (-21600 * A0 * pow(limit, 2) - 12592 * pow(A0, 2) * limit + 13888 * pow(A0, 3)) * pow(V0, 3),
      (-21600 * pow(A0, 2) * pow(limit, 2) - 2016 * pow(A0, 3) * limit + 5868 * pow(A0, 4)) * pow(err_dim, 2)
          + (129600 * A0 * pow(limit, 2) + 52848 * pow(A0, 2) * limit - 84096 * pow(A0, 3)) * pow(V0, 2) * err_dim
          + (-32400 * pow(limit, 2) - 44832 * A0 * limit + 73776 * pow(A0, 2)) * pow(V0, 4),
      (-259200 * A0 * pow(limit, 2) - 59616 * pow(A0, 2) * limit + 167760 * pow(A0, 3)) * V0 * pow(err_dim, 2)
          + (259200 * pow(limit, 2) + 270240 * A0 * limit - 599136 * pow(A0, 2)) * pow(V0, 3) * err_dim
          + (209664 * A0 - 64512 * limit) * pow(V0, 5),
      (172800 * A0 * pow(limit, 2) + 8640 * pow(A0, 2) * limit - 110080 * pow(A0, 3)) * pow(err_dim, 3)
          + (-777600 * pow(limit, 2) - 542592 * A0 * limit + 1807536 * pow(A0, 2)) * pow(V0, 2) * pow(err_dim, 2)
          + (517248 * limit - 2141568 * A0) * pow(V0, 4) * err_dim + 248832 * pow(V0, 6),
      (1036800 * pow(limit, 2) + 362880 * A0 * limit - 2399040 * pow(A0, 2)) * V0 * pow(err_dim, 3)
          + (8686080 * A0 - 1554048 * limit) * pow(V0, 3) * pow(err_dim, 2) - 3068928 * pow(V0, 5) * err_dim,
      (1180800 * pow(A0, 2) - 518400 * pow(limit, 2)) * pow(err_dim, 4)
          + (2073600 * limit - 17481600 * A0) * pow(V0, 2) * pow(err_dim, 3) + 15683328 * pow(V0, 4) * pow(err_dim, 2),
      (17452800 * A0 - 1036800 * limit) * V0 * pow(err_dim, 4) - 42508800 * pow(V0, 3) * pow(err_dim, 3),
      64454400 * pow(V0, 2) * pow(err_dim, 4) - 6912000 * A0 * pow(err_dim, 5), -51840000 * V0 * pow(err_dim, 5),
      17280000 * pow(err_dim, 6);

  return coeffs;
}

Vector3d MinimumJerk::computeCstCoeffX2(double dim, double limit)
{
  double A0 = A_0_dist[dim];
  double V0 = V_0[dim];
  double err_dim = err[dim];
  Eigen::VectorXd coeffs(13);
  coeffs << -25 * pow(A0, 4) * pow(limit, 2) + 4 * pow(A0, 5) * limit + 2 * pow(A0, 6),
      (-600 * pow(A0, 3) * pow(limit, 2) + 132 * pow(A0, 4) * limit + 84 * pow(A0, 5)) * V0,
      (1200 * pow(A0, 3) * pow(limit, 2) - 156 * pow(A0, 4) * limit - 168 * pow(A0, 5)) * err_dim
          + (-5400 * pow(A0, 2) * pow(limit, 2) + 1800 * pow(A0, 3) * limit + 1476 * pow(A0, 4)) * pow(V0, 2),
      (21600 * pow(A0, 2) * pow(limit, 2) - 4632 * pow(A0, 3) * limit - 5928 * pow(A0, 4)) * V0 * err_dim
          + (-21600 * A0 * pow(limit, 2) + 12592 * pow(A0, 2) * limit + 13888 * pow(A0, 3)) * pow(V0, 3),
      (-21600 * pow(A0, 2) * pow(limit, 2) + 2016 * pow(A0, 3) * limit + 5868 * pow(A0, 4)) * pow(err_dim, 2)
          + (129600 * A0 * pow(limit, 2) - 52848 * pow(A0, 2) * limit - 84096 * pow(A0, 3)) * pow(V0, 2) * err_dim
          + (-32400 * pow(limit, 2) + 44832 * A0 * limit + 73776 * pow(A0, 2)) * pow(V0, 4),
      (-259200 * A0 * pow(limit, 2) + 59616 * pow(A0, 2) * limit + 167760 * pow(A0, 3)) * V0 * pow(err_dim, 2)
          + (259200 * pow(limit, 2) - 270240 * A0 * limit - 599136 * pow(A0, 2)) * pow(V0, 3) * err_dim
          + (209664 * A0 + 64512 * limit) * pow(V0, 5),
      (172800 * A0 * pow(limit, 2) - 8640 * pow(A0, 2) * limit - 110080 * pow(A0, 3)) * pow(err_dim, 3)
          + (-777600 * pow(limit, 2) + 542592 * A0 * limit + 1807536 * pow(A0, 2)) * pow(V0, 2) * pow(err_dim, 2)
          + (-517248 * limit - 2141568 * A0) * pow(V0, 4) * err_dim + 248832 * pow(V0, 6),
      (1036800 * pow(limit, 2) - 362880 * A0 * limit - 2399040 * pow(A0, 2)) * V0 * pow(err_dim, 3)
          + (8686080 * A0 + 1554048 * limit) * pow(V0, 3) * pow(err_dim, 2) - 3068928 * pow(V0, 5) * err_dim,
      (1180800 * pow(A0, 2) - 518400 * pow(limit, 2)) * pow(err_dim, 4)
          + (-2073600 * limit - 17481600 * A0) * pow(V0, 2) * pow(err_dim, 3) + 15683328 * pow(V0, 4) * pow(err_dim, 2),
      (17452800 * A0 + 1036800 * limit) * V0 * pow(err_dim, 4) - 42508800 * pow(V0, 3) * pow(err_dim, 3),
      64454400 * pow(V0, 2) * pow(err_dim, 4) - 6912000 * A0 * pow(err_dim, 5), -51840000 * V0 * pow(err_dim, 5),
      17280000 * pow(err_dim, 6);

  return coeffs;
}

bool MinimumJerk::checkAccCstForDeltaT(double deltaT_target)
{
  double x1_new;
  double x2_new;
  double ax1_new;
  double ax2_new;
  bool allCstEnsured = true;
  auto newDiscriminant = computeDiscriminant(deltaT_target);
  for(size_t i = 0; i < 3; i++)
  {
    if(newDiscriminant[i] > 0)
    {
      x1_new = computeRootX1(i, deltaT_target);
      x2_new = computeRootX2(i, deltaT_target);
      ax1_new = computeAccAtRoot(i, deltaT_target, x1_new);
      ax2_new = computeAccAtRoot(i, deltaT_target, x2_new);

      if(x1_new > 0 and x1_new < 1)
      {
        if(abs(ax1_new) > (A_limit / 3.))
        {
          allCstEnsured = false;
          break;
        }
      }

      if(x2_new > 0 and x2_new < 1)
      {
        if(abs(ax2_new) > (A_limit / 3.))
        {
          allCstEnsured = false;
          break;
        }
      }
    }
  }

  return allCstEnsured;
}

double MinimumJerk::solveNewDelta(void)
{
  SetD roots;
  SetD roots_j0 = solveJerkCst();
  if(not roots_j0.empty())
  {
    roots.insert(roots_j0.begin(), roots_j0.end());
  }
  SetD roots_acc = solveAccCst();
  if(not roots_acc.empty())
  {
    roots.insert(roots_acc.begin(), roots_acc.end());
  }

  double deltaT_new = -1.0;
  for(auto it = roots.begin(); it != roots.end(); it++)
  {
    auto deltaT_tmp = *it;
    Vector3d j0_new = computeJerkAtZero(deltaT_tmp);
    if((j0_new.array().abs() > J_limit).any())
    {
      continue;
    }

    bool allCstEnsured = checkAccCstForDeltaT(deltaT_tmp);

    if(allCstEnsured)
    {
      deltaT_new = deltaT_tmp;
      break;
    }
  }

  if(deltaT_new >= deltaT)
  {
    assert(not(computeJerkAtZero(deltaT_new).array().abs() > J_limit).any());
    assert(checkAccCstForDeltaT(deltaT_new));
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
  for(size_t i = 0; i < 3; i++)
  {
    Vector4d coeffs;
    // Compute set of roots for upper limit
    coeffs << -J_limit, -9.0 * A_0[i], -36.0 * V_0[i], 60.0 * err[i];
    auto sub_root_set = RootFinder::solvePolynomial(coeffs, 0.0, INFINITY, 1e-6);
    auto valid_start = sub_root_set.upper_bound(deltaT);
    if(valid_start != sub_root_set.end())
    {
      roots.insert(valid_start, sub_root_set.end());
    }
    // Compute set of roots for lower limit
    coeffs << J_limit, -9.0 * A_0[i], -36.0 * V_0[i], 60.0 * err[i];
    sub_root_set = RootFinder::solvePolynomial(coeffs, 0.0, INFINITY, 1e-6);
    valid_start = sub_root_set.upper_bound(deltaT);
    if(valid_start != sub_root_set.end())
    {
      roots.insert(valid_start, sub_root_set.end());
    }
  }
  return roots;
}

SetD MinimumJerk::solveAccCst(void)
{
  SetD roots;
  auto discriminant = computeDiscriminant(deltaT);
  for(size_t i = 0; i < 3; i++)
  {
    if(discriminant[i] > 0)
    {
      auto sub_root_set_X1 = solveAccCstX1(i);
      if(!sub_root_set_X1.empty())
      {
        roots.insert(sub_root_set_X1.begin(), sub_root_set_X1.end());
      }
      auto sub_root_set_X2 = solveAccCstX2(i);
      if(!sub_root_set_X2.empty())
      {
        roots.insert(sub_root_set_X2.begin(), sub_root_set_X2.end());
      }
    }
  }
  return roots;
}

SetD MinimumJerk::solveAccCstX1(double dim)
{
  SetD roots;

  // Search roots for upper limit
  Vector3d coeffs = computeCstCoeffX1(dim, A_limit);
  auto sub_root_set = RootFinder::solvePolynomial(coeffs, 0.0, INFINITY, 1e-6);
  auto valid_start = sub_root_set.upper_bound(deltaT);
  if(valid_start != sub_root_set.end())
  {
    roots.insert(valid_start, sub_root_set.end());
  }

  // Search roots for lower limit
  coeffs = computeCstCoeffX1(dim, -A_limit);
  sub_root_set = RootFinder::solvePolynomial(coeffs, 0.0, INFINITY, 1e-6);
  valid_start = sub_root_set.upper_bound(deltaT);
  if(valid_start != sub_root_set.end())
  {
    roots.insert(valid_start, sub_root_set.end());
  }

  return roots;
}

SetD MinimumJerk::solveAccCstX2(double dim)
{
  SetD roots;

  // Search roots for upper limit
  Vector3d coeffs = computeCstCoeffX2(dim, A_limit);
  auto sub_root_set = RootFinder::solvePolynomial(coeffs, 0.0, INFINITY, 1e-6);
  auto valid_start = sub_root_set.upper_bound(deltaT);
  if(valid_start != sub_root_set.end())
  {
    roots.insert(valid_start, sub_root_set.end());
  }

  // Search roots for lower limit
  coeffs = computeCstCoeffX2(dim, -A_limit);
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
  logger.addLogEntry("MinimumJerk_error", [this]() { return err; });
  logger.addLogEntry("MinimumJerk_output_source", [this]() { return output_source; });
  logger.addLogEntry("MinimumJerk_delta_recompute_source", [this]() { return delta_recompute_source; });
  logger.addLogEntry("MinimumJerk_trajectory_time", [this]() { return deltaT; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_pos_vec", [this]() { return X_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_pos_norm", [this]() { return X_0.norm(); });
  logger.addLogEntry("MinimumJerk_trajectory_cur_vel_vec", [this]() { return V_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_vel_norm", [this]() { return V_0.norm(); });
  logger.addLogEntry("MinimumJerk_trajectory_cur_acc_vec", [this]() { return A_0; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_acc_norm", [this]() { return A_0.norm(); });
  logger.addLogEntry("MinimumJerk_trajectory_cur_dist_acc_vec", [this]() { return A_0_dist; });
  logger.addLogEntry("MinimumJerk_trajectory_cur_dist_acc_norm", [this]() { return A_0_dist.norm(); });
  logger.addLogEntry("MinimumJerk_trajectory_ref_norm_limit_acc", [this]() { return A_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_norm_limit_jerk", [this]() { return J_limit; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_pos", [this]() { return X_f; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_acc", [this]() { return target_acc; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_jerk", [this]() { return commanded_jerk; });
  logger.addLogEntry("MinimumJerk_trajectory_ref_jerk_norm", [this]() { return commanded_jerk.norm(); });
}
