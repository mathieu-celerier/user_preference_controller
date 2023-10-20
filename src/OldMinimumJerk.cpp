#include "./OldMinimumJerk.h"

OldMinimumJerk::OldMinimumJerk(double dt_, double mass_) : dt(dt_), mass(mass_), acc_lim(1)
{
  // Ri << -60., -36., -9., 60., -24., 3., 360., 192., 36., -360., 168., -24., -720., -360., -60., 720., -360., 60.;
  // Ri << 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., -60., -36., -9., 60., -24., 3.;
  Ri << -6., -3., -0.5, 6., -3., 0.5, 15., 8., 1.5, -15., 7., -1., -10., -6., -1.5, 10., -4., 0.5, 0., 0., 0.5, 0., 0.,
      0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.;

  vel_limit = 2.;

  Rs.setZero();
  RdeltaT.setZero();
  RdeltaT_d.setZero();
  initConditions.setZero();
  X_0.setZero();

  mc_rtc::log::info("Initialize Minimum-Jerk tracjectory generator with dt {}[s] for body of mass {}[kg]", dt, mass);
}

Eigen::Vector3d OldMinimumJerk::eval(void)
{
  return err;
}

void OldMinimumJerk::set_target(Eigen::Vector3d target, double duration)
{
  T = duration;
  t = 0;

  X_f = target;
  initConditions.row(3) = X_f;

  // mc_rtc::log::info("Ra = \n{}", Ra);
  // mc_rtc::log::info("Rb = \n{}", Rb);
  // mc_rtc::log::info("Initial conditions = \n{}", initConditions);
}

Eigen::Vector3d OldMinimumJerk::get_target_acc(void)
{
  auto target = acc_target_clipped;
  // std::cout << "==================== Return target acc ====================" << std::endl;
  // std::cout << fmt::format("Target acceleration = \n{}", target.transpose()) << std::endl;
  return target;
}

void OldMinimumJerk::set_force_limit(double limit_)
{
  acc_lim = limit_ / mass;
}

void OldMinimumJerk::update(Eigen::Vector3d pos_, Eigen::Vector3d vel_, Eigen::Vector3d acc_, bool discrete)
{
  // std::cout << "==================== Entering update function ====================" << std::endl;
  pos = pos_;
  vel = vel_;
  acc = acc_;
  err = X_f - pos;
  pos_norm = pos.norm();
  vel_norm = vel.norm();
  acc_norm = acc.norm();
  err_norm = err.norm();

  if(discrete)
  {
    discrete_acc();
  }
  else
  {
    continuous_acc();
  }
}

void OldMinimumJerk::continuous_acc(void)
{
  Eigen::Vector3d reference_acc = continuous_ref_acc();
  if(acc.transpose() * err < 0.0)
  {
    return;
  }
  else
  {
    if(vel.norm() >= vel_limit)
    {
      // For now let's not break and just stop accelerating
      reference_acc = Eigen::Vector3d::Zero();
    }
    else
    {
      auto vel_max = compute_vel_max();
      // if(vel_max > vel_limit)
      // {
      //   return;
      // }
      // else
      // {
      //   // Solve delta T such that vel_max = vel_limit
      // }
    }
  }

  acc_target_clipped = reference_acc;
}

Eigen::Vector3d OldMinimumJerk::continuous_ref_acc(void)
{
  Eigen::Vector3d reference_acc;
  // std::cout << fmt::format("Acc = \n{}", acc_) << std::endl;
  // std::cout << fmt::format("Commanded Acc = \n{}", c_acc_) << std::endl;

  // Define intial pose at the start of trajectory to compute non constrained reference trajectory
  if(t == 0)
  {
    X_0 = pos;
  }

  // Update current state
  state << pos.transpose(), vel.transpose(), acc.transpose();
  state_norm << pos_norm, vel_norm, acc_norm;

  // std::cout << "==================== Updating state variables ====================" << std::endl;
  // std::cout << fmt::format("Acc = \n{}", acc_.transpose()) << std::endl;
  // std::cout << fmt::format("Commanded Acc = \n{}", c_acc_.transpose()) << std::endl;
  // std::cout << fmt::format("State = \n{}", state) << std::endl;
  // std::cout << fmt::format("State commanded acc = \n{}", state_c_acc) << std::endl;

  // Update Matrices
  auto deltaT = T - t;
  RdeltaT(0, 0) = 1 / pow(deltaT, 3);
  RdeltaT(1, 1) = 1 / pow(deltaT, 4);
  RdeltaT(2, 2) = 1 / pow(deltaT, 5);
  Rs(0, 0) = 1;
  Rs(1, 1) = deltaT;
  Rs(2, 2) = pow(deltaT, 2);
  Eigen::MatrixXd coeffs = Ri * Rs * state + Ri * initConditions;
  // Eigen::VectorXd acc_root = compute_acc_roots(coeffs, deltaT);

  Rb = RdeltaT * Ri;
  Ra = Rb * Rs;

  // Compute state dot
  stateDot = Ra * state + Rb * initConditions;
  // std::cout << "==================== Computing state dot ====================" << std::endl;
  // std::cout << fmt::format("Acc = \n{}", acc_.transpose()) << std::endl;
  // std::cout << fmt::format("Commanded Acc = \n{}", c_acc_.transpose()) << std::endl;
  // std::cout << fmt::format("State = \n{}", state) << std::endl;
  // std::cout << fmt::format("State commanded acc = \n{}", state_c_acc) << std::endl;
  // std::cout << fmt::format("State dot = \n{}", stateDot) << std::endl;
  // std::cout << fmt::format("State dot commanded acc = \n{}", stateDotCAcc) << std::endl;

  // Update target acceleration
  Eigen::Vector3d jerk = stateDot.row(0);
  Eigen::Vector3d snap = stateDot.row(1);
  Eigen::Vector3d crackle = stateDot.row(2);
  acc_target = acc + jerk * dt + snap * dt * dt / 2.0 + crackle * dt * dt * dt / 6.0;
  std::cout << "Norm lim = " << acc_lim << std::endl;
  std::cout << "Norm = " << acc_target.norm() << std::endl;
  if(acc_lim < acc_target.norm())
  {
    acc_target_clipped = acc_target * (acc_lim / acc_target.norm());
  }
  else
  {
    acc_target_clipped = acc_target;
  }
  std::cout << "Clipped norm = " << acc_target_clipped.norm() << std::endl;

  // std::cout << "==================== Computing target acc ====================" << std::endl;
  // std::cout << fmt::format("Acc = \n{}", acc_.transpose()) << std::endl;
  // std::cout << fmt::format("Commanded Acc = \n{}", c_acc_.transpose()) << std::endl;
  // std::cout << fmt::format("State = \n{}", state) << std::endl;
  // std::cout << fmt::format("State commanded acc = \n{}", state_c_acc) << std::endl;
  // std::cout << fmt::format("State dot = \n{}", stateDot) << std::endl;
  // std::cout << fmt::format("State dot commanded acc = \n{}", stateDotCAcc) << std::endl;
  // std::cout << fmt::format("Jerk = \n{}", jerk.transpose()) << std::endl;
  // std::cout << fmt::format("Snap = \n{}", snap.transpose()) << std::endl;
  // std::cout << fmt::format("Crackle = \n{}", crackle.transpose()) << std::endl;
  // std::cout << fmt::format("Jerk commanded acc = \n{}", jerk_c_acc.transpose()) << std::endl;
  // std::cout << fmt::format("Snap commanded acc = \n{}", snap_c_acc.transpose()) << std::endl;
  // std::cout << fmt::format("Crackle commanded acc = \n{}", crackle_c_acc.transpose()) << std::endl;
  // std::cout << fmt::format("Target acceleration = \n{}", acc_target.transpose()) << std::endl;
  // std::cout << fmt::format("Target acceleration commanded acc = \n{}", c_acc_target.transpose()) << std::endl;

  ideal_pos = X_0 + (X_f - X_0) * (10 * pow(t / T, 3) - 15 * pow(t / T, 4) + 6 * pow(t / T, 5));

  t += dt;
  auto th = dt * 20;
  if(deltaT < th) t = T - th;

  return reference_acc;
}

Eigen::Vector3d OldMinimumJerk::discrete_ref_acc(void)
{
  Eigen::Vector3d reference_acc;
  auto deltaT = T - t;

  Rs(0, 0) = 1;
  Rs(1, 1) = deltaT;
  Rs(2, 2) = (deltaT * deltaT);
  auto rk = dt / deltaT;
  RdeltaT_d << 20.0 * pow(rk, 3), 12.0 * pow(rk, 2), 6. * rk, 2., 0., 0.;
  RdeltaT_d = (1 / pow(deltaT, 2)) * RdeltaT_d;

  Rb_d = RdeltaT_d * Ri_d;
  Ra_d = Rb_d * Rs;

  state << pos.transpose(), vel.transpose(), acc.transpose();

  acc_target << Ra_d * state + Rb_d * initConditions;

  t += dt;
  auto th = dt * 10;
  if(deltaT < th) t = T - th;

  return reference_acc;
}

Eigen::Vector3d OldMinimumJerk::compute_vel_max(void)
{

  return Eigen::Vector3d::Zero();
}

Eigen::VectorXd OldMinimumJerk::compute_acc_roots(Eigen::Vector6d eq_param, double deltaT)
{
  return Eigen::VectorXd::Zero(1);
}

void OldMinimumJerk::add_to_logger(mc_rtc::Logger & logger)
{
  logger.addLogEntry("OldMinimumJerk_deltaT", [this]() { return this->T - this->t; });
  logger.addLogEntry("OldMinimumJerk_ideal_pos", [this]() { return this->ideal_pos; });
  logger.addLogEntry("OldMinimumJerk_target_pos", [this]() { return this->X_f; });
  logger.addLogEntry("OldMinimumJerk_target_acceleration", [this]() { return this->acc_target; });
  logger.addLogEntry("OldMinimumJerk_target_acceleration_clipped", [this]() { return this->acc_target_clipped; });
  logger.addLogEntry("OldMinimumJerk_target_jerk", [this]() { return (Eigen::Vector3d)this->stateDot.row(0); });
  logger.addLogEntry("OldMinimumJerk_target_snap", [this]() { return (Eigen::Vector3d)this->stateDot.row(1); });
  logger.addLogEntry("OldMinimumJerk_target_crackle", [this]() { return (Eigen::Vector3d)this->stateDot.row(2); });
  logger.addLogEntry("OldMinimumJerk_current_err", [this]() { return this->err; });
  logger.addLogEntry("OldMinimumJerk_current_pos", [this]() { return this->pos; });
  logger.addLogEntry("OldMinimumJerk_current_vel", [this]() { return this->vel; });
  logger.addLogEntry("OldMinimumJerk_current_acc", [this]() { return this->acc; });
  logger.addLogEntry("OldMinimumJerk_error_norm", [this]() { return this->eval().norm(); });
}
