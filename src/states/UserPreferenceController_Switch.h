#pragma once

#include <mc_control/fsm/State.h>
#include "../MinimumJerk.h"

struct UserPreferenceController_Switch : mc_control::fsm::State
{
  void configure(const mc_rtc::Configuration & config) override;

  void start(mc_control::fsm::Controller & ctl) override;

  bool run(mc_control::fsm::Controller & ctl) override;

  void teardown(mc_control::fsm::Controller & ctl) override;

private:
  MinimumJerk* mj;

  bool targetIsFirst;
  Eigen::Vector3d first_target;
  Eigen::Vector3d second_target;
  Eigen::Vector3d current_target;

  Eigen::Vector3d vel;
  Eigen::Vector3d acc_deriv;
  Eigen::Vector3d acc;

  Eigen::Vector3d commanded_acc;
  Eigen::Vector3d commanded_vel;
  Eigen::Vector3d commanded_pos;

  Eigen::Vector3d eeAccel;

  void switch_target(mc_control::fsm::Controller & ctl_);

};
