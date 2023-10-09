#pragma once

#include <mc_control/fsm/Controller.h>
#include <mc_tasks/EndEffectorTask.h>

#include "api.h"

struct UserPreferenceController_DLLAPI UserPreferenceController : public mc_control::fsm::Controller
{
  UserPreferenceController(mc_rbdyn::RobotModulePtr rm, double dt, const mc_rtc::Configuration & config);

  bool run() override;

  void reset(const mc_control::ControllerResetData & reset_data) override;

  std::shared_ptr<mc_tasks::EndEffectorTask> eeTask;

private:
  mc_rtc::Configuration config_;
};
