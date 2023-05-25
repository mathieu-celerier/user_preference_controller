#include "UserPreferenceController_Initial.h"

#include "../UserPreferenceController.h"

void UserPreferenceController_Initial::configure(const mc_rtc::Configuration & config) {}

void UserPreferenceController_Initial::start(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<UserPreferenceController &>(ctl_);
}

bool UserPreferenceController_Initial::run(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<UserPreferenceController &>(ctl_);
  output("OK");
  return true;
}

void UserPreferenceController_Initial::teardown(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<UserPreferenceController &>(ctl_);
}

EXPORT_SINGLE_STATE("UserPreferenceController_Initial", UserPreferenceController_Initial)
