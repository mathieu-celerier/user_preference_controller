#include "UserPreferenceController.h"

UserPreferenceController::UserPreferenceController(mc_rbdyn::RobotModulePtr rm, double dt, const mc_rtc::Configuration & config)
: mc_control::fsm::Controller(rm, dt, config)
{
  dynamicsConstraint = mc_rtc::unique_ptr<mc_solver::DynamicsConstraint>(new mc_solver::DynamicsConstraint(robots(), 0, solver().dt(), {0.1, 0.01, 0.5}, 0.9, false, false));
  solver().addConstraintSet(dynamicsConstraint);

  eeTask = std::make_shared<mc_tasks::EndEffectorTask>(robot().frame("tool_frame"));

  getPostureTask("kinova")->weight(1);
  getPostureTask("kinova")->stiffness(1);

  datastore().make<std::string>("ControlMode", "Position");

  mc_rtc::log::success("UserPreferenceController init done ");
}

bool UserPreferenceController::run()
{
  auto ctrl_mode = datastore().get<std::string>("ControlMode");

  if (ctrl_mode.compare("Position") == 0)
  {
    return mc_control::fsm::Controller::run(mc_solver::FeedbackType::OpenLoop);
  }
  else
  {
    return mc_control::fsm::Controller::run(mc_solver::FeedbackType::ClosedLoopIntegrateReal);
  }
}

void UserPreferenceController::reset(const mc_control::ControllerResetData & reset_data)
{
  mc_control::fsm::Controller::reset(reset_data);
}
