#include "UserPreferenceController_Switch.h"

#include "../UserPreferenceController.h"

void UserPreferenceController_Switch::configure(const mc_rtc::Configuration & config) {}

void UserPreferenceController_Switch::start(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<UserPreferenceController &>(ctl_);
  auto & robot = ctl.robot("kinova");

  mj = new MinimumJerk(ctl.timeStep, robot.mb().bodies()[robot.bodyIndexByName("bracelet_link")].inertia().mass());
  mj->set_force_limit(4);
  first_target = robot.bodyPosW("tool_frame").translation() + Eigen::Vector3d(0.0,0.3,0.0);
  second_target = robot.bodyPosW("tool_frame").translation();
  current_target = first_target;
  ctl.eeTask->positionTask->position(first_target);
  mj->set_target(first_target,4);
  commanded_pos = robot.bodyPosW("tool_frame").translation();
  commanded_vel.setZero();
  targetIsFirst = true;

  ctl.eeTask->positionTask->stiffness(0);
  ctl.eeTask->positionTask->weight(10000);
  ctl.eeTask->orientationTask->weight(0);
  ctl.solver().addTask(ctl.eeTask);

  mj->add_to_logger(ctl.logger());

  eeAccel.setZero();

  ctl.logger().addLogEntry("body6d_kinova_tool_frame_position_curAccel", [this]() { return this->eeAccel; });
  ctl.logger().addLogEntry("commanded_vel", [this]() { return this->commanded_vel; });
  ctl.logger().addLogEntry("commanded_pos", [this]() { return this->commanded_pos; });
  ctl.logger().addLogEntry("acceleration_derivative", [this]() { return this->acc_deriv; });
  ctl.logger().addLogEntry("acceleration_body_transformed", [this]() { return this->acc; });

  ctl.gui()->addElement({"Controller","MinimumJerk"},
    mc_rtc::gui::NumberInput("Max force [N.m]", [this]() { return this->mj->get_force_limit(); }, [this](double l){ this->mj->set_force_limit(l); })
  );

  ctl.datastore().assign<std::string>("ControlMode","Torque");   
  
  mc_rtc::log::success("[UserPreferenceController] Switch state init done");
}

bool UserPreferenceController_Switch::run(mc_control::fsm::Controller & ctl_)
{
  // std::cout << "================================================================================" << std::endl;
  // std::cout << "=                             Running new run step                             =" << std::endl;
  // std::cout << "================================================================================" << std::endl;

  auto & ctl = static_cast<UserPreferenceController &>(ctl_);
  auto & robot = ctl.robot("kinova");

  auto bodyName = robot.frame("tool_frame").body();

  sva::PTransformd transform(robot.bodyPosW(bodyName));
  
  Eigen::Vector3d pos = transform.translation();
  acc_deriv = (robot.bodyVelW(bodyName).linear() - vel)/ctl.timeStep;
  vel = robot.bodyVelW(bodyName).linear();
  Eigen::Vector3d angVel = robot.bodyVelW(bodyName).angular();
  acc = transform.rotation().transpose() * robot.bodyAccB(bodyName).linear() + angVel.cross(vel);
  
  eeAccel = acc;
  // std::cout << "acc" << acc_deriv.transpose() << std::endl;
  // std::cout << "bodyAccW" << acc.transpose() << std::endl;
  // std::cout << "bodyAccB " << robot.bodyAccB(bodyName).linear().transpose() << std::endl;
  // std::cout << "Axe-Angle " << Eigen::AngleAxisd(transform.rotation().transpose()).axis().transpose() << " " << Eigen::AngleAxisd(transform.rotation().transpose()).angle() << std::endl;

  mj->update(pos, vel, acc);

  commanded_acc = mj->get_target_acc();

  // std::cout << "commanded acc " << commanded_acc.transpose() << std::endl;

  // commanded_vel += commanded_acc * ctl.timeStep;
  // commanded_pos += commanded_vel * ctl.timeStep;

  ctl.eeTask->positionTask->refAccel(commanded_acc);
  // ctl.eeTask->positionTask->refVel(commanded_vel);
  // ctl.eeTask->positionTask->position(commanded_pos);

  if (mj->eval().norm() < 0.001)
  {
    switch_target(ctl_);
  }

  // mc_rtc::log::error_and_throw<std::runtime_error>("");

  output("OK");
  return true;
}

void UserPreferenceController_Switch::teardown(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<UserPreferenceController &>(ctl_);
}

void UserPreferenceController_Switch::switch_target(mc_control::fsm::Controller & ctl_)
{
  auto & ctl = static_cast<UserPreferenceController &>(ctl_);

  if (targetIsFirst)
  {
    current_target = second_target;
    mj->set_target(second_target,4);
    commanded_pos = ctl.robot("kinova").bodyPosW("tool_frame").translation();
    commanded_vel.setZero();
    ctl.eeTask->positionTask->position(second_target);
    targetIsFirst = false;
  }
  else
  {
    current_target = first_target;
    mj->set_target(first_target,4);
    commanded_pos = ctl.robot("kinova").bodyPosW("tool_frame").translation();
    commanded_vel.setZero();
    ctl.eeTask->positionTask->position(first_target);
    targetIsFirst = true;
  }
  mc_rtc::log::info("[UserPreferenceController][Switch] Switching target");
}

EXPORT_SINGLE_STATE("UserPreferenceController_Switch", UserPreferenceController_Switch)