#include <stdio.h>

#include <mc_rtc/log/Logger.h>
#include <Eigen/Core>

class MinimumJerk
{
private:
  double dt;
  double T;
  double t;
  double mass;
  double acc_lim;
  double vel_limit;

  Eigen::Vector3d X_f;
  Eigen::Vector3d X_0;
  Eigen::Vector3d ideal_pos;
  Eigen::Vector3d err;
  Eigen::Vector3d pos;
  Eigen::Vector3d vel;
  Eigen::Vector3d acc;

  Eigen::Matrix3d state;

  Eigen::Vector3d acc_target;
  Eigen::Vector3d acc_target_clipped;
  Eigen::Vector3d commanded_acc;

  Eigen::Matrix3d stateDot;

  Eigen::MatrixXd Ri = Eigen::MatrixXd(3, 6);
  Eigen::MatrixXd Rs = Eigen::MatrixXd(6, 3);
  Eigen::Matrix3d RdeltaT = Eigen::Matrix3d();
  Eigen::Matrix3d Ra = Eigen::Matrix3d();
  Eigen::MatrixXd Rb = Eigen::MatrixXd(3, 6);
  Eigen::MatrixXd initConditions = Eigen::MatrixXd(6, 3);

  Eigen::Matrix6d Ri_d = Eigen::Matrix6d();
  Eigen::MatrixXd RdeltaT_d = Eigen::MatrixXd(1, 6);
  Eigen::MatrixXd Ra_d = Eigen::MatrixXd(1, 3);
  Eigen::MatrixXd Rb_d = Eigen::MatrixXd(1, 6);

public:
  MinimumJerk(double dt_, double mass_);

  Eigen::Vector3d eval(void);
  void update(Eigen::Vector3d pos_, Eigen::Vector3d vel_, Eigen::Vector3d acc_, bool discrete = false);
  void continuous_acc(void);
  Eigen::Vector3d continuous_ref_acc(void);
  void discrete_acc(void);
  Eigen::Vector3d discrete_ref_acc(void);
  Eigen::Vector3d compute_vel_max(void);

  void set_target(Eigen::Vector3d target, double duration);
  Eigen::Vector3d get_target_acc(void);
  void set_force_limit(double limit_);
  inline double get_force_limit(void)
  {
    return acc_lim * mass;
  };

  void add_to_logger(mc_rtc::Logger & logger);
};
