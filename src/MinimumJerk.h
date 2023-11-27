#include <stdio.h>

#include <mc_rtc/log/Logger.h>
#include <Eigen/Core>

#include <root_finder/root_finder.hpp>

#define EPSILON 1e-3

using namespace Eigen;

typedef Matrix<double, 6, 3> Matrix63d;
typedef Matrix<double, 3, 6> Matrix36d;
typedef std::set<double> SetD;

class MinimumJerk
{
private:
  // Control hyper-parameters
  double A_limit; // Maximum linear acceleration
  double J_limit; // Maximum linear jerk
  double Kd;

  // Control variables
  double dt;
  Vector3d err;
  double acc_mj;
  Vector3d target_acc;
  double commanded_jerk;
  std::string delta_recompute_source;
  std::string output_source;

  // Trajectory specific hyper-parameters
  Vector3d X_f;
  Vector3d X_0;
  Vector3d vel;
  double V_0;
  double A_0;
  double A_0_dist;
  double deltaT; // Remaining time for trajectory

  // Trajectory specific variables
  Vector3d X_s; // Current state [X,V,A]^T
  Vector6d u; // Input of state-space, here X_f
  Vector6d poly_coeffs; // Polynomial coefficients computed
  Matrix36d Rs; // Scaling Matrix | Normalized Polynomial -> Scaled Polynomial
  Matrix63d Rn; // Normalizing matrix | Scaled Polynomial -> Normalized Polynomial
  Matrix6d Ri; // Interpolation matrix
  Matrix36d Ris; // Interpolation matrix for state-space
  Matrix3d A; // State matrix of state-space (dX = [A]*X + B*u)
  Matrix3d B; // Input matrix of state-space (dX = A*X + [B]*u)

public:
  MinimumJerk(double dt_, double mass_);

  // Getters
  Vector3d getTargetAcceleration(void);
  Vector3d eval(void);

  // Setters
  void setTargetPos(Vector3d pos_);
  void setDuration(double T);
  void setJerkLimit(double limit);
  void setAccelerationLimit(double limit);
  void setDampingGain(double gain);

  // Member functions
  void update(Vector3d pos_, Vector3d vel_, Vector3d acc_, Vector3d acc_dist_);
  Vector3d decisionTreeControl(void);
  Vector3d computeTargetAcc(bool onlydDamp);
  double computeJerkAtZero(double deltaT_target);
  double solveNewDelta(void);
  SetD solveJerkCst(void);

  void add_to_logger(mc_rtc::Logger & logger);
};
