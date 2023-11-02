#include <stdio.h>

#include <mc_rtc/log/Logger.h>
#include <Eigen/Core>

#include <root_finder/root_finder.hpp>

using namespace Eigen;

typedef Matrix<double, 6, 3> Matrix63d;
typedef Matrix<double, 3, 6> Matrix36d;
typedef std::set<double> SetD;

class MinimumJerk
{
private:
  // Control hyper-parameters
  double V_limit; // Maximum linear velocity
  double A_limit; // Maximum linear acceleration
  double J_limit; // Maximum linear jerk

  // Control variables
  double dt;
  Vector3d err;
  Vector3d target_acc;
  Vector3d commanded_jerk;
  std::string output_source;

  // Trajectory specific hyper-parameters
  Vector3d X_f;
  Vector3d X_0;
  Vector3d V_0;
  Vector3d A_0;
  double deltaT; // Remaining time for trajectory

  // Trajectory specific variables
  Matrix3d X_s; // Current state [X,V,A]^T
  Matrix63d u; // Input of state-space, here X_f
  Matrix63d poly_coeffs; // Polynomial coefficients computed
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
  void setVelocityLimit(double limit);
  void setAccelerationLimit(double limit);
  void setJerkLimit(double limit);

  // Member functions
  void update(Vector3d pos_, Vector3d vel_, Vector3d acc_);
  Vector3d decisionTreeControl(void);
  Vector3d computeTargetAcc(void);
  Vector3d computeVelAtRoot(double t);
  Vector3d computeAccAtRoot(double t);
  Vector3d computeJerkAtRoot(double t);
  SetD computeRootsOfJerk(double deltaT);
  SetD solveJerkConstraint(void);
  SetD velTrajectoryExtremums(void);
  double solveNewDeltaT(void);

  void add_to_logger(mc_rtc::Logger & logger);
};
