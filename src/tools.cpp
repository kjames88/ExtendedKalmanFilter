#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse = VectorXd(4);  // 4d state
  rmse << 0, 0, 0, 0;
  if (estimations.size() > 0 &&
      estimations.size() == ground_truth.size()) {
    VectorXd sum = VectorXd(4);
    sum << 0, 0, 0, 0;
    for (int i=0; i < estimations.size(); ++i) {
      VectorXd d = estimations[i] - ground_truth[i];
      VectorXd d_square = d.array() * d.array();
      sum += d_square;
    }
    VectorXd mean = sum / estimations.size();
    rmse = mean.array().sqrt();
  }
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj = MatrixXd(3,4);
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];
  float d = (px * px) + (py * py);
  float d_12 = sqrt(d);
  float d_32 = d * d_12;
  float n = (vx * py) - (vy * px);

  if (fabs(d) > 0.0001) {
    Hj << px / d_12, py / d_12, 0, 0,
      -py / d, px / d, 0, 0,
      (py * n) / d_32, (-px * n) / d_32, px / d_12, py / d_12;
  } else {
    Hj << 0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0;
  }
  return Hj;
}

VectorXd Tools::PolarToCartesian(VectorXd const& polar) {
  float rho = polar[0];
  float phi = polar[1];
  float rho_dot = polar[2];
  float px = rho * cos(phi);
  float py = rho * sin(phi);
  float vx = rho_dot * cos(phi);
  float vy = rho_dot * sin(phi);
  VectorXd cartesian = VectorXd(4);
  cartesian << px, py, vx, vy;
  return cartesian;
}

VectorXd Tools::CartesianToPolar(VectorXd const& cartesian) {
  float px = cartesian[0];
  float py = cartesian[1];
  float vx = cartesian[2];
  float vy = cartesian[3];
  float rho = sqrt((px * px) + (py * py));
  
  double phi = atan2(py,px);
  
  float rho_dot = (rho > 0.0001) ? ((px * vx) + (py * vy)) / rho : 0.0;
  VectorXd polar = VectorXd(3);
  polar << rho, phi, rho_dot;
  return polar;
}
