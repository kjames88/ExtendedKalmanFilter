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
      
      std::cout << "    estimate " << std::endl << estimations[i] << std::endl << "    ground truth " << std::endl << ground_truth[i] << std::endl;
      
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
  /**
  TODO:
    * Calculate a Jacobian here.
  */
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
  std::cout << "phi is " <<  phi << std::endl;
  
  float rho_dot = ((px * vx) + (py * vy)) / rho;
  VectorXd polar = VectorXd(3);
  polar << rho, phi, rho_dot;
  return polar;
}
