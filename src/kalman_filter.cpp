#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // Cartesian coordinates only
  // Linear prediction step does not require Fj

  std::cout << "x_ = " << std::endl << x_ << std::endl;
  
  x_ = F_ * x_;

  std::cout << "x_' = " << std::endl << x_ << std::endl;
  std::cout << "Q_  = " << std::endl << Q_ << std::endl;
  
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Laser

  std::cout << "z = " << std::endl << z << std::endl;
  
  VectorXd y = z - H_ * x_;

  std::cout << "y = " << std::endl << y << std::endl;
  
  MatrixXd S = H_ * P_ * H_.transpose() + R_;

  std::cout << "S = " << std::endl << S << std::endl;
  
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  std::cout << "K = " << std::endl << K << std::endl;
  
  x_ = x_ + K * y;
  long dim = x_.size();
  MatrixXd I = MatrixXd::Identity(dim, dim);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(VectorXd const& z) {
  // Radar
  //   Convert cartesian coordinates to polar using h(x')
  VectorXd h_x = Tools::CartesianToPolar(x_);

  std::cout << "h_x = " << std::endl << h_x << std::endl;
  std::cout << "z = " << std::endl << z << std::endl;
  
  VectorXd y = z - h_x;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  long dim = x_.size();
  MatrixXd I = MatrixXd::Identity(dim, dim);
  P_ = (I - K * H_) * P_;  
}
