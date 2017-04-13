#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(int dim_state, int dim_pred_state) {
  x_ = VectorXd(dim_state);
  P_ = MatrixXd(dim_state, dim_state);
  F_ = MatrixXd(dim_state, dim_state);
  H_ = MatrixXd(dim_pred_state, dim_state);
  R_ = MatrixXd(dim_pred_state, dim_pred_state);
  Q_ = MatrixXd(dim_state, dim_state);
  I_ = MatrixXd::Identity(dim_state, dim_state);
}

void KalmanFilter::Predict() {
  // Cartesian coordinates only
  // Linear prediction step does not require Fj
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Laser
  VectorXd y = z - H_ * x_;
  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();
  x_ = x_ + K * y;
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(VectorXd const& z) {
  // Radar
  //   Convert cartesian coordinates to polar using h(x')
  VectorXd h_x = Tools::CartesianToPolar(x_);
  VectorXd y = z - h_x;
  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();
  x_ = x_ + K * y;
  P_ = (I_ - K * H_) * P_;  
}
