#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
 public:

  KalmanFilter();

  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param dim_state The number of elements in the state vector ([px, py, vx, vy])
   * @param dim_pred_state The number of elements in the predicted state ([px, py])
   */
  void Init(int dim_state, int dim_pred_state);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

  const Eigen::VectorXd getX() {
    return x_;
  }
  
  void setX(Eigen::VectorXd const& x) {
    x_ = x;
  }

  const Eigen::MatrixXd getP() {
    return P_;
  }
  
  void setP(Eigen::MatrixXd const& p) {
    P_ = p;
  }

  void setF(Eigen::MatrixXd const& f) {
    F_ = f;
  }

  void setQ(Eigen::MatrixXd const& q) {
    Q_ = q;
  }

  void setH(Eigen::MatrixXd const& h) {
    H_ = h;
  }

  void setR(Eigen::MatrixXd const& r) {
    R_ = r;
  }
  
 private:
  
  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transistion matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix
  Eigen::MatrixXd H_;

  // measurement covariance matrix
  Eigen::MatrixXd R_;

  Eigen::MatrixXd I_;
};

#endif /* KALMAN_FILTER_H_ */
