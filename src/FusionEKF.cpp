#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  
  H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;
  
  noise_ax_ = 9;
  noise_ay_ = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    ekf_.Init(4, 2);
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement

    cout << "EKF: " << measurement_pack.raw_measurements_ << endl;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
	 Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd meas = VectorXd(3);
      meas << measurement_pack.raw_measurements_[0],
	measurement_pack.raw_measurements_[1],
	measurement_pack.raw_measurements_[2];
      ekf_.setX(Tools::PolarToCartesian(meas));
    } else {
      /**
	 Initialize state.
      */
      VectorXd x(4);
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      ekf_.setX(x);
    }

    MatrixXd P(4, 4);
    P << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;  // process covariance matrix initialization
    ekf_.setP(P);
    
    previous_timestamp_ = measurement_pack.timestamp_;


    cout << "INIT DONE" << endl;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // Laser Measurements Part 3 shows 2D Laser KF
  
  // F = 1, 0, delta_t, 0,
  //     0, 1, 0, delta_t,
  //     0, 0, 1, 0,
  //     0, 0, 0, 1
  
  // Laser Measurements Part 3 defines sigma_ax^2 = noise_ax, sigma_ay^2 = noise_ay
  //   Q = delta_t^4/4 * noise_ax, 0, delta_t^3/2 * noise_ax, 0
  //       0, delta_t^4/4 * noise_ay, 0, delta_t^3/2 * noise_ay
  //       delta_t^3/2 * noise_ax, 0, delta_t^2 * noise_ax, 0
  //       0, delta_t^3/2 * noise_ay 0, delta_t^2 * noise_ay


  float dt = static_cast<float>(measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  // usec to sec
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = (dt * dt);
  float dt_3 = dt_2 * dt;
  float dt_4 = (dt_2 * dt_2);
  
  MatrixXd F = MatrixXd(4, 4);
  F << 1, 0, dt, 0,
    0, 1, 0, dt,
    0, 0, 1, 0,
    0, 0, 0, 1;
    
  MatrixXd Q = MatrixXd(4, 4);
  Q << (dt_4/4) * noise_ax_, 0, (dt_3/2) * noise_ax_, 0,
    0, (dt_4/4) * noise_ay_, 0, (dt_3/2) * noise_ay_,
    (dt_3/2) * noise_ax_, 0, dt_2 * noise_ax_, 0,
    0, (dt_3/2) * noise_ay_, 0, dt_2 * noise_ay_;
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /**
       Convert radar from polar to cartesian coordinates and initialize state.
    */
    Hj_ = Tools::CalculateJacobian(ekf_.getX());
    ekf_.setF(F);
    ekf_.setH(Hj_);
    ekf_.setR(R_radar_);
    ekf_.setQ(Q);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    /**
       Initialize state.
    */
    ekf_.setF(F);
    ekf_.setH(H_laser_);
    ekf_.setR(R_laser_);
    ekf_.setQ(Q);
  }

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    VectorXd z = VectorXd(3);
    z << measurement_pack.raw_measurements_[0],
      measurement_pack.raw_measurements_[1],
      measurement_pack.raw_measurements_[2];
    ekf_.UpdateEKF(z);
  } else {
    // Laser updates
    VectorXd z = VectorXd(2);
    z << measurement_pack.raw_measurements_[0],
      measurement_pack.raw_measurements_[1];
    ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.getX() << endl;
  cout << "P_ = " << ekf_.getP() << endl;
}
