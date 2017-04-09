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

// VectorXd FusionEKF::PolarToCartesian(VectorXd polar) {
//   float rho = polar[0];
//   float phi = polar[1];
//   float rho_dot = polar[2];
//   float px = rho * cos(phi);
//   float py = rho * sin(phi);
//   float vx = rho_dot * cos(phi);
//   float vy = rho_dot * sin(phi);
//   VectorXd cartesian = VectorXd(4);
//   cartesian << px, py, vx, vy;
//   return cartesian;
// }

// VectorXd FusionEKF::CartesianToPolar(VectorXd cartesian) {
//   float px = cartesian[0];
//   float py = cartesian[1];
//   float vx = cartesian[2];
//   float vy = cartesian[3];
//   float rho = sqrt((px * px) + (py * py));
//   float phi = atan(py/px);
//   float rho_dot = ((px * vx) + (py * vy)) / rho;
//   VectorXd polar = VectorXd(3);
//   polar << rho, phi, rho_dot;
//   return polar;
// }

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement

    cout << "EKF: " << measurement_pack.raw_measurements_ << endl;
    ekf_.x_ = VectorXd(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
	 Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd meas = VectorXd(3);
      meas << measurement_pack.raw_measurements_[0],
	measurement_pack.raw_measurements_[1],
	measurement_pack.raw_measurements_[2];
      ekf_.x_ = Tools::PolarToCartesian(meas);

      cout << "radar x_: " << ekf_.x_ << endl;
    } else {
      /**
	 Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;  // process covariance matrix initialization
    
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
   TODO:
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


  cout << "time in = " << measurement_pack.timestamp_ << " previous time = " << previous_timestamp_ << endl;
  
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
    float px = ekf_.x_[0];
    float py = ekf_.x_[1];
    float vx = ekf_.x_[2];
    float vy = ekf_.x_[3];
    float d = (px * px) + (py * py);
    float d_12 = sqrt(d);
    float d_32 = d * d_12;
    float n = (vx * py) - (vy * px);

    if (fabs(d) > 0.0001) {
      Hj_ << px / d_12, py / d_12, 0, 0,
	-py / d, px / d, 0, 0,
	(py * n) / d_32, (-px * n) / d_32, px / d_12, py / d_12;
    } else {
      Hj_ << 0, 0, 0, 0,
	0, 0, 0, 0,
	0, 0, 0, 0;
    }
    ekf_.F_ = F;
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.Q_ = Q;
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    /**
       Initialize state.
    */
    ekf_.F_ = F;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Q_ = Q;
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

    cout << "RADAR " << endl << measurement_pack.raw_measurements_ << endl;
    
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
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
