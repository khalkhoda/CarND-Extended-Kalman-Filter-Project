#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  is_second_meas_ready_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement matrix of Laser measurement
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

   ekf_.x_ = VectorXd(4);
   ekf_.x_ << 1, 1, 1, 1;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    // cout << "EKF: " << endl;

    if (is_second_meas_ready_ == true){
      // Initialize the speed too using a second measurement
      float x0_prev = ekf_.x_(0);
      float y0_prev = ekf_.x_(1);
      init_state_position(measurement_pack);
      cout << "timestamp: " << measurement_pack.timestamp_ <<endl;
      cout << "Second measurement position: " << ekf_.x_.transpose() << endl;
      float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
      previous_timestamp_ = measurement_pack.timestamp_;
      float x0_curr = ekf_.x_(0);
      float y0_curr = ekf_.x_(1);
      float vx0 = (x0_curr - x0_prev) / dt;
      float vy0 = (y0_curr - y0_prev) / dt;
      cout << "x0_prev" << x0_prev << endl;
      cout << "x0_curr" << x0_curr << endl;
      cout << "dt" << dt << endl;
      // Initialize velocity
      ekf_.x_(2) = vx0;
      ekf_.x_(3) = vy0;
      cout << "Second measurement velocity: " << ekf_.x_.transpose() << endl;
      // done initializing, no need to predict or update
      is_initialized_ = true;
      return;
    }

    init_state_position(measurement_pack);
    cout << "timestamp: " << measurement_pack.timestamp_ <<endl;
    cout << "First measurement: " << ekf_.x_.transpose() << endl;
    is_second_meas_ready_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 0.5, 0, 0, 0,
               0, 0.5, 0, 0,
               0, 0, 0.5, 0,
               0, 0, 0, 0.5;
    // done initializing using first meas, no need to predict or update
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // Micro second to seconds conversion
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;



  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ <<  1, 0, dt, 0,
              0, 1, 0, dt,
              0, 0, 1, 0,
              0, 0, 0, 1;
  float noise_ax = 9;
  float noise_ay = 9;

  // set the process covariance matrix Q
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    // Measurement Jacobian matrix
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    cout << "Radar: z = " << measurement_pack.raw_measurements_.transpose() << endl;


  } else {
    // TODO: Laser updates
     ekf_.H_ = H_laser_;
     ekf_.R_ = R_laser_;
     ekf_.Update(measurement_pack.raw_measurements_);
     cout << "Laser: z = " << measurement_pack.raw_measurements_.transpose() << endl;
  }

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}

void FusionEKF::init_state_position(const MeasurementPackage &measurement_pack){
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Convert radar from polar to cartesian coordinates
    //         and initialize state.
    float range = measurement_pack.raw_measurements_[0];
    float phi = measurement_pack.raw_measurements_[1];
    float px0 = range * cos(phi);
    float py0 = range * sin(phi);
    ekf_.x_(0) =  px0;
    ekf_.x_(1) =  py0;
    cout << "Init pos with Radar" << endl;
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // TODO: Initialize state.
    // set the state with the initial location and zero velocity
    ekf_.x_(0) = measurement_pack.raw_measurements_[0];
    ekf_.x_(1) = measurement_pack.raw_measurements_[1];
    cout << "Init pos with Laser" << endl;
  }
}
