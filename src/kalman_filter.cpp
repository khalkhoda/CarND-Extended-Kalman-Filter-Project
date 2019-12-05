#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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

// Time update (prediction)
void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
   // cout << "Kalman Prediction step" << endl;
   x_ = F_ * x_;
   MatrixXd Ft = F_.transpose();
   // cout << "P_" << P_ << endl;
   // cout << "F_" << F_ << endl;
   // cout << "Q_" << Q_ << endl;
   P_ = F_ * P_ * Ft + Q_;

}

// Laser measurement update
void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
   // cout << "Kalman Laser update step" << endl;
   VectorXd z_pred = H_ * x_;
   VectorXd y = z - z_pred;
   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_;
   MatrixXd S_inv = S.inverse();
   MatrixXd K = P_ * Ht * S_inv;
   x_ = x_ + K * y;
   long dim = x_.size();
   MatrixXd I = MatrixXd::Identity(dim, dim);
   P_ = (I - K * H_) * P_;
}

// Radar measurement update
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
   // std::cout << "Let us talk EKF Update step" << std::endl;
   VectorXd z_pred(3);
   /*State parameters for convenience*/
   float px = x_(0);
   float py = x_(1);
   float vx = x_(2);
   float vy = x_(3);
   float range = sqrt(px*px + py*py);
   z_pred(0) = range; // Predicted range
   z_pred(1) = atan2(py, px); // Predicted angle of arrival
   z_pred(2) = (px*vx + py*vy)/range; // Predicted range rate
   VectorXd y = z - z_pred;

   // Make sure that z(1) - z_pred(1) is in the range [-pi, pi]
   float angle = y(1);
   offsetInRange2pi(angle);
   y(1) = angle;
   // std::cout << "angle" << angle<< std::endl;

   // std::cout << "H_" << H_<< std::endl;
   MatrixXd Ht = H_.transpose(); // Here H_ is the Jacobian matrix
   // std::cout << "P_" << P_<< std::endl;
   MatrixXd S = H_ * P_ * Ht + R_;
   MatrixXd S_inv = S.inverse();
   MatrixXd K = P_ * Ht * S_inv;
   x_ = x_ + K * y;
   long dim = x_.size();
   MatrixXd I = MatrixXd::Identity(dim, dim);
   P_ = (I - K * H_) * P_;
}

void KalmanFilter::offsetInRange2pi(float &angle){
  while (angle < -M_PI){
    angle += 2*M_PI;
  }
  while (angle > M_PI){
    angle -= 2*M_PI;
  }
  return;
}
