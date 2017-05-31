#include "kalman_filter.h"
#include "tools.h"

# define PI 3.14159265358979323846

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
  /**
  ADDED:
    * predict the state
  */

  	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  ADDED:
    * update the state by using Kalman Filter equations
  */
  	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  ADDED:
    * update the state by using Extended Kalman Filter equations
  */

  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];

  double rho = sqrt(px * px + py * py);
  double phi = atan2(py, px);   // atan() leads to inaccurate values with negative x-values
  double rho_dot = (px * vx + py * vy) / (rho);
  VectorXd h = VectorXd(3);
  h << rho, phi, rho_dot;

    VectorXd y = z - h;
    MatrixXd Hj = H_;
    //VectorXd z_pred = Hj * x_;
	//VectorXd y = z - z_pred;
	MatrixXd Ht = Hj.transpose();
	MatrixXd S = Hj * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
    // normalizing

    // normalizing y(1)
	while (y(1) < -PI)
    y(1) += 2 * PI;
    while (y(1) > PI)
    y(1) -= 2 * PI;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
