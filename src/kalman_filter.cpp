#include "kalman_filter.h"
#include <iostream>
#include "math.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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

	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

	VectorXd y_ = z - H_ * x_;
	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = (H_ * PHt) + R_;
	MatrixXd Sinv = S.inverse();
	MatrixXd K = PHt * Sinv;
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);

	x_ += K * y_;
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

	float x = x_[0];
	float y = x_[1];
	float vx = x_[2];
	float vy = x_[3];

	float rho = sqrt(x*x + y * y);
	//ToDo: angle should be (-pi, pi)
	float theta = atan2(y, x);

	if (fabs(rho) < 0.001) {
		rho = 0.001;
	}

	float rho_dot = (x * vx + y * vy) / rho;

	VectorXd Z_pred = VectorXd(3);
	Z_pred << rho, theta, rho_dot;

	VectorXd y_ = z - Z_pred;
	float delta_theta = y_[1];
	
	while (delta_theta > M_PI)
		delta_theta -= 2* M_PI;
	while (delta_theta < -M_PI)
		delta_theta += 2* M_PI;

	y_[1] = delta_theta;

	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = (H_ * PHt) + R_;
	MatrixXd Sinv = S.inverse();
	MatrixXd K = PHt * Sinv;
	MatrixXd I = MatrixXd::Identity(4, 4);

	x_ += K * y_;
	P_ = (I - K * H_) * P_;
}
