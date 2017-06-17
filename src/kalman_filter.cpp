#include "kalman_filter.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <math.h>

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
  
  cnt = 0;
}

void KalmanFilter::Predict() {
  /**
   TODO:
   * predict the state
   */
  
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  //std::cout<<"x: \t"<<x_[0]<<"\t"<<x_[1]<<"\t"<<x_[2]<<"\t"<<x_[3]<<"\n";
  cnt++;
  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   TODO:
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
  //std::cout<<"x: \t"<<x_[0]<<"\t"<<x_[1]<<"\t"<<x_[2]<<"\t"<<x_[3]<<"\n";
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float distance = std::pow(std::pow(px,2)+std::pow(py,2),0.5);
  
  float angle = atan2f(py,px);
  float acc = ((px*vx)+(py*vy))/distance;
  
  VectorXd z_pred(3);
  z_pred << distance,angle,acc;
  VectorXd y = z - z_pred;
  
  while( y[1] < -M_PI || y[1] > M_PI){
    std::cout<<"Normalizing: "<<y[1]<<"\n";
    y[1] += 2 * M_PI * (y[1] < -M_PI?1:-1);
  }
  
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
//std::cout<<cnt<<":\t"<<(int)x_[0]<<" \t"<<(int)x_[1]<<" \t"<<(int)x_[2]<<" \t"<<(int)x_[3]<<" \t|\t"<<(int)distance<<" \t "<<(int)angle<<" \t"<<(int)acc<<"\n";
  x_ = x_ + (K * y);
  
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

