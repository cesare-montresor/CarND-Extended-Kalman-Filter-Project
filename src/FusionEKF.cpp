#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <sys/time.h>

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
  
  //H
  H_laser_ <<
  1,0,0,0,
  0,1,0,0;
  
  //measurement covariance matrix - laser
  R_laser_ <<
    0.0225, 0,
    0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ <<
    0.09, 0     , 0,
    0   , 0.0009, 0,
    0   , 0     , 0.09;
    
  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;
  


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::Init(const MeasurementPackage &measurement_pack){
  previous_timestamp_ = measurement_pack.timestamp_;
  VectorXd x = VectorXd(4);
  if ( measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    x << tools_.Polar2Cartesian(measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2]);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    x << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],4,2;
  }
  //cout<<"Init x:"<<x<<"\n";

  //state covariance matrix P
  MatrixXd P = MatrixXd(4, 4);
  P <<
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1000, 0,
  0, 0, 0, 1000;
  
  //the initial transition matrix F_
  MatrixXd F = MatrixXd(4, 4);
  F <<
  1, 0, 1, 0,
  0, 1, 0, 1,
  0, 0, 1, 0,
  0, 0, 0, 1;

  MatrixXd Q = MatrixXd(4, 4);
  Q <<
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0;

  ekf_.Init(x, P, F, H_laser_, R_laser_, Q);

  is_initialized_ = true;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  VectorXd measures = measurement_pack.raw_measurements_;
  MeasurementPackage::SensorType type = measurement_pack.sensor_type_;
  long long timestamp = measurement_pack.timestamp_;
  
  /*****************************************************************************
   *  Init
   ****************************************************************************/

  if (!is_initialized_) {
    Init(measurement_pack);
    return;
  }
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (timestamp - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = timestamp;
  
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  ekf_.Q_ <<
  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
  0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
  dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
  0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

   if ( type == MeasurementPackage::RADAR) {
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measures);
  } else if ( type == MeasurementPackage::LASER) {
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measures);
  }
    
    
    
  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
