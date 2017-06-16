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
  H_radar_ = MatrixXd(3, 4); //jacobian
  
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
  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::Init(const MeasurementPackage &measurement_pack){
  previous_timestamp_ = 0;
  VectorXd x = VectorXd(4);
  if ( measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    x << tools_.Polar2Cartesian(measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2]);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    x << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],0,0;
  }

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
  
  cout<<"sensor_type: "<<measurement_pack.sensor_type_<<"\n";
  cout<<"raw_measurements: "<<measurement_pack.raw_measurements_<<"\n";
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    Init(measurement_pack);
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

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<
  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
  0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
  dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
  0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if ( type == MeasurementPackage::RADAR) {
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    
  } else if ( type == MeasurementPackage::LASER) {
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
    
    
    
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}


/*
 
 #include "Dense"
 #include <iostream>
 #include "tracking.h"
 
 using namespace std;
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 using std::vector;
 
 Tracking::Tracking() {
	is_initialized_ = false;
	previous_timestamp_ = 0;
 
	//create a 4D state vector, we don't know yet the values of the x state
	kf_.x_ = VectorXd(4);
 
	//state covariance matrix P
	kf_.P_ = MatrixXd(4, 4);
	kf_.P_ << 1, 0, 0, 0,
 0, 1, 0, 0,
 0, 0, 1000, 0,
 0, 0, 0, 1000;
 
 
	//measurement covariance
	kf_.R_ = MatrixXd(2, 2);
	kf_.R_ << 0.0225, 0,
 0, 0.0225;
 
	//measurement matrix
	kf_.H_ = MatrixXd(2, 4);
	kf_.H_ << 1, 0, 0, 0,
 0, 1, 0, 0;
 
	//the initial transition matrix F_
	kf_.F_ = MatrixXd(4, 4);
	kf_.F_ << 1, 0, 1, 0,
 0, 1, 0, 1,
 0, 0, 1, 0,
 0, 0, 0, 1;
 
	//set the acceleration noise components
	noise_ax = 5;
	noise_ay = 5;
 
 }
 
 Tracking::~Tracking() {
 
 }
 
 // Process a single measurement
 void Tracking::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
	if (!is_initialized_) {
 //cout << "Kalman Filter Initialization " << endl;
 
 //set the state with the initial location and zero velocity
 kf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
 
 previous_timestamp_ = measurement_pack.timestamp_;
 is_initialized_ = true;
 return;
	}
 
	//compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;
 
	float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;
 
	//Modify the F matrix so that the time is integrated
	kf_.F_(0, 2) = dt;
	kf_.F_(1, 3) = dt;
 
	//set the process covariance matrix Q
	kf_.Q_ = MatrixXd(4, 4);
	kf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
 0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
 dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
 0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
 
	//predict
	kf_.Predict();
 
	//measurement update
	kf_.Update(measurement_pack.raw_measurements_);
 
	std::cout << "x_= " << kf_.x_ << std::endl;
	std::cout << "P_= " << kf_.P_ << std::endl;
 
 }

 
 */

