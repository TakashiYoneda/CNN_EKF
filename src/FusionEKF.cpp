#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

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
  // process and measurement noises
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;

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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      //vector for conversion ( rho, phi -> px, py )
      //px = x_state[0], py = x_state[1]
      VectorXd x_state(2);
      x_state(0) = measurement_pack.raw_measurements_[0] * cos (measurement_pack.raw_measurements_[1]);
      x_state(1) = measurement_pack.raw_measurements_[0] * sin (measurement_pack.raw_measurements_[1]);
      ekf_.x_ << x_state(0), x_state(1), 0, 0;

  		previous_timestamp_ = measurement_pack.timestamp_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
  		previous_timestamp_ = measurement_pack.timestamp_;
    }

    ekf_.P_ = MatrixXd(4, 4);


// 0.0983, 0.0852, 0.4088, 0.4698
    ekf_.P_ << 10, 0, 0, 0,
              0, 10, 0, 0,
              0, 0, 10, 0,
              0, 0, 0, 10;

/* 0.1153, 0.1249, 0.4886, 0.4269
    ekf_.P_ << 10, 0, 10, 0,
              0, 10, 0, 10,
              10, 0, 10, 0,
              0, 10, 0, 10;
*/

/* 0.0977, 0.0854, 0.4529, 0.4717
    ekf_.P_ << 100, 0, 0, 0,
              0, 100, 0, 0,
              0, 0, 100, 0,
              0, 0, 0, 100;
*/
/* 0.0973, 0.0855, 0.4656, 0.4723
    ekf_.P_ << 1000, 0, 0, 0,
              0, 1000, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
*/
/* 0.1376, 0.0863, 0.4801, 0.3997
    ekf_.P_ << 0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0;
*/
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

   //compute the time elapsed between the current and previous measurements
 	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
 	previous_timestamp_ = measurement_pack.timestamp_;

 	float dt_2 = dt * dt;
 	float dt_3 = dt_2 * dt;
 	float dt_4 = dt_3 * dt;

 	//Modify the F matrix so that the time is integrated
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ <<  1, 0, dt, 0,
              0, 1,  0, dt,
              0, 0,  1, 0,
              0, 0,  0, 1;

 	//set the process covariance matrix Q
//  float noise_ax = 9;
//	float noise_ay = 9;
  float noise_ax = 9;
	float noise_ay = 9;

 	ekf_.Q_ = MatrixXd(4, 4);
 	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
 			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
 			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
 			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

 	//predict
 	ekf_.Predict();


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /*
  TODO:
   * Use the sensor type to perform the update step.
   * Update the state and covariance matrices.
 */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);                   // 3dim
  } else {
    // Laser updates
    ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
    ekf_.Update(measurement_pack.raw_measurements_);                      // 4dim
  }

  // print the output
  cout << "Updated: x_ = " <<"\n"<< ekf_.x_ << endl;
  cout << "Updated: P_ = " <<"\n"<< ekf_.P_ << endl;
}
