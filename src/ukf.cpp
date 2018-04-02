#include "Eigen/Dense"
#include <iostream>
#include "ukf.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  x_ << 0., 0., 0., 0., 0.;

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  weights_ = VectorXd( 2*n_aug_+1 );
  weights_(0) = lambda_/( lambda_ + n_aug_ );
  for( int i=1; i<2*n_aug_+1; i++ )
     weights_(i) = 0.5/( n_aug_ + lambda_ );
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {

    if( meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {
      x_(0) =  meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  // Already initialized

  // Check the time delta
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; // in seconds
  time_us_ = meas_package.timestamp_;

  // Run prediction step
  Prediction(dt);

  // Run measurement step
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      cout << "Got LASER measurement" << endl;
      UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    cout << "Got RADAR measurement" << endl;
    UpdateRadar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  cout << "Started Prediction" << endl;

  // Augmented stated vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Augmented State covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;


  // Sigma point matrix for augmented state vector
  MatrixXd Xsig_aug = CreateSignmaPoints(x_aug, P_aug, n_aug_, 2 * n_aug_ + 1);

  // ------------------ Predict Sigma points -------------

  Xsig_pred_ = MatrixXd( n_x_, 2*n_aug_+1 );
  for (int i = 0; i< 2*n_aug_+1; i++) {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;
    //avoid division by zero
    if (fabs(yawd) > nearZero_) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // Sigma points to mean state conversion
  x_.fill(0.0);
  x_ = Xsig_pred_ * weights_;
  // State covariance matrix from sigma points
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  // PHEWWWW!! Done
  cout << "End Prediction" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd z = meas_package.raw_measurements_;
  int n_lidar = 2;

  // We need to transform sigma points to measurement space (px, py) for lidar
  MatrixXd Zsig_transformed = MatrixXd(n_lidar, 2 * n_aug_ + 1);
  for( int i = 0; i < 2*n_aug_ + 1; i++ )
  {
    Zsig_transformed(0,i) = Xsig_pred_(0,i);
    Zsig_transformed(1,i) = Xsig_pred_(1,i);
  }

  VectorXd z_pred = VectorXd(n_lidar);
  z_pred.fill(0.);
  for( int i = 0; i < 2*n_aug_ + 1; i++ ) {
    z_pred = z_pred + weights_(i)*Zsig_transformed.col(i);
  }

  //measurement covariance matrix
  MatrixXd S = MatrixXd(n_lidar, n_lidar);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig_transformed.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise. Measurement noise is just additive.
  MatrixXd R = MatrixXd(n_lidar, n_lidar);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S = S + R;

  // ---- KF Update -----

  //cross correlation Tc
   MatrixXd Tc = MatrixXd(n_x_, n_lidar);

   Tc.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd z_diff = Zsig_transformed.col(i) - z_pred;
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }

   MatrixXd K = Tc * S.inverse();
   VectorXd z_diff = z - z_pred;

   x_ = x_ + K * z_diff;
   P_ = P_ - K*S*K.transpose();

   //NIS
   NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
   cout << "================= NIS_lidar_: " << NIS_lidar_ << "================= " << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  VectorXd z = meas_package.raw_measurements_;
  int n_radar = 3;

  // We need to transform sigma points to measurement space for radar
  MatrixXd Zsig_transformed = MatrixXd(n_radar, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    Zsig_transformed(0, i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig_transformed(1, i) = atan2(p_y, p_x);
    Zsig_transformed(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
  }

  VectorXd z_pred = VectorXd(n_radar);
  z_pred.fill(0.);
  for( int i = 0; i < 2*n_aug_ + 1; i++ ) {
    z_pred = z_pred + weights_(i)*Zsig_transformed.col(i);
  }

  MatrixXd S = MatrixXd(n_radar, n_radar);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig_transformed.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise. Measurement noise is just additive.
   MatrixXd R = MatrixXd(n_radar, n_radar);
   R.fill(0.);
   R(0,0) = std_radr_*std_radr_;
   R(1,1) = std_radphi_*std_radphi_;
   R(2,2) = std_radrd_*std_radrd_;
   S = S + R;

   // ---- KF Update -----
   //cross correlation matrix
   MatrixXd Tc = MatrixXd(n_x_, n_radar);
   Tc.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd z_diff = Zsig_transformed.col(i) - z_pred;
      z_diff(1) = NormalizeAngle(z_diff(1));

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      x_diff(3) = NormalizeAngle(x_diff(3));

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }


   MatrixXd K = Tc * S.inverse();
   VectorXd z_diff = z - z_pred;

   // Normalize angles once more.
   NormalizeAngle(z_diff(1));

   x_ = x_ + K * z_diff;
   P_ = P_ - K*S*K.transpose();

   //NIS
   NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
   cout << "================= NIS_radar_: " << NIS_radar_ << "================= " << endl;
}

double UKF::NormalizeAngle(double phi) {
  // cout << "Normalize angle: " << phi << endl;
  long two_pi = 2 * M_PI;
  long to_add = (phi < 0 ? two_pi : -two_pi);
  while (fabs(phi) > M_PI) {
    phi += to_add;
  }
  // cout << "End Normalized angle: " << phi << endl;
  return phi;
}

MatrixXd UKF::CreateSignmaPoints(VectorXd state, MatrixXd P, int state_size, int num_sig) {
    MatrixXd Xsig = MatrixXd( state_size, num_sig );
    MatrixXd L = P.llt().matrixL();
    Xsig.col(0) = state;
    for (int i = 0; i < state_size; i++){
        Xsig.col( i + 1 ) = state + sqrt(lambda_ + state_size) * L.col(i);
        Xsig.col( i + 1 + state_size ) = state - sqrt(lambda_ + state_size) * L.col(i);
    }
    return Xsig;
}
