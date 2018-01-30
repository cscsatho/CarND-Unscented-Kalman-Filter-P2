#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const int UKF::n_z_las_ = 2;
const int UKF::n_z_rad_ = 3;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF(bool verbose, bool use_laser, bool use_radar)
 : verbose_(verbose),
   is_initialized_(false), // initially set to false, set to true in first call of ProcessMeasurement
   use_laser_(use_laser),  // if this is false, laser measurements will be ignored (except during init)
   use_radar_(use_radar),  // if this is false, radar measurements will be ignored (except during init)
   std_a_(1.8), // 30      // Process noise standard deviation longitudinal acceleration in m/s^2
   std_yawdd_(0.5), // 30  // Process noise standard deviation yaw acceleration in rad/s^2
   //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
   std_laspx_(0.15),       // Laser measurement noise standard deviation position1 in m
   std_laspy_(0.15),       // Laser measurement noise standard deviation position2 in m
   std_radr_(0.3),         // Radar measurement noise standard deviation radius in m
   std_radphi_(0.03),      // Radar measurement noise standard deviation angle in rad
   std_radrd_(0.3),        // Radar measurement noise standard deviation radius change in m/s
   //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
   n_x_(5),                // State dimension
   n_aug_(7),              // Augmented state dimension
   lambda_(3 - n_aug_),    // Sigma point spreading parameter
   x_(n_x_),               // initial state vector
   P_(n_x_, n_x_),         // initial covariance matrix
   Xsig_pred_(n_x_, 2 * n_aug_ + 1), // predicted sigma points matrix
   weights_(2 * n_aug_ + 1), // Weights of sigma points
   time_us_(0),            // time when the state is true, in us
   nis_las_(NAN), nis_rad_(NAN)
{
  // TODO: Hint: one or more values initialized above might be wildly off...

  // set weights
  weights_[0] = lambda_ / (lambda_ + n_aug_);
  weights_[1] = 1 / (2 * (lambda_ + n_aug_));
  for (int i = 2; i <= 2 * n_aug_; ++i) weights_[i] = weights_[1];
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void
UKF::processMeasurement(const MeasurementPackage& meas_package)
{
  if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) return;
  if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) return;

  const double dt = (meas_package.timestamp_ - time_us_) / 1.0E6; // dt - expressed in secs

  if (!is_initialized_ || fabs(dt) > 5.0)
  {
    if (verbose_) cout << "(UKF Init)" << endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Convert radar from polar to cartesian coordinates and initialize state
      // [0]=rho, [1]=theta, [2]=rho_dot
      x_ << cos(meas_package.raw_measurements_[1]) * meas_package.raw_measurements_[0],
            sin(meas_package.raw_measurements_[1]) * meas_package.raw_measurements_[0],
            0.0, 0.0, 0.0; // Note: ignoring initial velocity
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      // Initialize state
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0.0, 0.0, 0.0;
    }
    else
      assert(!"sensor_type_ not suppoerted.");

    P_ = MatrixXd::Identity(n_x_, n_x_);
    P_(2, 2) = 20;
    P_(3, 3) = 20;

    time_us_ = meas_package.timestamp_;

    // Done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  prediction(dt, meas_package.sensor_type_);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    nis_rad_ = updateLR<n_z_rad_, 'R'>(meas_package);
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    nis_las_ = updateLR<n_z_las_, 'L'>(meas_package);

  time_us_ = meas_package.timestamp_;

  if (verbose_)
    cout << "ts: " << time_us_ << ", dt: " << dt << endl
         << "x_:" << endl 
         << x_ << endl
         << "P_:" << endl
         << P_ << endl
         << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void
UKF::prediction(const double& delta_t, const MeasurementPackage::SensorType type)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // create augmented sigma point matrix
  { // scoping
    MatrixXd Xsig(n_aug_, 2 * n_aug_ + 1);
    setAugSigmaPoints(Xsig);
    predictSigmaPoints(Xsig, delta_t);
  }

  predictMeanAndCovariance(x_, P_);

  if (type == MeasurementPackage::LASER)
    predictLaserMeasurement(); 
  else if (type == MeasurementPackage::RADAR)
    predictRadarMeasurement();
}

/**
 * Updates the state and the state covariance matrix using a radar/lidar measurement.
 * @param {MeasurementPackage} meas_package
 */
template<const int n_z, const char lr>
double
UKF::updateLR(const MeasurementPackage& meas_package)
{
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculating cross correlation matrix
  Tc.fill(0.0);
  VectorXd z_pred_diff = VectorXd(n_z);
  VectorXd x_pred_diff = VectorXd(n_x_);
  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
  {
    z_pred_diff = Zsig_pred_.col(col) - z_pred_;
    if (lr == 'R')
    {
      while (z_pred_diff[1] >  M_PI) z_pred_diff[1] -= 2.*M_PI;
      while (z_pred_diff[1] < -M_PI) z_pred_diff[1] += 2.*M_PI;
    }

    x_pred_diff = Xsig_pred_.col(col) - x_;
    while (x_pred_diff[3] >  M_PI) x_pred_diff[3] -= 2.*M_PI;
    while (x_pred_diff[3] < -M_PI) x_pred_diff[3] += 2.*M_PI;

    Tc += x_pred_diff * z_pred_diff.transpose() * weights_(col);
  }

  // Calculating Kalman gain K
  MatrixXd K = Tc * S_pred_.inverse();

  z_pred_diff = meas_package.raw_measurements_ - z_pred_;
  if (lr == 'R')
  {
    while (z_pred_diff[1] >  M_PI) z_pred_diff[1] -= 2.*M_PI;
    while (z_pred_diff[1] < -M_PI) z_pred_diff[1] += 2.*M_PI;
  }

  // Updating state mean and covariance matrix
  x_ = x_ + K * (meas_package.raw_measurements_ - z_pred_);
  P_ = P_ - K * S_pred_ * K.transpose();

  double nis = z_pred_diff.transpose() * S_pred_.inverse() * z_pred_diff;

  if (verbose_)
    cout << "NIS{" << lr << "}: "
         << nis         << endl
         << "x_upd{" << lr << "}:" << endl
         << x_          << endl
         << "P_upd{" << lr << "}:" << endl
         << P_ << endl << endl;

  return nis;
}

void
UKF::setAugSigmaPoints(MatrixXd& Xsig_aug) const
{
  assert (Xsig_aug.rows() == n_aug_ && Xsig_aug.cols() == (2 * n_aug_ + 1));

  // caching this value
  static const double SqrtLNX = sqrt(lambda_ + n_aug_);

  // create augmented mean vector
  VectorXd x_aug(n_aug_);
  // set first column of sigma point matrix
  x_aug.head(n_x_) = x_;
  x_aug[n_x_]     = 0;
  x_aug[n_x_ + 1] = 0;

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd::Constant(n_aug_, n_aug_, 0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_    , n_x_)     = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // calculate square root of P
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  // assigning first column
  Xsig_aug.col(0) = x_aug;

  // set remaining sigma points
  for (int col = 0; col < n_aug_; ++col)
  {
    Xsig_aug.col(col + 1)          = x_aug + SqrtLNX * P_aug_sqrt.col(col);
    Xsig_aug.col(col + 1 + n_aug_) = x_aug - SqrtLNX * P_aug_sqrt.col(col);
  }
//  cout << "Xsig_aug:" << endl << Xsig_aug << endl;
}

void
UKF::predictSigmaPoints(const MatrixXd& Xsig_aug, const double& dt)
{
  assert (Xsig_pred_.rows() == n_x_ && Xsig_pred_.cols() == (2 * n_aug_ + 1));

  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
  {
    if (fabs(Xsig_aug(4, col)) > Tools::EPS) // psi_dot is not zero
    {
      Xsig_pred_(0, col) = Xsig_aug(0, col)
                         + Xsig_aug(2, col) / Xsig_aug(4, col) * (sin(Xsig_aug(3, col) + Xsig_aug(4, col) * dt) - sin(Xsig_aug(3, col)))
                         + 0.5 * dt*dt * cos(Xsig_aug(3, col)) * Xsig_aug(5, col);
      Xsig_pred_(1, col) = Xsig_aug(1, col)
                         + Xsig_aug(2, col) / Xsig_aug(4, col) * (-cos(Xsig_aug(3, col) + Xsig_aug(4, col) * dt) + cos(Xsig_aug(3, col)))
                         + 0.5 * dt*dt * sin(Xsig_aug(3, col)) * Xsig_aug(5, col);
      Xsig_pred_(2, col) = Xsig_aug(2, col)
                         + 0
                         + dt * Xsig_aug(5, col);
      Xsig_pred_(3, col) = Xsig_aug(3, col)
                         + Xsig_aug(4, col) * dt
                         + 0.5 * dt*dt * Xsig_aug(6, col);
      Xsig_pred_(4, col) = Xsig_aug(4, col)
                         + 0
                         + dt * Xsig_aug(6, col);
    }
    else
    {
      Xsig_pred_(0, col) = Xsig_aug(0, col) 
                         + Xsig_aug(2, col) * cos(Xsig_aug(3, col)) * dt
                         + 0.5 * dt*dt * cos(Xsig_aug(3, col)) * Xsig_aug(5, col);
      Xsig_pred_(1, col) = Xsig_aug(1, col)
                         + Xsig_aug(2, col) * sin(Xsig_aug(3, col)) * dt
                         + 0.5 * dt*dt * sin(Xsig_aug(3, col)) * Xsig_aug(5, col);
      Xsig_pred_(2, col) = Xsig_aug(2, col)
                         + 0
                         + dt * Xsig_aug(5, col);
      Xsig_pred_(3, col) = Xsig_aug(3, col)
                         + 0
                         + 0.5 * dt*dt * Xsig_aug(6, col);
      Xsig_pred_(4, col) = Xsig_aug(4, col)
                         + 0
                         + dt * Xsig_aug(6, col);
    }
  }
}

void
UKF::predictMeanAndCovariance(VectorXd& x_pred, MatrixXd& P_pred) const
{
  //create vector for predicted state
  x_pred = VectorXd::Constant(n_x_, 0.0);
  //create covariance matrix for prediction
  P_pred = MatrixXd::Constant(n_x_, n_x_, 0.0);

  // predict state mean
  for (int col = 0; col <= 2 * n_aug_; ++col)
    x_pred += Xsig_pred_.col(col) * weights_[col];
   
  VectorXd x_pred_diff = VectorXd::Constant(n_x_, 0.0);
  // predict state covariance matrix
  for (int col = 0; col <= 2 * n_aug_; ++col)
  {
    x_pred_diff = Xsig_pred_.col(col) - x_pred;

    while (x_pred_diff[3] >  M_PI) x_pred_diff[3] -= 2.*M_PI;
    while (x_pred_diff[3] < -M_PI) x_pred_diff[3] += 2.*M_PI;

    P_pred += x_pred_diff * weights_[col] * x_pred_diff.transpose();
  }

  if (verbose_)
    cout << "x_pred:" << endl
         << x_pred    << endl
         << "P_pred:" << endl
         << P_pred    << endl << endl;
}

void
UKF::predictRadarMeasurement()
{
  // create matrix for sigma points in measurement space
  Zsig_pred_ = MatrixXd(n_z_rad_, 2 * n_aug_ + 1);

  // transforming sigma points into measurement space
  VectorXd vecX(n_x_);
  VectorXd vecZ(n_z_rad_);
  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
  {
    vecX = Xsig_pred_.col(col);
    vecZ[0] = sqrt(vecX[0] * vecX[0] + vecX[1] * vecX[1]);
    vecZ[1] = std::atan2(vecX[1], vecX[0]);
    if (vecZ[0] > Tools::EPS)
      vecZ[2] = (vecX[0] * cos(vecX[3]) * vecX[2] + vecX[1] * sin(vecX[3]) * vecX[2]) / vecZ[0];
    else
      vecZ[2] = 0;
   
    Zsig_pred_.col(col) = vecZ;
  }

  // calculating mean predicted measurement
  z_pred_ = VectorXd::Constant(n_z_rad_, 0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
    z_pred_ += weights_[col] * Zsig_pred_.col(col);
 
  // calculating measurement covariance matrix S
  S_pred_ = MatrixXd::Constant(n_z_rad_, n_z_rad_, 0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
  {
    // state difference
    vecX = Zsig_pred_.col(col) - z_pred_;
    // angle normalization
    while (vecX[1] >  2.*M_PI) vecX[1] -= 2.*M_PI;
    while (vecX[1] < -2.*M_PI) vecX[1] += 2.*M_PI;
   
    S_pred_ += weights_[col] * vecX * vecX.transpose();
  }
  
  // adding measurement noise
  Eigen::Matrix3d R;
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;
  S_pred_ += R;

  if (verbose_)
    cout << "z_pred{R}: " << endl << z_pred_ << endl
         << "S_pred{R}: " << endl << S_pred_ << endl
         << endl;
}

void
UKF::predictLaserMeasurement()
{
  // create matrix for sigma points in measurement space
  Zsig_pred_ = MatrixXd(n_z_las_, 2 * n_aug_ + 1);

  // transforming sigma points into measurement space
  VectorXd vecX(n_x_);
//  VectorXd vecZ(n_z_las_);
  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
  {
    //Zsig_pred_(0, col) = Xsig_pred_(0, col);
    //Zsig_pred_(1, col) = Xsig_pred_(1, col);
    Zsig_pred_.col(col) = Xsig_pred_.col(col).head(n_z_las_);
  }

  // calculating mean predicted measurement
  z_pred_ = VectorXd::Constant(n_z_las_, 0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
    z_pred_ += weights_[col] * Zsig_pred_.col(col);
 
  // calculating measurement covariance matrix S
  S_pred_ = MatrixXd::Constant(n_z_las_, n_z_las_, 0.0);
  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
  {
    // state difference
    vecX = Zsig_pred_.col(col) - z_pred_;
    S_pred_ += weights_[col] * vecX * vecX.transpose();
  }
  
  // adding measurement noise
  Eigen::Matrix2d R;
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;
  S_pred_ += R;

  if (verbose_)
    cout << "z_pred{L}: " << endl << z_pred_ << endl
         << "S_pred{L}: " << endl << S_pred_ << endl
         << endl;
}
