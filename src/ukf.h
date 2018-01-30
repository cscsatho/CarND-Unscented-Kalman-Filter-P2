#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF
{
public:
  UKF(bool verbose, bool use_laser, bool use_radar);
  virtual ~UKF() = default;
  /**
   * processMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void processMeasurement(const MeasurementPackage& meas_package);
  /**
   * prediction predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void prediction(const double& delta_t, const MeasurementPackage::SensorType type);
  /**
   * updates the state and the state covariance matrix using a laser/radar measurement
   * @param meas_package The measurement at k+1
   */
  template<const int n_z, const char lr>
  double updateLR(const MeasurementPackage& meas_package);

private:
  void setAugSigmaPoints(MatrixXd& Xsig_aug) const;
  void predictSigmaPoints(const MatrixXd& Xsig_aug, const double& dt);
  void predictMeanAndCovariance(VectorXd& x_pred, MatrixXd& P_pred) const;
  void predictRadarMeasurement();
  void predictLaserMeasurement();

public:
  bool verbose_;
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;
  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;
  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

private:
  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_;
  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_;
  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_;
  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_;
  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_;
  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_;
  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_ ;

private:
  ///* State dimension
  const int n_x_;
  ///* Augmented state dimension
  const int n_aug_;
  ///* Measurement dimension - laser
  static const int n_z_las_;
  ///* Measurement dimension - radar
  static const int n_z_rad_;
  ///* Sigma point spreading parameter
  const double lambda_;

public: // TODO - create accessors
  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  ///* state covariance matrix
  MatrixXd P_;
  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;
  ///* Weights of sigma points
  VectorXd weights_;
  ///* time when the state is true, in us
  long time_us_;

  double nis_las_;
  double nis_rad_;

private:
  VectorXd z_pred_;
  MatrixXd S_pred_;
  MatrixXd Zsig_pred_;
};

#endif /* UKF_H */
