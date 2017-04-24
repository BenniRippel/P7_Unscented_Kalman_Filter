#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
// TODO: Interfaces der Funktionen stimmen nicht

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
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0, 0, 1, 0, 0 ;
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << std_laspx_*std_laspx_, 0, 0, std_laspy_*std_laspy_;

    //set state dimension
    int n_x_ = 5;
    //set augmented dimension
    int n_aug_ = 7;
    //define spreading parameter
    double lambda_ = 3 - n_aug_;
    //set vector for weights_
    VectorXd weights_ = VectorXd(2*n_aug_+1);
    //set weights_
    weights_.fill(0.5/(lambda_+n_aug_));
    weights_(0)=lambda_/(lambda_+n_aug_);
    // Predicted sigma points
    MatrixXd Xsig_pred_;




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

    // check for small measurements and adjust
    if (fabs(meas_package.raw_measurements_[0])<0.001 && fabs(meas_package.raw_measurements_[1])<0.001)
    {
        cout<<"Adjusting Measurement: Setting zeros to 0.1!"<<endl;
        meas_package.raw_measurements_[0]=0.1;
        meas_package.raw_measurements_[1]=0.1;
    }


    // INITIALISATION

    // PREDICTION
    // get delta T
    double delta_T = (meas_package.timestamp_ - time_us_) / 1000000.0;
    // only predict if delta_T is above threshold
    if (delta_T > 0.001)
    {
        UKF::Prediction(delta_T);
    }

    // UPDATE

    // set previous timestamp
    time_us_=meas_package.timestamp_;

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
    // create augmented sigma points
    UKF::AugmentedSigmaPoints(Xsig_pred_);
    // predict sigma points
    UKF::SigmaPointPrediction(Xsig_pred_, delta_t);
    // get mean and covariance
    UKF::PredictMeanAndCovariance(x_, P_);
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
    // update the state by using Kalman Filter equation
    Eigen::Vector2d z(meas_package.raw_measurements_[0], meas_package.raw_measurements_[1]);

    VectorXd z_pred = H_laser_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_laser_.transpose();
    MatrixXd PHt = P_ * Ht;
    MatrixXd S = H_laser_ * PHt + R_laser_;
    MatrixXd Si = S.inverse();
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_laser_) * P_;

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

    Eigen::Vector2d z(meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],
                      meas_package.raw_measurements_[2]);

    // Predict RADAR Measurements

    //

}


//---------------------------------------------------
// Functions aus der Vorlesung


void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    //create augmented mean state
    x_aug.head(n_x_) = x_;
    //create augmented covariance matrix
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_, n_x_) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;
    //create square root matrix
    MatrixXd Sqrt_P_aug = P_aug.llt().matrixL();
    //create augmented sigma points
    Xsig_aug.col(0) = x_aug;

    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * Sqrt_P_aug.col(i);
        Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * Sqrt_P_aug.col(i);
    }
    //write result
    *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, double delta_t) {
    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    for (int i=0; i<(2 * n_aug_ + 1); i++)
    {
        float x, y, v, yaw_a, yaw_r, n_a, n_y;
        x = Xsig_out(0, i);
        y = Xsig_out(1, i);
        v = Xsig_out(2, i);
        yaw_a = Xsig_out(3, i);
        yaw_r = Xsig_out(4, i);
        n_a = Xsig_out(5, i);
        n_y = Xsig_out(6, i);

        if (fabs(yaw_r)<0.001)
        {
            // formula for yaw_r==0
            Xsig_pred(0,i) = x + v*cos(yaw_a)*delta_t + 0.5*delta_t*delta_t*cos(yaw_a)*n_a;
            Xsig_pred(1,i) = y + v*sin(yaw_a)*delta_t + 0.5*delta_t*delta_t*sin(yaw_a)*n_a;
            Xsig_pred(3,i) = yaw_a + 0.5*delta_t*delta_t*n_y;
        }
        else
        {
            //standard formula
            Xsig_pred(0,i) = x + v/yaw_r*(sin(yaw_a+yaw_r*delta_t)-sin(yaw_a)) + 0.5*delta_t*delta_t*cos(yaw_a)*n_a;
            Xsig_pred(1,i) = y + v/yaw_r*(-1*cos(yaw_a+yaw_r*delta_t)+cos(yaw_a)) + 0.5*delta_t*delta_t*sin(yaw_a)*n_a;
            Xsig_pred(3,i) = yaw_a + yaw_r*delta_t +0.5*delta_t*delta_t*n_y;
        }
        // for both cases
        Xsig_pred(2,i) = v + delta_t*n_a;
        Xsig_pred(4,i) = yaw_r + delta_t*n_y;
    }
    //write result
    *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

    //create vector for predicted state
    VectorXd x = VectorXd(n_x_);
    //create covariance matrix for prediction
    MatrixXd P = MatrixXd(n_x_, n_x_);

    //set weights_
    weights_.fill(0.5/(lambda_+n_aug_));
    weights_(0)=lambda_/(lambda_+n_aug_);

    //predict state mean
    x.fill(0.0); // init with zeros
    for (int i=0; i<2*n_aug_+1; i++)
    {
        x = x+ weights_(i)*Xsig_pred_.col(i);
    }

    //predict state covariance matrix
    P.fill(0.0); // init with zeros
    for (int i=0; i<2*n_aug_+1; i++)
    {
        VectorXd vec = Xsig_pred_.col(i)-x;
        //angle normalization
        while (vec(3)> M_PI) vec(3)-=2.*M_PI;
        while (vec(3)<-M_PI) vec(3)+=2.*M_PI;

        P = P+ weights_(i)*vec*vec.transpose();
    }

    //write result
    *x_out = x;
    *P_out = P;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

    //transform sigma points into measurement space
    for (int i=0; i<2*n_aug_+1; i++)
    {
        double px, py, v, yaw;
        px = Xsig_pred_(0, i);
        py = Xsig_pred_(1, i);
        v= Xsig_pred_(2, i);
        yaw = Xsig_pred_(3, i);

        double radius = sqrt(px*px+py*py);
        Zsig(0, i) = radius;
        Zsig(1, i) = atan2(py, px);
        Zsig(2, i) = (px*v*cos(yaw)+py*v*sin(yaw))/radius;
    }
    //calculate mean predicted measurement
    z_pred.fill(0.0); //init with 0.0
    for (int i=0; i<2*n_aug_+1; i++)
    {
        z_pred = z_pred + weights_(i)*Zsig.col(i);
    }
    //calculate measurement covariance matrix S
    MatrixXd R_ = MatrixXd(n_z, n_z);
    R_(0,0) = std_radr_*std_radr_;
    R_(1,1) = std_radphi_*std_radphi_;
    R_(2,2) = std_radrd_*std_radrd_;

    S.fill(0.0);

    for (int i=0; i<2*n_aug_+1; i++)
    {
        VectorXd vec = Zsig.col(i)-z_pred;
        //angle normalization
        while (vec(1)> M_PI) vec(1)-=2.*M_PI;
        while (vec(1)<-M_PI) vec(1)+=2.*M_PI;
        S = S + weights_(i)*vec*vec.transpose();
    }
    S=S+R_;

    //write result
    *z_out = z_pred;
    *S_out = S;
}

void UKF::UpdateState() {
    int n_z=3;
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for(int i=0; i<2*n_aug_+1; i++)
    {
        VectorXd x_vec = Xsig_pred_.col(i)-x_;
        VectorXd z_vec = Zsig.col(i) - z_pred;

        //normalization
        while(x_vec(3)<M_PI){x_vec(3)+=2.0*M_PI;}
        while(x_vec(3)>M_PI){x_vec(3)-=2.0*M_PI;}
        while(z_vec(1)<M_PI){z_vec(1)+=2.0*M_PI;}
        while(z_vec(1)>M_PI){z_vec(1)-=2.0*M_PI;}

        Tc = Tc + weights_(i)*x_vec*z_vec.transpose();
    }
    //calculate Kalman gain K;
    MatrixXd K = MatrixXd(n_x_, n_z);
    K = Tc * S_.inverse();

    //update state mean and covariance matrix
    x_ = x_ + K*(z - z_pred);
    P_ = P_ - K*S_*K.transpose();
}

