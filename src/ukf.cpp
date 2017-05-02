#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */

UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = false;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    //set state dimension
    n_x_ = 5;
    //set augmented dimension
    n_aug_ = 7;
    //define spreading parameter
    lambda_ = 3 - n_aug_;

    // initial state vector
    x_ = VectorXd(n_x_);

    // initial covariance matrix
    P_ = MatrixXd::Identity(n_x_, n_x_);

    ///* predicted sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
    Xsig_pred_.fill(0);

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

    H_laser_ = MatrixXd(2, 5);
    H_laser_ << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 ;

    R_laser_ = MatrixXd(2, 2);
    R_laser_ << std_laspx_*std_laspx_, 0, 0, std_laspy_*std_laspy_;

    //set vector for weights_
    weights_ = VectorXd(2*n_aug_+1);
    //init weights_
    weights_.fill(0.5/(lambda_+n_aug_));
    weights_(0)=lambda_/(lambda_+n_aug_);

    // initialization
    is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  * Process incoming measurements:
   *    - initialize on first measurement
   *    - Prediction
   *    - Update dependend on measurement type
  */

    // INITIALISATION WITH FIRST MEASUREMENT
    if (!is_initialized_)
    {

        // check for small measurements and adjust
        if (fabs(meas_package.raw_measurements_[0])<0.001 && fabs(meas_package.raw_measurements_[1])<0.001)
        {
            std::cout<<"Adjusting Measurement: Setting zeros to 0.1!"<<std::endl;
            meas_package.raw_measurements_[0]=0.1;
            meas_package.raw_measurements_[1]=0.1;
        }


        if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            x_ << cos(meas_package.raw_measurements_[1])*meas_package.raw_measurements_[0],
                    sin(meas_package.raw_measurements_[1])*meas_package.raw_measurements_[0], 0, 0, 0;
            time_us_=meas_package.timestamp_;
            // done initializing, no need to predict or update
            is_initialized_ = true;


        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
            /**
            Initialize state.
            */
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
            time_us_=meas_package.timestamp_;
            // done initializing, no need to predict or update
            is_initialized_ = true;

        }

        std::cout<<"Init!"<<std::endl;
        std::cout<<"state: "<<x_<<std::endl;
        return;
    }
    // PREDICTION
    // get delta T
    double delta_T = (meas_package.timestamp_ - time_us_) / 1000000.0;
    // only predict if delta_T is above threshold
    if (delta_T > 0.001)
    {
        std::cout<<"Prediction Step"<<std::endl;
        UKF::Prediction(delta_T);
    }

    // UPDATE USING SENSOR MEASUREMENTS
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
        // for radar measurement
        std::cout<<"Radar Update"<<std::endl;
        UKF::UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {
        // for laser measurement
        std::cout<<"Laser Update"<<std::endl;
        UKF::UpdateLidar(meas_package);
    }

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
  * Predict sigma points, the state, and the state covariance matrix.
  */
    // create augmented sigma points
    MatrixXd Xsig_aug = UKF::AugmentedSigmaPoints();
    // predict sigma points
    UKF::SigmaPointPrediction(Xsig_aug, delta_t);
    // get mean and covariance
    UKF::PredictMeanAndCovariance();

    std::cout<<"Predicted state: "<<x_<<std::endl;
    //std::cout<<"Covariance: "<<P_<<std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    *Updates the state and the state covariance matrix using a lidar measurement
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

    std::cout<<"Update Lidar state: "<<x_<<std::endl;
    //std::cout<<"Covariance: "<<P_<<std::endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  *Updates the state and the state covariance matrix using a radar measurement
  */

    Eigen::Vector3d z(meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],
                      meas_package.raw_measurements_[2]);

    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

    // Predict RADAR Measurements
    UKF::PredictRadarMeasurement(z_pred, S, Zsig);
    // Update State
    UKF::UpdateState(z, z_pred, Zsig, S);
    std::cout<<"Update Radar state: "<<x_<<std::endl;

}

MatrixXd UKF::AugmentedSigmaPoints() {
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
    //return result
    return Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, const double &delta_t) {
    for (int i=0; i<(2 * n_aug_ + 1); i++)
    {
        float x, y, v, yaw_a, yaw_r, n_a, n_y;
        x = Xsig_aug(0, i);
        y = Xsig_aug(1, i);
        v = Xsig_aug(2, i);
        yaw_a = Xsig_aug(3, i);
        yaw_r = Xsig_aug(4, i);
        n_a = Xsig_aug(5, i);
        n_y = Xsig_aug(6, i);

        if (fabs(yaw_r)<0.001)
        {
            // formula for yaw_r==0
            Xsig_pred_(0,i) = x + v*cos(yaw_a)*delta_t + 0.5*delta_t*delta_t*cos(yaw_a)*n_a;
            Xsig_pred_(1,i) = y + v*sin(yaw_a)*delta_t + 0.5*delta_t*delta_t*sin(yaw_a)*n_a;
            Xsig_pred_(3,i) = yaw_a + 0.5*delta_t*delta_t*n_y;
        }
        else
        {
            //standard formula
            Xsig_pred_(0,i) = x + v/yaw_r*(sin(yaw_a+yaw_r*delta_t)-sin(yaw_a)) + 0.5*delta_t*delta_t*cos(yaw_a)*n_a;
            Xsig_pred_(1,i) = y + v/yaw_r*(-1*cos(yaw_a+yaw_r*delta_t)+cos(yaw_a)) + 0.5*delta_t*delta_t*sin(yaw_a)*n_a;
            Xsig_pred_(3,i) = yaw_a + yaw_r*delta_t +0.5*delta_t*delta_t*n_y;
        }
        // for both cases
        Xsig_pred_(2,i) = v + delta_t*n_a;
        Xsig_pred_(4,i) = yaw_r + delta_t*n_y;
    }
}

void UKF::PredictMeanAndCovariance() {
    //predict state mean
    x_.fill(0.0); // init with zeros
    for (int i=0; i<2*n_aug_+1; i++)
    {
        x_ = x_+ weights_(i)*Xsig_pred_.col(i);
    }

    //predict state covariance matrix
    P_.fill(0.0); // init with zeros
    for (int i=0; i<2*n_aug_+1; i++)
    {
        VectorXd vec = Xsig_pred_.col(i)-x_;
        //angle normalization
        while (vec(3)> M_PI) vec(3)-=2.*M_PI;
        while (vec(3)<-M_PI) vec(3)+=2.*M_PI;

        P_ = P_+ weights_(i)*vec*vec.transpose();
    }
}

void UKF::PredictRadarMeasurement(VectorXd &z_pred, MatrixXd &S, MatrixXd &Zsig) {
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
    MatrixXd R_ = MatrixXd(3,3);
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
}

void UKF::UpdateState(const VectorXd &z, const VectorXd &z_pred, const MatrixXd &Zsig, const MatrixXd &S) {
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
    K = Tc * S.inverse();

    //update state mean and covariance matrix
    x_ = x_ + K*(z - z_pred);
    P_ = P_ - (K*S*K.transpose());
}

