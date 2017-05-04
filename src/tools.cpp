#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
   Calculation of the RMSE.
   */
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    /*
    // check the validity of the inputs
    if ((estimations.size()==0) || (estimations.size()!=ground_truth.size())){
        std::cout<<"Error in CalculateRMSE(). Size of inputs are suspicious"<<std::endl;
    }
    */

    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        VectorXd res = estimations[i] - ground_truth[i];
        res = res.array()*res.array();
        rmse += res;
    }

    //calculate the mean
    rmse = rmse.array() * 1.0/estimations.size();
    //calculate the squared root
    rmse=rmse.array().sqrt();
    return rmse;
}
