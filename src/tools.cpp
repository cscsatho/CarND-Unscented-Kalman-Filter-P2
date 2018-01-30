#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

const float Tools::EPS = 1E-5;

VectorXd
Tools::CalculateRMSE(const vector<VectorXd>& estimations,
                     const vector<VectorXd>& ground_truth)
{
   VectorXd rmse(4);
   rmse << 0, 0, 0, 0;

   // check the validity of the following inputs:
   // the estimation vector size should not be zero
   if (estimations.empty())
   {
       std::cout << "estimations vector is empty" << std::endl;
       return rmse;
   }

   // the estimation vector size should equal ground truth vector size
   if (estimations.size() != ground_truth.size())
   {
       std::cout << "estimations vector and ground_truth vector sizes do not match" << std::endl;
       return rmse;
   }

   // accumulate squared residuals
   for (int i = 0; i < estimations.size(); ++i)
   {
       VectorXd diff = estimations[i] - ground_truth[i];
       diff = diff.array() * diff.array();
       rmse += diff;
   }

   // calculate the mean
   rmse /= estimations.size();

   // calculate the squared root
   rmse = rmse.array().sqrt();

   // return the result
   return rmse;
}
