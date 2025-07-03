#include "initial.hpp"
#include <vector>
#include <Eigen/Dense>


//double Thetap[11] = {0.0};

Eigen::Map<Eigen::Matrix<double,7,1>> q_map(q_twist);
Eigen::Map<Eigen::Matrix<double,7,1>> qdot_map(qdot_twist);


void initial(double &t, double &dt) {
   
    t = 66.8;
    dt = 0.01;
    w1 = a0;

    

}
