#include "differential_equations_dynamics.hpp"
#include "initial.hpp"    // u1…u12, l1…lv, ai の extern 宣言
#include <cmath>

using namespace std;

double f0 (const std::vector<double>& x) { 
    return 1.0; 
}

double f1 (const std::vector<double>& x) { 
    return u1 * std::cos(x[3] + x[4]); 
}

double f2 (const std::vector<double>& x) { 
    return u1 * std::sin(x[3] + x[4]); 
}


//theta1
double f3 (const std::vector<double>&x) { 
    return u1 * std::tan(x[4])/lv; 
}

//phi1
double f4 (const std::vector<double>&x) {  
    return u2; 
}

//phi2
double f5 (const std::vector<double>&x) {  
    return u3; 
}


double f6 (const std::vector<double>& x) { 
    return 1.0; 
}

double f7(const std::vector<double>& x) {
    double u1_dot =  0.0;
    double theta   = x[3];
    double theta_d = x_d[3];
    double phi1   = x[4];
    double phi1_d = x_d[4];
    return u1_dot * std::cos(theta + phi1)
         - u1    * std::sin(theta + phi1) * theta_d * phi1_d;
}

double f8(const std::vector<double>& x) {
    double u1_dot =  0.0;
    double theta   = x[3];
    double theta_d = x_d[3];
    double phi1   = x[4];
    double phi1_d = x_d[4];
    return u1_dot * std::sin(theta + phi1)
         + u1    * std::cos(theta + phi1) * theta_d * phi1_d;
}

double f9(const std::vector<double>& x) {
    double u1_dot = 0.0;
    double phi1     = x[4];
    double phi1_d   = x_d[4];
    double phi2     = x[5];
    double phi2_d   = x_d[5];
    return u1_dot * (-sin(phi2 - phi1)/(lv*cos(phi1)))
         + (u1/lv) * ((cos(phi2 - phi1)*(phi2_d - phi1_d))/cos(phi1) + (sin(phi2 - phi1)*sin(phi1) * phi1_d)/pow(cos(phi1),2));
}

double f10(const std::vector<double>& x) {
    // phi2_d の微分 = u2 の時間微分
    double u2_dot =  0.0;
    return u2_dot;
}

double f11(const std::vector<double>& x) {
    // phi2_d の微分 = u2 の時間微分
    double u3_dot =  0.0;
    return u3_dot;
}

const std::array<FunctionPtr, 6> fAll = {{
  f0, f1, f2, f3, f4, f5
}};
const std::array<FunctionPtr, 6> fdAll = {{
  f6, f7, f8, f9, f10, f11
}};
