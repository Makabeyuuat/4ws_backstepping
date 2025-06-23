#include "kinematics_solver.hpp"
#include "initial.hpp" 
#include "mathFunc.h"        // 数学関数のヘッダーファイル
#include <array>
#include <iostream>


using namespace std;

double KinematicsSolver::calc_alpha_1_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_alpha_1_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_alpha_2_1_()
{
double ret;
ret = -2*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*sr.Cs1 + sr.Cs*sr.d*Power(Tan(Thetap),2)*sr.Cs1 - (1 - sr.Cs*sr.d)*Power(Tan(Thetap),2)*sr.Cs1 + Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1) + (1 - sr.Cs*sr.d)*Tan(Thetap)*(-(Power(sr.Cs,2)*Power(Sec(Thetap),2)) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv) + Power(sr.Cs,2)*Power(Tan(Thetap),2) - Tan(Thetap)*sr.Cs1) - sr.d*Tan(Thetap)*sr.Cs2;

// for (size_t i = 0; i < x_old.size(); ++i) {
//   std::cout << "x_old["<< i << "] = " << x_old[i] << "\n";
// }
// std::cout << "Thetap=" << Thetap << ", thetaT=" << thetaT
//           << ", u1_act=" << u1_act << ", u2_act=" << u2_act << "\n";
return ret;
}

double KinematicsSolver::calc_alpha_2_2_()
{
    double ret;
    ret = (Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),3))/lv;
    return ret;
}

double KinematicsSolver::calc_aqd_1_()
{
double ret;
ret = Cos(Thetap + thetaT)*nu1 - Sin(Thetap + thetaT)*u1_act*(sr.Cs*calc_d_SX_d_t_1_1_() + calc_d_SX_d_t_3_1_());
return ret;
}

double KinematicsSolver::calc_aqd_2_()
{
double ret;
ret = nu1*Sin(Thetap + thetaT) + Cos(Thetap + thetaT)*u1_act*(sr.Cs*calc_d_SX_d_t_1_1_() + calc_d_SX_d_t_3_1_());
return ret;
}

double KinematicsSolver::calc_aqd_3_()
{
double ret;
ret = athetapd + asd*sr.Cs + sr.Cs1*Power(calc_d_SX_d_t_1_1_(),2);
return ret;
}

double KinematicsSolver::calc_aqd_4_()
{
double ret;
ret = nu2;
return ret;
}

double KinematicsSolver::calc_aqd_5_()
{
double ret;
ret = (Cos(q_map(3))*nu1 - Sin(q_map(3))*u1_act*u2_act)/wheelRadius;
return ret;
}

double KinematicsSolver::calc_aqd_6_()
{
double ret;
ret = nu1/wheelRadius;
return ret;
}

double KinematicsSolver::calc_Axi_1_1_()
{
double ret;
ret = Sin(q_map(3) + q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_1_2_()
{
double ret;
ret = -Cos(q_map(3) + q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_1_3_()
{
double ret;
ret = -(lv*Cos(q_map(3)));
return ret;
}

double KinematicsSolver::calc_Axi_1_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_1_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_1_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_2_1_()
{
double ret;
ret = Sin(q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_2_2_()
{
double ret;
ret = -Cos(q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_2_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_2_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_2_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_2_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_3_1_()
{
double ret;
ret = -Cos(q_map(3) + q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_3_2_()
{
double ret;
ret = -Sin(q_map(3) + q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_3_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_3_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_3_5_()
{
double ret;
ret = wheelRadius;
return ret;
}

double KinematicsSolver::calc_Axi_3_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_4_1_()
{
double ret;
ret = -Cos(q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_4_2_()
{
double ret;
ret = -Sin(q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_4_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_4_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_4_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_4_6_()
{
double ret;
ret = wheelRadius;
return ret;
}

double KinematicsSolver::calc_Cxi_1_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_1_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_1_3_()
{
double ret;
ret = -(lv*(m_b + 2*m_w)*Cos(q_map(2))*x_d[2]);
return ret;
}

double KinematicsSolver::calc_Cxi_1_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_1_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_1_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_2_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_2_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_2_3_()
{
double ret;
ret = -(lv*(m_b + 2*m_w)*Sin(q_map(2))*x_d[2]);
return ret;
}

double KinematicsSolver::calc_Cxi_2_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_2_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_2_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_3_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_3_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_3_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_3_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_3_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_3_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_4_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_4_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_4_3_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_Cxi_4_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_4_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_4_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_5_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_5_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_5_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_5_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_5_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_5_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_6_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_6_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_6_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_6_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_6_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_6_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_1_1_()
{
double ret;
ret = Cos(Thetap)/(1 - sr.Cs*sr.d);
return ret;
}

double KinematicsSolver::calc_SX_1_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_2_1_()
{
double ret;
ret = Sin(Thetap);
return ret;
}


double KinematicsSolver::calc_SX_2_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_3_1_()
{
double ret;
ret = -((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv;
return ret;
}

double KinematicsSolver::calc_SX_3_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_4_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_4_2_()
{
double ret;
ret = 1;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_1_1_()
{
double ret;
ret = -((Sin(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*u1_act)/(1 - sr.Cs*sr.d)) - (Cos(Thetap)*(-(sr.Cs*Sin(Thetap)*u1_act) - (Cos(Thetap)*sr.d*u1_act*sr.Cs1)/(1 - sr.Cs*sr.d)))/Power(1 - sr.Cs*sr.d,2);
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_1_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_2_1_()
{
double ret;
ret = Cos(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*u1_act;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_2_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_3_1_()
{
double ret;
ret = (sr.Cs*Sin(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*u1_act)/(1 - sr.Cs*sr.d) + (Power(Sec(x_old[4]),2)*u2_act)/lv - (Power(Cos(Thetap),2)*u1_act*sr.Cs1)/Power(1 - sr.Cs*sr.d,2) + (sr.Cs*Cos(Thetap)*(-(sr.Cs*Sin(Thetap)*u1_act) - (Cos(Thetap)*sr.d*u1_act*sr.Cs1)/(1 - sr.Cs*sr.d)))/Power(1 - sr.Cs*sr.d,2);
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_3_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_4_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_4_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Kxi_1_()
{
double ret;
ret = GRAV*(m_b + 2*m_w)*Sin(rho);
return ret;
}

double KinematicsSolver::calc_Kxi_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Kxi_3_()
{
double ret;
ret = -(GRAV*lv*(m_b + 2*m_w)*Sin(rho)*Sin(q_map(2)))/2.;
return ret;
}

double KinematicsSolver::calc_Kxi_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Kxi_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Kxi_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_1_1_()
{
double ret;
ret = m_b + 2*m_w;
return ret;
}

double KinematicsSolver::calc_Mxi_1_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_1_3_()
{
double ret;
ret = -(lv*(m_b + 2*m_w)*Sin(q_map(2)))/2.;
return ret;
}

double KinematicsSolver::calc_Mxi_1_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_1_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_1_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_2_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_2_2_()
{
double ret;
ret = m_b + 2*m_w;
return ret;
}

double KinematicsSolver::calc_Mxi_2_3_()
{
double ret;
ret = (lv*(m_b + 2*m_w)*Cos(q_map(2)))/2.;
return ret;
}

double KinematicsSolver::calc_Mxi_2_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_2_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_2_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_3_1_()
{
double ret;
ret = -(lv*(m_b + 2*m_w)*Sin(q_map(2)))/2.;
return ret;
}

double KinematicsSolver::calc_Mxi_3_2_()
{
double ret;
ret = (lv*(m_b + 2*m_w)*Cos(q_map(2)))/2.;
return ret;
}

double KinematicsSolver::calc_Mxi_3_3_()
{
double ret;
ret = I_theta + (Power(lv,2)*(m_b + 4*m_w))/4.;
return ret;
}

double KinematicsSolver::calc_Mxi_3_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_3_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_3_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_4_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_4_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_4_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_4_4_()
{
double ret;
ret = I_phi;
return ret;
}

double KinematicsSolver::calc_Mxi_4_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_4_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_5_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_5_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_5_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_5_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_5_5_()
{
double ret;
ret = I_psif;
return ret;
}

double KinematicsSolver::calc_Mxi_5_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_6_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_6_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_6_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_6_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_6_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_6_6_()
{
double ret;
ret = I_psir;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_1_()
{
double ret;
ret = 2*Power(sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Power(sr.Cs1,2) + 2*sr.d*Power(Tan(Thetap),2)*Power(sr.Cs1,2) - 4*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*sr.Cs1*(-((sr.Cs*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) - sr.d*Sec(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*sr.Cs1*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1) + (1 - sr.Cs*sr.d)*Sec(Thetap)*(-((sr.Cs*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d))*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1) - sr.d*Tan(Thetap)*sr.Cs1*(-(Power(sr.Cs,2)*Power(Sec(Thetap),2)) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv) + Power(sr.Cs,2)*Power(Tan(Thetap),2) - Tan(Thetap)*sr.Cs1) - 2*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*sr.Cs2 + sr.Cs*sr.d*Power(Tan(Thetap),2)*sr.Cs2 - (1 - sr.Cs*sr.d)*Power(Tan(Thetap),2)*sr.Cs2 + Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*((-2*sr.Cs*Cos(Thetap)*Power(sr.d,2)*Power(sr.Cs1,2))/Power(1 - sr.Cs*sr.d,3) - (2*Cos(Thetap)*sr.d*Power(sr.Cs1,2))/Power(1 - sr.Cs*sr.d,2) - (sr.Cs*Cos(Thetap)*sr.d*sr.Cs2)/Power(1 - sr.Cs*sr.d,2) - (Cos(Thetap)*sr.Cs2)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*(sr.Cs*sr.d*Power(Sec(Thetap),2)*Tan(Thetap)*sr.Cs1 - (1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)*sr.Cs1 - 6*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap)*sr.Cs1 + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*Tan(Thetap)*(-((sr.Cs*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) - sr.d*Power(Sec(Thetap),2)*sr.Cs2) + (1 - sr.Cs*sr.d)*Tan(Thetap)*(-2*sr.Cs*Power(Sec(Thetap),2)*sr.Cs1 + 2*sr.Cs*sr.d*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*sr.Cs1 - 2*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*sr.Cs1 + 2*sr.Cs*Power(Tan(Thetap),2)*sr.Cs1 - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) - Tan(Thetap)*sr.Cs2) - sr.d*Tan(Thetap)*sr.Cs3;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_2_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*((2*Power(sr.Cs,3)*Power(Sec(Thetap),2))/(1 - sr.Cs*sr.d) + 2*Power(sr.Cs,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv))*Tan(Thetap) + (2*Power(sr.Cs,2)*sr.d*Power(Sec(Thetap),2)*sr.Cs1)/(1 - sr.Cs*sr.d) + 2*sr.Cs*sr.d*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*sr.Cs1 - 2*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*sr.Cs1 + 2*sr.Cs*Power(Tan(Thetap),2)*sr.Cs1 + Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*((-2*Power(sr.Cs,2)*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,3) - (2*sr.Cs*Cos(Thetap)*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*(-2*Power(sr.Cs,2)*Power(Sec(Thetap),2)*Tan(Thetap) - 6*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - Power(Sec(Thetap),2)*sr.Cs1) - (Power(sr.Cs,2)*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1))/(1 - sr.Cs*sr.d) - sr.Cs*Sec(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1) - sr.Cs*Tan(Thetap)*(-(Power(sr.Cs,2)*Power(Sec(Thetap),2)) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv) + Power(sr.Cs,2)*Power(Tan(Thetap),2) - Tan(Thetap)*sr.Cs1) - Tan(Thetap)*sr.Cs2;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_3_()
{
double ret;
ret = -2*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)*sr.Cs1 - 6*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap)*sr.Cs1 + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*Tan(Thetap)*(-((sr.Cs*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Tan(Thetap)*(-2*Power(sr.Cs,2)*Power(Sec(Thetap),2)*Tan(Thetap) - 6*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - Power(Sec(Thetap),2)*sr.Cs1) + sr.Cs*Tan(Thetap)*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1) + (1 - sr.Cs*sr.d)*Sec(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap)*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1) + Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*((sr.Cs*sr.d*Sin(Thetap)*sr.Cs1)/Power(1 - sr.Cs*sr.d,2) + (Sin(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*(-(Power(sr.Cs,2)*Power(Sec(Thetap),2)) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv) + Power(sr.Cs,2)*Power(Tan(Thetap),2) - Tan(Thetap)*sr.Cs1) + (1 - sr.Cs*sr.d)*Sec(Thetap)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),4)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),5)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv) + sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Power(Tan(Thetap),2) + 9*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Power(Tan(Thetap),2) - 2*sr.d*Power(Sec(Thetap),2)*Tan(Thetap)*sr.Cs1) - sr.d*Power(Sec(Thetap),2)*sr.Cs2;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_4_()
{
double ret;
ret = (-2*sr.Cs*Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),3)*Tan(Thetap))/lv + (3*Power(1 - sr.Cs*sr.d,3)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),4)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap))/lv - (2*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),3)*sr.Cs1)/lv + ((1 - sr.Cs*sr.d)*Power(Sec(x_old[4]),2)*Sec(Thetap)*(-(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1))/lv;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_1_()
{
double ret;
ret = (-2*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),3)*sr.Cs1)/lv;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_2_()
{
double ret;
ret = (-2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),3))/lv;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_3_()
{
double ret;
ret = (3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),3)*Tan(Thetap))/lv;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_4_()
{
double ret;
ret = (2*Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),3)*Tan(x_old[4]))/lv;
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_1_()
{
double ret;
ret = (Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2);
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_2_()
{
double ret;
ret = (sr.Cs*Cos(Thetap))/Power(1 - sr.Cs*sr.d,2);
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_3_()
{
double ret;
ret = -(Sin(Thetap)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_t_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_t_2_()
{
double ret;
ret = u2_act*(P21*calc_pd_Z2_pd_X_1_4_() + P22*calc_pd_Z2_pd_X_2_4_() + P23*calc_pd_Z2_pd_X_3_4_()) + u1_act*((P21*Tan(x_old[4])*calc_pd_Z2_pd_X_1_3_())/lv + (P22*Tan(x_old[4])*calc_pd_Z2_pd_X_2_3_())/lv + (P23*Tan(x_old[4])*calc_pd_Z2_pd_X_3_3_())/lv + (sr.Cs*Cos(Thetap)*(P21*calc_pd_Z2_pd_X_1_3_() + P22*calc_pd_Z2_pd_X_2_3_() + P23*calc_pd_Z2_pd_X_3_3_()))/(-1 + sr.Cs*sr.d) + P21*Sin(Thetap)*calc_pd_Z2_pd_X_1_2_() + P22*Sin(Thetap)*calc_pd_Z2_pd_X_2_2_() + P23*Sin(Thetap)*calc_pd_Z2_pd_X_3_2_() + (P21*Cos(Thetap)*calc_pd_Z2_pd_X_1_1_())/(1 - sr.Cs*sr.d) + (P22*Cos(Thetap)*calc_pd_Z2_pd_X_2_1_())/(1 - sr.Cs*sr.d) + (P23*Cos(Thetap)*calc_pd_Z2_pd_X_3_1_())/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_1_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_1_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_1_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_1_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_1_()
{
double ret;
ret = P21*calc_pd_Z2_pd_X_1_1_() + P22*calc_pd_Z2_pd_X_2_1_() + P23*calc_pd_Z2_pd_X_3_1_();
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_2_()
{
double ret;
ret = P21*calc_pd_Z2_pd_X_1_2_() + P22*calc_pd_Z2_pd_X_2_2_() + P23*calc_pd_Z2_pd_X_3_2_();
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_3_()
{
double ret;
ret = P21*calc_pd_Z2_pd_X_1_3_() + P22*calc_pd_Z2_pd_X_2_3_() + P23*calc_pd_Z2_pd_X_3_3_();
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_4_()
{
double ret;
ret = P21*calc_pd_Z2_pd_X_1_4_() + P22*calc_pd_Z2_pd_X_2_4_() + P23*calc_pd_Z2_pd_X_3_4_();
return ret;
}



//Z21,22,23とその偏微分の計算
double KinematicsSolver::calc_Z_2_1_()
{
double ret;
ret = Power(1 - sr.Cs*sr.d,2)*Power(Sec(thetaT),3)*(-((sr.Cs*Cos(thetaT))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv) - sr.Cs*(1 - sr.Cs*sr.d)*Power(Tan(thetaT),2) - sr.d*Tan(thetaT)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_Z_2_2_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Tan(thetaT);
return ret;
}

double KinematicsSolver::calc_Z_2_3_()
{
double ret;
ret = sr.d;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_1_()
{
double ret;
ret = -2*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*sr.Cs1 + sr.Cs*sr.d*Power(Tan(Thetap),2)*sr.Cs1 - (1 - sr.Cs*sr.d)*Power(Tan(Thetap),2)*sr.Cs1 + Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) - sr.d*Tan(Thetap)*sr.Cs2;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_2_()
{
double ret;
ret = -(Power(sr.Cs,2)*Power(Sec(Thetap),2)) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv) + Power(sr.Cs,2)*Power(Tan(Thetap),2) - Tan(Thetap)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_3_()
{
double ret;
ret = -(sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(Thetap),2)*Tan(Thetap)) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(Thetap),3)*(-((sr.Cs*Cos(Thetap))/(1 - sr.Cs*sr.d)) + Tan(x_old[4])/lv)*Tan(Thetap) - sr.d*Power(Sec(Thetap),2)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_4_()
{
double ret;
ret = (Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4]),2)*Power(Sec(Thetap),3))/lv;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_1_()
{
double ret;
ret = -(sr.d*Tan(Thetap)*sr.Cs1);
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_2_()
{
double ret;
ret = -(sr.Cs*Tan(Thetap));
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_3_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Power(Sec(Thetap),2);
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_3_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_3_2_()
{
double ret;
ret = 1;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_3_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_3_4_()
{
double ret;
ret = 0;
return ret;
}

//目標速度の時間微分
double KinematicsSolver::calc_pd_ud_pd_t_1_()
{
double ret;
ret = -((w1*(calc_d_SX_d_t_4_1_()*calc_pd_G11_pd_X_4_() + calc_d_SX_d_t_3_1_()*calc_pd_G11_pd_X_3_() + calc_d_SX_d_t_2_1_()*calc_pd_G11_pd_X_2_() + calc_d_SX_d_t_1_1_() * calc_pd_G11_pd_X_1_()))/Power(calc_SX_1_1_(),2));
return ret;
}

double KinematicsSolver::calc_pd_ud_pd_t_2_()
{
double ret;
ret = ((w1*calc_alpha_2_1_() - w2)*(calc_d_SX_d_t_4_1_()*calc_pd_alpha2_pd_X_2_4_() + calc_d_SX_d_t_3_1_()*calc_pd_alpha2_pd_X_2_3_() + calc_d_SX_d_t_2_1_()*calc_pd_alpha2_pd_X_2_2_() + calc_d_SX_d_t_1_1_()*calc_pd_alpha2_pd_X_2_1_()) + calc_alpha_2_2_()*(calc_d_SX_d_t_4_1_()*(-(w1*calc_pd_alpha2_pd_X_1_4_()) 
    + calc_pd_W_pd_X_2_4_()) + calc_d_SX_d_t_3_1_()*(-(w1*calc_pd_alpha2_pd_X_1_3_()) + calc_pd_W_pd_X_2_3_()) - w1*calc_d_SX_d_t_2_1_()*calc_pd_alpha2_pd_X_1_2_() + calc_d_SX_d_t_2_1_()*calc_pd_W_pd_X_2_2_() - w1*calc_d_SX_d_t_1_1_()*calc_pd_alpha2_pd_X_1_1_() + calc_d_SX_d_t_1_1_()*calc_pd_W_pd_X_2_1_()))/Power(calc_alpha_2_2_(),2);
return ret;
}



