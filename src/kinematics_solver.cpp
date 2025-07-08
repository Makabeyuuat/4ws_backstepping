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

double KinematicsSolver::calc_alpha_1_3_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_alpha_2_1_()
{
double ret;
ret = Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv) - sr.Cs*(1 - sr.Cs*sr.d)*Power(Tan(x_old[4] + Thetap),2) - sr.d*Tan(x_old[4] + Thetap)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_alpha_2_2_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2);
return ret;
}

double KinematicsSolver::calc_alpha_2_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_alpha_3_1_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*(-(Power(sr.Cs,2)/(1 - sr.Cs*sr.d)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv))*Tan(x_old[4] + Thetap) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) - sr.d*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*sr.Cs1 + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_alpha_3_2_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap);
return ret;
}

double KinematicsSolver::calc_alpha_3_3_()
{
double ret;
ret = (Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap))/lv;
return ret;
}

//Z21,22,23とその偏微分の計算
double KinematicsSolver::calc_Z_2_1_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Tan(x_old[4] + Thetap);
return ret;
}

double KinematicsSolver::calc_Z_2_2_()
{
double ret;
ret = sr.d;
return ret;
}

double KinematicsSolver::calc_Z_3_1_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv);
return ret;
}

double KinematicsSolver::calc_Z_3_2_()
{
double ret;
ret = Thetap;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_1_()
{
double ret;
ret = -(sr.d*Tan(x_old[4] + Thetap)*sr.Cs1);
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_2_()
{
double ret;
ret = -(sr.Cs*Tan(x_old[4] + Thetap));
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_3_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2);
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_4_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2);
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_2_()
{
double ret;
ret = 1;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_2_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_1_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv);
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_2_()
{
double ret;
ret = -(Power(sr.Cs,2)/(1 - sr.Cs*sr.d)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv);
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_3_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap);
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_4_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap);
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_5_()
{
double ret;
ret = (Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap))/lv;
return ret;
}


double KinematicsSolver::calc_pd_Z3_pd_X_2_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_2_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_2_3_()
{
double ret;
ret = 1;
return ret;
}


double KinematicsSolver::calc_pd_Z3_pd_X_2_4_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_pd_Z3_pd_X_2_5_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_aqd_1_()
{
double ret;
ret = Cos(q_map(3) + Thetap + thetaT)*nu1 - Sin(q_map(3) + Thetap + thetaT)*u1_act*(u2_act + sr.Cs*calc_SX_1_1_()*u1_act + calc_SX_3_1_()*u1_act);
return ret;
}

double KinematicsSolver::calc_aqd_2_()
{
double ret;
ret = nu1*Sin(q_map(3) + Thetap + thetaT) + Cos(q_map(3) + Thetap + thetaT)*u1_act*(u2_act + sr.Cs*calc_SX_1_1_()*u1_act + calc_SX_3_1_()*u1_act);
return ret;
}

double KinematicsSolver::calc_aqd_3_()
{
double ret;
ret = athetapd + asd*sr.Cs + sr.Cs1*Power(calc_SX_1_1_()*u1_act,2);
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
ret = nu1/wheelRadius;
return ret;
}

double KinematicsSolver::calc_aqd_6_()
{
double ret;
ret = nu3;
return ret;
}

double KinematicsSolver::calc_aqd_7_()
{
double ret;
ret = (Cos(q_map(5) - q_map(3))*nu1 + Sin(q_map(5) - q_map(3))*u1_act*(u2_act - u3_act))/wheelRadius;
return ret;
}

double KinematicsSolver::calc_Axi_1_1_()
{
double ret;
ret = Sin(q_map(5) + q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_1_2_()
{
double ret;
ret = -Cos(q_map(5) + q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_1_3_()
{
double ret;
ret = -(lv*Cos(q_map(5)));
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

double KinematicsSolver::calc_Axi_1_7_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_2_1_()
{
double ret;
ret = Sin(q_map(3) + q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_2_2_()
{
double ret;
ret = -Cos(q_map(3) + q_map(2));
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

double KinematicsSolver::calc_Axi_2_7_()
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

double KinematicsSolver::calc_Axi_3_7_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_4_1_()
{
double ret;
ret = -Cos(q_map(5) + q_map(2));
return ret;
}

double KinematicsSolver::calc_Axi_4_2_()
{
double ret;
ret = -Sin(q_map(5) + q_map(2));
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
ret = 0;
return ret;
}

double KinematicsSolver::calc_Axi_4_7_()
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
ret = -(lv*(m_b + 2*m_w)*Cos(q_map(2))*qdot_map(2));
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

double KinematicsSolver::calc_Cxi_1_7_()
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
ret = -(lv*(m_b + 2*m_w)*Sin(q_map(2))*qdot_map(2));
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

double KinematicsSolver::calc_Cxi_2_7_()
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

double KinematicsSolver::calc_Cxi_3_7_()
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

double KinematicsSolver::calc_Cxi_4_7_()
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

double KinematicsSolver::calc_Cxi_5_7_()
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

double KinematicsSolver::calc_Cxi_6_7_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_7_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_7_2_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_Cxi_7_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_7_4_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_7_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_7_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Cxi_7_7_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_1_1_()
{
double ret;
ret = Cos(x_old[4] + Thetap)/(1 - sr.Cs*sr.d);
return ret;
}

double KinematicsSolver::calc_SX_1_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_1_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_2_1_()
{
double ret;
ret = Sin(x_old[4] + Thetap);
return ret;
}


double KinematicsSolver::calc_SX_2_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_2_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_3_1_()
{
double ret;
ret = -((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv;
return ret;
}

double KinematicsSolver::calc_SX_3_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_3_3_()
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

double KinematicsSolver::calc_SX_4_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_5_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_5_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_SX_5_3_()
{
double ret;
ret =1;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_1_1_()
{
double ret;
ret = -((Sin(x_old[4] + Thetap)*((-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*u1_act + u2_act))/(1 - sr.Cs*sr.d)) - (Cos(x_old[4] + Thetap)*(-(sr.Cs*Sin(x_old[4] + Thetap)*u1_act) - (Cos(x_old[4] + Thetap)*sr.d*u1_act*sr.Cs1)/(1 - sr.Cs*sr.d)))/Power(1 - sr.Cs*sr.d,2);
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_1_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_1_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_2_1_()
{
double ret;
ret = Cos(x_old[4] + Thetap)*((-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*u1_act + u2_act);
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_2_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_2_3_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_d_SX_d_t_3_1_()
{
double ret;
ret = -((Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4])*u2_act)/lv) + (sr.Cs*Sin(x_old[4] + Thetap)*((-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*u1_act + u2_act))/(1 - sr.Cs*sr.d) - (Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4])*(-((-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*u1_act) + u2_act - u3_act))/lv - (Power(Cos(x_old[4] + Thetap),2)*u1_act*sr.Cs1)/Power(1 - sr.Cs*sr.d,2) + (sr.Cs*Cos(x_old[4] + Thetap)*(-(sr.Cs*Sin(x_old[4] + Thetap)*u1_act) - (Cos(x_old[4] + Thetap)*sr.d*u1_act*sr.Cs1)/(1 - sr.Cs*sr.d)))/Power(1 - sr.Cs*sr.d,2);
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_3_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_3_3_()
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

double KinematicsSolver::calc_d_SX_d_t_4_3_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_5_1_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_5_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_5_3_()
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

double KinematicsSolver::calc_Kxi_7_()
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

double KinematicsSolver::calc_Mxi_1_7_()
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

double KinematicsSolver::calc_Mxi_2_7_()
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

double KinematicsSolver::calc_Mxi_3_7_()
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
ret = I_phiR;
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

double KinematicsSolver::calc_Mxi_4_7_()
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
ret = I_varphiR;
return ret;
}

double KinematicsSolver::calc_Mxi_5_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_5_7_()
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
ret = I_phiF;
return ret;
}

double KinematicsSolver::calc_Mxi_6_7_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_Mxi_7_1_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_Mxi_7_2_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_7_3_()
{
double ret;
ret = 0;
return ret;
}       

double KinematicsSolver::calc_Mxi_7_4_()
{
double ret;
ret = 0;             
return ret;
}

double KinematicsSolver::calc_Mxi_7_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_7_6_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_Mxi_7_7_()
{
double ret;
ret = I_varphiF;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_1_()
{
double ret;
ret = -2*sr.d*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*sr.Cs1 + sr.Cs*sr.d*Power(Tan(x_old[4] + Thetap),2)*sr.Cs1 - (1 - sr.Cs*sr.d)*Power(Tan(x_old[4] + Thetap),2)*sr.Cs1 + Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) - sr.d*Tan(x_old[4] + Thetap)*sr.Cs2;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_2_()
{
double ret;
ret = -(Power(sr.Cs,2)*Power(Sec(x_old[4] + Thetap),2)) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv) + Power(sr.Cs,2)*Power(Tan(x_old[4] + Thetap),2) - Tan(x_old[4] + Thetap)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_3_()
{
double ret;
ret = Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4] + Thetap),3)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2)*Tan(x_old[4] + Thetap) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap) - sr.d*Power(Sec(x_old[4] + Thetap),2)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_4_()
{
double ret;
ret = Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4] + Thetap),3)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv) - 2*sr.Cs*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2)*Tan(x_old[4] + Thetap) + 3*Power(1 - sr.Cs*sr.d,2)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap) - sr.d*Power(Sec(x_old[4] + Thetap),2)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_5_()
{
double ret;
ret = (Cos(x_old[4] - x_old[5] - Thetap)*Power(1 - sr.Cs*sr.d,2)*Sec(x_old[4])*Power(Sec(x_old[4] + Thetap),3))/lv;
return ret;
}


double KinematicsSolver::calc_pd_alpha2_pd_X_2_1_()
{
double ret;
ret = -(sr.d*Power(Sec(x_old[4] + Thetap),2)*sr.Cs1);
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_2_()
{
double ret;
ret = -(sr.Cs*Power(Sec(x_old[4] + Thetap),2));
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_3_()
{
double ret;
ret = 2*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2)*Tan(x_old[4] + Thetap);
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_4_()
{
double ret;
ret = 2*(1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2)*Tan(x_old[4] + Thetap);
return ret;
}


double KinematicsSolver::calc_pd_alpha2_pd_X_2_5_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_pd_alpha2_pd_X_3_1_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_pd_alpha2_pd_X_3_2_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_pd_alpha2_pd_X_3_3_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_pd_alpha2_pd_X_3_4_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_pd_alpha2_pd_X_3_5_()
{
double ret;
ret = 0;
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_1_1_()
{
double ret;
ret = -(sr.d*(-(Power(sr.Cs,2)/(1 - sr.Cs*sr.d)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv))*Tan(x_old[4] + Thetap)*sr.Cs1) - sr.d*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap))*sr.Cs1 + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap))*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) - 2*sr.d*Sec(x_old[4] + Thetap)*sr.Cs1*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Tan(x_old[4] + Thetap)*(-((Power(sr.Cs,2)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (2*sr.Cs*sr.Cs1)/(1 - sr.Cs*sr.d) - Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*sr.Cs1 - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d))) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*(-(sr.d*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d))*sr.Cs1) - sr.d*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)*sr.Cs1 + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((sr.Cs*sr.d*Sin(x_old[4] + Thetap)*sr.Cs1)/Power(1 - sr.Cs*sr.d,2) + (Sin(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d))) - sr.d*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*sr.Cs2 + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((-2*sr.Cs*Cos(x_old[4] + Thetap)*Power(sr.d,2)*Power(sr.Cs1,2))/Power(1 - sr.Cs*sr.d,3) - (2*Cos(x_old[4] + Thetap)*sr.d*Power(sr.Cs1,2))/Power(1 - sr.Cs*sr.d,2) - (sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs2)/Power(1 - sr.Cs*sr.d,2) - (Cos(x_old[4] + Thetap)*sr.Cs2)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_1_2_()
{
double ret;
ret = -(sr.Cs*(-(Power(sr.Cs,2)/(1 - sr.Cs*sr.d)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv))*Tan(x_old[4] + Thetap)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*(-(sr.Cs*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d))) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) - (Power(sr.Cs,2)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)))/(1 - sr.Cs*sr.d) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) + (Power(sr.Cs,2)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2) - Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*sr.Cs1 + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((-2*Power(sr.Cs,2)*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,3) - (2*sr.Cs*Cos(x_old[4] + Thetap)*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_1_3_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2)*(-(Power(sr.Cs,2)/(1 - sr.Cs*sr.d)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)) + (1 - sr.Cs*sr.d)*Tan(x_old[4] + Thetap)*(-(sr.Cs*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d))) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d))*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*((1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) + (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv) + 2*(1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d))*Tan(x_old[4] + Thetap) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Power(Tan(x_old[4] + Thetap),2)) - sr.d*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d))*sr.Cs1 - sr.d*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)*sr.Cs1 + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((sr.Cs*sr.d*Sin(x_old[4] + Thetap)*sr.Cs1)/Power(1 - sr.Cs*sr.d,2) + (Sin(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_1_4_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),2)*(-(Power(sr.Cs,2)/(1 - sr.Cs*sr.d)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)) + (1 - sr.Cs*sr.d)*Tan(x_old[4] + Thetap)*(-(sr.Cs*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*((1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv + (Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4])*Tan(x_old[4]))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d))*Tan(x_old[4] + Thetap) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv)*Tan(x_old[4] + Thetap) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Power(Tan(x_old[4] + Thetap),2)) - sr.d*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv)*sr.Cs1 - sr.d*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)*sr.Cs1 + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((sr.Cs*sr.d*Sin(x_old[4] + Thetap)*sr.Cs1)/Power(1 - sr.Cs*sr.d,2) + (Sin(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_1_5_()
{
double ret;
ret = -((sr.Cs*Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*(((1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Sin(x_old[4] - x_old[5] - Thetap))/lv + (Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap))/lv) + (Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*((1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)))/lv - (Cos(x_old[4] - x_old[5] - Thetap)*sr.d*Sec(x_old[4])*Sec(x_old[4] + Thetap)*sr.Cs1)/lv;
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_2_1_()
{
double ret;
ret = -(sr.d*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv)*sr.Cs1) - sr.d*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap)*sr.Cs1 + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2)) - (Cos(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d)) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((sr.Cs*sr.d*Sin(x_old[4] + Thetap)*sr.Cs1)/Power(1 - sr.Cs*sr.d,2) + (Sin(x_old[4] + Thetap)*sr.Cs1)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_2_2_()
{
double ret;
ret = -(sr.Cs*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv)) - sr.Cs*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Tan(x_old[4] + Thetap);
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_2_3_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv + (Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4])*Tan(x_old[4]))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d))*Tan(x_old[4] + Thetap) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv)*Tan(x_old[4] + Thetap) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Power(Tan(x_old[4] + Thetap),2);
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_2_4_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Power(Sec(x_old[4] + Thetap),3)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) + (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv - (Power(Sec(x_old[4]),3)*Sin(x_old[4] - x_old[5] - Thetap))/lv - (2*Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4])*Tan(x_old[4]))/lv - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Power(Tan(x_old[4]),2))/lv) + 2*(1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4]))/lv) + (sr.Cs*Sin(x_old[4] + Thetap))/(1 - sr.Cs*sr.d) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap)*Tan(x_old[4]))/lv)*Tan(x_old[4] + Thetap) + (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((sr.Cs*Cos(x_old[4] + Thetap))/(1 - sr.Cs*sr.d)) - (Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv)*Power(Tan(x_old[4] + Thetap),2);
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_2_5_()
{
double ret;
ret = (1 - sr.Cs*sr.d)*Sec(x_old[4] + Thetap)*(-((Sec(x_old[4])*Sin(x_old[4] - x_old[5] - Thetap))/lv) + (Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4])*Tan(x_old[4]))/lv) + (Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap))/lv;
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_3_1_()
{
double ret;
ret = -((Cos(x_old[4] - x_old[5] - Thetap)*sr.d*Sec(x_old[4])*Sec(x_old[4] + Thetap)*sr.Cs1)/lv);
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_3_2_()
{
double ret;
ret = -((sr.Cs*Cos(x_old[4] - x_old[5] - Thetap)*Sec(x_old[4])*Sec(x_old[4] + Thetap))/lv);
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_3_3_()
{
double ret;
ret = ((1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Sin(x_old[4] - x_old[5] - Thetap))/lv + (Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap))/lv;
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_3_4_()
{
double ret;
ret = -(((1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Sin(x_old[4] - x_old[5] - Thetap))/lv) + (Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Tan(x_old[4]))/lv + (Cos(x_old[4] - x_old[5] - Thetap)*(1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Tan(x_old[4] + Thetap))/lv;
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_3_5_()
{
double ret;
ret = ((1 - sr.Cs*sr.d)*Sec(x_old[4])*Sec(x_old[4] + Thetap)*Sin(x_old[4] - x_old[5] - Thetap))/lv;
return ret;
}




double KinematicsSolver::calc_pd_G11_pd_X_1_()
{
double ret;
ret = (Cos(x_old[4] + Thetap)*sr.d*sr.Cs1)/Power(1 - sr.Cs*sr.d,2);
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_2_()
{
double ret;
ret = (sr.Cs*Cos(x_old[4] + Thetap))/Power(1 - sr.Cs*sr.d,2);
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_3_()
{
double ret;
ret = -(Sin(x_old[4] + Thetap)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_4_()
{
double ret;
ret = -(Sin(x_old[4] + Thetap)/(1 - sr.Cs*sr.d));
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_5_()
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
ret = ((k1*k2*z22 - k1*k2*d0d - ddd0d - k1*dd0d - k2*dd0d)*a0_dot + a0*(k1*k2*dd0d + dddd0d + (k1 + k2)*ddd0d))/Power(a0,2);
return ret;
}

double KinematicsSolver::calc_pd_W_pd_t_3_()
{
double ret;
ret = ((-ddthetap1d - (k3 + k4)*dthetap1d + k3*k4*(z32 - thetap1d))*a0_dot + a0*(dddthetap1d + (k3 + k4)*ddthetap1d + k3*k4*dthetap1d))/Power(a0,2);
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


double KinematicsSolver::calc_pd_W_pd_X_1_5_()
{
double ret;
ret = 0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_1_()
{
double ret;
ret = -((k1 + k2)*calc_pd_Z2_pd_X_1_1_()) - (k1*k2*calc_pd_Z2_pd_X_2_1_())/a0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_2_()
{
double ret;
ret = -((k1 + k2)*calc_pd_Z2_pd_X_1_2_()) - (k1*k2*calc_pd_Z2_pd_X_2_2_())/a0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_3_()
{
double ret;
ret = -((k1 + k2)*calc_pd_Z2_pd_X_1_3_()) - (k1*k2*calc_pd_Z2_pd_X_2_3_())/a0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_4_()
{
double ret;
ret = -((k1 + k2)*calc_pd_Z2_pd_X_1_4_()) - (k1*k2*calc_pd_Z2_pd_X_2_4_())/a0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_2_5_()
{
double ret;
ret = -((k1 + k2)*calc_pd_Z2_pd_X_1_5_()) - (k1*k2*calc_pd_Z2_pd_X_2_5_())/a0;
return ret;
}


double KinematicsSolver::calc_pd_W_pd_X_3_1_()
{
double ret;
ret = -((k3 + k4)*calc_pd_Z3_pd_X_1_1_()) - (k3*k4*calc_pd_Z3_pd_X_2_1_())/a0;
return ret;
}


double KinematicsSolver::calc_pd_W_pd_X_3_2_()
{
double ret;
ret = -((k3 + k4)*calc_pd_Z3_pd_X_1_2_()) - (k3*k4*calc_pd_Z3_pd_X_2_2_())/a0;
return ret;
}


double KinematicsSolver::calc_pd_W_pd_X_3_3_()
{
double ret;
ret = -((k3 + k4)*calc_pd_Z3_pd_X_1_3_()) - (k3*k4*calc_pd_Z3_pd_X_2_3_())/a0;
return ret;
}

double KinematicsSolver::calc_pd_W_pd_X_3_4_()
{
double ret;
ret = -((k3 + k4)*calc_pd_Z3_pd_X_1_4_()) - (k3*k4*calc_pd_Z3_pd_X_2_4_())/a0;
return ret;
}


double KinematicsSolver::calc_pd_W_pd_X_3_5_()
{
double ret;
ret = -((k3 + k4)*calc_pd_Z3_pd_X_1_5_()) - (k3*k4*calc_pd_Z3_pd_X_2_5_())/a0;
return ret;
}




//目標速度の時間微分
double KinematicsSolver::calc_pd_ud_pd_t_1_()
{
double ret;
ret = -((w1*(u3_act*calc_pd_G11_pd_X_5_() + u2_act*calc_pd_G11_pd_X_4_() + calc_SX_3_1_()*u1_act*calc_pd_G11_pd_X_3_() + calc_SX_2_1_()*u1_act*calc_pd_G11_pd_X_2_() + calc_SX_1_1_()*u1_act*calc_pd_G11_pd_X_1_()))/Power(calc_SX_1_1_(),2));
return ret;
}

double KinematicsSolver::calc_pd_ud_pd_t_2_()
{
double ret;
ret = (-((calc_alpha_3_3_()*(w1*calc_alpha_2_1_() - w2) + calc_alpha_2_3_()*(-(w1*calc_alpha_3_1_()) + w3))*(-(calc_alpha_3_3_()*(u3_act*calc_pd_alpha2_pd_X_2_5_() + u2_act*calc_pd_alpha2_pd_X_2_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha2_pd_X_2_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha2_pd_X_2_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha2_pd_X_2_1_())) + calc_alpha_3_2_()*(u3_act*calc_pd_alpha2_pd_X_3_5_() + u2_act*calc_pd_alpha2_pd_X_3_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha2_pd_X_3_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha2_pd_X_3_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha2_pd_X_3_1_()) + calc_alpha_2_3_()*(u3_act*calc_pd_alpha3_pd_X_2_5_() + u2_act*calc_pd_alpha3_pd_X_2_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha3_pd_X_2_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha3_pd_X_2_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha3_pd_X_2_1_()) - calc_alpha_2_2_()*(u3_act*calc_pd_alpha3_pd_X_3_5_() + u2_act*calc_pd_alpha3_pd_X_3_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha3_pd_X_3_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha3_pd_X_3_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha3_pd_X_3_1_()))) + (calc_alpha_2_3_()*calc_alpha_3_2_() - calc_alpha_2_2_()*calc_alpha_3_3_())*((-(w1*calc_alpha_3_1_()) + w3)*(u3_act*calc_pd_alpha2_pd_X_3_5_() + u2_act*calc_pd_alpha2_pd_X_3_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha2_pd_X_3_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha2_pd_X_3_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha2_pd_X_3_1_()) + (w1*calc_alpha_2_1_() - w2)*(u3_act*calc_pd_alpha3_pd_X_3_5_() + u2_act*calc_pd_alpha3_pd_X_3_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha3_pd_X_3_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha3_pd_X_3_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha3_pd_X_3_1_()) + calc_alpha_3_3_()*(u3_act*(w1*calc_pd_alpha2_pd_X_1_5_() - calc_pd_W_pd_X_2_5_()) + u2_act*(w1*calc_pd_alpha2_pd_X_1_4_() - calc_pd_W_pd_X_2_4_()) + w1*calc_SX_3_1_()*u1_act*calc_pd_alpha2_pd_X_1_3_() - calc_SX_3_1_()*u1_act*calc_pd_W_pd_X_2_3_() + w1*calc_SX_2_1_()*u1_act*calc_pd_alpha2_pd_X_1_2_() - calc_SX_2_1_()*u1_act*calc_pd_W_pd_X_2_2_() + w1*calc_SX_1_1_()*u1_act*calc_pd_alpha2_pd_X_1_1_() - calc_SX_1_1_()*u1_act*calc_pd_W_pd_X_2_1_()) + calc_alpha_2_3_()*(u3_act*calc_pd_W_pd_X_3_5_() + u2_act*calc_pd_W_pd_X_3_4_() + calc_SX_3_1_()*u1_act*calc_pd_W_pd_X_3_3_() + calc_SX_2_1_()*u1_act*calc_pd_W_pd_X_3_2_() - w1*(u3_act*calc_pd_alpha3_pd_X_1_5_() + u2_act*calc_pd_alpha3_pd_X_1_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha3_pd_X_1_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha3_pd_X_1_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha3_pd_X_1_1_()) + calc_SX_1_1_()*u1_act*calc_pd_W_pd_X_3_1_())))/Power(calc_alpha_2_3_()*calc_alpha_3_2_() - calc_alpha_2_2_()*calc_alpha_3_3_(),2);
return ret;
}

double KinematicsSolver::calc_pd_ud_pd_t_3_()
{
double ret;
ret = (-((calc_alpha_3_2_()*(-(w1*calc_alpha_2_1_()) + w2) + calc_alpha_2_2_()*(w1*calc_alpha_3_1_() - w3))*(-(calc_alpha_3_3_()*(u3_act*calc_pd_alpha2_pd_X_2_5_() + u2_act*calc_pd_alpha2_pd_X_2_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha2_pd_X_2_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha2_pd_X_2_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha2_pd_X_2_1_())) + calc_alpha_3_2_()*(u3_act*calc_pd_alpha2_pd_X_3_5_() + u2_act*calc_pd_alpha2_pd_X_3_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha2_pd_X_3_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha2_pd_X_3_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha2_pd_X_3_1_()) + calc_alpha_2_3_()*(u3_act*calc_pd_alpha3_pd_X_2_5_() + u2_act*calc_pd_alpha3_pd_X_2_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha3_pd_X_2_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha3_pd_X_2_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha3_pd_X_2_1_()) - calc_alpha_2_2_()*(u3_act*calc_pd_alpha3_pd_X_3_5_() + u2_act*calc_pd_alpha3_pd_X_3_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha3_pd_X_3_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha3_pd_X_3_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha3_pd_X_3_1_()))) + (calc_alpha_2_3_()*calc_alpha_3_2_() - calc_alpha_2_2_()*calc_alpha_3_3_())*((w1*calc_alpha_3_1_() - w3)*(u3_act*calc_pd_alpha2_pd_X_2_5_() + u2_act*calc_pd_alpha2_pd_X_2_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha2_pd_X_2_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha2_pd_X_2_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha2_pd_X_2_1_()) + (-(w1*calc_alpha_2_1_()) + w2)*(u3_act*calc_pd_alpha3_pd_X_2_5_() + u2_act*calc_pd_alpha3_pd_X_2_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha3_pd_X_2_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha3_pd_X_2_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha3_pd_X_2_1_()) + calc_alpha_3_2_()*(u3_act*calc_pd_W_pd_X_2_5_() + u2_act*calc_pd_W_pd_X_2_4_() + calc_SX_3_1_()*u1_act*calc_pd_W_pd_X_2_3_() + calc_SX_2_1_()*u1_act*calc_pd_W_pd_X_2_2_() - w1*(u3_act*calc_pd_alpha2_pd_X_1_5_() + u2_act*calc_pd_alpha2_pd_X_1_4_() + calc_SX_3_1_()*u1_act*calc_pd_alpha2_pd_X_1_3_() + calc_SX_2_1_()*u1_act*calc_pd_alpha2_pd_X_1_2_() + calc_SX_1_1_()*u1_act*calc_pd_alpha2_pd_X_1_1_()) + calc_SX_1_1_()*u1_act*calc_pd_W_pd_X_2_1_()) + calc_alpha_2_2_()*(u3_act*(w1*calc_pd_alpha3_pd_X_1_5_() - calc_pd_W_pd_X_3_5_()) + u2_act*(w1*calc_pd_alpha3_pd_X_1_4_() - calc_pd_W_pd_X_3_4_()) + w1*calc_SX_3_1_()*u1_act*calc_pd_alpha3_pd_X_1_3_() - calc_SX_3_1_()*u1_act*calc_pd_W_pd_X_3_3_() + w1*calc_SX_2_1_()*u1_act*calc_pd_alpha3_pd_X_1_2_() - calc_SX_2_1_()*u1_act*calc_pd_W_pd_X_3_2_() + w1*calc_SX_1_1_()*u1_act*calc_pd_alpha3_pd_X_1_1_() - calc_SX_1_1_()*u1_act*calc_pd_W_pd_X_3_1_())))/Power(calc_alpha_2_3_()*calc_alpha_3_2_() - calc_alpha_2_2_()*calc_alpha_3_3_(),2);
return ret;
}



