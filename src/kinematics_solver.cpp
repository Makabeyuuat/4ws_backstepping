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
ret = Power(1 - c(s(t))*d(t),2)*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv) - c(s(t))*(1 - c(s(t))*d(t))*Power(Tan(phi1(t) + thetap(t)),2) - d(t)*Tan(phi1(t) + thetap(t))*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_alpha_2_2_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2);
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
ret = (1 - c(s(t))*d(t))*(-(Power(c(s(t)),2)/(1 - c(s(t))*d(t))) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) - d(t)*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*sr.Cs1 + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)));
return ret;
}

double KinematicsSolver::calc_alpha_3_2_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t));
return ret;
}

double KinematicsSolver::calc_alpha_3_3_()
{
double ret;
ret = -((Cos(phi1(t) - phi2(t))*(1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t)))/lv);
return ret;
}

//Z21,22,23とその偏微分の計算
double KinematicsSolver::calc_Z_2_1_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Tan(phi1(t) + thetap(t));
return ret;
}

double KinematicsSolver::calc_Z_2_2_()
{
double ret;
ret = d(t);
return ret;
}

double KinematicsSolver::calc_Z_3_1_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv);
return ret;
}

double KinematicsSolver::calc_Z_3_2_()
{
double ret;
ret = thetap(t);
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_1_()
{
ret = -(d(t)*Tan(phi1(t) + thetap(t))*sr.Cs1);
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_2_()
{
double ret;
ret = -(c(s(t))*Tan(phi1(t) + thetap(t)));
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_3_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2);
return ret;
}

double KinematicsSolver::calc_pd_Z2_pd_X_1_4_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2);
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
ret = -(d(t)*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*sr.Cs1) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)));
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_2_()
{
double ret;
ret = -(Power(c(s(t)),2)/(1 - c(s(t))*d(t))) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv);
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_3_()
{
double ret;
ret = c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t));
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_4_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t));
return ret;
}

double KinematicsSolver::calc_pd_Z3_pd_X_1_5_()
{
double ret;
ret = -((Cos(phi1(t) - phi2(t))*(1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t)))/lv);
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
ret = Cos(phiR(t) + thetap(t) + thetat(t))*nu1(t) - Sin(phiR(t) + thetap(t) + thetat(t))*uact1(t)*(uact2(t) + c(s(t))*(calc_SX_1_1_()*u1) + (calc_SX_3_1_()*u1));
return ret;
}

double KinematicsSolver::calc_aqd_2_()
{
double ret;
ret = nu1(t)*Sin(phiR(t) + thetap(t) + thetat(t)) + Cos(phiR(t) + thetap(t) + thetat(t))*uact1(t)*(uact2(t) + c(s(t))*(calc_SX_1_1_()*u1) + (calc_SX_3_1_()*u1));
return ret;
}

double KinematicsSolver::calc_aqd_3_()
{
double ret;
ret = athetapd + asd*c(s(t)) + sr.Cs1*Power((calc_SX_1_1_()*u1),2);
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
ret = (Cos(phiR(t))*nu1(t) - Sin(phiR(t))*uact1(t)*uact2(t))/RADIUS;
return ret;
}

double KinematicsSolver::calc_aqd_6_()
{
double ret;
ret = nu3(t);
return ret;
}

double KinematicsSolver::calc_aqd_7_()
{
double ret;
ret = (Cos(phiF(t))*nu1(t) - Sin(phiF(t))*uact1(t)*uact3(t))/RADIUS;
return ret;
}

double KinematicsSolver::calc_Axi_1_1_()
{
double ret;
ret = Sin(phiF(t) + theta0);
return ret;
}

double KinematicsSolver::calc_Axi_1_2_()
{
double ret;
ret = -Cos(phiF(t) + theta0);
return ret;
}

double KinematicsSolver::calc_Axi_1_3_()
{
double ret;
ret = -(LENGTH*Cos(phi2(t)));
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
ret = Sin(phiR(t) + theta0);
return ret;
}

double KinematicsSolver::calc_Axi_2_2_()
{
double ret;
ret = -Cos(phiR(t) + theta0);
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
ret = -Cos(phiR(t) + theta0);
return ret;
}

double KinematicsSolver::calc_Axi_3_2_()
{
double ret;
ret = -Sin(phiR(t) + theta0);
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
ret = -Cos(phiF(t) + theta0);
return ret;
}

double KinematicsSolver::calc_Axi_4_2_()
{
double ret;
ret = -Sin(phiF(t) + theta0);
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
ret = RADIUS;
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
ret = -(LENGTH*(MASSCAR + 2*MASSWHEEL)*Cos(theta0)*x_d[3]);
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
ret = -(LENGTH*(MASSCAR + 2*MASSWHEEL)*Sin(theta0)*x_d[3]);
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
ret = Cos(phi1(t) + thetap(t))/(1 - c(s(t))*d(t));
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
ret = Sin(phi1(t) + thetap(t));
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
ret = -((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv;
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
ret = 0;
return ret;
}

double KinematicsSolver::calc_d_SX_d_t_1_1_()
{
double ret;
ret = -((Sin(phi1(t) + thetap(t))*((-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*u1(t) + u2(t)))/(1 - c(s(t))*d(t))) - (Cos(phi1(t) + thetap(t))*(-(c(s(t))*Sin(phi1(t) + thetap(t))*u1(t)) - (Cos(phi1(t) + thetap(t))*d(t)*u1(t)*sr.Cs1)/(1 - c(s(t))*d(t))))/Power(1 - c(s(t))*d(t),2);
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
ret = Cos(phi1(t) + thetap(t))*((-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*u1(t) + u2(t));
return ret
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
ret = (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t))*u2(t))/lv + (c(s(t))*Sin(phi1(t) + thetap(t))*((-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*u1(t) + u2(t)))/(1 - c(s(t))*d(t)) + (Cos(phi1(t) - phi2(t))*Sec(phi1(t))*(u2(t) - u3(t)))/lv - (Power(Cos(phi1(t) + thetap(t)),2)*u1(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2) + (c(s(t))*Cos(phi1(t) + thetap(t))*(-(c(s(t))*Sin(phi1(t) + thetap(t))*u1(t)) - (Cos(phi1(t) + thetap(t))*d(t)*u1(t)*sr.Cs1)/(1 - c(s(t))*d(t))))/Power(1 - c(s(t))*d(t),2);
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
ret = IvarphiR;;
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
ret = IphiF;
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
ret = IvarphiF;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_1_()
{
double ret;
ret = -2*d(t)*(1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*sr.Cs1 + c(s(t))*d(t)*Power(Tan(phi1(t) + thetap(t)),2)*sr.Cs1 - (1 - c(s(t))*d(t))*Power(Tan(phi1(t) + thetap(t)),2)*sr.Cs1 + Power(1 - c(s(t))*d(t),2)*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t))) - d(t)*Tan(phi1(t) + thetap(t))*sr.Cs2;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_2_()
{
double ret;
ret = -(Power(c(s(t)),2)*Power(Sec(phi1(t) + thetap(t)),2)) - 2*c(s(t))*(1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv) + Power(c(s(t)),2)*Power(Tan(phi1(t) + thetap(t)),2) - Tan(phi1(t) + thetap(t))*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_3_()
{
double ret;
ret = -(c(s(t))*(1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2)*Tan(phi1(t) + thetap(t))) + 3*Power(1 - c(s(t))*d(t),2)*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t)) - d(t)*Power(Sec(phi1(t) + thetap(t)),2)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_4_()
{
double ret;
ret = Power(1 - c(s(t))*d(t),2)*Power(Sec(phi1(t) + thetap(t)),3)*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv) - 2*c(s(t))*(1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2)*Tan(phi1(t) + thetap(t)) + 3*Power(1 - c(s(t))*d(t),2)*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t)) - d(t)*Power(Sec(phi1(t) + thetap(t)),2)*sr.Cs1;
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_1_5_()
{
double ret;
ret = -((Cos(phi1(t) - phi2(t))*Power(1 - c(s(t))*d(t),2)*Sec(phi1(t))*Power(Sec(phi1(t) + thetap(t)),3))/lv);
return ret;
}


double KinematicsSolver::calc_pd_alpha2_pd_X_2_1_()
{
double ret;
ret = -(d(t)*Power(Sec(phi1(t) + thetap(t)),2)*sr.Cs1);
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_2_()
{
double ret;
ret = -(c(s(t))*Power(Sec(phi1(t) + thetap(t)),2));
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_3_()
{
double ret;
ret = 2*(1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2)*Tan(phi1(t) + thetap(t));
return ret;
}

double KinematicsSolver::calc_pd_alpha2_pd_X_2_4_()
{
double ret;
ret = 2*(1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2)*Tan(phi1(t) + thetap(t));
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
ret = -(d(t)*(-(Power(c(s(t)),2)/(1 - c(s(t))*d(t))) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv))*Tan(phi1(t) + thetap(t))*sr.Cs1) - d(t)*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t)))*sr.Cs1 + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t)))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t))) - 2*d(t)*Sec(phi1(t) + thetap(t))*sr.Cs1*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t))) + (1 - c(s(t))*d(t))*Tan(phi1(t) + thetap(t))*(-((Power(c(s(t)),2)*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (2*c(s(t))*sr.Cs1)/(1 - c(s(t))*d(t)) - Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*sr.Cs1 - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*(Tan(phi1(t) + thetap(t))*sr.Cs1 - d(t)*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))*sr.Cs1 + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)))) - d(t)*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*sr.Cs2 + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((-2*c(s(t))*Cos(phi1(t) + thetap(t))*Power(d(t),2)*Power(sr.Cs1,2))/Power(1 - c(s(t))*d(t),3) - (2*Cos(phi1(t) + thetap(t))*d(t)*Power(sr.Cs1,2))/Power(1 - c(s(t))*d(t),2) - (c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs2)/Power(1 - c(s(t))*d(t),2) - (Cos(phi1(t) + thetap(t))*sr.Cs2)/(1 - c(s(t))*d(t)));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_1_2_()
{
double ret;
ret = -(c(s(t))*(-(Power(c(s(t)),2)/(1 - c(s(t))*d(t))) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv))*Tan(phi1(t) + thetap(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*(-((Power(c(s(t)),2)*Tan(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) - (Power(c(s(t)),2)*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))))/(1 - c(s(t))*d(t)) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) + (Power(c(s(t)),2)*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2) - Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*sr.Cs1 + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((-2*Power(c(s(t)),2)*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),3) - (2*c(s(t))*Cos(phi1(t) + thetap(t))*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_1_3_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2)*(-(Power(c(s(t)),2)/(1 - c(s(t))*d(t))) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)) + (1 - c(s(t))*d(t))*Tan(phi1(t) + thetap(t))*(-((Power(c(s(t)),2)*Tan(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) + c(s(t))*Tan(phi1(t) + thetap(t))*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*(c(s(t))*Power(Sec(phi1(t) + thetap(t)),2) + (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv) + c(s(t))*Power(Tan(phi1(t) + thetap(t)),2) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Power(Tan(phi1(t) + thetap(t)),2)) - (c(s(t))*d(t)*Tan(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)) - d(t)*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))*sr.Cs1 + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((c(s(t))*d(t)*Sin(phi1(t) + thetap(t))*sr.Cs1)/Power(1 - c(s(t))*d(t),2) + (Sin(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_1_4_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),2)*(-(Power(c(s(t)),2)/(1 - c(s(t))*d(t))) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)) + (1 - c(s(t))*d(t))*Tan(phi1(t) + thetap(t))*(-(c(s(t))*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv)) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv)*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*(c(s(t))*Power(Sec(phi1(t) + thetap(t)),2) + (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv)*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Power(Tan(phi1(t) + thetap(t)),2)) - d(t)*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv)*sr.Cs1 - d(t)*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))*sr.Cs1 + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((c(s(t))*d(t)*Sin(phi1(t) + thetap(t))*sr.Cs1)/Power(1 - c(s(t))*d(t),2) + (Sin(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_1_5_()
{
double ret;
ret = (c(s(t))*Cos(phi1(t) - phi2(t))*(1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t) + thetap(t)))/lv - (Cos(phi1(t) - phi2(t))*Power(1 - c(s(t))*d(t),2)*Sec(phi1(t))*Power(Sec(phi1(t) + thetap(t)),2)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t)))/lv - (Cos(phi1(t) - phi2(t))*(1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*(c(s(t))*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))))/lv + (Cos(phi1(t) - phi2(t))*d(t)*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*sr.Cs1)/lv;
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_2_1_()
{
double ret;
ret = -(d(t)*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv)*sr.Cs1) - d(t)*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t))*sr.Cs1 + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2)) - (Cos(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t))) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((c(s(t))*d(t)*Sin(phi1(t) + thetap(t))*sr.Cs1)/Power(1 - c(s(t))*d(t),2) + (Sin(phi1(t) + thetap(t))*sr.Cs1)/(1 - c(s(t))*d(t)));
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_2_2_()
{
double ret;
ret = -(c(s(t))*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv)) - c(s(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Tan(phi1(t) + thetap(t));
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_2_3_()
{
double ret;
ret = c(s(t)) + (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv)*Tan(phi1(t) + thetap(t)) + c(s(t))*Power(Tan(phi1(t) + thetap(t)),2) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Power(Tan(phi1(t) + thetap(t)),2);
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_2_4_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Power(Sec(phi1(t) + thetap(t)),3)*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) - (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv + (Power(Sec(phi1(t)),3)*Sin(phi1(t) - phi2(t)))/lv + (2*Cos(phi1(t) - phi2(t))*Sec(phi1(t))*Tan(phi1(t)))/lv + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Power(Tan(phi1(t)),2))/lv) + 2*(1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((Cos(phi1(t) - phi2(t))*Sec(phi1(t)))/lv + (c(s(t))*Sin(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t)) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t))*Tan(phi1(t)))/lv)*Tan(phi1(t) + thetap(t)) + (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*(-((c(s(t))*Cos(phi1(t) + thetap(t)))/(1 - c(s(t))*d(t))) + (Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv)*Power(Tan(phi1(t) + thetap(t)),2);
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_2_5_()
{
double ret;
ret = (1 - c(s(t))*d(t))*Sec(phi1(t) + thetap(t))*((Sec(phi1(t))*Sin(phi1(t) - phi2(t)))/lv - (Cos(phi1(t) - phi2(t))*Sec(phi1(t))*Tan(phi1(t)))/lv) - (Cos(phi1(t) - phi2(t))*(1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t) + thetap(t)))/lv;
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_3_1_()
{
double ret;
ret = (Cos(phi1(t) - phi2(t))*d(t)*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*sr.Cs1)/lv;
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_3_2_()
{
double ret;
ret = (c(s(t))*Cos(phi1(t) - phi2(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t)))/lv;
return ret;
}


double KinematicsSolver::calc_pd_alpha3_pd_X_3_3_()
{
double ret;
ret = -((Cos(phi1(t) - phi2(t))*(1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t) + thetap(t)))/lv);
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_3_4_()
{
double ret;
ret = ((1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*Sin(phi1(t) - phi2(t)))/lv - (Cos(phi1(t) - phi2(t))*(1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t)))/lv - (Cos(phi1(t) - phi2(t))*(1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*Tan(phi1(t) + thetap(t)))/lv;
return ret;
}

double KinematicsSolver::calc_pd_alpha3_pd_X_3_5_()
{
double ret;
ret = -(((1 - c(s(t))*d(t))*Sec(phi1(t))*Sec(phi1(t) + thetap(t))*Sin(phi1(t) - phi2(t)))/lv);
return ret;
}




double KinematicsSolver::calc_pd_G11_pd_X_1_()
{
double ret;
ret = (Cos(phi1(t) + thetap(t))*d(t)*sr.Cs1)/Power(1 - c(s(t))*d(t),2);
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_2_()
{
double ret;
ret = (c(s(t))*Cos(phi1(t) + thetap(t)))/Power(1 - c(s(t))*d(t),2);
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_3_()
{
double ret;
ret = -(Sin(phi1(t) + thetap(t))/(1 - c(s(t))*d(t)));
return ret;
}

double KinematicsSolver::calc_pd_G11_pd_X_4_()
{
double ret;
ret = -(Sin(phi1(t) + thetap(t))/(1 - c(s(t))*d(t)));
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
ret = ((k1*k2*Z22 - k1*k2*Pattern(d,Blank(d))(t) - Pattern(ddotd,Blank(d))(t) - k1*Pattern(dotd,Blank(d))(t) - k2*Pattern(dotd,Blank(d))(t))*a0_dot + a0*(k1*k2*Derivative(1)(Pattern(d,Blank(d)))(t) + Derivative(1)(Pattern(ddotd,Blank(d)))(t) + (k1 + k2)*Derivative(1)(Pattern(dotd,Blank(d)))(t)))/Power(a0,2);
return ret;
}

double KinematicsSolver::calc_pd_W_pd_t_3_()
{
double ret;
ret = ((-Pattern(ddotthetap,Blank(d))(t) - (k3 + k4)*Pattern(dotthetap,Blank(d))(t) + k3*k4*(Z32 - Pattern(thetap,Blank(d))(t)))*a0_dot + a0*(Derivative(1)(Pattern(ddotthetap,Blank(d)))(t) + (k3 + k4)*Derivative(1)(Pattern(dotthetap,Blank(d)))(t) + k3*k4*Derivative(1)(Pattern(thetap,Blank(d)))(t)))/Power(a0,2);
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
ret = -((W1*(x_d[5]*calc_pd_G11_pd_X_5_() + x_d[4]*calc_pd_G11_pd_X_4_() + (calc_SX_3_1_()*u1)*calc_pd_G11_pd_X_3_() + (calc_SX_2_1_()*u1)*calc_pd_G11_pd_X_2_() + (calc_SX_1_1_()*u1)*calc_pd_G11_pd_X_1_()))/Power(calc_SX_1_1_(),2));
return ret;
}

double KinematicsSolver::calc_pd_ud_pd_t_2_()
{
double ret;
ret = (-((calc_alpha_3_3_()*(W1*calc_alpha_2_1_() - w2) + calc_alpha_2_3_()*(-(W1*calc_alpha_3_1_()) + w3))*(-(calc_alpha_3_3_()*(x_d[5]*calc_pd_alpha2_pd_X_2_5_() + x_d[4]*calc_pd_alpha2_pd_X_2_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha2_pd_X_2_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha2_pd_X_2_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha2_pd_X_2_1_())) + calc_alpha_3_2_()*(x_d[5]*calc_pd_alpha2_pd_X_3_5_() + x_d[4]*calc_pd_alpha2_pd_X_3_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha2_pd_X_3_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha2_pd_X_3_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha2_pd_X_3_1_()) + calc_alpha_2_3_()*(x_d[5]*calc_pd_alpha3_pd_X_2_5_() + x_d[4]*calc_pd_alpha3_pd_X_2_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha3_pd_X_2_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha3_pd_X_2_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha3_pd_X_2_1_()) - calc_alpha_2_2_()*(x_d[5]*calc_pd_alpha3_pd_X_3_5_() 
        + x_d[4]*calc_pd_alpha3_pd_X_3_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha3_pd_X_3_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha3_pd_X_3_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha3_pd_X_3_1_()))) + (calc_alpha_2_3_()*calc_alpha_3_2_() - calc_alpha_2_2_()*calc_alpha_3_3_())*((-(W1*calc_alpha_3_1_()) + w3)*(x_d[5]*calc_pd_alpha2_pd_X_3_5_() + x_d[4]*calc_pd_alpha2_pd_X_3_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha2_pd_X_3_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha2_pd_X_3_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha2_pd_X_3_1_()) + (W1*calc_alpha_2_1_() - w2)*(x_d[5]*calc_pd_alpha3_pd_X_3_5_() + x_d[4]*calc_pd_alpha3_pd_X_3_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha3_pd_X_3_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha3_pd_X_3_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha3_pd_X_3_1_()) + calc_alpha_3_3_()*(x_d[5]*(W1*calc_pd_alpha2_pd_X_1_5_() - calc_pd_W_pd_X_2_5_()) + x_d[4]*(W1*calc_pd_alpha2_pd_X_1_4_() 
        - calc_pd_W_pd_X_2_4_()) + W1*(calc_SX_3_1_()*u1)*calc_pd_alpha2_pd_X_1_3_() - (calc_SX_3_1_()*u1)*calc_pd_W_pd_X_2_3_() + W1*(calc_SX_2_1_()*u1)*calc_pd_alpha2_pd_X_1_2_() - (calc_SX_2_1_()*u1)*calc_pd_W_pd_X_2_2_() + W1*(calc_SX_1_1_()*u1)*calc_pd_alpha2_pd_X_1_1_() - (calc_SX_1_1_()*u1)*calc_pd_W_pd_X_2_1_()) + calc_alpha_2_3_()*(x_d[5]*calc_pd_W_pd_X_3_5_() + x_d[4]*calc_pd_W_pd_X_3_4_() + (calc_SX_3_1_()*u1)*calc_pd_W_pd_X_3_3_() + (calc_SX_2_1_()*u1)*calc_pd_W_pd_X_3_2_() - W1*(x_d[5]*calc_pd_alpha3_pd_X_1_5_() + x_d[4]*calc_pd_alpha3_pd_X_1_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha3_pd_X_1_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha3_pd_X_1_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha3_pd_X_1_1_()) + (calc_SX_1_1_()*u1)*calc_pd_W_pd_X_3_1_())))/Power(calc_alpha_2_3_()*calc_alpha_3_2_() - calc_alpha_2_2_()*calc_alpha_3_3_(),2);
return ret;
}

double KinematicsSolver::calc_pd_ud_pd_t_3_()
{
double ret;
ret = (-((calc_alpha_3_2_()*(-(W1*calc_alpha_2_1_()) + w2) + calc_alpha_2_2_()*(W1*calc_alpha_3_1_() - w3))*(-(calc_alpha_3_3_()*(x_d[5]*calc_pd_alpha2_pd_X_2_5_() + x_d[4]*calc_pd_alpha2_pd_X_2_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha2_pd_X_2_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha2_pd_X_2_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha2_pd_X_2_1_())) + calc_alpha_3_2_()*(x_d[5]*calc_pd_alpha2_pd_X_3_5_() + x_d[4]*calc_pd_alpha2_pd_X_3_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha2_pd_X_3_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha2_pd_X_3_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha2_pd_X_3_1_()) + calc_alpha_2_3_()*(x_d[5]*calc_pd_alpha3_pd_X_2_5_() + x_d[4]*calc_pd_alpha3_pd_X_2_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha3_pd_X_2_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha3_pd_X_2_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha3_pd_X_2_1_()) - calc_alpha_2_2_()*(x_d[5]*calc_pd_alpha3_pd_X_3_5_() + x_d[4]*calc_pd_alpha3_pd_X_3_4_() 
        + (calc_SX_3_1_()*u1)*calc_pd_alpha3_pd_X_3_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha3_pd_X_3_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha3_pd_X_3_1_()))) + (calc_alpha_2_3_()*calc_alpha_3_2_() - calc_alpha_2_2_()*calc_alpha_3_3_())*((W1*calc_alpha_3_1_() - w3)*(x_d[5]*calc_pd_alpha2_pd_X_2_5_() + x_d[4]*calc_pd_alpha2_pd_X_2_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha2_pd_X_2_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha2_pd_X_2_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha2_pd_X_2_1_()) + (-(W1*calc_alpha_2_1_()) + w2)*(x_d[5]*calc_pd_alpha3_pd_X_2_5_() + x_d[4]*calc_pd_alpha3_pd_X_2_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha3_pd_X_2_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha3_pd_X_2_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha3_pd_X_2_1_()) + calc_alpha_3_2_()*(x_d[5]*calc_pd_W_pd_X_2_5_() + x_d[4]*calc_pd_W_pd_X_2_4_() + (calc_SX_3_1_()*u1)*calc_pd_W_pd_X_2_3_() + (calc_SX_2_1_()*u1)*calc_pd_W_pd_X_2_2_() - W1*(x_d[5]*calc_pd_alpha2_pd_X_1_5_() 
        + x_d[4]*calc_pd_alpha2_pd_X_1_4_() + (calc_SX_3_1_()*u1)*calc_pd_alpha2_pd_X_1_3_() + (calc_SX_2_1_()*u1)*calc_pd_alpha2_pd_X_1_2_() + (calc_SX_1_1_()*u1)*calc_pd_alpha2_pd_X_1_1_()) + (calc_SX_1_1_()*u1)*calc_pd_W_pd_X_2_1_()) + calc_alpha_2_2_()*(x_d[5]*(W1*calc_pd_alpha3_pd_X_1_5_() - calc_pd_W_pd_X_3_5_()) + x_d[4]*(W1*calc_pd_alpha3_pd_X_1_4_() - calc_pd_W_pd_X_3_4_()) + W1*(calc_SX_3_1_()*u1)*calc_pd_alpha3_pd_X_1_3_() - (calc_SX_3_1_()*u1)*calc_pd_W_pd_X_3_3_() + W1*(calc_SX_2_1_()*u1)*calc_pd_alpha3_pd_X_1_2_() - (calc_SX_2_1_()*u1)*calc_pd_W_pd_X_3_2_() + W1*(calc_SX_1_1_()*u1)*calc_pd_alpha3_pd_X_1_1_() - (calc_SX_1_1_()*u1)*calc_pd_W_pd_X_3_1_())))/Power(calc_alpha_2_3_()*calc_alpha_3_2_() - calc_alpha_2_2_()*calc_alpha_3_3_(),2);
return ret;
}



