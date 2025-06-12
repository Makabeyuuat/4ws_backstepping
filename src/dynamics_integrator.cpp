#include "dynamics_integrator.hpp"
#include "initial.hpp"
#include "mathFunc.h"		// 数学関数のヘッダーファイル
#include "Bezier.h"			// Bezier 曲線の関数
#include "vehicle.hpp"       // Vehicle クラスの宣言
#include "callback.hpp"     // コールバック関数の宣言
#include "initial.hpp"		// 初期値設定のヘッダーファイル
#include "getInputValue.hpp"
#include "pidTau.hpp"
#include "differential_equations_dynamics.hpp"
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/SVD> 

using namespace std;
DynamicsIntegrator::DynamicsIntegrator(double m_b,
                                       double I_theta,
                                       double lv,
                                       double g,
                                       double rho,
                                       const PIDGains& drive_gains,
                                       const PIDGains& steer_gains,
                                       double dt)
    : m_b(m_b), I_theta_(I_theta), lv_(lv), g_(g), rho_(rho), dt_(dt),
      drive_pid_(drive_gains, dt), steer_pid_(steer_gains, dt) {}



Eigen::Matrix<double,4,1> DynamicsIntegrator::computeXAlpha(
    const std::vector<double> x,
    const std::vector<double> x_d,
    double u1,
    double u2)
{
    // 実際の前進速度と前輪操舵角速度
    double u1_act = xdot * std::cos(theta) + ydot * std::sin(theta);
    double u2_act = x_d[4];

    // 結果格納用ベクトル
    Eigen::Matrix<double,4,1> Xalpha;

    //状態変数ベクトル
    Eigen::Matrix<double,4,1> Sx;
    Sx<<
        Cos(Thetap)/(1 - c(s(t))*d(t)), 
        Sin(Thetap),  
        -((c(s(t))*Cos(Thetap))/(1 - c(s(t))*d(t))) + Tan(phi(t))/LENGTH, 
        1.0;

    Eigen::Matrix<double,4,1> dSx;
    dSx <<
        -((Sin(Thetap)*(-((c(s(t))*Cos(Thetap))/(1 - c(s(t))*d(t))) + Tan(phi(t))/LENGTH)*uact1(t))/(1 - c(s(t))*d(t))) - (Cos(Thetap)*(-(c(s(t))*Sin(Thetap)*uact1(t)) - (Cos(Thetap)*d(t)*uact1(t)*Derivative(1)(c)(s(t)))/(1 - c(s(t))*d(t))))/Power(1 - c(s(t))*d(t),2),
        Cos(Thetap)*(-((c(s(t))*Cos(Thetap))/(1 - c(s(t))*d(t))) + Tan(phi(t))/LENGTH)*uact1(t),
         (c(s(t))*Sin(Thetap)*(-((c(s(t))*Cos(Thetap))/(1 - c(s(t))*d(t))) + Tan(phi(t))/LENGTH)*uact1(t))/(1 - c(s(t))*d(t)) + (Power(Sec(phi(t)),2)*uact2(t))/LENGTH - (Power(Cos(Thetap),2)*uact1(t)*Derivative(1)(c)(s(t)))/Power(1 - c(s(t))*d(t),2) + (c(s(t))*Cos(Thetap)*(-(c(s(t))*Sin(Thetap)*uact1(t)) - (Cos(Thetap)*d(t)*uact1(t)*Derivative(1)(c)(s(t)))/(1 - c(s(t))*d(t))))/Power(1 - c(s(t))*d(t),2),
         0.0;

  

    
    // 速度誤差
    double r_b1 = u1 - u1_act;
    double r_b2 = u2 - phi_dot;

    double pdud_pdt1 = -((W1*(Derivative(1)(phi)(t)*Derivative(List(0,0,0,1))(G11)(List(s(t),d(t),Thetap,phi(t))) + Derivative(1)(thetap)(t)*Derivative(List(0,0,1,0))(G11)(List(s(t),d(t),Thetap,phi(t))) + Derivative(1)(d)(t)*Derivative(List(0,1,0,0))(G11)(List(s(t),d(t),Thetap,phi(t))) + Derivative(1)(s)(t)*Derivative(List(1,0,0,0))(G11)(List(s(t),d(t),Thetap,phi(t)))))/Power(G11(List(s(t),d(t),Thetap,phi(t))),2))
    double pdud_pdt2 = ((W1*alpha21(List(s(t),d(t),Thetap,phi(t))) - W2(List(s(t),d(t),Thetap,phi(t))))*(Derivative(1)(phi)(t)*Derivative(List(0,0,0,1))(alpha22)(List(s(t),d(t),Thetap,phi(t))) + Derivative(1)(thetap)(t)*Derivative(List(0,0,1,0))(alpha22)(List(s(t),d(t),Thetap,phi(t))) + Derivative(1)(d)(t)*Derivative(List(0,1,0,0))(alpha22)(List(s(t),d(t),Thetap,phi(t))) + Derivative(1)(s)(t)*Derivative(List(1,0,0,0))(alpha22)(List(s(t),d(t),Thetap,phi(t)))) + alpha22(List(s(t),d(t),Thetap,phi(t)))*(Derivative(1)(phi)(t)*(-(W1*Derivative(List(0,0,0,1))(alpha21)(List(s(t),d(t),Thetap,phi(t)))) + Derivative(List(0,0,0,1))(W2)(List(s(t),d(t),Thetap,phi(t)))) + Derivative(1)(thetap)(t)*(-(W1*Derivative(List(0,0,1,0))(alpha21)(List(s(t),d(t),Thetap,phi(t)))) + Derivative(List(0,0,1,0))(W2)(List(s(t),d(t),Thetap,phi(t)))) - W1*Derivative(1)(d)(t)*Derivative(List(0,1,0,0))(alpha21)(List(s(t),d(t),Thetap,phi(t))) + Derivative(1)(d)(t)*Derivative(List(0,1,0,0))(W2)(List(s(t),d(t),Thetap,phi(t))) - W1*Derivative(1)(s)(t)*Derivative(List(1,0,0,0))(alpha21)(List(s(t),d(t),Thetap,phi(t))) + Derivative(1)(s)(t)*Derivative(List(1,0,0,0))(W2)(List(s(t),d(t),Thetap,phi(t)))))/Power(alpha22(List(s(t),d(t),Thetap,phi(t))),2)

    //偏差ベクトル
    Eigen::Matrix<double,2,1> r_b;
    r_b <<r_b1, r_b2

    //ゲイン
    Eigen::Matrix<double,2,1> C;
    C << 
        1.0, 1.0,
        1.0, 1.0;
        

    //目標加速度νを計算
    nu1 = -C * r_b + pdud_pdt1;
    nu2 = -C * r_b + pdud_pdt2;


    // 目標加速度を各成分に割り当て
    Xalpha(0) = dSx(0) * u1act + Sx(0) * nu1;   // alpha_sd
    Xalpha(1) = dSx(1) * u1act + Sx(1) * nu1;   // alpha_dd
    Xalpha(2) = dSx(2) * u1act + Sx(2) * nu1;   // alpha_theta             
    Xalpha(3) = dSx(3) * u1act + Sx(3) * nu2;   // alpha_phi

    return Xalpha;
}



Eigen::Matrix<double,6,1> DynamicsIntegrator::computeAlpha(
    const Eigen::Matrix<double,6,1>& q,
    const Eigen::Matrix<double,6,1>& qdot,
    double u1,
    double u2)
{
    //状態変数ベクトルの目標加速度 α_d の取得
    Eigen::Matrix<double,4,1> Xalpha = computeXAlpha(x_old, x_dd, u1, u2);
    
     
    double theta   = q(2);            // 車体ヨー角
    double xdot    = qdot(0), ydot   = qdot(1);
    double phi_dot = qdot(3);

    // 実際の前進速度（路面接触点速度の x,y 成分投影）
    double u1_act = xdot * std::cos(theta) + ydot * std::sin(theta);
    // 速度誤差
    double r_b1 = u1 - u1_act;
    double r_b2 = u2 - phi_dot;

    // ゲイン設定(要調整)
    const double Kp1 = 1.0, Kp2 = 1.0;
    // 目標加速度を各成分に割り当て
    alpha(0) = Kp1 * e1;   // α_x
    alpha(1) = 0.0;        // α_y (通常は運動学モデル由来)
    alpha(2) = 0.0;        // α_θ
    alpha(3) = Kp2 * e2;   // α_φ (ステア角加速度)
    // 前輪／後輪回転角加速度
    alpha(4) = alpha(0) / wheelRadius;  // α_ψf
    alpha(5) = alpha(0) / wheelRadius;  // α_ψr

    return alpha;
}


void DynamicsIntegrator::step(
    const Eigen::Matrix<double,6,1>& q,
    const Eigen::Matrix<double,6,1>& qdot,
    double u1,
    double u2,
    Eigen::Matrix<double,6,1>& qdd)
{
    // --- 1) 状態分解 ---
    double x      = q(0);
    double y      = q(1);
    double theta = q(2);
    double phi    = q(3);
    double xdot   = qdot(0);
    double ydot   = qdot(1);
    double thetadot = qdot(2);
    double phidot = qdot(3);

    // --- 2) 目標加速度 α_d の取得 ---
    Eigen::Matrix<double,6,1> alpha = computeAlpha(q, qdot, u1, u2);
    // (x,y,θ) 部分だけ抜き出し
    Eigen::Vector3d alpha3 = alpha.head<3>();

    // --- 3) 拘束 λ の導出 (最小ノルム解) ---
    // 拘束式 f1,f2,f3,f4 の時間微分から得る 3×4 行列 A3
    Eigen::Matrix<double,3,4> A3;
    A3 <<
      sin(theta+phi), sin(theta), cos(theta+phi), cos(theta),
     -cos(theta+phi),-cos(theta), sin(theta+phi), sin(theta),
      lv_/2 * cos(phi), 0.0,              0.0,            0.0;

    // 運動方程式上位3行: 慣性 M3, コリオリ C3, 重力 K3
    Eigen::Matrix3d M3;
    M3 <<
      m_b + 2*m_w_r, 0, -m_b*lv_/2*sin(theta),
      0, m_b + 2*m_w_r,  m_b*lv_/2*cos(theta),
     -m_b*lv_/2*sin(theta), m_b*lv_/2*cos(theta), I_theta_ + m_b*lv_*lv_/4;

    Eigen::Vector3d C3;
    C3 <<
      -m_b*lv_*thetadot*cos(theta)*thetadot,
      -m_b*lv_*thetadot*sin(theta)*thetadot,
      -m_b*lv_*thetadot*(xdot*cos(theta)+ydot*sin(theta));

    Eigen::Vector3d K3;
    K3 << (m_b+2*m_w_r)*g_*sin(rho_), 0.0, 0.0;
    // 右辺 = -(M3*α3 + C3 + K3)
    
    Eigen::Vector3d rhs = - (M3*alpha3 + C3 + K3);
    // 擬似逆行列(最小ノルム)による λ(4×1)
    Eigen::Matrix4x3d At = A3.transpose();
    Eigen::Matrix3d AAt = A3 * At;
    Eigen::Vector4d lambda = At * AAt.inverse() * rhs;

    // --- 4) 駆動力 Q の計算 ---
    Q_phi   = J_phi   * alpha(3);          // ステアトルク
    Q_psi_f = I_psif  * alpha(4) + wheelRadius * lambda(2); // 前輪トルク
    Q_psi_r = I_psir  * alpha(5) + wheelRadius * lambda(3); // 後輪トルク

    // --- 5) 全自由度動力学の解 ---
    Eigen::Matrix<double,6,6> Mxi = Eigen::Matrix<double,6,6>::Zero();
    Mxi.block<3,3>(0,0) = M3;
    Mxi(3,3) = J_phi;
    Mxi(4,4) = I_psif;
    Mxi(5,5) = I_psir;

    Eigen::Matrix<double,6,6> Cxi = Eigen::Matrix<double,6,6>::Zero();
    Cxi.block<3,3>(0,0) <<
      0, 0, -m_b*lv_*thetadot*cos(theta),
      0, 0, -m_b*lv_*thetadot*sin(theta),
      -m_b*lv_*thetadot*cos(theta), -m_b*lv_*thetadot*sin(theta), 0;

    Eigen::Matrix<double,6,1> Kxi = Eigen::Matrix<double,6,1>::Zero();
    Kxi(0) = (m_b+2*m_w_r)*g_*sin(rho_);

    Eigen::Matrix<double,6,4> Axi = Eigen::Matrix<double,6,4>::Zero();
    Axi.block<3,4>(0,0) = A3;

    Eigen::Matrix<double,6,1> Qfull = Eigen::Matrix<double,6,1>::Zero();
    Qfull(3) = Q_phi;    Qfull(4) = Q_psi_f;    Qfull(5) = Q_psi_r;

    Eigen::Matrix<double,6,1> rhs6 = Qfull + Axi*lambda - Cxi*qdot - Kxi;
    // 逆行列で q̈ を求める
    qdd = Mxi.inverse() * rhs6;

    //積分用の配列に代入
    x_dd[0] = 0.0;
    x_dd[1] = qdd(0);
    x_dd[2] = qdd(1);
    x_dd[3] = qdd(2);
    x_dd[4] = qdd(3);
}




// void DynamicsIntegrator::step(const Eigen::Vector3d& q,
//                               const Eigen::Vector3d& qdot,
//                               double& phi,
//                               double& phidot,
//                               double u1,
//                               double u2)
// {
//     // 状態展開
//     double x        = q(0);
//     double y        = q(1);
//     double theta    = q(2);
//     double xdot     = qdot(0);
//     double ydot     = qdot(1);
//     double thetadot = qdot(2);

//     //PID制御でtauを計算
//     double u1_act = xdot * cos(theta) + ydot * sin(theta);
//     double tau1   = drive_pid_.compute(u1, u1_act);
//     Tau1 = tau1;
//     double tau2   = steer_pid_.compute(u2, phidot);
//     Tau2 = tau2;
//     //駆動力
//     Eigen::Vector3d Qc;
//     Qc << tau1 * std::cos(theta),
//           tau1 * std::sin(theta),
//           0.0;
//     double Qphi = tau2;

//     //各行列を定義
//     // 質量行列(3x3)
//     Eigen::Matrix3d Mxi;
//     Mxi <<  m_b, 0.0, -(lv_ *m_b * sin(theta))/2.0,
//             0.0, m_b, (lv_ *m_b * cos(theta))/2.0,
//             -(lv_ *m_b * sin(theta))/2.0, (lv_ *m_b * cos(theta))/2.0, (2*I_theta_ + m_b*((pow(lv_,2)*pow(cos(theta),2))/2.0 + (pow(lv_,2)*pow(sin(theta),2))/2.0))/2.0;

//     // コリオリ行列(3x3)
//     Eigen::Matrix3d Cxi;
//     Cxi <<  0.0, 0.0, -(lv_*m_b*cos(theta)*thetadot),
//             0.0, 0.0, -(lv_*m_b*sin(theta)*thetadot),
//             0.0, 0.0, 0.0;

//     // 重力ベクトル(3x1)
//     Eigen::Vector3d Kxi;
//     Kxi << GRAV*m_b*sin(rho_), 0.0, -(GRAV*lv_*m_b*sin(rho_)*sin(theta))/2.0;

//     // 拘束行列(3x2)
//     Eigen::Matrix<double,3,2> Axi;
//     Axi <<  sin(theta + phi), sin(theta),
//             -cos(theta + phi), -cos(theta),
//             -lv_ * cos(phi), 0.0;

//     // ヤコビ行列(2x3)
//     Eigen::Matrix<double,2,3> J;
//     J <<  sin(theta + phi), -cos(theta + phi), -lv_ * cos(phi),
//             sin(theta), -cos(theta), 0.0;

//     // 拘束付き運動方程式行列 H (5x5)
//     Eigen::Matrix<double,5,5> H;
//     H.setZero();
//     H.block<3,3>(0,0) = Mxi;
//     H.block<3,2>(0,3) = Axi;
//     H.block<2,3>(3,0) = J;
//     // H.block<2,2>(3,3) = Zero

//     // 右辺ベクトル b (5x1)
//     Eigen::Matrix<double,5,1> b;
//     b.setZero();
//     b.block<3,1>(0,0) = Qc - Cxi * qdot - Kxi;
//     // 非ホロを微分した値dのddq以外の項を代入
//     b(3) = -xdot*(thetadot+phidot)*cos(theta+phi) - ydot*(thetadot+phidot)*sin(theta+phi)-lv_*thetadot*phidot*sin(phi);
//     b(4) = -xdot*thetadot*cos(theta) + ydot*thetadot*sin(theta);

//     // Hの疑似逆行列を計算
//     Eigen::MatrixXd Hxi_pinv = (H.transpose() * H).inverse()* H.transpose();
//     Eigen::VectorXd sol= Hxi_pinv * b;
    
//     Eigen::Vector3d qdd = sol.block<3,1>(0,0);
//     Eigen::Vector2d lambda = sol.block<2,1>(3,0); 
//     double phidd = Qphi;

//     //積分用の配列に代入
//     x_dd[0] = 0.0;
//     x_dd[1] = qdd(0);
//     x_dd[2] = qdd(1);
//     x_dd[3] = qdd(2);
//     x_dd[4] = phidd;
// }