#include "getInputValue_dynamics.hpp"
#include "initial.hpp"
#include "differential_equations_dynamics.hpp"
#include "mathFunc.h"
#include "kinematics_solver.hpp"
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/SVD> 
#include <Eigen/Core>   
#include <iostream>

using namespace std;

// コンストラクタ
getInputValue::getInputValue(double h)
    : h(h),
    k(Dim+1, std::vector<double>(4,0.0)),
    r(Dim+1, std::vector<double>(4,0.0)),
    q(Dim+1, std::vector<double>(4,0.0)),
    x(3,   std::vector<double>(Dim+1,0.0)),
    rearOmega{0.0, 0.0},
    rearTorque{0.0, 0.0},
    frontTorque{0.0, 0.0},
    fAllVec(fAll.begin(), fAll.end()),
    fdAllVec(fdAll.begin(), fdAll.end()),
    kinematics_solver_()
   
{}

array<double,2> getInputValue::computeRearWheelOmegas(double speed, double steeringAngle) {
    const double W = 0.05;            // トレッド幅[m]
    array<double,2> omegas;

    if (fabs(steeringAngle) < 1e-6) {
        double omega = speed / wheelRadius;
        omegas[0] = omega;
        omegas[1] = omega;
        return omegas;
    }
    double absPhi = fabs(steeringAngle);
    double R = lv / tan(absPhi);
    double R_in  = R - W/2.0;
    double R_out = R + W/2.0;
    double v_in  = speed * (R_in  / R);
    double v_out = speed * (R_out / R);
    double omega_in  = v_in  / wheelRadius;
    double omega_out = v_out / wheelRadius;

    if (steeringAngle > 0) {
        // 左折: 左が内輪
        omegas[0] = omega_in;
        omegas[1] = omega_out;
    } else {
        // 右折: 右が内輪
        omegas[0] = omega_out;
        omegas[1] = omega_in;
    }
    return omegas;
}



// 入力がアクスルごとのトルク Q の場合
std::array<double,2> getInputValue::computeFrontWheelTorque(
    double Qf,
    double steeringAngleFront,
    double steeringAngleRear)
{
    const double Wf = 0.1;  // 前輪トレッド幅 [m]
    std::array<double,2> torques;

    double tan_diff = std::tan(steeringAngleFront) - std::tan(steeringAngleRear);
    if (std::fabs(tan_diff) < 1e-9) {
        torques[0] = Qf * 0.5;
        torques[1] = Qf * 0.5;
        return torques;
    }

    // 1) リアアクスル基準での回転中心半径
    double R_rear_center = lv / tan_diff;

    // 2) 前輪アクスルまで平行移動
    double Rf_center = R_rear_center - lv;

    // 3) 内輪／外輪の絶対半径
    double Rf_abs   = std::abs(Rf_center);
    double Rf_inner = Rf_abs - Wf/2.0;
    double Rf_outer = Rf_abs + Wf/2.0;

    // 4) 内外で逆比（パワー均等）にトルクを配分
    double sum = Rf_inner + Rf_outer;
    double Tin = Qf * (Rf_outer / sum);
    double Tou = Qf * (Rf_inner / sum);

    // 5) 旋回方向に応じて左右に割り当て
    if (tan_diff > 0) {
        // 左折：左が内輪
        torques[0] = Tin;  // 左
        torques[1] = Tou;  // 右
    } else {
        // 右折：右が内輪
        torques[0] = Tou;  // 左
        torques[1] = Tin;  // 右
    }
    return torques;
}

std::array<double,2> getInputValue::computeRearWheelTorque(
    double Qr,
    double steeringAngleFront,
    double steeringAngleRear)
{
    const double Wr = 0.1;  // 後輪トレッド幅 [m]
    std::array<double,2> torques;

    double tan_diff = std::tan(steeringAngleFront) - std::tan(steeringAngleRear);
    if (std::fabs(tan_diff) < 1e-9) {
        torques[0] = Qr * 0.5;
        torques[1] = Qr * 0.5;
        return torques;
    }

    // 1) リアアクスル基準での回転中心半径
    double Rr_center = lv / tan_diff;

    // 2) 内輪／外輪の絶対半径
    double Rr_abs   = std::abs(Rr_center);
    double Rr_inner = Rr_abs - Wr/2.0;
    double Rr_outer = Rr_abs + Wr/2.0;

    // 3) パワー均等配分
    double sum = Rr_inner + Rr_outer;
    double Tin = Qr * (Rr_outer / sum);
    double Tou = Qr * (Rr_inner / sum);

    // 4) 左右アサイン
    if (tan_diff > 0) {
        // 左折
        torques[0] = Tin;
        torques[1] = Tou;
    } else {
        // 右折
        torques[0] = Tou;
        torques[1] = Tin;
    }
    return torques;
}

void getInputValue::rungeKutta(std::vector<double>& x_old, std::vector<double>& x_d) {
    int n = static_cast<int>(x_old.size());  // n = Dim+1
    //qdを積分
    // 第1段階
        for (int i = 0; i < n; i++) {
            double dx = x_d[i];
            k[i][0] = h * dx;
            r[i][0] = (k[i][0] - 2.0 * q[i][3]) / 2.0;
            x[0][i] = x_old[i] + r[i][0];
            q[i][0] = q[i][3] + 3.0 * r[i][0] - k[i][0] / 2.0;
        }   
        // 第2段階
        for (int i = 0; i < n; i++) {
            double dx = fAllVec[i](x[0]);
            k[i][1] = h * dx;
            r[i][1] = (1.0 - sqrt(0.5)) * (k[i][1] - q[i][0]);
            x[1][i] = x[0][i] + r[i][1];
            q[i][1] = q[i][0] + 3.0 * r[i][1] - (1.0 - sqrt(0.5)) * k[i][1];
        }
        // 第3段階
        for (int i = 0; i < n; i++) {
            double dx = fAllVec[i](x[0]);
            k[i][2] = h * dx;
            r[i][2] = (1.0 + sqrt(0.5)) * (k[i][2] - q[i][1]);
            x[2][i] = x[1][i] + r[i][2];
            q[i][2] = q[i][1] + 3.0 * r[i][2] - (1.0 + sqrt(0.5)) * k[i][2];
        }
        // 第4段階
        for (int i = 0; i < n; i++) {
            double dx = fAllVec[i](x[0]);
            k[i][3] = h * dx;
            r[i][3] = (k[i][3] - 2.0 * q[i][2]) / 6.0;
            x_new[i] = x[2][i] + r[i][3];
            q[i][3] = q[i][2] + 3.0 * r[i][3] - k[i][3] / 2.0;

            
        }
        // ... 更新結果を x_new に格納
        for (int i = 1; i < n; i++) {
            x_old[i] = x_new[i];
        }
        x_old[0] = x_old[0] + h;  // 時間の更新

        // dynamic_v = x_d[1]*cos(x_old[3]) + x_d[2]*sin(x_old[3]);

        // // 後輪左右の角速度をクラスメンバに格納
        // rearOmega = computeRearWheelOmegas(dynamic_v, x_old[4]);
        // omega_rear[0] = rearOmega[0];  // 左後輪
        // omega_rear[1] = rearOmega[1];  // 右後輪

         // 後輪左右のトルクをクラスメンバに格納
        // rearTorque = computeRearWheelTorque(Q_varphiR, x_old[5], x_old[4]);
        // frontTorque = computeFrontWheelTorque(Q_varphiF, x_old[5], x_old[4]);
        
        // torque_rear[0] = rearTorque[0];  // 左後輪
        // torque_rear[1] = rearTorque[1];  // 右後輪

        // torque_front[0] = frontTorque[0];  // 左後輪
        // torque_front[1] = frontTorque[1];  // 右後輪
    
}

void getInputValue::ddrungeKutta(std::vector<double>& x_d, std::vector<double>& x_dd) {
    std::vector<double> x_newd = std::vector<double>(Dim +1, 0.0);
    int n = static_cast<int>(x_dd.size()); 
    //qddを積分
    // 第1段階
        for (int i = 0; i < n; i++) {
            double ddx = x_dd[i];
            k[i][0] = h * ddx;
            r[i][0] = (k[i][0] - 2.0 * q[i][3]) / 2.0;
            x[0][i] = x_d[i] + r[i][0];
            q[i][0] = q[i][3] + 3.0 * r[i][0] - k[i][0] / 2.0;
        }   
        // 第2段階
        for (int i = 0; i < n; i++) {
            double ddx = fdAllVec[i](x[0]);
            k[i][1] = h * ddx;
            r[i][1] = (1.0 - sqrt(0.5)) * (k[i][1] - q[i][0]);
            x[1][i] = x[0][i] + r[i][1];
            q[i][1] = q[i][0] + 3.0 * r[i][1] - (1.0 - sqrt(0.5)) * k[i][1];
        }
        // 第3段階
        for (int i = 0; i < n; i++) {
            double ddx = fdAllVec[i](x[0]);
            k[i][2] = h * ddx;
            r[i][2] = (1.0 + sqrt(0.5)) * (k[i][2] - q[i][1]);
            x[2][i] = x[1][i] + r[i][2];
            q[i][2] = q[i][1] + 3.0 * r[i][2] - (1.0 + sqrt(0.5)) * k[i][2];
        }
        // 第4段階
        for (int i = 0; i < n; i++) {
            double ddx = fdAllVec[i](x[0]);
            k[i][3] = h * ddx;
            r[i][3] = (k[i][3] - 2.0 * q[i][2]) / 6.0;
            x_newd[i] = x[2][i] + r[i][3];
            q[i][3] = q[i][2] + 3.0 * r[i][3] - k[i][3] / 2.0;
            // ... 更新結果を x_new に格納 
            x_d[i] = x_newd[i];
        }
}

void getInputValue::getU(std::vector<double>& x_old, int sr_j) {
    // --- 制御入力の計算 ---
    // 各内部関数を呼び出して制御入力を計算
    thetaT = atan2(dRdq[sr_j][1], dRdq[sr_j][0]);
    Thetap = x_old[3] - thetaT;

    U1(x_old, sr_j);
    U2(x_old, sr_j);
    U3(x_old, sr_j);


}

void getInputValue::getXInput(std::vector<double>& x_old, std::vector<double>& x_input){
	
}


//制御入力計算用内部関数
void getInputValue::U1(const std::vector<double>& x_old, int sr_j) {


	u1 = ((1 - sr.d * sr.Cs) / cos(Thetap + x_old[4])) * w1;


}


void getInputValue::U2(const std::vector<double>& x_old, int sr_j) {
	z21 = kinematics_solver_.Z_funcs[0]();
	z22 = kinematics_solver_.Z_funcs[1]();
    alpha21 = kinematics_solver_.alpha_funcs[3]();
    alpha22 = kinematics_solver_.alpha_funcs[4]();

	//経路追従
	d0d = 0.0;
	dd0d = 0.0;
	ddd0d = 0.0;

	//重心を経路に対して周期的に変化させる
	/*d0d = -0.7 *  sin(2* PAI * x[0] / t_max);
	printf("%lf\n", d0d);
	dd0d =  -0.7 * (2* PAI / t_max) * cos(2 * PAI * x[0] / t_max);
	ddd0d = 0.7 * (2 * PAI / t_max) * (2 * PAI / t_max) * sin(2 * PAI * x[0] / t_max);*/

	w2 = ddd0d / a0 + (k1 + k2) * ((dd0d / a0) - z21) + k1 * k2 * ((d0d / a0) - z22 / a0);

	u2 = (1 / alpha22) * (w2 - alpha21 * w1);
}


void getInputValue::U3(const std::vector<double>& x_old, int sr_j) {
    z31 = kinematics_solver_.Z_funcs[2]();
	z32 = kinematics_solver_.Z_funcs[3]();
    alpha31 = kinematics_solver_.alpha_funcs[6]();
    alpha32 = kinematics_solver_.alpha_funcs[7]();
    alpha33 = kinematics_solver_.alpha_funcs[8]();

    thetap1d = 0.0;
	dthetap1d = 0.0;
	ddthetap1d = 0.0;

	w3 = ddthetap1d / a0 + (k3 + k4) * ((dthetap1d / a0) - z31) + k3 * k4 * ((thetap1d / a0) - z32 / a0);

    double w31 = (k3 + k4) * ((dthetap1d / a0) - z31); 
    double w32 = k3 * k4 * ((thetap1d / a0) - z32 / a0);

	u3 = (1 / alpha33) * (w3 - (alpha31 * w1 + alpha32 * u2));

    std::cout << "w31=" << w31 << ", w32=" << w32 << "\n";
    std::cout << "z21=" << z21 << ", z22=" << z22 << "\n";
    std::cout << "z31=" << z31 << ", z32=" << z32 << "\n";
    std::cout << "alpha21=" << alpha21 << ", alpha22=" << alpha22 << ", w2=" << w2 << "\n";
    std::cout << "alpha31=" << alpha31 << ", alpha32=" << alpha32 << ", alpha33=" << alpha33 << ", w3=" << w3 << "\n";
    std::cout << "u1 =" <<u1 << ", u2 =" <<u2 <<", u3 =" <<u3 << "\n\n";
}