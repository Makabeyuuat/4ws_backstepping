#pragma once

#include <Eigen/Dense>
#include "pidTau.hpp"

/**
 * @brief DynamicsIntegrator: 制約付き動力学＋操舵ダイナミクス統合クラス
 */
class DynamicsIntegrator {
public:

    //新規メソッド: バックステッピングから得られる目標加速度を計算する
    Eigen::Matrix<double,6,1> computeAlpha(
        const Eigen::Matrix<double,6,1>& q, // 状態ベクトル (x,y,θ,φ,ψ_f,ψ_r)
        const Eigen::Matrix<double,6,1>& qdot,     // 速度ベクトル         
        double u1,                                   // 前輪/後輪目標前進速度
        double u2                                    // ステア角速度指令
    );

    Eigen::Matrix<double,4,1> DynamicsIntegrator::computeXAlpha(
    std::vector<double> x_d,
    std::vector<double> x_dd,
    double u1,
    double u2)

    DynamicsIntegrator(double m_b,
                       double I_theta,
                       double lv,
                       double g,
                       double rho,
                       const PIDGains& drive_gains,
                       const PIDGains& steer_gains,
                       double dt);

    // @brief 1ステップの動力学計算 & 数値積分
    // @param q      [x, y, θ]
    // @param qdot   [ẋ, ẏ, θ̇]
    // @param phi    操舵角 φ
    // @param phi_dot 操舵角速度 φ̇
    // @param u1     目標前進速度 u1
    // @param u2     目標操舵速度 u2

    void step(const Eigen::Matrix<double,6,1>& q,
              const Eigen::Matrix<double,6,1>& qdot,
              double u1,
              double u2
              const Eigen::Matrix<double,6,1>& qddot,);

private:
    double m_b_, I_theta_, lv_, g_, rho_, dt_;
    PIDTorque drive_pid_;
    PIDTorque steer_pid_;
};