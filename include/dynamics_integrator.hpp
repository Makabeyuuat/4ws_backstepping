#pragma once

#include "kinematics_solver.hpp"
#include <Eigen/Dense>
#include "pidTau.hpp"
#include <vector>


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

    Eigen::Matrix<double,4,1> computeXAlpha(
    std::vector<double> x_d,
    std::vector<double> x_dd,
    double u1,
    double u2);

     Eigen::Matrix<double,6,1> computeAlpha(const Eigen::Matrix<double,6,1>& q,
                                            const Eigen::Matrix<double,6,1>& qdot,
                                            double u1,
                                            double u2,
                                            const Eigen::Matrix<double,6,1>& qddot);

    DynamicsIntegrator(double m_b,
                       double I_theta,
                       double lv,
                       double g,
                       double rho,
                       const PIDGains& drive_gains,
                       const PIDGains& steer_gains,
                       double dt);


    void step(const Eigen::Matrix<double,6,1>& q,
              const Eigen::Matrix<double,6,1>& qdot,
              double u1,
              double u2,
              const Eigen::Matrix<double,6,1>& qddot);

private:
        KinematicsSolver kinematics_solver_;
};