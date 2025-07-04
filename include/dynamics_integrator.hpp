#pragma once

#include "kinematics_solver.hpp"
#include <Eigen/Dense>
#include "pidTau.hpp"
#include <vector>
#include "getInputValue_dynamics.hpp"

/**
 * @brief DynamicsIntegrator: 制約付き動力学＋操舵ダイナミクス統合クラス
 */
class DynamicsIntegrator {
public:

    //新規メソッド: バックステッピングから得られる目標加速度を計算する
    Eigen::Matrix<double,7,1> computeAlpha(
        const Eigen::Matrix<double,7,1>& q, // 状態ベクトル (x,y,θ,φ,ψ_f,ψ_r)
        const Eigen::Matrix<double,7,1>& qdot,     // 速度ベクトル         
        double u1,                                   // 前輪/後輪目標前進速度
        double u2);

    Eigen::Matrix<double,5,1> computeXAlpha(
        std::vector<double> x_d,
        std::vector<double> x_dd,
        double u1,
        double u2,
        double u3);

     Eigen::Matrix<double,7,1> computeAlpha(const Eigen::Matrix<double,7,1>& q,
                                            const Eigen::Matrix<double,7,1>& qdot,
                                            double u1,
                                            double u2,
                                            const Eigen::Matrix<double,7,1>& qddot);

    DynamicsIntegrator(double m_b,
                       double I_theta,
                       double lv,
                       double g,
                       double rho,
                       const PIDGains& drive_gains,
                       const PIDGains& steer_gains,
                       double dt,
                       getInputValue& inputValue);


    void step(const Eigen::Matrix<double,7,1>& q,
              const Eigen::Matrix<double,7,1>& qdot,
              double u1,
              double u2);

private:
        KinematicsSolver kinematics_solver_;
        getInputValue& inputValue_ref_;
};