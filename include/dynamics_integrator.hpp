#pragma once

#include "kinematics_solver.hpp"
#include <Eigen/Dense>
#include "pidTau.hpp"
#include <vector>
#include "getInputValue_dynamics.hpp"

// === 追加: Ni連動補償で使うパラメータ/結果 ===
struct NiCompParams {
  // 幾何・物理
  double m;      // 総質量 [kg]
  double g;      // 重力加速度 [m/s^2]
  double L;      // ホイールベース [m]
  double t;      // トレッド [m]
  double h;      // 重心高 [m]
  double r;      // 車輪半径 [m]
  // 重心から各軸/各側までの半距離（一般式用）
  double lF;     // 前軸まで [m]（lF + lR = L）
  double lR;     // 後軸まで [m]
  double wL;     // 左側まで [m]（wL + wR = t）
  double wR;     // 右側まで [m]

  // 摩擦・補償係数
  double mu_t;     // 接地（縦）摩擦（飽和用）
  double mu_r;     // 転がり抵抗係数
  double kc_phi;   // 駆動クーロン補償係数 [m] (= 名目friction/N0)
  double kc_varphi;// ステアクーロン補償係数 [m]
  double b_phi;    // 駆動粘性 [N·m·s/rad]
  double b_varphi; // ステア粘性 [N·m·s/rad]
  double eps;      // 滑らか符号の微小量

  // 任意：実加速度（重力除去後）を既に持っていれば入れる
  double ax_dyn{0.0}; // 車体x前後[m/s^2]
  double ay_dyn{0.0}; // 車体y左右[m/s^2]
};

struct NiCompTorques {
  // この補正トルクを Q_* に「加算」してください
  double d_phiR{0}, d_phiF{0};
  double d_varphiR{0}, d_varphiF{0};

  // デバッグ/飽和用
  double N_FL{0}, N_FR{0}, N_RL{0}, N_RR{0};
  double tau_max_FL{0}, tau_max_FR{0}, tau_max_RL{0}, tau_max_RR{0};
};

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

        // === 追加: Ni連動補償の計算メソッド ===
        NiCompTorques computeCompensationTorques(
            const Eigen::Matrix<double,7,1>& q,
            const Eigen::Matrix<double,7,1>& qdot,
            double alpha,   // 斜面角（世界x軸に沿って上り）
            double psi,     // 登り方向に対する車体ヨーずれ
            const NiCompParams& P);
};