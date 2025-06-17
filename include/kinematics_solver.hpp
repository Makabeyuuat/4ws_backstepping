#pragma once

#include <vector>
#include <functional>

class KinematicsSolver {
public:
    // メソッド宣言
    double calc_alpha_1_1_();    double calc_alpha_1_2_();
    double calc_alpha_2_1_();    double calc_alpha_2_2_();

    // ALPHA2の偏微分
    double calc_pd_alpha2_pd_X_1_1_(); double calc_pd_alpha2_pd_X_1_2_(); double calc_pd_alpha2_pd_X_1_3_(); double calc_pd_alpha2_pd_X_1_4_();
    double calc_pd_alpha2_pd_X_2_1_(); double calc_pd_alpha2_pd_X_2_2_(); double calc_pd_alpha2_pd_X_2_3_(); double calc_pd_alpha2_pd_X_2_4_();


    // 一般化座標の目標加速度
    double calc_aqd_1_(); double calc_aqd_2_(); double calc_aqd_3_();
    double calc_aqd_4_(); double calc_aqd_5_(); double calc_aqd_6_();

    // Axiの計算
    double calc_Axi_1_1_(); double calc_Axi_1_2_(); double calc_Axi_1_3_(); double calc_Axi_1_4_(); double calc_Axi_1_5_(); double calc_Axi_1_6_();
    double calc_Axi_2_1_(); double calc_Axi_2_2_(); double calc_Axi_2_3_(); double calc_Axi_2_4_(); double calc_Axi_2_5_(); double calc_Axi_2_6_();
    double calc_Axi_3_1_(); double calc_Axi_3_2_(); double calc_Axi_3_3_(); double calc_Axi_3_4_(); double calc_Axi_3_5_(); double calc_Axi_3_6_();
    double calc_Axi_4_1_(); double calc_Axi_4_2_(); double calc_Axi_4_3_(); double calc_Axi_4_4_(); double calc_Axi_4_5_(); double calc_Axi_4_6_();

    // Cxiの計算
    double calc_Cxi_1_1_(); double calc_Cxi_1_2_(); double calc_Cxi_1_3_(); double calc_Cxi_1_4_(); double calc_Cxi_1_5_(); double calc_Cxi_1_6_();
    double calc_Cxi_2_1_(); double calc_Cxi_2_2_(); double calc_Cxi_2_3_(); double calc_Cxi_2_4_(); double calc_Cxi_2_5_(); double calc_Cxi_2_6_();
    double calc_Cxi_3_1_(); double calc_Cxi_3_2_(); double calc_Cxi_3_3_(); double calc_Cxi_3_4_(); double calc_Cxi_3_5_(); double calc_Cxi_3_6_();
    double calc_Cxi_4_1_(); double calc_Cxi_4_2_(); double calc_Cxi_4_3_(); double calc_Cxi_4_4_(); double calc_Cxi_4_5_(); double calc_Cxi_4_6_();
    double calc_Cxi_5_1_(); double calc_Cxi_5_2_(); double calc_Cxi_5_3_(); double calc_Cxi_5_4_(); double calc_Cxi_5_5_(); double calc_Cxi_5_6_();
    double calc_Cxi_6_1_(); double calc_Cxi_6_2_(); double calc_Cxi_6_3_(); double calc_Cxi_6_4_(); double calc_Cxi_6_5_(); double calc_Cxi_6_6_();

    // Mxiの計算
    double calc_Mxi_1_1_(); double calc_Mxi_1_2_(); double calc_Mxi_1_3_(); double calc_Mxi_1_4_(); double calc_Mxi_1_5_(); double calc_Mxi_1_6_();
    double calc_Mxi_2_1_(); double calc_Mxi_2_2_(); double calc_Mxi_2_3_(); double calc_Mxi_2_4_(); double calc_Mxi_2_5_(); double calc_Mxi_2_6_();
    double calc_Mxi_3_1_(); double calc_Mxi_3_2_(); double calc_Mxi_3_3_(); double calc_Mxi_3_4_(); double calc_Mxi_3_5_(); double calc_Mxi_3_6_();
    double calc_Mxi_4_1_(); double calc_Mxi_4_2_(); double calc_Mxi_4_3_(); double calc_Mxi_4_4_(); double calc_Mxi_4_5_(); double calc_Mxi_4_6_();
    double calc_Mxi_5_1_(); double calc_Mxi_5_2_(); double calc_Mxi_5_3_(); double calc_Mxi_5_4_(); double calc_Mxi_5_5_(); double calc_Mxi_5_6_();
    double calc_Mxi_6_1_(); double calc_Mxi_6_2_(); double calc_Mxi_6_3_(); double calc_Mxi_6_4_(); double calc_Mxi_6_5_(); double calc_Mxi_6_6_();

    // Kxiの計算
    double calc_Kxi_1_(); double calc_Kxi_2_(); double calc_Kxi_3_(); double calc_Kxi_4_(); double calc_Kxi_5_(); double calc_Kxi_6_();

    
    // G11の偏微分
    double calc_pd_G11_pd_X_1_(); double calc_pd_G11_pd_X_2_(); double calc_pd_G11_pd_X_3_(); double calc_pd_G11_pd_X_4_();

    // chainedform Wの時間微分
    double calc_pd_W_pd_t_1_(); double calc_pd_W_pd_t_2_();

    // chainedform Wの偏微分
    double calc_pd_W_pd_X_1_1_(); double calc_pd_W_pd_X_1_2_(); double calc_pd_W_pd_X_1_3_(); double calc_pd_W_pd_X_1_4_();
    double calc_pd_W_pd_X_2_1_(); double calc_pd_W_pd_X_2_2_(); double calc_pd_W_pd_X_2_3_(); double calc_pd_W_pd_X_2_4_();

    // chainedform Z
    // double calc_Z_1_1_(); double calc_Z_1_2_(); double calc_Z_1_3_(); double calc_Z_2_1_(); double calc_Z_2_2_(); double calc_Z_2_3_();

    // chainedform Zの偏微分
    double calc_pd_Z2_pd_X_1_1_(); double calc_pd_Z2_pd_X_1_2_(); double calc_pd_Z2_pd_X_1_3_(); double calc_pd_Z2_pd_X_1_4_();
    double calc_pd_Z2_pd_X_2_1_(); double calc_pd_Z2_pd_X_2_2_(); double calc_pd_Z2_pd_X_2_3_(); double calc_pd_Z2_pd_X_2_4_();
    double calc_pd_Z2_pd_X_3_1_(); double calc_pd_Z2_pd_X_3_2_(); double calc_pd_Z2_pd_X_3_3_(); double calc_pd_Z2_pd_X_3_4_();

    //状態変数ベクトル
    double calc_SX_1_1_(); double calc_SX_1_2_(); double calc_SX_2_1_(); double calc_SX_2_2_();
    double calc_SX_3_1_(); double calc_SX_3_2_(); double calc_SX_4_1_(); double calc_SX_4_2_();
    // 状態変数ベクトルの時間微分
    double calc_d_SX_d_t_1_1_(); double calc_d_SX_d_t_1_2_(); double calc_d_SX_d_t_2_1_(); double calc_d_SX_d_t_2_2_();
    double calc_d_SX_d_t_3_1_(); double calc_d_SX_d_t_3_2_(); double calc_d_SX_d_t_4_1_(); double calc_d_SX_d_t_4_2_();

    //目標速度の時間微分
    double calc_pd_ud_pd_t_1_();
    double calc_pd_ud_pd_t_2_();
    // ------------------------------------------------------
    // 関数ポインタ（std::function）グループ
    using Func = std::function<double()>;
    std::vector<Func> alpha_funcs;
    std::vector<Func> pd_alpha2_funcs;
    std::vector<Func> aqd_funcs;
    std::vector<Func> Axi_funcs;
    std::vector<Func> Cxi_funcs;
    std::vector<Func> Mxi_funcs;
    std::vector<Func> Kxi_funcs;
    std::vector<Func> pd_G11_funcs;
    std::vector<Func> pd_Wt_funcs;
    std::vector<Func> pd_WX_funcs;
    std::vector<Func> Z_funcs;
    std::vector<Func> pd_Z2_funcs;
    std::vector<Func> SX_funcs;
    std::vector<Func> dSXdt_funcs;
    std::vector<Func> pdud_funcs;

    // コンストラクタで各ベクターを初期化
    KinematicsSolver() {
        alpha_funcs = {
            [this](){ return calc_alpha_1_1_(); },
            [this](){ return calc_alpha_1_2_(); },
            [this](){ return calc_alpha_2_1_(); },
            [this](){ return calc_alpha_2_2_(); }
        };
        pd_alpha2_funcs = {
            [this](){ return calc_pd_alpha2_pd_X_1_1_(); },
            [this](){ return calc_pd_alpha2_pd_X_1_2_(); },
            [this](){ return calc_pd_alpha2_pd_X_1_3_(); },
            [this](){ return calc_pd_alpha2_pd_X_1_4_(); },
            [this](){ return calc_pd_alpha2_pd_X_2_1_(); },
            [this](){ return calc_pd_alpha2_pd_X_2_2_(); },
            [this](){ return calc_pd_alpha2_pd_X_2_3_(); },
            [this](){ return calc_pd_alpha2_pd_X_2_4_(); }
        };
        aqd_funcs = {
            [this](){ return calc_aqd_1_(); },
            [this](){ return calc_aqd_2_(); },
            [this](){ return calc_aqd_3_(); },
            [this](){ return calc_aqd_4_(); },
            [this](){ return calc_aqd_5_(); },
            [this](){ return calc_aqd_6_(); }
        };
        Axi_funcs = {
            [this](){ return calc_Axi_1_1_(); }, [this](){ return calc_Axi_1_2_(); }, [this](){ return calc_Axi_1_3_(); }, [this](){ return calc_Axi_1_4_(); }, [this](){ return calc_Axi_1_5_(); }, [this](){ return calc_Axi_1_6_(); },
            [this](){ return calc_Axi_2_1_(); }, [this](){ return calc_Axi_2_2_(); }, [this](){ return calc_Axi_2_3_(); }, [this](){ return calc_Axi_2_4_(); }, [this](){ return calc_Axi_2_5_(); }, [this](){ return calc_Axi_2_6_(); },
            [this](){ return calc_Axi_3_1_(); }, [this](){ return calc_Axi_3_2_(); }, [this](){ return calc_Axi_3_3_(); }, [this](){ return calc_Axi_3_4_(); }, [this](){ return calc_Axi_3_5_(); }, [this](){ return calc_Axi_3_6_(); },
            [this](){ return calc_Axi_4_1_(); }, [this](){ return calc_Axi_4_2_(); }, [this](){ return calc_Axi_4_3_(); }, [this](){ return calc_Axi_4_4_(); }, [this](){ return calc_Axi_4_5_(); }, [this](){ return calc_Axi_4_6_(); }
        };
        Cxi_funcs = {
            [this](){ return calc_Cxi_1_1_(); }, [this](){ return calc_Cxi_1_2_(); }, [this](){ return calc_Cxi_1_3_(); }, [this](){ return calc_Cxi_1_4_(); }, [this](){ return calc_Cxi_1_5_(); }, [this](){ return calc_Cxi_1_6_(); },
            [this](){ return calc_Cxi_2_1_(); }, [this](){ return calc_Cxi_2_2_(); }, [this](){ return calc_Cxi_2_3_(); }, [this](){ return calc_Cxi_2_4_(); }, [this](){ return calc_Cxi_2_5_(); }, [this](){ return calc_Cxi_2_6_(); },
            [this](){ return calc_Cxi_3_1_(); }, [this](){ return calc_Cxi_3_2_(); }, [this](){ return calc_Cxi_3_3_(); }, [this](){ return calc_Cxi_3_4_(); }, [this](){ return calc_Cxi_3_5_(); }, [this](){ return calc_Cxi_3_6_(); },
            [this](){ return calc_Cxi_4_1_(); }, [this](){ return calc_Cxi_4_2_(); }, [this](){ return calc_Cxi_4_3_(); }, [this](){ return calc_Cxi_4_4_(); }, [this](){ return calc_Cxi_4_5_(); }, [this](){ return calc_Cxi_4_6_(); },
            [this](){ return calc_Cxi_5_1_(); }, [this](){ return calc_Cxi_5_2_(); }, [this](){ return calc_Cxi_5_3_(); }, [this](){ return calc_Cxi_5_4_(); }, [this](){ return calc_Cxi_5_5_(); }, [this](){ return calc_Cxi_5_6_(); },
            [this](){ return calc_Cxi_6_1_(); }, [this](){ return calc_Cxi_6_2_(); }, [this](){ return calc_Cxi_6_3_(); }, [this](){ return calc_Cxi_6_4_(); }, [this](){ return calc_Cxi_6_5_(); }, [this](){ return calc_Cxi_6_6_(); }
        };
        Mxi_funcs = {
            [this](){ return calc_Mxi_1_1_(); }, [this](){ return calc_Mxi_1_2_(); }, [this](){ return calc_Mxi_1_3_(); }, [this](){ return calc_Mxi_1_4_(); }, [this](){ return calc_Mxi_1_5_(); }, [this](){ return calc_Mxi_1_6_(); },
            [this](){ return calc_Mxi_2_1_(); }, [this](){ return calc_Mxi_2_2_(); }, [this](){ return calc_Mxi_2_3_(); }, [this](){ return calc_Mxi_2_4_(); }, [this](){ return calc_Mxi_2_5_(); }, [this](){ return calc_Mxi_2_6_(); },
            [this](){ return calc_Mxi_3_1_(); }, [this](){ return calc_Mxi_3_2_(); }, [this](){ return calc_Mxi_3_3_(); }, [this](){ return calc_Mxi_3_4_(); }, [this](){ return calc_Mxi_3_5_(); }, [this](){ return calc_Mxi_3_6_(); },
            [this](){ return calc_Mxi_4_1_(); }, [this](){ return calc_Mxi_4_2_(); }, [this](){ return calc_Mxi_4_3_(); }, [this](){ return calc_Mxi_4_4_(); }, [this](){ return calc_Mxi_4_5_(); }, [this](){ return calc_Mxi_4_6_(); },
            [this](){ return calc_Mxi_5_1_(); }, [this](){ return calc_Mxi_5_2_(); }, [this](){ return calc_Mxi_5_3_(); }, [this](){ return calc_Mxi_5_4_(); }, [this](){ return calc_Mxi_5_5_(); }, [this](){ return calc_Mxi_5_6_(); },
            [this](){ return calc_Mxi_6_1_(); }, [this](){ return calc_Mxi_6_2_(); }, [this](){ return calc_Mxi_6_3_(); }, [this](){ return calc_Mxi_6_4_(); }, [this](){ return calc_Mxi_6_5_(); }, [this](){ return calc_Mxi_6_6_(); }
        };
        Kxi_funcs = {
            [this](){ return calc_Kxi_1_(); },
            [this](){ return calc_Kxi_2_(); },
            [this](){ return calc_Kxi_3_(); },
            [this](){ return calc_Kxi_4_(); },
            [this](){ return calc_Kxi_5_(); },
            [this](){ return calc_Kxi_6_(); }
        };
        pd_G11_funcs = {
            [this](){ return calc_pd_G11_pd_X_1_(); },
            [this](){ return calc_pd_G11_pd_X_2_(); },
            [this](){ return calc_pd_G11_pd_X_3_(); },
            [this](){ return calc_pd_G11_pd_X_4_(); }
        };
        pd_Wt_funcs = {
            [this](){ return calc_pd_W_pd_t_1_(); },
            [this](){ return calc_pd_W_pd_t_2_(); }
        };
        pd_WX_funcs = {
            [this](){ return calc_pd_W_pd_X_1_1_(); },
            [this](){ return calc_pd_W_pd_X_1_2_(); },
            [this](){ return calc_pd_W_pd_X_1_3_(); },
            [this](){ return calc_pd_W_pd_X_1_4_(); },
            [this](){ return calc_pd_W_pd_X_2_1_(); },
            [this](){ return calc_pd_W_pd_X_2_2_(); },
            [this](){ return calc_pd_W_pd_X_2_3_(); },
            [this](){ return calc_pd_W_pd_X_2_4_(); }
        };
        // Z_funcs = {
        //     [this](){ return calc_Z_1_1_(); },
        //     [this](){ return calc_Z_1_2_(); },
        //     [this](){ return calc_Z_1_3_(); },
        //     [this](){ return calc_Z_2_1_(); },
        //     [this](){ return calc_Z_2_2_(); },
        //     [this](){ return calc_Z_2_3_(); }double calc_Z_1_1_(); double calc_Z_1_2_(); double calc_Z_1_3_(); double calc_Z_2_1_(); double calc_Z_2_2_(); double calc_Z_2_3_();

        // };
        pd_Z2_funcs = {
            [this](){ return calc_pd_Z2_pd_X_1_1_(); },
            [this](){ return calc_pd_Z2_pd_X_1_2_(); },
            [this](){ return calc_pd_Z2_pd_X_1_3_(); },
            [this](){ return calc_pd_Z2_pd_X_1_4_(); },
            [this](){ return calc_pd_Z2_pd_X_2_1_(); },
            [this](){ return calc_pd_Z2_pd_X_2_2_(); },
            [this](){ return calc_pd_Z2_pd_X_2_3_(); },
            [this](){ return calc_pd_Z2_pd_X_2_4_(); },
            [this](){ return calc_pd_Z2_pd_X_3_1_(); },
            [this](){ return calc_pd_Z2_pd_X_3_2_(); },
            [this](){ return calc_pd_Z2_pd_X_3_3_(); },
            [this](){ return calc_pd_Z2_pd_X_3_4_(); }
        };
        SX_funcs = {
            [this](){ return calc_SX_1_1_(); },
            [this](){ return calc_SX_1_2_(); },
            [this](){ return calc_SX_2_1_(); },
            [this](){ return calc_SX_2_2_(); },
            [this](){ return calc_SX_3_1_(); },
            [this](){ return calc_SX_3_2_(); },
            [this](){ return calc_SX_4_1_(); },
            [this](){ return calc_SX_4_2_(); }
        };

        dSXdt_funcs = {
            [this](){ return calc_d_SX_d_t_1_1_(); },
            [this](){ return calc_d_SX_d_t_1_2_(); },
            [this](){ return calc_d_SX_d_t_2_1_(); },
            [this](){ return calc_d_SX_d_t_2_2_(); },
            [this](){ return calc_d_SX_d_t_3_1_(); },
            [this](){ return calc_d_SX_d_t_3_2_(); },
            [this](){ return calc_d_SX_d_t_4_1_(); },
            [this](){ return calc_d_SX_d_t_4_2_(); }
        };
        pdud_funcs = {
            [this](){ return calc_pd_ud_pd_t_1_(); },
            [this](){ return calc_pd_ud_pd_t_2_(); }
        };
        
    }
};