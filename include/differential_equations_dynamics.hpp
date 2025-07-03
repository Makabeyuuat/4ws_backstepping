#ifndef DIFFERENTIAL_EQUATIONS_HPP
#define DIFFERENTIAL_EQUATIONS_HPP

#include <vector>
#include <cmath>
#include <array>

// 状態ベクトルの次元
extern const int DIM;

// 各状態の微分方程式のシグネチャ
double f0 (const std::vector<double>& x);
double f1 (const std::vector<double>& x);
double f2 (const std::vector<double>& x);
double f3 (const std::vector<double>& x);
double f4 (const std::vector<double>& x);
double f5 (const std::vector<double>& x);
double f6 (const std::vector<double>& x);
double f7 (const std::vector<double>& x);
double f8 (const std::vector<double>& x);
double f9 (const std::vector<double>& x);
double f10 (const std::vector<double>& x);
double f11 (const std::vector<double>& x);

// 関数ポインタ型エイリアス（getInputValue からも使えます）
using FunctionPtr = double(*)(const std::vector<double>&);
extern const std::array<FunctionPtr, /*DIM+1=*/6> fAll;
extern const std::array<FunctionPtr, /*DIM+1=*/6> fdAll;

#endif // DIFFERENTIAL_EQUATIONS_HPP
