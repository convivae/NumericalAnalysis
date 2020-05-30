//
// 只有一些函数
// Created by convi on 2020/5/12.
//

#ifndef NUMERICALANALYSIS_FUN_H
#define NUMERICALANALYSIS_FUN_H

#include "../Type/Type.h"

//二分法的迭代函数
double fun_dichotomy(double x) {
    return x + sin(x) - 1;
}

//简单迭代法的迭代函数
double fun_simple(double x) {
    return 1.0 / (1.0 + x);
}

//Steffensen 法迭代函数
double fun_Steffensen(double x) {
    return (2 - pow(M_E, x)) / 10;
}

//Newton 法迭代函数
double fun_Newton(double x) {
    return pow(1 + x * x, 1 / 3);
}

double fun_Newton_1(double x) {
    return 1 + cos(x);
}

// 割线法
double fun_secant(double x) {
    return x + sin(x) - 1;
}

//方程组的简单迭代法
double fun_equations_simple1(double x, double y) {
    return 1.2 - pow(M_E, -2 * y);
}

double fun_equations_simple2(double x, double y) {
    return 0.5 * (1.97 - pow(M_E, -x));
}

//方程组的 Newton 迭代法
double fun_equations_Newton1(double x, double y) {
    auto a = 1 - cos(x + y);
    auto b = -cos(x + y);
    auto c = -sin(x + y);
    auto d = 1 - sin(x + y);
    auto e = x - sin(x + y) - 1.2;
    auto f = y + cos(x + y) - 0.5;

    return x + (f / d - e / b) / (a / b - c / d);
}

double fun_equations_Newton2(double x, double y) {
    auto a = 1 - cos(x + y);
    auto b = -cos(x + y);
    auto c = -sin(x + y);
    auto d = 1 - sin(x + y);
    auto e = x - sin(x + y) - 1.2;
    auto f = y + cos(x + y) - 0.5;

    return y + (f / c - e / a) / (b / a - d / c);
}

// Romberg 积分法的原函数
double fun_Romberg(double x) {
    return pow(M_E, -x * x);
}

// Gauss_Legendre 积分法的原函数
double fun_Gauss_Legendre(double x) {
    return 1 / (x+2);
}

// Gauss_Lagueree 积分法的原函数
double fun_Gauss_Lagueree(double x) {
    return pow(M_E, -x) * sqrt(x);
}

// Gauss_Hermite 积分法的原函数
double fun_Gauss_Hermite(double x) {
    return pow(M_E, -x * x) * cos(x);
}

// Gauss_Chebyshev 积分法的原函数
double fun_Gauss_Chebyshev(double x) {
    return sqrt((2 + x) / (1.0 - x * x));
}

#endif //NUMERICALANALYSIS_FUN_H
