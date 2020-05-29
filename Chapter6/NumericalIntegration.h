//
// Created by convi on 2020/5/29.
//

#ifndef NUMERICALANALYSIS_NUMERICALINTEGRATION_H
#define NUMERICALANALYSIS_NUMERICALINTEGRATION_H

#include "../lib/OrthogonalPolynomial.h"
#include "../Type/Type.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdarg>

using namespace std;

class NumericalIntegration {
private:
    //定义一种 pFun 的的函数指针，这种函数以一个 double 为参数并返回 double 类型
    typedef double (*pFun)(double x);

    static void draw_dividing_line();

    /**
     * 幂
     * @param a
     * @param b
     * @return a ^ b
     */
    int my_pow(int a, int b);

    /**
     * 阶乘
     * @param n
     * @return
     */
    double my_fac(int n);

public:
    /**
     * Romberg 积分法
     * @param fun 原函数
     * @param up 积分上限
     * @param down 积分下限
     * @param epsilon 迭代精度控制
     */
    void
    Romberg_numerical_integration_method(pFun fun, double down, double up, double epsilon = 1e-5, int max_steps = 200);


    /**
     * Gauss_Legendre 积分法 [-1,1]
     * @param fun 原函数
     * @param n n 点
     * @param ... 求积节点: x_i (共有 n 个)(课本 P170)
     */
    void Gauss_Legendre_integration_method(pFun fun, int n, ...);

    /**
     * Gauss_Laguerre 积分法 [0,inf]
     * @param fun 原函数
     * @param n n 点
     * @param ... 求积节点: x_i (共有 n 个)(课本 P171)
     */
    void Gauss_Laguerre_integration_method(pFun fun, int n, ...);

    /**
     * Gauss_Hermite 积分法 [-inf,inf]
     * @param fun 原函数
     * @param n n 点
     * @param ... 求积节点: x_i (共有 n 个)(课本 P172)
     */
    void Gauss_Hermite_integration_method(pFun fun, int n, ...);

    /**
     * Gauss_Chebyshev 积分法 [-1,1]
     * @param fun 原函数
     * @param n n 点
     * @param ... 求积节点: x_i (共有 n 个)
     */
    void Gauss_Chebyshev_integration_method(pFun fun, int n, ...);
};


#endif //NUMERICALANALYSIS_NUMERICALINTEGRATION_H
