//
// 数值分析
// 第四章 非线性方程与非线性方程组的迭代解法
// Created by convi on 2020/4/30.
//

#ifndef NUMERICALANALYSIS_NONLINEAREQUATIONS_H
#define NUMERICALANALYSIS_NONLINEAREQUATIONS_H

#include "../Type/Type.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

namespace convivae {
    class NonlinearEquations {
    private:
        //定义一种 pFun 的的函数指针，这种函数以一个 double 为参数并返回 double 类型
        typedef f8 (*pFun)(f8 x);

        //定义一种 ppFun 的的函数指针，这种函数以两个 double 为参数并返回 double 类型
        typedef f8 (*ppFun)(f8 x, f8 y);

        static void draw_dividing_line();

    public:
        /**
         * 二分法（对分法）
         * 只能求单根和奇数重根
         * @param fun 迭代函数的函数指针
         * @param left 左区间
         * @param right 右区间
         * @param epsilon 迭代精度
         * @param max_steps 最大迭代次数
         */
        void dichotomy_method(pFun fun, double left, double right, double epsilon = 1e-6, bool show_details = false,
                              int max_steps = 2000) const;

        /**
         * 简单迭代法
         * @param fun 函数指针
         * @param x0 迭代初始值
         * @param eta 迭代精度
         * @param show_details
         * @param max_steps
         */
        void simple_iteration_method(pFun fun, double x0, double eta = 1e-6, bool show_details = false,
                                     int max_steps = 2000) const;


        /**
         * Steffensen 迭代法
         * @param fun 函数指针
         * @param x0 迭代初始值
         * @param eta 迭代精度
         * @param show_details
         * @param max_steps
         */
        void Steffensen_iteration_method(pFun fun, double x0, double eta = 1e-6, bool show_details = false,
                                         int max_steps = 2000) const;


        /**
         * Newton 迭代法
         * 使用了牛顿下山法来解决初始值选定后不收敛的问题，因此无需给定初始值
         * @param fun 迭代函数
         * @param fun_1 迭代函数的一阶导数
         * @param x0
         * @param eta
         * @param show_details
         * @param down_hill 是否采用下山法，若选定的初始值不收敛可以选择开启
         * @param max_steps
         */
        void Newton_iteration_method(pFun fun, pFun fun_1, double eta = 1e-6, double x0 = 0, bool show_details = false,
                                     bool down_hill = false, int max_steps = 2000) const;

        /**
         * 割线法
         * 把 Newton 法中的导数替换为增量比，收敛速度低于 Newton 法，但高于一阶
         * @param fun
         * @param x0 两个初始值
         * @param x1 两个初始值
         * @param eta
         * @param show_details
         * @param max_steps
         */
        void secant_iteration_method(pFun fun, double x0, double x1, double eta = 1e-6, bool show_details = false,
                                     int max_steps = 2000) const;

        /**
         * 单点割线法
         * 在割线法中用固定点 (x0, f(x0)) 代替 (x_(k-1), f(x_(k-1)))，即 (x0, f(x0)) 永远是割线上的一点
         * @param fun
         * @param x0 两个初始值
         * @param x1 两个初始值
         * @param eta
         * @param show_details
         * @param max_steps
         */
        void secant_with_single_point_iteration_method(pFun fun, double x0, double x1, double eta = 1e-6,
                                                       bool show_details = false, int max_steps = 2000) const;

        /**
         * 非线性方程组的简单迭代法（仅针对两个未知数）
         * 完整写法比较麻烦，这里只是针对课后练习题的算法
         * @param fun1
         * @param fun2
         * @param x1_0 x1 的初始值
         * @param x2_0 x2 的初始值
         * @param eta
         * @param show_details
         * @param max_steps
         */
        void simple_iteration_equations_method(ppFun fun1, ppFun fun2, double x1_0, double x2_0, double eta = 1e-6,
                                               bool show_details = false, int max_steps = 2000) const;

        /**
         * 非线性方程组的 Newton 迭代法（仅针对两个未知数）
         * @param fun1
         * @param fun2
         * @param x1_0
         * @param x2_0
         * @param eta
         * @param show_details
         * @param max_steps
         */
        void Newton_iteration_equations_method(ppFun fun1, ppFun fun2, double x1_0, double x2_0, double eta = 1e-6,
                                               bool show_details = false, int max_steps = 2000) const;
    };
}

#endif //NUMERICALANALYSIS_NONLINEAREQUATIONS_H
