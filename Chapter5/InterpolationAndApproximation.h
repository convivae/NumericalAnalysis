//
// 数值分析
// 第五章 插值与逼近
// Created by convi on 2020/5/12.
//

#ifndef NUMERICALANALYSIS_INTERPOLATIONANDAPPROXIMATION_H
#define NUMERICALANALYSIS_INTERPOLATIONANDAPPROXIMATION_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

#include "../Type/Type.h"
#include "../lib/Polymerization.h"

namespace convivae {
    class InterpolationAndApproximation {
    private:
        typedef std::vector<std::vector<double>> mat_type;
        typedef std::vector<double> vec_type;

        i4 _n{};
        bool _initialized;
        std::string _rdr;
        vector<f8> X, Y;

        /**
         * 计算差商表 f[x0,x1,...,xk]
         *
         * f[x0,x1] = (f(x1) - f(x0)) / (x1-x0)
         * f[x0,x1,...xk] = (f[x0,x1,...x_(k-2),xk]-f[x0,x1,...x_(k-2),xk])/(xk - x_(k-1))
         *
         * @param start_index 开始下标
         * @param before_end_index 倒数第二个下标，start_index == end_index, before_index = -1
         * @param end_index 最后一个下标
         * @param table 差商表（只记录 f(x0), f(x0,x1), f(x0,x1,...xk)）
         * @return
         */
        double
        differential_quotient(int start_index, int before_end_index, int end_index, vector<vector<double>> &table,
                              vector<vector<bool>> &flag);

        static void draw_dividing_line();

        static void print_vec(const string &s, vector<double> v, int start_pos = 0);

        static void print_mat(const string &s, vector<vector<double>> v);

        void read_data();

    public:
        explicit InterpolationAndApproximation(std::string filename)
                : _rdr(std::move(filename)), _initialized(false) {}

        /**
         * Lagrange 插值法
         * @param n 次数
         * @return Lagrange 插值多项式
         */
        Polymerization Lagrange_polynomial(int n);

        /**
         * Newton 插值法
         * @param n 次数
         * @return Newton 插值多项式
         */
        Polymerization Newton_polynomial(int n);

        /**
         * 正交多项式位于 lib/OrthogonalPolynomial.h
         */

        /**
         * 三次样条函数的三弯矩法
         * 第一种边界条件（给定两边界节点的二阶导数）
         * @param f2_0 f''(x_0)
         * @param f2_n f''(x_n)
         * @param simplify 是否显示简化的多项式
         * @return
         */
        void three_moment_method_bound1(double f2_0, double f2_n, bool simplify = false);

        /**
         * 三次样条函数的三弯矩法
         * 第二种边界条件（给定两边界节点的一阶导数）
         * @param f2_0 f'(x_0)
         * @param f2_n f'(x_n)
         * @param simplify 是否显示简化的多项式
         * @return
         */
        void three_moment_method_bound2(double f1_0, double f1_n, bool simplify = false);

        /**
         * 三次样条函数的三弯矩法
         * 第三种边界条件（s'(x_0+) = s'(x_n-), s''(x_0+) = s''(x_n-)）
         * @param simplify 是否显示简化的多项式
         * @return
         */
        void three_moment_method_bound3(bool simplify = false);
    };
}


#endif //NUMERICALANALYSIS_INTERPOLATIONANDAPPROXIMATION_H
