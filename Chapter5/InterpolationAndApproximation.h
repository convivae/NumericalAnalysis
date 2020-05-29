//
// 数值分析
// 第五章 插值与逼近
// Created by convi on 2020/5/12.
//

#ifndef NUMERICALANALYSIS_INTERPOLATIONANDAPPROXIMATION_H
#define NUMERICALANALYSIS_INTERPOLATIONANDAPPROXIMATION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

#include "../Type/Type.h"
#include "../lib/Polymerization.h"

namespace convivae {
    class InterpolationAndApproximation {
    private:
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
    };
}


#endif //NUMERICALANALYSIS_INTERPOLATIONANDAPPROXIMATION_H
