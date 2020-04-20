//
// 数值分析
// 第三章 矩阵特征值与特征向量的计算
// Created by convi on 2020/4/5.
//

#ifndef NUMERICALANALYSIS_MATRIXEIGENVALUESANDEIGENVECTORS_H
#define NUMERICALANALYSIS_MATRIXEIGENVALUESANDEIGENVECTORS_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <iomanip>

#include "../Type/Type.h"
#include "../lib/MatrixOperation.h"

namespace convivae {

    class MatrixEigenvaluesAndEigenvectors {
    private:
        typedef std::vector<std::vector<f8>> mat_type;
        typedef std::vector<f8> vec_type;

        std::string _rdr; //文件名
        mat_type _original_mat; //读入的初始矩阵
        bool _initialized; //如果未初始化，就读入矩阵
        i4 _dimension{}; //方阵的维数 n，在读入时初始化

        void read_mat();

        static void draw_dividing_line();

        static void print_vec(const std::string &info, vec_type const &a);

    public:
        explicit MatrixEigenvaluesAndEigenvectors(std::string filename)
                : _rdr(std::move(filename)), _initialized(false) {}

        /**
         * 幂法
         * 用于计算矩阵主特征值 (矩阵按模最大的特征值)及对应特征向量的迭代方法
         * @param epsilon 允许误差
         * @param norm 范数的选择，只能选二范数和无穷范数，默认是无穷范数
         * @param show_details
         * @param max_steps 非稀疏矩阵收敛过慢
         */
        void
        power_method(double epsilon = 1e-5, bool show_details = false, norm norm = infinite_norm, int max_steps = 2000);


        /**
         * 反幂法
         * 1. 计算矩阵按模最小的特征值与其对应的特征向量
         * @param epsilon 允许误差
         * @param show_details
         * @param norm 范数的选择，只能选二范数和无穷范数，默认是二范数
         * @param max_steps 非稀疏矩阵收敛过慢
         */
        void inverse_power_method(double epsilon = 1e-5, bool show_details = false, norm norm = two_norm,
                                  int max_steps = 2000);


        /**
         * 反幂法
         * 2. 求一个给定近似特征值对应的特征向量
         * @param eigenvalue 给定近似特征值
         * @param epsilon 允许误差
         * @param show_details
         * @param norm 范数的选择，只能选二范数和无穷范数，默认是二范数
         * @param max_steps 非稀疏矩阵收敛过慢
         */
        void
        inverse_power_method(double eigenvalue, double epsilon = 1e-5, bool show_details = false, norm norm = two_norm,
                             int max_steps = 2000);

        /**
         * Jacobi 方法
         * 只适用于实对称方阵
         * 求出所有特征值和特征向量
         *
         */
        void Jacobi_method()
    };
}


#endif //NUMERICAL ANALYSIS_MATRIXEIGENVALUESANDEIGENVECTORS_H
