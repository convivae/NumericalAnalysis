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

        const f8 PI = atan(1.0) * 4.0;

        void read_mat();

        static void draw_dividing_line();

        static void print_vec(const std::string &info, vec_type const &a);

        static void print_mat(const std::string &info, mat_type const &a);

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
        void inverse_power_method(norm norm = two_norm, double epsilon = 1e-5, bool show_details = false,
                                  int max_steps = 2000);


        /**
         * 反幂法
         * 2. 求一个给定近似特征值对应的特征向量
         * @param eigenvalue 给定近似特征值
         * @param epsilon 允许误差
         * @param show_details
         * @param norm 范数的选择，只能选二范数和无穷范数，默认是二范数
         * @param inner_use 是否是由其他函数调用（若是，只输出最终结果，不输出多余信息）
         * @param max_steps 非稀疏矩阵收敛过慢
         */
        void
        inverse_power_method(double eigenvalue, double epsilon = 1e-5, bool show_details = false,
                             norm norm = two_norm, bool inner_use = false, int max_steps = 2000);


        /**
         * Jacobi 方法
         * 只适用于实对称方阵，求出所有特征值和特征向量
         * @param epsilon
         * @param show_details
         */
        void Jacobi_method(double epsilon = 1e-5, bool show_details = false);


        /**
         * 利用 Householder 矩阵（镜面映射矩阵）对矩阵 A 做 QR 分解
         * A = QR
         * @param A n*n 实矩阵
         * @param Q 返回值，正交矩阵
         * @param R 返回值，上三角矩阵
         *
         */
        void QR_decomposition_Householder(mat_type A, mat_type &Q, mat_type &R, bool show_QR = false);


        /**
         * QR 方法求矩阵的特征值（对应的特征向量可利用反幂法得到）
         * 主要用来求上 Hessberg 矩阵的全部特征值（矩阵A将会转化为上 Hessberg 矩阵再做计算）
         * QR 分解部分是利用 Householder 矩阵（镜面映射矩阵）对矩阵 A 做 QR 分解
         * 没做任何优化
         * @param epsilon
         * @param show_details
         * @param max_steps
         */
        void QR_method_normal(double epsilon = 1e-5, bool show_details = false, int max_steps = 2000);

        /**
         * QR 方法求矩阵的特征值（对应的特征向量可利用反幂法得到）
         * 主要用来求上 Hessberg 矩阵的全部特征值（矩阵A将会转化为上 Hessberg 矩阵再做计算）
         * QR 分解部分是利用平面旋转矩阵对矩阵 A 做 QR 分解
         * 下三角部分会收敛到 0，上三角部分可能不收敛，主对角线元素收敛到特征值
         * 没做任何优化
         * @param epsilon
         * @param show_details
         */
        void
        QR_method_Transformed_Rotational_Matrix(double epsilon = 1e-5, bool show_details = false, int max_steps = 2000);

        /**
         * 带原点位移的 QR 分解，加快迭代的收敛速度
         * QR 方法求矩阵的特征值（对应的特征向量可利用反幂法得到）
         * 主要用来求上 Hessberg 矩阵的全部特征值（矩阵A将会转化为上 Hessberg 矩阵再做计算）
         * 再利用平面旋转矩阵对矩阵 A 做 QR 分解
         * @param epsilon
         * @param show_details
         * @param max_steps
         */
        void QR_method_with_one_step_shifted(double epsilon = 1e-5, bool show_details = false, int max_steps = 2000);

        /**
         * 带双步位移的 QR 分解，加快迭代的收敛速度
         * QR 方法求矩阵的特征值（对应的特征向量可利用反幂法得到）
         * 主要用来求上 Hessberg 矩阵的全部特征值（矩阵A将会转化为上 Hessberg 矩阵再做计算）
         * 再利用平面旋转矩阵对矩阵 A 做 QR 分解
         * @param epsilon
         * @param show_details
         * @param max_steps
         */
        void QR_method_with_two_steps_shifted(double epsilon = 1e-5, bool show_details = false, int max_steps = 2000);

    private:

        /**
         * 用 Householder 方法把矩阵 A 化为上 Hessenberg 矩阵（拟上三角矩阵）
         * 把 A 变成上 Hessenberg 矩阵（拟上三角矩阵）的目是减少 QR 方法的计算量。原因是：
         * 1. 对拟上三角矩阵作 QR 分解时，Q 一定是拟上三角矩阵
         * 2. RQ (＝A_(k＋1) )的乘积为拟上三角矩阵。
         * @param a
         * @param show_P 显示 A_(hess) = P'AP
         * @return
         */
        mat_type Hessenberg_matrix(mat_type A) const;
    };
}


#endif //NUMERICAL ANALYSIS_MATRIXEIGENVALUESANDEIGENVECTORS_H
