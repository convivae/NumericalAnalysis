//
// 数值分析
// 第二章 线性方程组的解法
// Created by convi on 2020/3/21.
//

#ifndef NUMERICALANALYSIS_LINEAREQUATIONS_H
#define NUMERICALANALYSIS_LINEAREQUATIONS_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <iomanip>

#include "../Type/Type.h"

namespace convivae {

    template<typename T>
    class LinearEquations {
    private:
        typedef std::vector<std::vector<T>> mat_type;
        typedef std::vector<T> vec_type;

        std::string _rdr; //文件名
        mat_type _original_mat; //读入的初始矩阵（增广矩阵）
        vec_type _b_mat;    // b
        bool _initialized; //如果未初始化，就读入矩阵
        i4 _dimension{}; //增广矩阵的维数 n，在读入时初始化

        void read_mat();

        static void print_mat(const std::string &info, mat_type const &a);

        static void print_vec(const std::string &info, const std::string &name, vec_type const &a);

        static void print_vec(const std::string &info, vec_type const &a);

        static void draw_dividing_line();

        static bool is_convergent(vec_type const &a, vec_type const &b, f8 epsilon);

    public:
        explicit LinearEquations(std::string filename)
                : _rdr(std::move(filename)), _initialized(false) {}

        /**
         * Pivot element 主元 在对矩阵做某种算法时,首先进行的部分元素
         * partial pivoting就是只在当前进行变换的列中选择主元，只需要进行行交换
         * 列主元素 Gauss 消去法
         * @param show_details
         */
        void Gaussian_elimination_with_partial_pivoting_method(bool show_details = false);

        /**
         * 高斯若当消去法
         * @param show_details
         */
        void Gauss_Jordan_elimination(bool show_details = false);

        /**
         * 三角分解法：(A = LU)
         * Doolittle 分解法
         * @param show_details
         */
        void LU_factorization_Doolittle(bool show_details = false);

        /**
         * 三角分解法：(A = LU)
         * Crout 分解法
         * @param show_details
         */
        void LU_factorization_Crout(bool show_details = false);

        /**
         * 三角分解法：(A = LU)
         * 选主元的 Doolittle 分解法
         * @param show_details
         */
        void LU_factorization_Doolittle_with_partial_pivoting_method(bool show_details = false);


        /**
         * Jacobi 迭代法
         * @param epsilon 收敛条件
         * @param max_steps 最大迭代次数
         * @param show_details 显示每一步的细节
         */
        void Jacobi_iteration_method(double epsilon = 1e-5, bool show_details = false, int max_steps = 200);


        /**
         * Gauss Seidel 迭代法
         * @param epsilon
         * @param max_steps
         * @param show_details
         */
        void Gauss_Seidel_method(double epsilon = 1e-5, bool show_details = false, int max_steps = 200);


        /**
         * SOR 迭代法
         * @param omega 松弛因子
         * @param epsilon
         * @param show_details
         * @param max_steps
         */
        void Successive_Over_Relaxation_method(double omega = 1.25, double epsilon = 1e-5, bool show_details = false,
                                               int max_steps = 200);

    };


    /**
     * SOR 迭代法
     * @tparam T
     * @param omega 松弛因子
     * @param epsilon
     * @param show_details
     * @param max_steps
     */
    template<typename T>
    void LinearEquations<T>::Successive_Over_Relaxation_method(double omega, double epsilon, bool show_details,
                                                               int max_steps) {
        std::cout << "SOR 迭代法" << std::endl;
        if (!_initialized)
            read_mat();

        i4 n = _dimension;
        mat_type a(_original_mat);
        vec_type b(_b_mat);
        vec_type x0(n), x(n);  //初始值设为全零

        if (show_details)
            print_vec("start:", x);

        i4 k = max_steps;
        while (k--) {
            for (auto i = 0; i < n; ++i) {
                auto sum = 0.0;
                for (auto j = 0; j < i; ++j)
                    sum -= a[i][j] * x[j];
                for (auto j = i + 1; j < n; ++j)
                    sum -= a[i][j] * x0[j];
                sum = (sum + b[i]) / a[i][i];
                x[i] = omega * sum - (omega - 1) * x0[i];
            }

            if (show_details)
                print_vec("next iteration:", x);

            if (is_convergent(x0, x, epsilon))
                break;

            for (auto i = 0; i < n; ++i)
                x0[i] = x[i];
        }
        if (k <= 0) {
            std::cout << "已经进行" << max_steps - k << " 次迭代, 无法收敛到 " << epsilon << std::endl;
            print_vec("最终结果 x：", "x", x);
        } else {
            std::cout << "迭代次数为: " << max_steps - k << std::endl;
            print_vec("计算可得 x：", "x", x);
        }

        draw_dividing_line();
    }

    /**
     * Gauss Seidel 迭代法
     * @tparam T
     * @param epsilon
     * @param max_steps
     * @param show_details
     */
    template<typename T>
    void LinearEquations<T>::Gauss_Seidel_method(double epsilon, bool show_details, int max_steps) {
        std::cout << "Gauss Seidel 迭代法" << std::endl;
        if (!_initialized)
            read_mat();

        i4 n = _dimension;
        mat_type a(_original_mat);
        vec_type b(_b_mat);
        vec_type x0(n), x(n);  //初始值设为全零

        if (show_details)
            print_vec("start:", x);

        i4 k = max_steps;
        while (k--) {
            for (auto i = 0; i < n; ++i) {
                auto sum = 0.0;
                for (auto j = 0; j < i; ++j)
                    sum += a[i][j] * x[j];
                for (auto j = i + 1; j < n; ++j)
                    sum += a[i][j] * x0[j];

                x[i] = (b[i] - sum) / a[i][i];
            }

            if (show_details)
                print_vec("next iteration:", x);

            if (is_convergent(x0, x, epsilon))
                break;

            for (auto i = 0; i < n; ++i)
                x0[i] = x[i];
        }
        if (k <= 0) {
            std::cout << "已经进行" << max_steps - k << " 次迭代, 无法收敛到 " << epsilon << std::endl;
            print_vec("最终结果 x：", "x", x);
        } else {
            std::cout << "迭代次数为: " << max_steps - k << std::endl;
            print_vec("计算可得 x：", "x", x);
        }

        draw_dividing_line();
    }

    /**
     * Jacobi 迭代法
     * @tparam T
     * @param epsilon
     * @param max_steps
     * @param show_details
     */
    template<typename T>
    void LinearEquations<T>::Jacobi_iteration_method(double epsilon, bool show_details, int max_steps) {
        std::cout << "Jacobi 迭代法" << std::endl;
        if (!_initialized)
            read_mat();

        i4 n = _dimension;
        mat_type a(_original_mat);
        vec_type b(_b_mat);
        vec_type x0(n), x(n);  //初始值设为全零

        if (show_details)
            print_vec("start:", x);

        i4 k = max_steps;
        while (k--) {
            for (auto i = 0; i < n; ++i) {
                auto sum = 0.0;
                for (auto j = 0; j < n; ++j) {
                    if (j == i)
                        continue;
                    sum += a[i][j] * x0[j];
                }
                x[i] = (b[i] - sum) / a[i][i];
            }

            if (show_details)
                print_vec("next iteration:", x);

            if (is_convergent(x0, x, epsilon))
                break;

            for (auto i = 0; i < n; ++i)
                x0[i] = x[i];
        }
        if (k <= 0) {
            std::cout << "已经进行" << max_steps - k << " 次迭代, 无法收敛到 " << epsilon << std::endl;
            print_vec("最终结果 x：", "x", x);
        } else {
            std::cout << "迭代次数为: " << max_steps - k << std::endl;
            print_vec("计算可得 x：", "x", x);
        }

        draw_dividing_line();
    }

    /**
     * 选主元的 Doolittle 分解法
     * @tparam T
     */
    template<typename T>
    void LinearEquations<T>::LU_factorization_Doolittle_with_partial_pivoting_method(bool show_details) {
        std::cout << "三角分解法之选主元的 Doolittle 分解法" << std::endl;
        if (!_initialized)
            read_mat();

        i4 n = _dimension;
        mat_type a(_original_mat);
        vec_type b(_b_mat);

        vec_type x(n), y(n), M(n), s(n);          //中间变量

        if (show_details)
            print_mat("start:", a);

        // (1) 做分解 QA = LU
        for (auto k = 0; k < n; ++k) {
            auto max_num = 0.0;
            for (auto i = k; i < n; ++i) {
                auto sum = 0.0;
                for (auto t = 0; t < k; ++t) {
                    sum += a[i][t] * a[t][k];
                }
                s[i] = a[i][k] - sum;
                auto tmp = s[i];
                if (tmp < 0)
                    tmp *= -1;

                if (tmp > max_num) {
                    max_num = tmp;
                    M[k] = i;
                }
            }

            if (M[k] != k) {
                a[k].swap(a[M[k]]);
                std::swap(s[k], s[M[k]]);
            }

            if (show_details)
                print_mat("换行:", a);

            a[k][k] = s[k];
            for (auto j = k + 1; j < n; ++j) {
                auto sum = 0.0;
                for (auto t = 0; t < k; ++t) {
                    sum += a[k][t] * a[t][j];
                }
                a[k][j] = a[k][j] - sum;
            }

            for (auto i = k + 1; i < n; ++i) {
                a[i][k] = s[i] / a[k][k];
            }
            if (show_details)
                print_mat("消元:", a);
        }
        print_mat("LU 的合矩阵为：", a);

        // (2) 求 Qb
        for (auto k = 0; k < n - 1; ++k) {
            std::swap(b[k], b[M[k]]);
        }
        print_vec("Qb为：", "Qb", b);

        //(3) 求解 Ly=Qb,Ux=y
        for (auto i = 0; i < n; ++i) {
            auto sum = 0.0;
            for (auto t = 0; t < i; ++t) {
                sum += a[i][t] * y[t];
            }
            y[i] = b[i] - sum;
        }
        print_vec("计算可得 y：", "y", y);

        for (auto i = n - 1; i >= 0; --i) {
            auto sum = 0.0;
            for (auto t = i + 1; t < n; ++t) {
                sum += a[i][t] * x[t];
            }
            x[i] = (y[i] - sum) / a[i][i];
        }
        print_vec("计算可得 x：", "x", x);
        draw_dividing_line();
    }

    /**
     * Crout 分解法
     * @tparam T
     */
    template<typename T>
    void LinearEquations<T>::LU_factorization_Crout(bool show_details) {
        std::cout << "三角分解法之 Crout 分解法" << std::endl;
        if (!_initialized)
            read_mat();
        i4 n = _dimension;
        mat_type a(_original_mat);
        vec_type b(_b_mat);
        mat_type combine_LU(n);    //同时存储 LU 的矩阵
        vec_type x(n), y(n);          //中间变量

        for (auto i = 0; i < n; i++)
            combine_LU[i].resize(n);

        // 计算 LU 矩阵
        for (auto k = 0; k < n; ++k) {
            for (auto i = k; i < n; ++i) {
                auto sum = 0.0;
                for (auto t = 0; t < k; ++t) {
                    sum += combine_LU[i][t] * combine_LU[t][k];
                }
                combine_LU[i][k] = a[i][k] - sum;
            }

            if (show_details)
                print_mat("消元:", combine_LU);

            for (auto j = k + 1; j < n; ++j) {
                auto sum = 0.0;
                for (auto t = 0; t < k; ++t) {
                    sum += combine_LU[k][t] * combine_LU[t][j];
                }
                combine_LU[k][j] = (a[k][j] - sum) / combine_LU[k][k];
            }
            if (show_details)
                print_mat("消元:", combine_LU);
        }
        print_mat("LU 的合矩阵为：", combine_LU);

        //计算 y、x
        for (auto i = 0; i < n; ++i) {
            auto sum = 0.0;
            for (auto t = 0; t < i; ++t) {
                sum += combine_LU[i][t] * y[t];
            }
            y[i] = (b[i] - sum) / combine_LU[i][i];
        }
        print_vec("计算可得 y：", "y", y);

        for (auto i = n - 1; i >= 0; --i) {
            auto sum = 0.0;
            for (auto t = i + 1; t < n; ++t) {
                sum += combine_LU[i][t] * x[t];
            }
            x[i] = y[i] - sum;
        }
        print_vec("计算可得 x：", "x", x);

        draw_dividing_line();
    }

    /**
     * Doolittle 分解法
     * 输入为增广矩阵
     * @tparam T
     */
    template<typename T>
    void LinearEquations<T>::LU_factorization_Doolittle(bool show_details) {
        std::cout << "三角分解法之 Doolittle 分解法" << std::endl;
        if (!_initialized)
            read_mat();
        i4 n = _dimension;
        mat_type a(_original_mat);
        vec_type b(_b_mat);
        mat_type combine_LU(n);    //同时存储 LU 的矩阵
        vec_type x(n), y(n);          //中间变量

        for (auto i = 0; i < n; i++)
            combine_LU[i].resize(n);

        // 计算 LU 矩阵
        for (auto k = 0; k < n; ++k) {
            for (auto j = k; j < n; ++j) {
                auto sum = 0.0;
                for (auto t = 0; t < k; ++t) {
                    sum += combine_LU[k][t] * combine_LU[t][j];
                }
                combine_LU[k][j] = a[k][j] - sum;
            }

            if (show_details)
                print_mat("消元:", combine_LU);

            for (auto i = k + 1; i < n; ++i) {
                auto sum = 0.0;
                for (auto t = 0; t < k; ++t) {
                    sum += combine_LU[i][t] * combine_LU[t][k];
                }
                combine_LU[i][k] = (a[i][k] - sum) / combine_LU[k][k];
            }
            if (show_details)
                print_mat("消元:", combine_LU);
        }
        print_mat("LU 的合矩阵为：", combine_LU);

        //计算 y、x
        for (auto i = 0; i < n; ++i) {
            auto sum = 0.0;
            for (auto t = 0; t < i; ++t) {
                sum += combine_LU[i][t] * y[t];
            }
            y[i] = b[i] - sum;
        }
        print_vec("计算可得 y：", "y", y);

        for (auto i = n - 1; i >= 0; --i) {
            auto sum = 0.0;
            for (auto t = i + 1; t < n; ++t) {
                sum += combine_LU[i][t] * x[t];
            }
            x[i] = (y[i] - sum) / combine_LU[i][i];
        }
        print_vec("计算可得 x：", "x", x);

        draw_dividing_line();
    }

    /**
     * 高斯若当消去法
     * 输入为增广矩阵
     * 输出系数矩阵的逆矩阵以及线性方程的解
     * @tparam T
     */
    template<typename T>
    void LinearEquations<T>::Gauss_Jordan_elimination(bool show_details) {
        std::cout << "高斯若当消去法" << std::endl;
        if (!_initialized)
            read_mat();

        i4 n = _dimension; //n 行，n+1 列
        i4 ik; //行号

        // 构造 (a, I)
        mat_type a(n);
        for (auto i = 0; i < n; ++i) {
            a[i].resize(2 * n);
            for (auto j = 0; j < n; ++j) {
                a[i][j] = _original_mat[i][j];
            }
            for (auto j = n; j < 2 * n; ++j) {
                a[i][j] = i == j - n ? 1.0 : 0.0;
            }
        }

        if (show_details)
            print_mat("start:", a);

        vec_type b(_b_mat);

        // 消元
        for (auto k = 0; k < n; ++k) {
            auto max_num = a[k][k] > 0 ? a[k][k] : -a[k][k];
            ik = k;

            // 找到绝对值最大的行号交换
            for (auto i = k; i < n; ++i) {
                auto tmp = a[i][k] > 0 ? a[i][k] : -a[i][k];
                if (tmp > max_num) {
                    max_num = tmp;
                    ik = i;
                }
            }
            if (ik != k) {
                a[k].swap(a[ik]);
            }
            if (show_details)
                print_mat("换行:", a);

            // 令 a[k][k]为1
            for (auto i = 2 * n - 1; i >= k; --i)
                a[k][i] /= a[k][k];

            //消去其下方
            for (auto i = k + 1; i < n; ++i) {
                auto m = a[i][k];
                for (auto j = k; j < 2 * n; ++j) {
                    a[i][j] -= m * a[k][j];
                }
            }

            //消去其上方
            for (auto i = k - 1; i >= 0; --i) {
                auto m = a[i][k];
                for (auto j = k; j < 2 * n; ++j) {
                    a[i][j] -= m * a[k][j];
                }
            }
            if (show_details)
                print_mat("消元:", a);
        }
        // 输出消元后的结果
        print_mat("消元后的增广矩阵为：", a);

        // 得到逆矩阵
        mat_type a_inverse(n);
        for (auto i = 0; i < n; ++i) {
            a_inverse[i].resize(n);
            for (auto j = 0; j < n; ++j) {
                a_inverse[i][j] = a[i][j + n];
            }
        }

        print_mat("系数矩阵的逆矩阵为：", a_inverse);


        //回代
        vec_type x(n);
        for (auto k = 0; k < n; ++k) {
            x[k] = 0;
            for (auto i = 0; i < n; ++i) {
                x[k] += a_inverse[k][i] * b[i];
            }
        }

        print_vec("求得：", "x", x);
        draw_dividing_line();
    }

    /**
     * 列主元素 Gauss 消去法
     * 输入为增广矩阵
     * @tparam T
     */
    template<typename T>
    void LinearEquations<T>::Gaussian_elimination_with_partial_pivoting_method(bool show_details) {
        std::cout << "列主元素 Gauss 消去法" << std::endl;
        if (!_initialized)
            read_mat();

        mat_type a(_original_mat);
        i4 ik; //行号
        i4 n = _dimension; //n 行，n+1 列

        if (show_details)
            print_mat("start:", a);

        for (auto k = 0; k < n - 1; ++k) {
            auto max_num = a[k][k] > 0 ? a[k][k] : -a[k][k];
            ik = k;

            // 找到绝对值最大的行号交换
            for (auto i = k; i < n; ++i) {
                auto tmp = a[i][k] > 0 ? a[i][k] : -a[i][k];
                if (tmp > max_num) {
                    max_num = tmp;
                    ik = i;
                }
            }
            if (ik != k) {
                a[k].swap(a[ik]);
            }

            if (show_details)
                print_mat("换行:", a);

            //消元
            for (auto i = k + 1; i < n; ++i) {
                auto m = a[i][k] / a[k][k];
                for (auto j = k; j < n + 1; ++j) {
                    a[i][j] -= m * a[k][j];
                }
            }
            if (show_details)
                print_mat("消元:", a);
        }
        // 输出消元后的结果
        print_mat("消元后的增广矩阵为：", a);

        //回代
        vec_type x(n);
        x[n - 1] = a[n - 1][n] / a[n - 1][n - 1];
        for (auto k = n - 2; k >= 0; --k) {
            auto sum = 0.0;
            for (auto j = k + 1; j < n; ++j) {
                sum += a[k][j] * x[j];
            }

            x[k] = (a[k][n] - sum) / a[k][k];
        }

        print_vec("回代后，得：", "x", x);
        draw_dividing_line();
    }

    template<typename T>
    void LinearEquations<T>::read_mat() {
        std::ifstream in(_rdr);
        for (std::string s; getline(in, s);) {
            vec_type tmp;

            std::istringstream sin(s);
            for (T ia; sin >> ia;) {
                tmp.push_back(ia);
            }
            _original_mat.push_back(tmp);
        }
        in.close();
        _initialized = true;
        _dimension = _original_mat.size();

        _b_mat.resize(_dimension);
        for (auto i = 0; i < _dimension; ++i) {
            _b_mat[i] = _original_mat[i][_dimension];
        }
    }

    template<typename T>
    void LinearEquations<T>::print_mat(const std::string &info, const mat_type &a) {
        std::cout << info << std::endl;
        for (auto i : a) {
            for (auto j : i) {
                std::cout << std::left << std::setw(14) << j;
            }
            std::cout << std::endl;
        }
    }

    template<typename T>
    void LinearEquations<T>::print_vec(const std::string &info, const std::string &name,
                                       const LinearEquations::vec_type &a) {
        std::cout << info << std::endl;
        for (auto i = 0; i < a.size(); ++i)
            std::cout << name << i + 1 << " = " << a[i] << std::endl;
    }

    template<typename T>
    void LinearEquations<T>::print_vec(const std::string &info, const LinearEquations::vec_type &a) {
        std::cout << info << std::endl;
        for (auto i : a)
            std::cout << i << " ";
        std::cout << std::endl;
    }

    template<typename T>
    void LinearEquations<T>::draw_dividing_line() {
        std::cout << "------------------------------------------------------" << std::endl;
    }

    template<typename T>
    bool LinearEquations<T>::is_convergent(const LinearEquations::vec_type &a, const LinearEquations::vec_type &b,
                                           f8 epsilon) {
        vec_type x(a.size());
        for (auto i = 0; i < a.size(); ++i) {
            auto tmp = a[i] - b[i];
            x[i] = tmp > 0 ? tmp : -tmp;
        }

        for (auto i : x) {
            if (i > epsilon)
                return false;
        }
        return true;
    }
}

#endif //NUMERICALANALYSIS_LINEAREQUATIONS_H
