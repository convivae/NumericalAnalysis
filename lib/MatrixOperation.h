//
// Created by convi on 2020/4/13.
//

#ifndef NUMERICALANALYSIS_MATRIXOPERATION_H
#define NUMERICALANALYSIS_MATRIXOPERATION_H


#include <vector>
#include <cmath>
#include <iostream>
#include "../Type/Type.h"

namespace convivae {

    template<typename T>
    class MatrixOperation {
    private:
        typedef std::vector<std::vector<T>> mat_type;
        typedef std::vector<T> vec_type;
    public:
        /**
         * 矩阵乘以矩阵
         * @param a
         * @param b
         * @return a*b
         */
        static mat_type matrix_multiply(mat_type const &a, mat_type const &b);

        /**
         * 矩阵数乘
         * @param a
         * @param b
         * @return  矩阵
         */
        static mat_type matrix_multiply(mat_type const &a, f8 const &b);

        /**
         * 矩阵乘以列向量
         * @param a
         * @param b
         * @return
         */
        static vec_type matrix_multiply(mat_type const &a, vec_type const &b);

        /**
         * 行向量乘以列向量
         * @param a
         * @param b
         * @return 一个数
         */
        static T matrix_multiply(vec_type const &a, vec_type const &b);


        /**
         * 矩阵除法
         * @param a
         * @param b
         * @return
         */
        static mat_type matrix_div(mat_type const &a, f8 const &b);

        /**
         * 矩阵求逆
         * @param m
         * @return
         */
        static mat_type matrix_inverse(mat_type const &m);


        /**
         * 矩阵转置
         * @param m
         * @return
         */
        static mat_type matrix_transposition(mat_type const &m);
    };

    template<typename T>
    inline std::vector<T> MatrixOperation<T>::matrix_multiply(mat_type const &a, vec_type const &b) {
        if (a[0].size() != b.size()) {
            std::cerr << "This matrix and vector can't multiply." << std::endl;
            exit(-1);
        }

        vec_type res(a.size());
        for (auto i = 0; i < a.size(); ++i) {
            for (auto j = 0; j < a[i].size(); ++j) {
                res[i] += a[i][j] * b[j];
            }
        }

        return res;
    }

    template<typename T>
    inline T MatrixOperation<T>::matrix_multiply(vec_type const &a, vec_type const &b) {
        if (a.size() != b.size()) {
            std::cerr << "These two vectors can't multiply." << std::endl;
            exit(-1);
        }

        T res(0);
        for (auto i = 0; i < a.size(); ++i)
            res += a[i] * b[i];
        return res;
    }

    template<typename T>
    inline std::vector<std::vector<T>> MatrixOperation<T>::matrix_multiply(mat_type const &a, mat_type const &b) {
        int row = a.size(), column = b[0].size(), mid = b.size();
        if (a[0].size() != mid) {
            std::cerr << "These tow matrix can't multiply." << std::endl;
            exit(-1);
        }
        mat_type res(row);
        for (auto i : res)
            i.resize(column);

        for (auto i = 0; i < row; ++i) {
            for (auto j = 0; j < column; ++j) {
                for (auto k = 0; k < mid; ++k) {
                    res[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return res;
    }

    template<typename T>
    inline std::vector<std::vector<T>> MatrixOperation<T>::matrix_multiply(mat_type const &a, f8 const &b) {
        mat_type res(a);
        for (auto i = 0; i < res.size(); ++i) {
            for (auto j = 0; j < res[i].size(); ++j) {
                res[i][j] *= b;
            }
        }
        return res;
    }


    template<typename T>
    inline std::vector<std::vector<T>> MatrixOperation<T>::matrix_div(mat_type const &a, f8 const &b) {
        mat_type res(a);
        for (auto i = 0; i < res.size(); ++i) {
            for (auto j = 0; j < res[i].size(); ++j) {
                res[i][j] /= b;
            }
        }
        return res;
    }

    template<typename T>
    inline std::vector<std::vector<T>> MatrixOperation<T>::matrix_inverse(mat_type const &m) {
        i4 n = m.size(); //n 阶方阵

        if (m[0].size() != m.size()) {
            std::cerr << "This matrix can't inverse." << std::endl;
            exit(-1);
        }


        // 构造 (a, I)
        mat_type a(n);
        for (auto i = 0; i < n; ++i) {
            a[i].resize(2 * n);
            for (auto j = 0; j < n; ++j) {
                a[i][j] = m[i][j];
            }
            for (auto j = n; j < 2 * n; ++j) {
                a[i][j] = i == j - n ? 1.0 : 0.0;
            }
        }

        // 消元
        for (auto k = 0; k < n; ++k) {
            auto max_num = a[k][k] > 0 ? a[k][k] : -a[k][k];
            i4 ik = k;

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

            // 令 a[k][k]为1
            for (auto i = 2 * n - 1; i >= k; --i)
                a[k][i] /= a[k][k];

            //消去其下方
            for (auto i = k + 1; i < n; ++i) {
                auto m1 = a[i][k];
                for (auto j = k; j < 2 * n; ++j) {
                    a[i][j] -= m1 * a[k][j];
                }
            }

            //消去其上方
            for (auto i = k - 1; i >= 0; --i) {
                auto m1 = a[i][k];
                for (auto j = k; j < 2 * n; ++j) {
                    a[i][j] -= m1 * a[k][j];
                }
            }
        }

        // 得到逆矩阵
        mat_type res(n);
        for (auto i = 0; i < n; ++i) {
            res[i].resize(n);
            for (auto j = 0; j < n; ++j) {
                res[i][j] = a[i][j + n];
            }
        }

        return res;
    }

    template<typename T>
    inline std::vector<std::vector<T>> MatrixOperation<T>::matrix_transposition(mat_type const &m) {
        mat_type res(m[0].size());
        for (auto i : res) {
            i.resize(m.size());
        }

        for (auto i = 0; i < res.size(); ++i) {
            for (auto j = 0; j < res[0].size(); ++j) {
                res[i][j] = m[j][i];
            }
        }

        return res;
    }


}


#endif //NUMERICAL ANALYSIS_MATRIX OPERATION_H
