//
// Created by convi on 2020/4/5.
//

#include "MatrixEigenvaluesAndEigenvectors.h"

void convivae::MatrixEigenvaluesAndEigenvectors::read_mat() {
    std::ifstream in(_rdr);
    for (std::string s; getline(in, s);) {
        vec_type tmp;

        std::istringstream sin(s);
        for (f8 ia; sin >> ia;) {
            tmp.push_back(ia);
        }
        _original_mat.push_back(tmp);
    }
    in.close();
    _initialized = true;
    _dimension = _original_mat.size();
}

void convivae::MatrixEigenvaluesAndEigenvectors::draw_dividing_line() {
    std::cout << "------------------------------------------------------" << std::endl;
}

void convivae::MatrixEigenvaluesAndEigenvectors::print_vec(const std::string &info,
                                                           const convivae::MatrixEigenvaluesAndEigenvectors::vec_type &a) {
    auto n = a.size();
    std::cout << info << "\n(";
    for (auto i = 0; i < n - 1; ++i)
        std::cout << a[i] << ",";
    std::cout << a[n - 1] << ")" << std::endl;
}

void
convivae::MatrixEigenvaluesAndEigenvectors::power_method(double epsilon, bool show_details, norm norm, int max_steps) {
    std::cout << "幂法求主特征值" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    mat_type a(_original_mat);

    //选取 v0 != 0
    vec_type v0(n);
    v0[n - 1] = 1.0;
    vec_type v(v0);

    //迭代
    auto k = 0;
    auto h = 0.0, h0 = 1.0, beta = 0.0, beta0 = 0.0;

    if (norm == two_norm) {
        for (auto i : v0) {
            h += i * i;
        }
        h = sqrt(h);
    } else {
        for (auto i : v0) {
            const auto tmp = fabs(i);
            if (tmp > h)
                h = tmp;
        }
    }

    while (k++ < max_steps) {
        h0 = h;

        if (show_details)
            print_vec("迭代向量 u_k:", v0);

        //规范化
        for (auto &i : v0)
            i /= fabs(h);

        if (show_details)
            print_vec("规范化后的向量 y_k:", v0);

        // 迭代
        for (auto i = 0; i < n; ++i)
            v[i] = v0[i];

        v0 = MatrixOperation<f8>::matrix_multiply(a, v0);

        h = 0.0;
        if (norm == two_norm) {
            for (auto i : v0) {
                h += i * i;
            }
            h = sqrt(h);
        } else {
            for (auto i : v0) {
                const auto tmp = fabs(i);
                if (tmp > h)
                    h = tmp;
            }
        }

        //结束条件
        beta0 = beta;
        auto sign = h0 > 0.0 ? 1.0 : -1.0;
        beta = h * sign;
        if (fabs(beta) <= epsilon || fabs(beta - beta0) / fabs(beta) <= epsilon) {
            break;
        }
    }

    if (show_details)
        print_vec("迭代向量 u_k:", v0);


    std::cout << "迭代次数:" << k << std::endl;
    std::cout << "主特征值:" << beta << std::endl;
    std::cout << "对应的特征向量: (";
    for (auto i = 0; i < n - 1; ++i)
        std::cout << v[i] << ",";
    std::cout << v[n - 1] << ")" << std::endl;
    draw_dividing_line();
}


void
convivae::MatrixEigenvaluesAndEigenvectors::inverse_power_method(double epsilon, bool show_details, convivae::norm norm,
                                                                 int max_steps) {
    std::cout << "反幂法求按模最小的特征值" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    mat_type a(_original_mat);
    a = MatrixOperation<f8>::matrix_inverse(a);

    //选取 v0 != 0
    vec_type v0(n, 1.0);
    //v0[n - 1] = 1.0;
    vec_type v(v0);

    //迭代
    auto k = 0;
    auto h = 0.0, h0 = 1.0, beta = 0.0, beta0 = 0.0;

    if (norm == two_norm) {
        for (auto i : v0) {
            h += i * i;
        }
        h = sqrt(h);
    } else {
        for (auto i : v0) {
            const auto tmp = fabs(i);
            if (tmp > h)
                h = tmp;
        }
    }

    while (k++ < max_steps) {
        h0 = h;

        if (show_details)
            print_vec("迭代向量 u_k:", v0);

        //规范化
        for (auto &i : v0)
            i /= fabs(h);

        if (show_details)
            print_vec("规范化后的向量 y_k:", v0);

        // 迭代
        for (auto i = 0; i < n; ++i)
            v[i] = v0[i];

        v0 = MatrixOperation<f8>::matrix_multiply(a, v0);

        h = 0.0;
        if (norm == two_norm) {
            for (auto i : v0) {
                h += i * i;
            }
            h = sqrt(h);
        } else {
            for (auto i : v0) {
                const auto tmp = fabs(i);
                if (tmp > h)
                    h = tmp;
            }
        }

        //结束条件
        beta0 = beta;
        beta = 1.0 / MatrixOperation<f8>::matrix_multiply(v, v0);
        if (fabs(beta) <= epsilon || fabs(beta - beta0) / fabs(beta) <= epsilon) {
            break;
        }
    }

    if (show_details)
        print_vec("迭代向量 u_k:", v0);


    std::cout << "迭代次数:" << k << std::endl;
    std::cout << "按模最小的特征值:" << beta << std::endl;
    std::cout << "对应的特征向量: (";
    for (auto i = 0; i < n - 1; ++i)
        std::cout << v[i] << ",";
    std::cout << v[n - 1] << ")" << std::endl;
    draw_dividing_line();
}

void
convivae::MatrixEigenvaluesAndEigenvectors::inverse_power_method(double eigenvalue, double epsilon, bool show_details,
                                                                 convivae::norm norm, int max_steps) {
    std::cout << "反幂法求给定特征值对应的特征向量" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    mat_type a(_original_mat);
    for (auto i = 0; i < n; ++i) {
        a[i][i] -= eigenvalue;
    }

    a = MatrixOperation<f8>::matrix_inverse(a);

    //选取 v0 != 0
    vec_type v0(n, 1.0);
    //v0[n - 1] = 1.0;
    vec_type v(v0);

    //迭代
    auto k = 0;
    auto h = 0.0, h0 = 1.0, beta = 0.0, beta0 = 0.0;

    if (norm == two_norm) {
        for (auto i : v0) {
            h += i * i;
        }
        h = sqrt(h);
    } else {
        for (auto i : v0) {
            const auto tmp = fabs(i);
            if (tmp > h)
                h = tmp;
        }
    }

    while (k++ < max_steps) {
        h0 = h;

        if (show_details)
            print_vec("迭代向量 u_k:", v0);

        //规范化
        for (auto &i : v0)
            i /= fabs(h);

        if (show_details)
            print_vec("规范化后的向量 y_k:", v0);

        // 迭代
        for (auto i = 0; i < n; ++i)
            v[i] = v0[i];

        v0 = MatrixOperation<f8>::matrix_multiply(a, v0);

        h = 0.0;
        if (norm == two_norm) {
            for (auto i : v0) {
                h += i * i;
            }
            h = sqrt(h);
        } else {
            for (auto i : v0) {
                const auto tmp = fabs(i);
                if (tmp > h)
                    h = tmp;
            }
        }

        //结束条件
        beta0 = beta;
        beta = 1.0 / MatrixOperation<f8>::matrix_multiply(v, v0);
        if (fabs(beta) <= epsilon || fabs(beta - beta0) / fabs(beta) <= epsilon) {
            break;
        }
    }

    if (show_details)
        print_vec("迭代向量 u_k:", v0);


    std::cout << "迭代次数:" << k << std::endl;
    std::cout << "特征值的精确值:" << eigenvalue + beta << std::endl;
    std::cout << "对应的特征向量: (";
    for (auto i = 0; i < n - 1; ++i)
        std::cout << v[i] << ",";
    std::cout << v[n - 1] << ")" << std::endl;
    draw_dividing_line();
}


