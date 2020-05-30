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
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

void convivae::MatrixEigenvaluesAndEigenvectors::print_vec(const std::string &info,
                                                           const convivae::MatrixEigenvaluesAndEigenvectors::vec_type &
                                                           a) {
    auto n = a.size();
    std::cout << info << "\n(";
    for (auto i = 0; i < n - 1; ++i)
        std::cout << a[i] << ",";
    std::cout << a[n - 1] << ")" << std::endl;
}

void convivae::MatrixEigenvaluesAndEigenvectors::print_mat(const std::string &info,
                                                           const convivae::MatrixEigenvaluesAndEigenvectors::mat_type &
                                                           a) {
    std::cout << info << std::endl;
    for (auto i : a) {
        for (auto j : i) {
            std::cout << std::left << std::setw(14) << j;
        }
        std::cout << std::endl;
    }
}

void
convivae::MatrixEigenvaluesAndEigenvectors::power_method(double epsilon, bool show_details, norm norm, int max_steps) {
    std::cout << "幂法求主特征值" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    mat_type a(_original_mat);

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

        auto tmp = MatrixOperation<f8>::matrix_multiply(a, MatrixOperation<f8>::matrix_transposition(v0));
        for (auto i = 0; i < tmp.size(); ++i)
            v0[i] = tmp[i][0];

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

void convivae::MatrixEigenvaluesAndEigenvectors::inverse_power_method(norm norm, double epsilon, bool show_details,
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

        auto tmp = MatrixOperation<f8>::matrix_multiply(a, MatrixOperation<f8>::matrix_transposition(v0));
        for (auto i = 0; i < tmp.size(); ++i)
            v0[i] = tmp[i][0];

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
        beta = 1.0 / MatrixOperation<f8>::matrix_multiply(v, MatrixOperation<f8>::matrix_transposition(v0));
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
                                                                 norm norm, bool inner_use, int max_steps) {
    if (!inner_use)
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

        if (show_details && !inner_use)
            print_vec("迭代向量 u_k:", v0);

        //规范化
        for (auto &i : v0)
            i /= fabs(h);

        if (show_details && !inner_use)
            print_vec("规范化后的向量 y_k:", v0);

        // 迭代
        for (auto i = 0; i < n; ++i)
            v[i] = v0[i];

        auto tmp = MatrixOperation<f8>::matrix_multiply(a, MatrixOperation<f8>::matrix_transposition(v0));
        for (auto i = 0; i < tmp.size(); ++i)
            v0[i] = tmp[i][0];

        h = 0.0;
        if (norm == two_norm) {
            for (auto i : v0) {
                h += i * i;
            }
            h = sqrt(h);
        } else {
            for (auto i : v0) {
                if (fabs(i) > h)
                    h = fabs(i);
            }
        }

        //结束条件
        beta0 = beta;
        beta = 1.0 / MatrixOperation<f8>::matrix_multiply(v, MatrixOperation<f8>::matrix_transposition(v0));
        if (k > 1 && (fabs(beta) <= epsilon || fabs(beta - beta0) / fabs(beta) <= epsilon)) {
            break;
        }
    }

    if (show_details && !inner_use)
        print_vec("迭代向量 u_k:", v0);


    if (!inner_use) {
        std::cout << "迭代次数:" << k << std::endl;
        std::cout << "特征值的给定值:" << eigenvalue << std::endl;
        std::cout << "特征值的精确值:" << eigenvalue + beta << std::endl;
        std::cout << "对应的特征向量: (";
    }

    for (auto i = 0; i < n - 1; ++i)
        std::cout << v[i] << ",";
    std::cout << v[n - 1] << ")" << std::endl;

    if (!inner_use)
        draw_dividing_line();
}


void convivae::MatrixEigenvaluesAndEigenvectors::Jacobi_method(double epsilon, bool show_details) {
    std::cout << "Jacobi 方法求实对称矩阵所有特征值和特征向量" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    mat_type a(_original_mat), res_U = MatrixOperation<f8>::identity_matrix(n);

    if (!MatrixOperation<f8>::is_symmetric_matrix(a)) {
        std::cout << "不是实对称矩阵，无法使用 Jacobi 方法求特征值" << std::endl;
        draw_dividing_line();
        return;
    }

    auto count = 0;
    while (true) {
        count++;
        auto phi = 0.0, sin_phi = 0.0, cos_phi = 0.0; //角度PHI
        // 1. 选取非主对角线元素的主元素（按模最大的元素）
        auto apq = 0.0;
        auto p = -1, q = -1;
        for (auto i = 1; i < n; ++i) {
            for (auto j = 0; j < i; ++j) {
                if (fabs(apq) < fabs(a[i][j])) {
                    apq = a[i][j];
                    p = i;
                    q = j;
                }
            }
        }

        if (fabs(apq) < epsilon)
            break;

        if (show_details)
            std::cout << "选主元: A[" << p << "][" << q << "] = " << apq << std::endl;

        if (a[p][p] == a[q][q]) {
            auto sign = apq > 0.0 ? 1.0 : -1.0;
            phi = PI / 4.0 * sign;
            sin_phi = sin(phi);
            cos_phi = cos(phi);
        } else {
            auto c = (a[p][p] - a[q][q]) / (2.0 * a[p][q]);
            auto sign = c > 0.0 ? 1.0 : -1.0;
            auto t = sign / (fabs(c) + sqrt(c * c + 1));
            cos_phi = 1.0 / sqrt(1 + t * t);
            sin_phi = t / sqrt(1 + t * t);
        }

        if (show_details)
            std::cout << "sin(phi)=" << sin_phi << " cos(phi)=" << cos_phi << std::endl;

        //2. 构造平面旋转矩阵
        auto U = MatrixOperation<f8>::identity_matrix(n);
        U[p][p] = cos_phi;
        U[q][q] = cos_phi;
        U[p][q] = -sin_phi;
        U[q][p] = sin_phi;

        if (show_details)
            print_mat("平面旋转矩阵 U =", U);

        res_U = MatrixOperation<f8>::matrix_multiply(res_U, U);

        //3. 计算 a=U'aU
        auto UT = MatrixOperation<f8>::matrix_transposition(U);
        a = MatrixOperation<f8>::matrix_multiply(UT, a);
        a = MatrixOperation<f8>::matrix_multiply(a, U);
        if (show_details)
            print_mat("A = U'AU =", a);
    }

    std::cout << "迭代次数: " << count << std::endl;
    for (auto j = 0; j < n; ++j) {
        std::cout << "特征值: " << a[j][j] << "  对应的特征向量: (";
        for (auto i = 0; i < n - 1; ++i)
            std::cout << res_U[i][j] << ",";
        std::cout << res_U[n - 1][j] << ")" << std::endl;
    }
    draw_dividing_line();
}

void convivae::MatrixEigenvaluesAndEigenvectors::cal_Householder_matrix_Hessenberg(bool show_details) {
    std::cout << "计算初等反射矩阵（Householder 矩阵）使 A 化为拟上三角矩阵" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    mat_type A(_original_mat);
    mat_type I = MatrixOperation<f8>::identity_matrix(n);
    mat_type H(I);

    if (show_details)
        print_mat("start: A =", A);

    vec_type s(n);
    for (auto i = 0; i < n - 2; ++i) {
        //1. s
        s[i] = 0.0;
        for (auto j = i + 1; j < n; ++j)
            s[j] = A[j][i];

        if (show_details) {
            std::cout << "s = (";
            for (auto index = 0; index < n - 1; ++index)
                std::cout << s[index] << ",";
            std::cout << s[n - 1] << ")" << std::endl;
        }

        //2. c
        auto c = MatrixOperation<f8>::matrix_multiply(s, MatrixOperation<f8>::matrix_transposition(s));
        auto sign = s[i + 1] > 0.0 ? -1.0 : 1.0;
        c = sqrt(c) * sign;

        if (show_details)
            std::cout << "c = " << c << std::endl;

        //3. u
        vec_type u(s);
        u[i + 1] -= c;

        if (show_details) {
            std::cout << "u = (";
            for (auto index = 0; index < n - 1; ++index)
                std::cout << u[index] << ",";
            std::cout << u[n - 1] << ")" << std::endl;
        }

        //4. H
        auto u_value = MatrixOperation<f8>::matrix_multiply(u, MatrixOperation<f8>::matrix_transposition(u));
        H = MatrixOperation<f8>::matrix_multiply(MatrixOperation<f8>::matrix_transposition(u), u);
        H = MatrixOperation<f8>::matrix_multiply(H, 2.0 / u_value);
        H = MatrixOperation<f8>::matrix_sub(I, H);

        A = MatrixOperation<f8>::matrix_multiply(H, A);
        A = MatrixOperation<f8>::matrix_multiply(A, H);

        if (show_details) {
            std::cout << "\n第" << i + 1 << "次迭代计算可得: ";
            print_mat("H =", H);
            print_mat("A =", A);
        }
    }

    std::cout << "\n最终结果为: ";
    print_mat("H =", H);
    print_mat("A =", A);

    draw_dividing_line();
}

void convivae::MatrixEigenvaluesAndEigenvectors::cal_Householder_matrix_normal(bool show_details) {
    std::cout << "计算初等反射矩阵（Householder 矩阵）使 A 化为上三角矩阵" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    mat_type A(_original_mat);
    mat_type I = MatrixOperation<f8>::identity_matrix(n);
    mat_type H(I), Q(I);

    if (show_details)
        print_mat("start: A =", A);

    vec_type s(n);
    for (auto i = 0; i < n - 1; ++i) {
        //1. s
        for (auto j = i; j < n; ++j)
            s[j] = A[j][i];

        if (show_details) {
            std::cout << "s = (";
            for (auto index = 0; index < n - 1; ++index)
                std::cout << s[index] << ",";
            std::cout << s[n - 1] << ")" << std::endl;
        }

        //2. c
        auto c = MatrixOperation<f8>::matrix_multiply(s, MatrixOperation<f8>::matrix_transposition(s));
        auto sign = s[i] > 0.0 ? -1.0 : 1.0;
        c = sqrt(c) * sign;

        if (show_details)
            std::cout << "c = " << c << std::endl;

        //3. u
        vec_type u(s);
        u[i] -= c;

        if (show_details) {
            std::cout << "u = (";
            for (auto index = 0; index < n - 1; ++index)
                std::cout << u[index] << ",";
            std::cout << u[n - 1] << ")" << std::endl;
        }

        //4. H
        auto u_value = MatrixOperation<f8>::matrix_multiply(u, MatrixOperation<f8>::matrix_transposition(u));
        H = MatrixOperation<f8>::matrix_multiply(MatrixOperation<f8>::matrix_transposition(u), u);
        H = MatrixOperation<f8>::matrix_multiply(H, 2.0 / u_value);
        H = MatrixOperation<f8>::matrix_sub(I, H);

        Q = MatrixOperation<f8>::matrix_multiply(Q, H);

        A = MatrixOperation<f8>::matrix_multiply(H, A);

        if (show_details) {
            std::cout << "\n第" << i + 1 << "次迭代计算可得: ";
            print_mat("H =", H);
            print_mat("A =", A);
        }

        s[i] = 0.0;
    }

    std::cout << "\n最终结果为: ";
    print_mat("H =", H);
    print_mat("Q = ", Q);
    print_mat("R = A_n =", A);
    print_mat("验证：A = Q * R =", MatrixOperation<f8>::matrix_multiply(Q, A));

    draw_dividing_line();
}

void convivae::MatrixEigenvaluesAndEigenvectors::QR_decomposition_Householder(
        convivae::MatrixEigenvaluesAndEigenvectors::mat_type A, convivae::MatrixEigenvaluesAndEigenvectors::mat_type &Q,
        convivae::MatrixEigenvaluesAndEigenvectors::mat_type &R, bool show_QR) {
    if (show_QR) {
        std::cout << "QR 分解" << std::endl;
        if (A.empty()) {
            if (!_initialized)
                read_mat();
            A = mat_type(_original_mat);
        }
        print_mat("A=:", A);
    }


    auto n = A.size();
    Q = MatrixOperation<f8>::identity_matrix(n);

    for (auto r = 0; r < n - 1; ++r) {
        // 1.
        for (auto i = r + 1; i < n; ++i) {
            if (A[i][r] != 0.0)
                break;
            if (i == n - 1) {
                continue;
            }
        }

        // 2.
        auto d = 0.0;
        for (auto i = r; i < n; ++i) {
            d += A[i][r] * A[i][r];
        }
        d = sqrt(d);

        auto c = A[r][r] > 0.0 ? -d : d;
        auto h = c * c - c * A[r][r];

        // 3.
        vec_type u(n);
        for (auto i = 0; i < r; ++i) {
            u[i] = 0.0;
        }
        u[r] = A[r][r] - c;
        for (auto i = r + 1; i < n; ++i) {
            u[i] = A[i][r];
        }

        // 4.
        auto ut = MatrixOperation<f8>::matrix_transposition(u);
        auto omega = MatrixOperation<f8>::matrix_multiply(Q, ut);
        auto tmp = MatrixOperation<f8>::matrix_multiply(omega, u);
        tmp = MatrixOperation<f8>::matrix_multiply(tmp, 1.0 / h);
        Q = MatrixOperation<f8>::matrix_sub(Q, tmp);

        auto AT = MatrixOperation<f8>::matrix_transposition(A);
        tmp = MatrixOperation<f8>::matrix_multiply(AT, ut);
        tmp = MatrixOperation<f8>::matrix_multiply(tmp, 1.0 / h);

        tmp = MatrixOperation<f8>::matrix_transposition(tmp);
        tmp = MatrixOperation<f8>::matrix_multiply(ut, tmp);
        A = MatrixOperation<f8>::matrix_sub(A, tmp);
    }
    R = A;
    if (show_QR) {
        print_mat("Q:", Q);
        print_mat("R:", R);
        print_mat("验算 Q*R=:", MatrixOperation<f8>::matrix_multiply(Q, R));
        draw_dividing_line();
    }
}

convivae::MatrixEigenvaluesAndEigenvectors::mat_type
convivae::MatrixEigenvaluesAndEigenvectors::Hessenberg_matrix(
        convivae::MatrixEigenvaluesAndEigenvectors::mat_type A) const {
    auto n = A.size();

    for (auto r = 0; r < n - 2; ++r) {
        // 1.
        auto k = r + 2;
        for (; k < n; ++k) {
            if (A[k][r] != 0)
                break;
        }
        if (k >= n)
            continue;

        // 2.
        auto d = 0.0;
        for (auto i = r + 1; i < n; ++i) {
            d += A[i][r] * A[i][r];
        }
        d = sqrt(d);

        auto c = A[r + 1][r] > 0.0 ? -d : d;
        auto h = c * c - c * A[r + 1][r];

        // 3.
        vec_type u(n);
        for (auto i = 0; i < r + 1; ++i) {
            u[i] = 0.0;
        }
        u[r + 1] = A[r + 1][r] - c;
        for (auto i = r + 2; i < n; ++i) {
            u[i] = A[i][r];
        }

        // 4.
        auto ut = MatrixOperation<f8>::matrix_transposition(u);
        auto AT = MatrixOperation<f8>::matrix_transposition(A);

        auto p = MatrixOperation<f8>::matrix_multiply(AT, ut);
        p = MatrixOperation<f8>::matrix_multiply(p, 1.0 / h);

        auto q = MatrixOperation<f8>::matrix_multiply(A, ut);
        q = MatrixOperation<f8>::matrix_multiply(q, 1.0 / h);

        vec_type pt(n);
        for (auto i = 0; i < n; ++i)
            pt[i] = p[i][0];
        auto t = MatrixOperation<f8>::matrix_multiply(pt, ut) / h;
        auto omega = MatrixOperation<f8>::matrix_multiply(ut, t);
        omega = MatrixOperation<f8>::matrix_sub(q, omega);

        A = MatrixOperation<f8>::matrix_sub(A, MatrixOperation<f8>::matrix_multiply(omega, u));
        A = MatrixOperation<f8>::matrix_sub(A, MatrixOperation<f8>::matrix_multiply(ut, pt));
    }

    return A;
}

void convivae::MatrixEigenvaluesAndEigenvectors::QR_method_normal(double epsilon, bool show_details, int max_steps) {
    std::cout << "QR 方法（QR分解由 Householder 矩阵求出）" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    if (show_details)
        print_mat("start:", _original_mat);

    auto A = Hessenberg_matrix(_original_mat);

    if (show_details)
        print_mat("化为拟上三角矩阵:", A);

    mat_type Q, R;

    auto k = 0;
    while (k++ < max_steps) {
        if (show_details)
            std::cout << "第" << k << "次迭代:" << std::endl;

        QR_decomposition_Householder(A, Q, R);

        if (show_details) {
            std::cout << "利用 Householder 方法进行 QR 分解" << std::endl;
            print_mat("Q:", Q);
            print_mat("R:", R);
        }

        A = MatrixOperation<f8>::matrix_multiply(R, Q);

        if (show_details)
            print_mat("A = R * Q", A);

        //停止条件
        auto i = 1;
        for (; i < n; ++i) {
            if (fabs(A[i][i - 1]) > epsilon)
                break;
        }
        if (i >= n)
            break;
    }

    std::cout << "迭代次数: " << k << std::endl;
    for (auto j = 0; j < n; ++j) {
        std::cout << "特征值: " << A[j][j] << "  采用反幂法求对应的特征向量: (";
        inverse_power_method(A[j][j], 1e-5, false, norm::two_norm, true);
    }
    draw_dividing_line();
}

void
convivae::MatrixEigenvaluesAndEigenvectors::QR_method_Transformed_Rotational_Matrix(double epsilon, bool show_details,
                                                                                    int max_steps) {
    std::cout << "QR 方法（QR分解由平面旋转矩阵求出）" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    if (show_details)
        print_mat("start:", _original_mat);

    auto a = Hessenberg_matrix(_original_mat);

    if (show_details)
        print_mat("化为拟上三角矩阵:", a);

    mat_type Q = MatrixOperation<f8>::identity_matrix(n);

    auto count = 0;
    while (count++ < max_steps) {
        if (show_details)
            std::cout << "第" << count << "次迭代:" << std::endl;

        for (auto i = 0; i < n - 1; ++i) {
            auto h1 = a[i][i], h2 = a[i + 1][i];
            auto r = sqrt(h1 * h1 + h2 * h2);
            auto cos_phi = h1 / r, sin_phi = -h2 / r;

            if (show_details)
                std::cout << "cos_phi = " << cos_phi << " sin_phi = " << sin_phi << std::endl;

            //2. 构造平面旋转矩阵
            auto U = MatrixOperation<f8>::identity_matrix(n);
            U[i][i] = cos_phi;
            U[i + 1][i + 1] = cos_phi;
            U[i][i + 1] = -sin_phi;
            U[i + 1][i] = sin_phi;

            if (show_details)
                print_mat("平面旋转矩阵: U = ", U);


            a = MatrixOperation<f8>::matrix_multiply(U, a);
            auto UT = MatrixOperation<f8>::matrix_transposition(U);
            Q = MatrixOperation<f8>::matrix_multiply(Q, UT);
            if (show_details) {
                print_mat("R_0 = a, R = U * a = ", a);
                print_mat("Q0 = I, Q = Q * U' = ", Q);
            }

        }

        a = MatrixOperation<f8>::matrix_multiply(a, Q);
        Q = MatrixOperation<f8>::identity_matrix(n);

        if (show_details) {
            print_mat("a = ", a);
            std::cout << std::endl;
        }

        auto i = 1;
        for (; i < n; ++i) {
            if (fabs(a[i][i - 1]) >= epsilon)
                break;
        }
        if (i >= n)
            break;
    }

    std::cout << "迭代次数: " << count << std::endl;
    for (auto j = 0; j < n; ++j) {
        std::cout << "特征值: " << a[j][j] << "  采用反幂法求对应的特征向量: (";
        inverse_power_method(a[j][j], 1e-5, false, norm::two_norm, true);
    }
    draw_dividing_line();
}

void convivae::MatrixEigenvaluesAndEigenvectors::QR_method_with_one_step_shifted(double epsilon, bool show_details,
                                                                                 int max_steps) {
    std::cout << "带原点位移的 QR 方法（QR分解由平面旋转矩阵求出）" << std::endl;
    if (!_initialized)
        read_mat();

    if (show_details)
        print_mat("start:", _original_mat);

    auto a = Hessenberg_matrix(_original_mat);
    i4 n = _dimension;

    if (show_details)
        print_mat("化为拟上三角矩阵:", a);

    mat_type I = MatrixOperation<f8>::identity_matrix(n);
    vec_type res;

    //加速因子取 a[n-1][n-1]
    auto m = n - 1, count = 0;

    while (count++ < max_steps) {
        auto s = a[m][m];
        mat_type Q = MatrixOperation<f8>::identity_matrix(n);

        if (show_details)
            std::cout << "第" << count << "次迭代:" << std::endl;

        //A - s * I
        a = MatrixOperation<f8>::matrix_sub(a, MatrixOperation<f8>::matrix_multiply(I, s));

        for (auto i = 0; i < n - 1; ++i) {
            auto h1 = a[i][i], h2 = a[i + 1][i];
            auto r = sqrt(h1 * h1 + h2 * h2);
            auto cos_phi = h1 / r, sin_phi = -h2 / r;

            if (show_details)
                std::cout << "cos_phi = " << cos_phi << " sin_phi = " << sin_phi << std::endl;

            //2. 构造平面旋转矩阵
            auto U = MatrixOperation<f8>::identity_matrix(n);
            U[i][i] = cos_phi;
            U[i + 1][i + 1] = cos_phi;
            U[i][i + 1] = -sin_phi;
            U[i + 1][i] = sin_phi;

            if (show_details)
                print_mat("平面旋转矩阵: U = ", U);

            a = MatrixOperation<f8>::matrix_multiply(U, a);
            auto UT = MatrixOperation<f8>::matrix_transposition(U);
            Q = MatrixOperation<f8>::matrix_multiply(Q, UT);
            if (show_details) {
                print_mat("R_0 = a, R = U * a = ", a);
                print_mat("Q0 = I, Q = Q * U' = ", Q);
            }
        }

        a = MatrixOperation<f8>::matrix_multiply(a, Q);
        a = MatrixOperation<f8>::matrix_add(a, MatrixOperation<f8>::matrix_multiply(I, s));

        if (show_details)
            print_mat("A = ", a);

        if (show_details) {
            print_mat("a = RQ + sI = ", a);
            std::cout << std::endl;
        }

        if (fabs(a[m][m - 1]) < epsilon) {
            res.push_back(a[m][m]);
            m--;
        }

        if (m == 0) {
            res.push_back(a[m][m]);
            break;
        }
    }

    std::cout << "迭代次数: " << count << std::endl;
    for (auto i : res) {
        std::cout << "特征值: " << i << "  采用反幂法求对应的特征向量: (";
        inverse_power_method(i, 1e-5, false, norm::two_norm, true);
    }
    draw_dividing_line();
}

void convivae::MatrixEigenvaluesAndEigenvectors::QR_method_with_two_steps_shifted(double epsilon, bool show_details,
                                                                                  int max_steps) {
    std::cout << "带双步位移的 QR 方法（QR分解由平面旋转矩阵求出）" << std::endl;
    if (!_initialized)
        read_mat();

    i4 n = _dimension;
    if (show_details)
        print_mat("start:", _original_mat);

    auto A = Hessenberg_matrix(_original_mat);
    auto I = MatrixOperation<f8>::identity_matrix(n);

    //存储特征值
    vec_type res;

    if (show_details)
        print_mat("化为拟上三角矩阵:", A);

    //2.
    auto k = 0, m = n - 1;
    f8 s1, s2;

    //3.
    while (k++ < max_steps) {
        if (show_details)
            std::cout << "第" << k << "次迭代" << std::endl;
        if (fabs(A[m][m - 1]) < epsilon) {
            res.push_back(A[m][m]);
            m--;
        }
            //5.
        else {
            auto b = A[m - 1][m - 1] + A[m][m];
            auto c = A[m - 1][m - 1] * A[m][m] - A[m][m - 1] * A[m - 1][m];
            auto delta = b * b - 4 * c;
            if (delta < 0 && (m == 1 || m > 1 && fabs(A[m - 1][m - 2]) < epsilon)) {
                std::cout << "有一对复特征值: " << std::endl;
                std::cout << -b / 2.0 << "+" << sqrt(-delta) / 2.0 << "i, " << -b / 2.0 << "-" << sqrt(-delta) / 2.0
                          << "i" << std::endl;
                break;
            }

            s1 = -b / 2.0 + sqrt(delta) / 2.0;
            s2 = -b / 2.0 - sqrt(delta) / 2.0;

            //6.
            if (m == 1) {
                res.push_back(s1);
                res.push_back(s2);
                break; //转 11
            }

            //7.
            if (fabs(A[m - 1][m - 2]) < epsilon) {
                res.push_back(s1);
                res.push_back(s2);
                m -= 2;
            } else {
                //9.
                auto s = A[m - 1][m - 1] + A[m][m];
                auto t = A[m - 1][m - 1] * A[m][m] - A[m][m - 1] * A[m - 1][m];

                auto M = MatrixOperation<f8>::matrix_multiply(A, A);
                M = MatrixOperation<f8>::matrix_sub(M, MatrixOperation<f8>::matrix_multiply(A, s));
                M = MatrixOperation<f8>::matrix_add(M, MatrixOperation<f8>::matrix_multiply(I, t));

                if (show_details)
                    print_mat("B0 = M = ", M);

                mat_type B(M);

                for (auto r = 0; r < n - 1; ++r) {
                    auto k = r + 1;
                    for (; k < n; ++k) {
                        if (A[k][r] != 0)
                            break;
                    }
                    if (k >= n)
                        continue;

                    // 9.2
                    auto d = 0.0;
                    for (auto i = r; i < n; ++i) {
                        d += B[i][r] * B[i][r];
                    }
                    d = sqrt(d);

                    auto c = B[r][r] > 0.0 ? -d : d;
                    auto h = c * c - c * B[r][r];

                    // 9.3
                    vec_type u(n);
                    for (auto i = 0; i < r; ++i) {
                        u[i] = 0.0;
                    }
                    u[r] = B[r][r] - c;
                    for (auto i = r + 1; i < n; ++i) {
                        u[i] = B[i][r];
                    }

                    // 9.4
                    auto ut = MatrixOperation<f8>::matrix_transposition(u);
                    auto BT = MatrixOperation<f8>::matrix_transposition(B);
                    auto AT = MatrixOperation<f8>::matrix_transposition(A);

                    auto v = MatrixOperation<f8>::matrix_multiply(BT, ut);
                    v = MatrixOperation<f8>::matrix_multiply(v, 1.0 / h);
                    vec_type vt(n);
                    for (auto i = 0; i < n; ++i)
                        vt[i] = v[i][0];

                    B = MatrixOperation<f8>::matrix_sub(B, MatrixOperation<f8>::matrix_multiply(ut, vt));

                    auto p = MatrixOperation<f8>::matrix_multiply(AT, ut);
                    p = MatrixOperation<f8>::matrix_multiply(p, 1.0 / h);
                    vec_type pt(n);
                    for (auto i = 0; i < n; ++i)
                        pt[i] = p[i][0];

                    auto q = MatrixOperation<f8>::matrix_multiply(A, ut);
                    q = MatrixOperation<f8>::matrix_multiply(q, 1.0 / h);


                    t = MatrixOperation<f8>::matrix_multiply(pt, ut) / h;
                    auto omega = MatrixOperation<f8>::matrix_multiply(ut, t);
                    omega = MatrixOperation<f8>::matrix_sub(q, omega);

                    A = MatrixOperation<f8>::matrix_sub(A, MatrixOperation<f8>::matrix_multiply(omega, u));
                    A = MatrixOperation<f8>::matrix_sub(A, MatrixOperation<f8>::matrix_multiply(ut, pt));

                }
                if (show_details)
                    print_mat("A = ", A);
                continue;
            }
        }
        //4.
        if (m <= 0) {
            if (m == 0)
                res.push_back(A[0][0]);
            break;  //转 11
        }
    }

    if (k >= max_steps) {
        std::cout << "超过最大迭代次数，未得到全部特征值" << std::endl;
    }

    std::cout << "迭代次数: " << k << std::endl;
    for (auto i : res) {
        std::cout << "特征值: " << i << "  采用反幂法求对应的特征向量: (";
        inverse_power_method(i, 1e-5, false, norm::two_norm, true);
    }
    draw_dividing_line();
}


