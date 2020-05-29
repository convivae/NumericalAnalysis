//
// Created by convi on 2020/5/12.
//

#include "InterpolationAndApproximation.h"

void convivae::InterpolationAndApproximation::draw_dividing_line() {
    cout << "----------------------------------------------------------------------------------" << endl;
}

void convivae::InterpolationAndApproximation::print_vec(const string &s, vector<double> v, int start_pos) {
    cout << s << ":" << endl;
    for (auto i = start_pos; i < v.size(); ++i) {
        cout << s << i << "=" << v[i] << "\t";
    }
    cout << endl;
}

void convivae::InterpolationAndApproximation::print_mat(const string &s, vector<vector<double>> v) {
    cout << s << ":" << endl;
    for (const auto &i : v) {
        for (auto j :i) {
            printf("%.9f ", j);
        }
        cout << endl;
    }
}

void convivae::InterpolationAndApproximation::read_data() {
    if (_initialized)
        return;

    ifstream in(_rdr);
    string s;

    getline(in, s);
    istringstream sx(s);
    for (f8 i; sx >> i;) {
        this->X.push_back(i);
    }

    getline(in, s);
    istringstream sy(s);
    for (f8 i; sy >> i;) {
        this->Y.push_back(i);
    }

    in.close();

    if (this->X.size() != this->Y.size()) {
        cerr << "data error" << endl;
        exit(-1);
    }

    this->_initialized = true;
    this->_n = this->X.size();
}

Polymerization convivae::InterpolationAndApproximation::Lagrange_polynomial(int n) {
    cout << "Lagrange 插值法" << endl;
    read_data();
    if (n > _n) {
        cerr << "error" << endl;
        exit(-1);
    }
    Polymerization res("0");
    string s;
    stringstream ss;
    for (auto k = 0; k <= n; ++k) {
        Polymerization tmp("1");
        for (auto j = 0; j <= n; ++j) {
            if (j == k)
                continue;
            double a = 1 / (X[k] - X[j]);
            double b = X[j] == 0 ? 0.0 : -X[j] / (X[k] - X[j]);
            ss << a << "x";
            b < 0.0 ? ss << b : ss << "+" << b;
            ss >> s;
            ss.clear();
            ss.str("");
            tmp *= Polymerization(s);
        }
        res += tmp * Y[k];
    }

    cout << res.to_string() << endl;
    draw_dividing_line();
    return res;
}

Polymerization convivae::InterpolationAndApproximation::Newton_polynomial(int n) {
    cout << "Newton 插值法" << endl;
    read_data();
    if (n > _n) {
        cerr << "error" << endl;
        exit(-1);
    }

    vector<vector<double>> table(n + 1, vector<double>(n + 2, 0));
    vector<vector<bool>> flag(n + 1, vector<bool>(n + 2, false));
    for (auto i = 0; i <= n; ++i) {
        table[i][0] = X[i];
        table[i][1] = Y[i];
        flag[i][0] = true;
        flag[i][1] = true;
    }

    for (auto i = 1; i <= n; ++i) {
        for (auto j = 2; j <= i + 1; ++j) {
            table[i][j] = differential_quotient(i - j + 1, i - 1, i, table, flag);
        }
    }

    cout << "table:" << endl;
    cout << "xi\tf[xi]\t";
    for (auto i = 0; i < n; ++i) {
        cout << "No." << i + 1 << "\t";
    }
    cout << endl;

    for (auto i = 0; i <= n; ++i) {
        for (auto j = 0; j < i + 2; ++j) {
            cout << table[i][j] << "\t";
        }
        cout << endl;
    }

    cout << "The Newton polynomial is:" << endl;
    Polymerization res(to_string(table[0][1]));
    Polymerization p("1");

    for (auto i = 1; i <= n; ++i) {
        auto s = to_string(-table[i - 1][0]);
        if (-table[i - 1][0] >= 0)
            s = string("x+").append(s);
        else
            s = string("x").append(s);
        p *= Polymerization(s);
        res += p * table[i][i + 1];

    }

    cout << res.to_string() << endl;
    draw_dividing_line();
    return res;
}

double
convivae::InterpolationAndApproximation::differential_quotient(int start_index, int before_end_index, int end_index,
                                                               vector<vector<double>> &table,
                                                               vector<vector<bool>> &flag) {
    if (start_index == end_index) { // f(xk) 零阶
        table[start_index][1] = Y[start_index];
        flag[start_index][1] = true;
        return Y[start_index];
    }
    if (before_end_index == start_index) { //f(x_start, x_end) 一阶
        double res = (Y[end_index] - Y[start_index]) / (X[end_index] - X[start_index]);
        if (end_index - before_end_index == 1) {
            table[end_index][2] = res;
            flag[end_index][2] = true;
        }
        return res;
    }

    // f(x0, x1, ... , x_(k-1), xk) k阶
    double res;
    int row = end_index, column = before_end_index - start_index + 2;
    if (end_index - before_end_index == 1 && flag[row][column]) {
        res = table[row][column];
    } else {
        double a = differential_quotient(start_index, before_end_index - 1, end_index, table, flag);
        double b = differential_quotient(start_index, before_end_index - 1, before_end_index, table, flag);
        res = (a - b) / (X[end_index] - X[before_end_index]);
        if (end_index - before_end_index == 1) {
            table[row][column] = res;
            flag[row][column] = true;
        }
    }

    return res;
}

void convivae::InterpolationAndApproximation::three_moment_method_bound1(double f2_0, double f2_n, bool simplify) {
    cout << "三弯矩法求三次样条插值函数_第一种边界条件" << endl;
    read_data();
    int n = X.size() - 1;

    vec_type h(n + 1, U4_INF), alpha(n, U4_INF), beta(n + 1, U4_INF), gama(n + 1, U4_INF);
    for (auto i = 1; i <= n; ++i) {
        h[i] = X[i] - X[i - 1];
    }

    alpha[0] = 0.0;
    gama[n] = 0.0;
    beta[0] = 2 * f2_0;
    beta[n] = 2 * f2_n;
    for (auto i = 1; i <= n - 1; i++) {
        alpha[i] = h[i + 1] / (h[i] + h[i + 1]);
        gama[i] = 1 - alpha[i];
        beta[i] = 6.0 / (h[i] + h[i + 1]) * ((Y[i + 1] - Y[i]) / h[i + 1] - (Y[i] - Y[i - 1]) / h[i]);
    }

    print_vec("alpha", alpha);
    print_vec("beta", beta);
    print_vec("gama", gama, 1);

    //构造增广矩阵
    mat_type mat(n + 1, vec_type(n + 2));
    for (auto i = 0; i < n + 1; i++) {
        mat[i][i] = 2.0;
        mat[i][n + 1] = beta[i];
        if (i != n)
            mat[i][i + 1] = alpha[i];
        if (i != 0)
            mat[i][i - 1] = gama[i];
    }

    print_mat("增广矩阵", mat);


    n++;
    mat_type a(mat);
    vec_type b(beta);
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

        for (auto i = k + 1; i < n; ++i) {
            auto sum = 0.0;
            for (auto t = 0; t < k; ++t) {
                sum += combine_LU[i][t] * combine_LU[t][k];
            }
            combine_LU[i][k] = (a[i][k] - sum) / combine_LU[k][k];
        }
    }

    //计算 y、x
    for (auto i = 0; i < n; ++i) {
        auto sum = 0.0;
        for (auto t = 0; t < i; ++t) {
            sum += combine_LU[i][t] * y[t];
        }
        y[i] = b[i] - sum;
    }

    for (auto i = n - 1; i >= 0; --i) {
        auto sum = 0.0;
        for (auto t = i + 1; t < n; ++t) {
            sum += combine_LU[i][t] * x[t];
        }
        x[i] = (y[i] - sum) / combine_LU[i][i];
    }
    print_vec("M", x);
    vec_type M(x);

    //s(x)
    n--;
    string s;
    for (auto i = 1; i <= n; i++) {
        Polymerization p("0");
        stringstream ss;
        auto tmp = M[i - 1] / (6 * h[i]);
        ss << tmp << "(" << X[i] << "-x)^3";
        Polymerization p_tmp1(to_string(X[i]) + "-x");
        p += p_tmp1 * p_tmp1 * p_tmp1 * tmp;

        tmp = M[i] / (6 * h[i]);
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp2(to_string(X[i - 1]) + "-x");
        p += p_tmp2 * p_tmp2 * p_tmp2 * -tmp;

        ss << tmp << "(" << "x";
        if (X[i - 1] >= 0)
            ss << "-";
        ss << X[i - 1] << ")^3";

        tmp = Y[i - 1] / h[i] - M[i - 1] * h[i] / 6;
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp3(to_string(X[i]) + "-x");
        p += p_tmp3 * tmp;

        ss << tmp << "(" << X[i] << "-x)";

        tmp = Y[i] / h[i] - M[i] * h[i] / 6;
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp4(to_string(X[i - 1]) + "-x");
        p += p_tmp4 * -tmp;

        ss << tmp << "(" << "x";
        if (X[i - 1] >= 0)
            ss << "-";
        ss << X[i - 1] << ")";

        ss >> s;
        cout << s << " \t" << X[i - 1] << "=<x<=" << X[i];
        if (simplify)
            cout << " simplify: " << p.to_string();
        cout << endl;
    }
    draw_dividing_line();
}

void convivae::InterpolationAndApproximation::three_moment_method_bound2(double f1_0, double f1_n, bool simplify) {
    cout << "三弯矩法求三次样条插值函数_第二种边界条件" << endl;
    read_data();
    int n = X.size() - 1;

    vec_type h(n + 1, U4_INF), alpha(n, U4_INF), beta(n + 1, U4_INF), gama(n + 1, U4_INF);
    for (auto i = 1; i <= n; ++i) {
        h[i] = X[i] - X[i - 1];
    }

    alpha[0] = 1.0;
    gama[n] = 1.0;
    beta[0] = 6 / h[1] * ((Y[1] - Y[0]) / h[1] - f1_0);
    beta[n] = 6 / h[n] * (f1_n - (Y[n] - Y[n - 1]) / h[n]);

    for (auto i = 1; i <= n - 1; i++) {
        alpha[i] = h[i + 1] / (h[i] + h[i + 1]);
        gama[i] = 1 - alpha[i];
        beta[i] = 6.0 / (h[i] + h[i + 1]) * ((Y[i + 1] - Y[i]) / h[i + 1] - (Y[i] - Y[i - 1]) / h[i]);
    }

    print_vec("alpha", alpha);
    print_vec("beta", beta);
    print_vec("gama", gama, 1);

    //构造增广矩阵
    mat_type mat(n + 1, vec_type(n + 2));
    for (auto i = 0; i < n + 1; i++) {
        mat[i][i] = 2.0;
        mat[i][n + 1] = beta[i];
        if (i != n)
            mat[i][i + 1] = alpha[i];
        if (i != 0)
            mat[i][i - 1] = gama[i];
    }

    print_mat("增广矩阵", mat);


    n++;
    mat_type a(mat);
    vec_type b(beta);
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

        for (auto i = k + 1; i < n; ++i) {
            auto sum = 0.0;
            for (auto t = 0; t < k; ++t) {
                sum += combine_LU[i][t] * combine_LU[t][k];
            }
            combine_LU[i][k] = (a[i][k] - sum) / combine_LU[k][k];
        }
    }

    //计算 y、x
    for (auto i = 0; i < n; ++i) {
        auto sum = 0.0;
        for (auto t = 0; t < i; ++t) {
            sum += combine_LU[i][t] * y[t];
        }
        y[i] = b[i] - sum;
    }

    for (auto i = n - 1; i >= 0; --i) {
        auto sum = 0.0;
        for (auto t = i + 1; t < n; ++t) {
            sum += combine_LU[i][t] * x[t];
        }
        x[i] = (y[i] - sum) / combine_LU[i][i];
    }
    print_vec("M", x);
    vec_type M(x);

    //s(x)
    n--;
    string s;
    for (auto i = 1; i <= n; i++) {
        Polymerization p("0");
        stringstream ss;
        auto tmp = M[i - 1] / (6 * h[i]);
        ss << tmp << "(" << X[i] << "-x)^3";
        Polymerization p_tmp1(to_string(X[i]) + "-x");
        p += p_tmp1 * p_tmp1 * p_tmp1 * tmp;

        tmp = M[i] / (6 * h[i]);
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp2(to_string(X[i - 1]) + "-x");
        p += p_tmp2 * p_tmp2 * p_tmp2 * -tmp;

        ss << tmp << "(" << "x";
        if (X[i - 1] >= 0)
            ss << "-";
        ss << X[i - 1] << ")^3";

        tmp = Y[i - 1] / h[i] - M[i - 1] * h[i] / 6;
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp3(to_string(X[i]) + "-x");
        p += p_tmp3 * tmp;

        ss << tmp << "(" << X[i] << "-x)";

        tmp = Y[i] / h[i] - M[i] * h[i] / 6;
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp4(to_string(X[i - 1]) + "-x");
        p += p_tmp4 * -tmp;

        ss << tmp << "(" << "x";
        if (X[i - 1] >= 0)
            ss << "-";
        ss << X[i - 1] << ")";

        ss >> s;
        cout << s << " \t" << X[i - 1] << "=<x<=" << X[i];
        if (simplify)
            cout << " simplify: " << p.to_string();
        cout << endl;
    }
    draw_dividing_line();
}

void convivae::InterpolationAndApproximation::three_moment_method_bound3(bool simplify) {
    cout << "三弯矩法求三次样条插值函数_第三种边界条件" << endl;
    read_data();
    int n = X.size() - 1;

    vec_type h(n + 1, U4_INF), alpha(n + 1, U4_INF), beta(n + 1, U4_INF), gama(n + 1, U4_INF);
    for (auto i = 1; i <= n; ++i) {
        h[i] = X[i] - X[i - 1];
    }

    alpha[n] = h[1] / (h[1] + h[n]);
    gama[n] = 1.0 - alpha[n];
    beta[n] = 6 / (h[1] + h[n]) * ((Y[1] - Y[0]) / h[1] - (Y[n] - Y[n - 1]) / h[n]);

    for (auto i = 1; i <= n - 1; i++) {
        alpha[i] = h[i + 1] / (h[i] + h[i + 1]);
        gama[i] = 1 - alpha[i];
        beta[i] = 6.0 / (h[i] + h[i + 1]) * ((Y[i + 1] - Y[i]) / h[i + 1] - (Y[i] - Y[i - 1]) / h[i]);
    }

    print_vec("alpha", alpha, 1);
    print_vec("beta", beta, 1);
    print_vec("gama", gama, 1);

    //构造增广矩阵
    mat_type mat(n, vec_type(n + 1));
    for (auto i = 0; i < n; i++) {
        mat[i][i] = 2.0;
        mat[i][n] = beta[i + 1];
        if (i != n - 1)
            mat[i][i + 1] = alpha[i + 1];
        if (i != 0)
            mat[i][i - 1] = gama[i + 1];
    }
    mat[0][n - 1] = gama[1];
    mat[n - 1][0] = alpha[n];

    print_mat("增广矩阵", mat);


    mat_type a(mat);
    vec_type b(beta);
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

        for (auto i = k + 1; i < n; ++i) {
            auto sum = 0.0;
            for (auto t = 0; t < k; ++t) {
                sum += combine_LU[i][t] * combine_LU[t][k];
            }
            combine_LU[i][k] = (a[i][k] - sum) / combine_LU[k][k];
        }
    }

    //计算 y、x
    for (auto i = 0; i < n; ++i) {
        auto sum = 0.0;
        for (auto t = 0; t < i; ++t) {
            sum += combine_LU[i][t] * y[t];
        }
        y[i] = b[i] - sum;
    }

    for (auto i = n - 1; i >= 0; --i) {
        auto sum = 0.0;
        for (auto t = i + 1; t < n; ++t) {
            sum += combine_LU[i][t] * x[t];
        }
        x[i] = (y[i] - sum) / combine_LU[i][i];
    }
    print_vec("M", x);
    vec_type M(x);

    //s(x)
    string s;
    for (auto i = 1; i <= n; i++) {
        Polymerization p("0");
        stringstream ss;
        auto tmp = M[i - 1] / (6 * h[i]);
        ss << tmp << "(" << X[i] << "-x)^3";
        Polymerization p_tmp1(to_string(X[i]) + "-x");
        p += p_tmp1 * p_tmp1 * p_tmp1 * tmp;

        tmp = M[i] / (6 * h[i]);
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp2(to_string(X[i - 1]) + "-x");
        p += p_tmp2 * p_tmp2 * p_tmp2 * -tmp;

        ss << tmp << "(" << "x";
        if (X[i - 1] >= 0)
            ss << "-";
        ss << X[i - 1] << ")^3";

        tmp = Y[i - 1] / h[i] - M[i - 1] * h[i] / 6;
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp3(to_string(X[i]) + "-x");
        p += p_tmp3 * tmp;

        ss << tmp << "(" << X[i] << "-x)";

        tmp = Y[i] / h[i] - M[i] * h[i] / 6;
        if (tmp >= 0)
            ss << "+";

        Polymerization p_tmp4(to_string(X[i - 1]) + "-x");
        p += p_tmp4 * -tmp;

        ss << tmp << "(" << "x";
        if (X[i - 1] >= 0)
            ss << "-";
        ss << X[i - 1] << ")";

        ss >> s;
        cout << s << " \t" << X[i - 1] << "=<x<=" << X[i];
        if (simplify)
            cout << " simplify: " << p.to_string();
        cout << endl;
    }
    draw_dividing_line();
}







