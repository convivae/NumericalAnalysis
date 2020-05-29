//
// Created by convi on 2020/5/29.
//

#include "NumericalIntegration.h"

void NumericalIntegration::draw_dividing_line() {
    cout << "-----------------------------------------------------" << endl;
}

int NumericalIntegration::my_pow(int a, int b) {
    int res = 1;
    for (auto i = 0; i < b; ++i) {
        res *= a;
    }
    return res;
}

double NumericalIntegration::my_fac(int n) {
    double res = 1.0;
    for (auto i = 2; i <= n; ++i) {
        res *= i;
    }
    return res;
}

void NumericalIntegration::Romberg_numerical_integration_method(NumericalIntegration::pFun fun, double down, double up,
                                                                double epsilon, int max_steps) {
    cout << "Romberg 积分法" << endl;
    if (up < down) {
        cerr << "error! up less than down!" << endl;
        return;
    }
    vector<vector<double>> T(4);
    for (auto &i : T)
        i.resize(max_steps);
    int steps = 0;

    T[0][0] = (up - down) / 2.0 * (fun(up) + fun(down));
    for (auto m = 1; m < max_steps; ++m) {
        auto h = (up - down) / pow(2, m);

        auto tmp = 0.0;
        int tmp_up = my_pow(2, m - 1);
        for (auto i = 1; i <= tmp_up; ++i) {
            tmp += fun(down + (2.0 * i - 1.0) * h);
        }

        T[0][m] = 0.5 * T[0][m - 1] + h * tmp;

        int tmp_end = m < 3 ? m : 3;
        for (auto j = 1; j <= tmp_end; ++j) {
            T[j][m - j] = (my_pow(4, j) * T[j - 1][m - j + 1] - T[j - 1][m - j]) / (my_pow(4, j) - 1);
        }
        if (abs(T[3][1]) > 1e-5) {
            steps = m;
            if (abs(T[3][m - 3] - T[3][m - 4]) / abs(T[3][m - 3]) <= epsilon)
                break;
        }
    }

    cout << "由 Romberg 积分公式(P162)，计算结果为:" << endl;

    for (auto i = 0; i < 4; ++i) {
        printf("%8s", "T_m^");
        printf("(%d)\t", i);
    }

    cout << endl;

    for (auto i = 0; i < steps + 1; i++) {
        for (auto j = 0; j < i + 1; ++j) {
            if (j > 3)
                break;
            printf("%.9f\t", T[j][i - j]);
        }
        cout << endl;
    }


    draw_dividing_line();
}


void NumericalIntegration::Gauss_Legendre_integration_method(NumericalIntegration::pFun fun, int n, ...) {
    cout << "利用 Gauss_Legendre 求积公式计算" << endl;
    va_list ap;
    va_start(ap, n);
    auto L = OrthogonalPolynomial::Legendre_orthogonal_polynomial(n);
    double res = 0.0;

    for (auto i = 1; i <= n; ++i) {
        double x = va_arg(ap, double);
        double tmp = L.derivative().value_at_point(x);
        double A = 2.0 / ((1 - x * x) * tmp * tmp);
        cout << "A" << i << "=";
        printf("%.9f\t", A);
        cout << "x" << i << "=";
        printf("%.9f", x);
        cout << endl;
        res += A * fun(x);
    }
    va_end(ap);

    cout << "结果为:" << endl;
    printf("%.9f", res);
    cout << endl;
    draw_dividing_line();
}

void NumericalIntegration::Gauss_Laguerre_integration_method(NumericalIntegration::pFun fun, int n, ...) {
    cout << "利用 Gauss_Laguerre 求积公式计算" << endl;
    va_list ap;
    va_start(ap, n);
    auto L = OrthogonalPolynomial::Laguerre_orthogonal_polynomial(n);
    auto res = 0.0;

    for (auto i = 1; i <= n; ++i) {
        double x = va_arg(ap, double);
        double a = my_fac(n);
        double b = L.derivative().value_at_point(x);

        auto A = a * a / (x * b * b);
        cout << "A" << i << "=";
        printf("%.9f\t", A);
        cout << "x" << i << "=";
        printf("%.9f", x);
        cout << endl;
        res += A * fun(x);
    }
    va_end(ap);

    cout << "结果为:" << endl;
    printf("%.9f", res);
    cout << endl;
    draw_dividing_line();
}

void NumericalIntegration::Gauss_Hermite_integration_method(NumericalIntegration::pFun fun, int n, ...) {
    cout << "利用 Gauss_Hermite 求积公式计算" << endl;
    va_list ap;
    va_start(ap, n);
    auto L = OrthogonalPolynomial::Hermite_orthogonal_polynomial(n);
    auto res = 0.0;

    for (auto i = 1; i <= n; ++i) {
        double x = va_arg(ap, double);
        double a = pow(2,n+1) * my_fac(n) * sqrt(M_PI);
        double b = L.derivative().value_at_point(x);

        auto A = a / (b * b);
        cout << "A" << i << "=";
        printf("%.9f\t", A);
        cout << "x" << i << "=";
        printf("%.9f", x);
        cout << endl;
        res += A * fun(x);
    }
    va_end(ap);

    cout << "结果为:" << endl;
    printf("%.9f", res);
    cout << endl;
    draw_dividing_line();
}

void NumericalIntegration::Gauss_Chebyshev_integration_method(NumericalIntegration::pFun fun, int n, ...) {
    cout << "利用 Gauss_Chebyshev 求积公式计算" << endl;
    va_list ap;
    va_start(ap, n);
    auto L = OrthogonalPolynomial::Chebyshev_orthogonal_polynomial(n);
    auto res = 0.0;
    auto A = M_PI / n;

    for (auto i = 1; i <= n; ++i) {
        auto x = va_arg(ap, double);
        cout << "A" << i << "=" ;
        printf("%.9f\t", A);
        cout << "x" << i << "=";
        printf("%.9f", x);
        cout << endl;
        res += A * fun(x);
    }
    va_end(ap);

    cout << "结果为:" << endl;
    printf("%.9f", res);
    cout << endl;
    draw_dividing_line();
}


