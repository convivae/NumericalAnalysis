//
// Created by convi on 2020/4/30.
//

#include "NonlinearEquations.h"

void convivae::NonlinearEquations::draw_dividing_line() {
    std::cout << "------------------------------------------------" << std::endl;
}

void convivae::NonlinearEquations::dichotomy_method(const pFun fun, double left, double right, double epsilon,
                                                    bool show_details, int max_steps) {
    std::cout << "二分法" << std::endl;
    auto k = 0;
    auto mid = 0.0;
    while (k++ < max_steps) {
        mid = (left + right) / 2;

        if (right - left < epsilon) {
            if (show_details)
                std::cout << "right - left = " << right - left << " < " << epsilon << std::endl;
            break;
        }

        if (show_details)
            std::cout << "第" << k << "次迭代:\nmid = " << mid;

        auto a = fun(left), b = fun(right), c = fun(mid);

        if (fabs(c) < 1e-10) {
            if (show_details)
                std::cout << "fabs(f(mid)) = " << fabs(c) << " < 1e-10" << std::endl;
            break;
        }

        if (show_details) {
            std::cout << ", f(left) = " << a << ", f(mid) = " << c << ", f(left) * f(mid)";
            a * c > 0 ? std::cout << " > 0" << std::endl : std::cout << " < 0" << std::endl;
        }


        if (a * c < 0)
            right = mid;
        else
            left = mid;

        if (show_details)
            std::cout << "在 (" << left << ", " << right << ") 之间" << std::endl;
    }

    std::cout << "迭代次数: " << k - 1 << std::endl;
    std::cout << "结果: x = " << mid << std::endl;
    std::cout << "验算: fun(x) = " << fun(mid) << std::endl;
    draw_dividing_line();
}

void
convivae::NonlinearEquations::simple_iteration_method(convivae::NonlinearEquations::pFun fun, double x0, double eta,
                                                      bool show_details, int max_steps) {
    std::cout << "简单迭代法" << std::endl;
    auto k = 0;
    auto xk = x0;

    if (show_details)
        std::cout << "x0 = " << x0 << std::endl;

    while (k++ < max_steps) {
        xk = fun(x0);

        if (show_details)
            std::cout << "x" << k << " = " << xk << std::endl;

        if (fabs(xk - x0) / fabs(xk) <= eta) {
            if (show_details)
                std::cout << "|x" << k << " - x" << k - 1 << "| / |x" << k << "| = " << fabs(xk - x0) / fabs(xk)
                          << " < " << eta << std::endl;
            break;
        }

        x0 = xk;
    }

    std::cout << "迭代次数: " << k << std::endl;
    std::cout << "结果: x = " << xk << std::endl;

    draw_dividing_line();
}

void
convivae::NonlinearEquations::Steffensen_iteration_method(convivae::NonlinearEquations::pFun fun, double x0, double eta,
                                                          bool show_details, int max_steps) {
    std::cout << "Steffensen 迭代法" << std::endl;
    auto k = 0;
    auto xk = x0;

    if (show_details)
        std::cout << "x0 = " << x0 << std::endl;

    while (k++ < max_steps) {
        auto y = fun(x0), z = fun(y);
        xk = x0 - (y - x0) * (y - x0) / (z - 2 * y + x0);

        if (show_details)
            std::cout << "x" << k << " = " << xk << std::endl;

        if (fabs(xk - x0) / fabs(xk) <= eta) {
            if (show_details)
                std::cout << "|x" << k << " - x" << k - 1 << "| / |x" << k << "| = " << fabs(xk - x0) / fabs(xk)
                          << " < " << eta << std::endl;
            break;
        }

        x0 = xk;
    }

    std::cout << "迭代次数: " << k << std::endl;
    std::cout << "结果: x = " << xk << std::endl;
    draw_dividing_line();
}

void convivae::NonlinearEquations::Newton_iteration_method(convivae::NonlinearEquations::pFun fun,
                                                           convivae::NonlinearEquations::pFun fun_1, double eta,
                                                           double x0, bool show_details, bool down_hill,
                                                           int max_steps) {
    // 下山法确定初始值
    if (down_hill) {
        auto lambda = 1.0;
        auto y = x0 - lambda * fun(x0) / fun_1(x0);
        while (fabs(fun(y)) > fabs(fun(x0))) {
            lambda /= 2;
            y = x0 - lambda * fun(x0) / fun_1(x0);
        }
        x0 = y;
    }

    std::cout << "Newton 迭代法" << std::endl;
    auto k = 0;
    auto xk = x0;

    if (show_details)
        std::cout << "x0 = " << x0 << std::endl;

    while (k++ < max_steps) {
        xk = x0 - fun(x0) / fun_1(x0);

        if (show_details)
            std::cout << "x" << k << " = " << xk << std::endl;

        if (fabs(xk - x0) / fabs(xk) <= eta) {
            if (show_details)
                std::cout << "|x" << k << " - x" << k - 1 << "| / |x" << k << "| = " << fabs(xk - x0) / fabs(xk)
                          << " < " << eta << std::endl;
            break;
        }

        x0 = xk;
    }

    std::cout << "迭代次数: " << k << std::endl;
    std::cout << "结果: x = " << xk << std::endl;
    draw_dividing_line();
}

void convivae::NonlinearEquations::secant_iteration_method(convivae::NonlinearEquations::pFun fun, double x0, double x1,
                                                           double eta, bool show_details, int max_steps) {
    std::cout << "割线法" << std::endl;
    auto k = 1;
    auto xk = x0;

    if (show_details)
        std::cout << "x0 = " << x0 << ", x1 = " << x1 << std::endl;

    while (k++ < max_steps) {
        xk = x1 - (fun(x1) * (x1 - x0)) / (fun(x1) - fun(x0));

        if (show_details)
            std::cout << "x" << k << " = " << xk << std::endl;

        if (fabs(xk - x1) / fabs(xk) <= eta) {
            if (show_details)
                std::cout << "|x" << k << " - x" << k - 1 << "| / |x" << k << "| = " << fabs(xk - x1) / fabs(xk)
                          << " < " << eta << std::endl;
            break;
        }

        x0 = x1;
        x1 = xk;
    }

    std::cout << "迭代次数: " << k - 1 << std::endl;
    std::cout << "结果: x = " << xk << std::endl;
    draw_dividing_line();
}

void convivae::NonlinearEquations::secant_with_single_point_iteration_method(convivae::NonlinearEquations::pFun fun,
                                                                             double x0, double x1, double eta,
                                                                             bool show_details, int max_steps) {
    std::cout << "单点割线法" << std::endl;
    auto k = 1;
    auto xk = x0, f_x0 = fun(x0);

    if (show_details)
        std::cout << "x0 = " << x0 << ", x1 = " << x1 << std::endl;

    while (k++ < max_steps) {
        xk = x1 - (fun(x1) * (x1 - x0)) / (fun(x1) - f_x0);

        if (show_details)
            std::cout << "x" << k << " = " << xk << std::endl;

        if (fabs(xk - x1) / fabs(xk) <= eta) {
            if (show_details)
                std::cout << "|x" << k << " - x" << k - 1 << "| / |x" << k << "| = " << fabs(xk - x1) / fabs(xk)
                          << " < " << eta << std::endl;
            break;
        }

        x1 = xk;
    }

    std::cout << "迭代次数: " << k - 1 << std::endl;
    std::cout << "结果: x = " << xk << std::endl;
    draw_dividing_line();
}

void convivae::NonlinearEquations::simple_iteration_equations_method(convivae::NonlinearEquations::ppFun fun1,
                                                                     convivae::NonlinearEquations::ppFun fun2,
                                                                     double x1_0, double x2_0, double eta,
                                                                     bool show_details, int max_steps) {
    std::cout << "非线性方程组的简单迭代法" << std::endl;
    auto k = 0;

    if (show_details)
        std::cout << "x1_0 = " << x1_0 << ", x2_0 = " << x2_0 << std::endl;

    auto x1 = 0.0, x2 = 0.0;

    while (k++ < max_steps) {
        x1 = fun1(x1_0, x2_0);
        x2 = fun2(x1_0, x2_0);

        if (show_details)
            std::cout << "x1_" << k << " = " << x1 << ", x2_" << k << " = " << x2 << std::endl;

        if (fabs(x1 - x1_0) / fabs(x1) <= eta && fabs(x2 - x2_0) / fabs(x2) <= eta) {
            if (show_details) {
                std::cout << "|x1_" << k << " - x1_" << k - 1 << "| / |x1_" << k << "| = " << fabs(x1 - x1_0) / fabs(x1)
                          << " < " << eta << std::endl;
                std::cout << "|x2_" << k << " - x2_" << k - 1 << "| / |x2_" << k << "| = " << fabs(x2 - x2_0) / fabs(x2)
                          << " < " << eta << std::endl;
                break;
            }
        }

        x1_0 = x1;
        x2_0 = x2;
    }

    std::cout << "迭代次数: " << k << std::endl;
    std::cout << "结果: x1 = " << x1 << ", x2 = " << x2 << std::endl;
    draw_dividing_line();

}

void convivae::NonlinearEquations::Newton_iteration_equations_method(convivae::NonlinearEquations::ppFun fun1,
                                                                     convivae::NonlinearEquations::ppFun fun2,
                                                                     double x1_0, double x2_0, double eta,
                                                                     bool show_details, int max_steps) {
    std::cout << "非线性方程组的 Newton 迭代法" << std::endl;
    auto k = 0;

    if (show_details)
        std::cout << "x1_0 = " << x1_0 << ", x2_0 = " << x2_0 << std::endl;

    auto x1 = 0.0, x2 = 0.0;

    while (k++ < max_steps) {
        x1 = fun1(x1_0, x2_0);
        x2 = fun2(x1_0, x2_0);

        if (show_details)
            std::cout << "x1_" << k << " = " << x1 << ", x2_" << k << " = " << x2 << std::endl;

        if (fabs(x1 - x1_0) / fabs(x1) <= eta && fabs(x2 - x2_0) / fabs(x2) <= eta) {
            if (show_details) {
                std::cout << "|x1_" << k << " - x1_" << k - 1 << "| / |x1_" << k << "| = " << fabs(x1 - x1_0) / fabs(x1)
                          << " < " << eta << std::endl;
                std::cout << "|x2_" << k << " - x2_" << k - 1 << "| / |x2_" << k << "| = " << fabs(x2 - x2_0) / fabs(x2)
                          << " < " << eta << std::endl;
                break;
            }
        }

        x1_0 = x1;
        x2_0 = x2;
    }

    std::cout << "迭代次数: " << k << std::endl;
    std::cout << "结果: x1 = " << x1 << ", x2 = " << x2 << std::endl;
    draw_dividing_line();
}
