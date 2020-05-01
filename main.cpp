#include "Chapter2/LinearEquations.h"
#include "Chapter3/MatrixEigenvaluesAndEigenvectors.h"
#include "Chapter4/NonlinearEquations.h"

constexpr auto _e = 2.718281828459;

using namespace convivae;

//二分法的迭代函数
double fun_dichotomy(double x) {
    return x + sin(x) - 1;
}

//简单迭代法的迭代函数
double fun_simple(double x) {
    return 2 * acos(sqrt((pow(_e, -x) + 1) / 2));
}

//Steffensen 法迭代函数
double fun_Steffensen(double x) {
    return (2 - pow(_e, x)) / 10;
}

//Newton 法迭代函数
double fun_Newton(double x) {
    return x + sin(x) - 1;
}

double fun_Newton_1(double x) {
    return 1 + cos(x);
}

// 割线法
double fun_secant(double x) {
    return x + sin(x) - 1;
}

int main(int argc, char *argv[], char *env[]) {
    /**
     * 第二章 解线性方程组
     */
    //    auto test02 = LinearEquations<f8>("../2.txt");
    //
    //    // 高斯消元
    //    test02.Gaussian_elimination_with_partial_pivoting_method();
    //    test02.Gauss_Jordan_elimination();
    //
    //    // 三角分解法
    //    test02.LU_factorization_Doolittle(false);
    //    test02.LU_factorization_Crout(false);
    //    test02.LU_factorization_Doolittle_with_partial_pivoting_method(false);
    //
    //    // 迭代法
    //    test02.Jacobi_iteration_method(1e-4, false);
    //    test02.Gauss_Seidel_method(1e-4, false);
    //    test02.Successive_Over_Relaxation_method(1.25, 1e-5, false);

    /**
     * 第三章 求特征值和特征向量
     */
    //    auto test03 = MatrixEigenvaluesAndEigenvectors("../3.txt");
    //
    //    // 幂法与反幂法
    //    std::vector<std::vector<f8>> A, Q, R;
    //    test03.QR_decomposition_Householder(A, Q, R, true);
    //
    //    test03.power_method(1e-5);
    //    test03.inverse_power_method(norm::two_norm, 1e-5, false);
    //    test03.inverse_power_method(20.9681);
    //    test03.Jacobi_method(1e-5, false);
    //
    //    test03.QR_method_normal(1e-5, false);
    //    test03.QR_method_Transformed_Rotational_Matrix(1e-5, false);
    //    test03.QR_method_with_one_step_shifted(1e-5, false);
    //    test03.QR_method_with_two_steps_shifted(1e-5, false);

    /**
     * 第四章 非线性方程与非线性方程组的迭代解法
     */
    auto *test04 = new NonlinearEquations();
    //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(9);

    //传递函数指针
    test04->dichotomy_method(fun_dichotomy, 0, 1.6, 1e-3, false);
    test04->simple_iteration_method(fun_simple, 1, 1e-6, false);
    test04->Steffensen_iteration_method(fun_Steffensen, 0, 1e-6, false);
    test04->Newton_iteration_method(fun_Newton, fun_Newton_1, 1e-6, 0.6, false, true);
    test04->secant_iteration_method(fun_secant, 0, 1.6, 1e-6, false);
    test04->secant_with_single_point_iteration_method(fun_secant, 0, 1.6, 1e-6, false);
}
