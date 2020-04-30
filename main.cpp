#include "Chapter2/LinearEquations.h"
#include "Chapter3/MatrixEigenvaluesAndEigenvectors.h"

using namespace convivae;

int main(int argc, char *argv[], char *env[]) {

    /**
     * 第二章 解线性方程组
     */
    auto test02 = LinearEquations<f8>("../2.txt");

    // 高斯消元
    test02.Gaussian_elimination_with_partial_pivoting_method();
    test02.Gauss_Jordan_elimination();

    // 三角分解法
    test02.LU_factorization_Doolittle(false);
    test02.LU_factorization_Crout(false);
    test02.LU_factorization_Doolittle_with_partial_pivoting_method(false);

    // 迭代法
    test02.Jacobi_iteration_method(1e-4, false);
    test02.Gauss_Seidel_method(1e-4, false);
    test02.Successive_Over_Relaxation_method(1.25, 1e-5, false);

    /**
     * 第三章 求特征值和特征向量
     */
    auto test03 = MatrixEigenvaluesAndEigenvectors("../3.txt");

    // 幂法与反幂法
    std::vector<std::vector<f8>> A, Q, R;
    test03.QR_decomposition_Householder(A, Q, R, true);

    test03.power_method(1e-5);
    test03.inverse_power_method(norm::two_norm, 1e-5, false);
    test03.inverse_power_method(20.9681);
    test03.Jacobi_method(1e-5, false);

    test03.QR_method_normal(1e-5, false);
    test03.QR_method_Transformed_Rotational_Matrix(1e-5, false);
    test03.QR_method_with_one_step_shifted(1e-5, false);
    test03.QR_method_with_two_steps_shifted(1e-5, false);

}