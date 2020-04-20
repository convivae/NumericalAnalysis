#include "Chapter2/LinearEquations.h"
#include "Chapter3/MatrixEigenvaluesAndEigenvectors.h"

using namespace convivae;

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
//    test02.Jacobi_iteration_method(1e-4, true);
//    test02.Gauss_Seidel_method(1e-4, true);
//    test02.Successive_Over_Relaxation_method(1.25, 1e-5, true);

    /**
     * 第三章 求特征值和特征向量
     */
    auto test03 = MatrixEigenvaluesAndEigenvectors("../3.txt");
    test03.power_method(1e-5);
    test03.inverse_power_method(1e-5, false, norm::two_norm);
    test03.inverse_power_method(1.2679, 1e-5, false, norm::infinite_norm);


}