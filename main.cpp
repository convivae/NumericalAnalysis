#include "include/Chapter2/LinearEquations.h"

using namespace chapter2;

int main(int argc, char* argv[], char* env[])
{
    auto test = LinearEquations<f8>("../in.txt");

    // 高斯消元
    test.Gaussian_elimination_with_partial_pivoting_method(false);
    test.Gauss_Jordan_elimination(false);

    // 三角分解法
    test.LU_factorization_Doolittle(false);
    test.LU_factorization_Crout(false);
    test.LU_factorization_Doolittle_with_partial_pivoting_method(false);

    // 迭代法
    test.Jacobi_iteration_method(5e-5);
    test.Gauss_Seidel_method(5e-5);
    test.Successive_Over_Relaxation_method(1.25, 5e-5, true);
}