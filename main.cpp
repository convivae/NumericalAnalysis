#include "Chapter2/LinearEquations.h"

using namespace chapter2;

int main() {
    auto test = LinearEquations<f8>("../in.txt");

    // 高斯消元
    test.Gaussian_elimination_with_partial_pivoting_method(true);
    test.Gauss_Jordan_elimination();

    // 三角分解法
    test.LU_factorization_Doolittle(true);
    test.LU_factorization_Crout();
    test.LU_factorization_Doolittle_with_partial_pivoting_method();
}