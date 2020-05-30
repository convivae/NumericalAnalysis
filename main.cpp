#include "lib/fun.h"
#include "lib/OrthogonalPolynomial.h"
#include "Chapter2/LinearEquations.h"
#include "Chapter3/MatrixEigenvaluesAndEigenvectors.h"
#include "Chapter4/NonlinearEquations.h"
#include "Chapter5/InterpolationAndApproximation.h"
#include "Chapter6/NumericalIntegration.h"

using namespace std;
using namespace convivae;

int main(int argc, char* argv[], char* env[])
{
	/**
	 * 第二章 解线性方程组    
	 */
	auto test02 = LinearEquations<f8>("../2.txt");

	// 高斯消元
	test02.Gaussian_elimination_with_partial_pivoting_method(false);
	test02.Gauss_Jordan_elimination();

	// 三角分解法
	test02.LU_factorization_Doolittle(false);
	test02.LU_factorization_Crout(false);
	test02.LU_factorization_Doolittle_with_partial_pivoting_method(false);

	// 迭代法
	test02.Jacobi_iteration_method(1e-4, false);
	test02.Gauss_Seidel_method(1e-4, false);
	test02.Successive_Over_Relaxation_method(1.25, 1e-5, false);
	//
	/**
	 * 第三章 求特征值和特征向量
	 */
	auto test03 = MatrixEigenvaluesAndEigenvectors("../3.txt");
	// QR 分解
	// vector<vector<f8>> A, Q, R;
	// test03.QR_decomposition_Householder(A, Q, R, true);

	// 幂法与反幂法
	test03.power_method(1e-5, true, infinite_norm, 2);
	test03.inverse_power_method(norm::two_norm, 1e-5, false);
	test03.inverse_power_method(20.9681);

	//Jacobi 法
	test03.Jacobi_method(1e-5, false);

	// QR 方法
	test03.QR_method_normal(1e-5, 1);
	test03.QR_method_Transformed_Rotational_Matrix(1e-5, 1);
	test03.QR_method_with_one_step_shifted(1e-5, 1);
	test03.QR_method_with_two_steps_shifted(1e-5, false);

	/**
	 * 第四章 非线性方程与非线性方程组的迭代解法
	 */
	auto* test04 = new NonlinearEquations();
	//cout << setiosflags(ios::fixed) << setprecision(9);

	// 以下均为传递函数指针
	// 非线性方程的迭代解法
	test04->dichotomy_method(fun_dichotomy, 0, 1.6, 1e-3, false);
	test04->simple_iteration_method(fun_simple, 1, 1e-6, true);
	test04->Steffensen_iteration_method(fun_Steffensen, 0, 1e-6, false);
	test04->Newton_iteration_method(fun_Newton, fun_Newton_1, 1e-6, 1.4, false, false);
	test04->secant_iteration_method(fun_secant, 0, 1.6, 1e-6, false);
	test04->secant_with_single_point_iteration_method(fun_secant, 0, 1.6, 1e-6, false);

	// 方程组的迭代解法
	test04->simple_iteration_equations_method(fun_equations_simple1, fun_equations_simple2, 1, 1, 5e-4, false);
	test04->Newton_iteration_equations_method(fun_equations_Newton1, fun_equations_Newton2, 1, 1, 1e-4, false);
	//
	/**
	 * 第五章 插值与逼近
	 */
	auto test05 = new InterpolationAndApproximation("../5.txt");
	test05->Lagrange_polynomial(3);
	test05->Newton_polynomial(3);
	test05->three_moment_method_bound1(0, 0, true);
	test05->three_moment_method_bound2(-1, 2);
	test05->three_moment_method_bound3(false);


	cout << OrthogonalPolynomial::Legendre_orthogonal_polynomial(3).to_string() << endl;
	cout << OrthogonalPolynomial::Chebyshev_orthogonal_polynomial(3).to_string() << endl;
	cout << OrthogonalPolynomial::Laguerre_orthogonal_polynomial(3).to_string() << endl;
	cout << OrthogonalPolynomial::Hermite_orthogonal_polynomial(3).to_string() << endl;

	/**
	 * 第六章 数值积分
	 */
	auto test06 = new NumericalIntegration();
	test06->Romberg_numerical_integration_method(fun_Romberg, 0, 1, 1e-4);
	test06->Gauss_Legendre_integration_method(fun_Gauss_Legendre, 4, -0.8611363116, 0.8611363116, -0.3399810436,
	                                          0.3399810436);
	test06->Gauss_Legendre_integration_method(fun_Gauss_Legendre, 3, 0, -0.7745966692, 0.7745966692);

	test06->Gauss_Legendre_integration_method(fun_Gauss_Legendre, 5, 0, -0.5384693101, 0.5384693101,
	                                          -0.9061798459, 0.9061798459);


	test06->Gauss_Laguerre_integration_method(fun_Gauss_Lagueree, 3, 0.4157745568, 2.2942803603, 6.2899450829);
	test06->Gauss_Hermite_integration_method(fun_Gauss_Hermite, 5, 0.0, -2.0201828705, 2.0201828705, -0.9585724646,
	                                         0.9585724646);
	test06->Gauss_Chebyshev_integration_method(fun_Gauss_Chebyshev, 3, -0.866025403, 0.0, 0.866025403);
}
