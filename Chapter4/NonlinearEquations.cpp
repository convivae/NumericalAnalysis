//
// Created by convi on 2020/4/30.
//

#include "NonlinearEquations.h"

void convivae::NonlinearEquations::draw_dividing_line()
{
	std::cout << "------------------------------------------------" << std::endl;
}

void convivae::NonlinearEquations::Dichotomy(const pFun fun, double left, double right, double epsilon,
                                             bool show_details, int max_steps)
{
	auto k = 0;
	auto mid = 0.0;
	while (k++ < max_steps) {
		mid = (left + right) / 2;

		if (show_details)
			std::cout << "第" << k << "次迭代:\nmid = " << mid;

		if (right - left < epsilon) {
			if (show_details)
				std::cout << ", right - left < " << epsilon << std::endl;
			break;
		}


		auto a = fun(left), b = fun(right), c = fun(mid);

		if (fabs(c) < 1e-10) {
			if (show_details)
				std::cout << "fabs(f(mid)) = " << fabs(c) << " < 1e-10" << std::endl;
			break;
		}

		if (show_details)
			std::cout << ", a * c = " << a * c << std::endl;

		if (a * c < 0)
			right = mid;
		else
			left = mid;

		if (show_details)
			std::cout << "left = " << left << ", right = " << right << std::endl;
	}

	std::cout << "迭代次数: " << k << std::endl;
	std::cout << "结果: " << mid << std::endl;
	draw_dividing_line();
}
