//
// 数值分析
// 第四章 非线性方程与非线性方程组的迭代解法
// Created by convi on 2020/4/30.
//

#ifndef NUMERICALANALYSIS_NONLINEAREQUATIONS_H
#define NUMERICALANALYSIS_NONLINEAREQUATIONS_H

#include "../Type/Type.h"
#include <iostream>
#include <vector>

namespace convivae {
	class NonlinearEquations {
	private:
		//定义一种 pFun 的的函数指针，这种函数以一个 double 为参数并返回 double 类型
		typedef f8 (*pFun)(f8 x);

		// 函数指针，分别指向迭代函数、迭代函数的一阶导数
		//pFun fun_, fun_1_;

		static void draw_dividing_line();

	public:
		// NonlinearEquations(const pFun fun, const pFun fun_1)
		//  : fun_(fun),fun_1_(fun_1){}

		/**
		 * 二分法（对分法）
		 * 只能求单根和奇数重根
		 * @param fun 迭代函数的函数指针
		 * @param left 左区间
		 * @param right 右区间
		 * @param epsilon 迭代精度
		 * @param max_steps 最大迭代次数
		 */
		void Dichotomy(const pFun fun, double left, double right, double epsilon = 1e-5, bool show_details = false,
		               int max_steps = 2000);


	};
}

#endif //NUMERICALANALYSIS_NONLINEAREQUATIONS_H
