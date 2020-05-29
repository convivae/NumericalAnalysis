//
// Created by convi on 2020/5/29.
// 各种正交多项式
//

#ifndef NUMERICALANALYSIS_ORTHOGONALPOLYNOMIAL_H
#define NUMERICALANALYSIS_ORTHOGONALPOLYNOMIAL_H


#include "Polymerization.h"
#include <string>

using namespace std;

class OrthogonalPolynomial {
public:
    /**
     * 正交多项式
     * Legendre 多项式（勒让德多项式）
     * @param n L_n(x)
     * @return
     */
    static Polymerization Legendre_orthogonal_polynomial(int n);

    /**
     * 正交多项式
     * Chebyshev 多项式（切比雪夫多项式）
     * @param n L_n(x)
     * @return
     */
    static Polymerization Chebyshev_orthogonal_polynomial(int n);

    /**
     * 正交多项式
     * Laguerre 多项式(拉盖尔多项式)
     * @param n L_n(x)
     * @return
     */
    static Polymerization Laguerre_orthogonal_polynomial(int n);

    /**
     * 正交多项式
     * Hermite 多项式(埃尔米特多项式)
     * @param n L_n(x)
     * @return
     */
    static Polymerization Hermite_orthogonal_polynomial(int n);
};


#endif //NUMERICALANALYSIS_ORTHOGONALPOLYNOMIAL_H
