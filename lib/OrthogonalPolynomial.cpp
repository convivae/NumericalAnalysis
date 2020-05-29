//
// Created by convi on 2020/5/29.
//

#include "OrthogonalPolynomial.h"

Polymerization OrthogonalPolynomial::Legendre_orthogonal_polynomial(int n) {
    Polymerization L0("1"), L1("x"), L("0");
    if (n == 0)
        return L0;
    if (n == 1)
        return L1;
    for (auto i = 1; i < n; ++i) {
        Polymerization a = Polymerization("x") * L1 * ((2.0 * i + 1.0) / (i + 1.0));
        Polymerization b = L0 * (-1.0 * i / (i + 1.0));
        L = a + b;
        L0 = L1;
        L1 = L;
    }

    return L;}

Polymerization OrthogonalPolynomial::Chebyshev_orthogonal_polynomial(int n) {
    Polymerization L0("1"), L1("x"), L("0");
    if (n == 0)
        return L0;
    if (n == 1)
        return L1;
    for (auto i = 1; i < n; ++i) {
        L = Polymerization("2x") * L1 - L0;
        L0 = L1;
        L1 = L;
    }

    return L;
}

Polymerization OrthogonalPolynomial::Laguerre_orthogonal_polynomial(int n) {
    Polymerization L0("1"), L1("1-x"), L("0");
    if (n == 0)
        return L0;
    if (n == 1)
        return L1;
    for (auto i = 1; i < n; ++i) {
        string s = to_string(2.0 * i + 1) + "-x";
        L = Polymerization(s) * L1 - L0 * (1.0 * i * i);
        L0 = L1;
        L1 = L;
    }

    return L;
}

Polymerization OrthogonalPolynomial::Hermite_orthogonal_polynomial(int n) {
    Polymerization L0("1"), L1("2x"), L("0");
    if (n == 0)
        return L0;
    if (n == 1)
        return L1;
    for (auto i = 1; i < n; ++i) {
        L = Polymerization("2x") * L1 - L0 * (2.0 * i);
        L0 = L1;
        L1 = L;
    }

    return L;
}
