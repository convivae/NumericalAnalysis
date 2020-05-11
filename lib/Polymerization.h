//
// 多项式
// Created by convi on 2020/5/12.
//

#ifndef NUMERICALANALYSIS_POLYMERIZATION_H
#define NUMERICALANALYSIS_POLYMERIZATION_H


#include <map>
#include <utility>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

using std::map;
using std::string;
using std::stringstream;

namespace convivae {
    class Polymerization {
    private:
        typedef map<int, double> poly_type;

        poly_type _polynome;

        void read_from_string(const string &s);

        explicit Polymerization(map<int, double> p) : _polynome(std::move(p)) {}

    public:
        explicit Polymerization(const string &s) {
            read_from_string(s);
        }

        static string to_string(const poly_type &p);

        Polymerization operator+(const Polymerization &p);

        Polymerization operator-(const Polymerization &p);

        Polymerization operator*(const Polymerization &p);

        // 求一阶导数
        Polymerization derivative();

        // 二阶导数
        Polymerization second_derivative();

        // 三阶导数
        Polymerization third_derivative();
    };

}


#endif //NUMERICALANALYSIS_POLYMERIZATION_H
