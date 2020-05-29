//
// Created by convi on 2020/5/12.
//

#include "Polymerization.h"

void Polymerization::read_from_string(const string &s) {
    std::vector<string> poly;

    size_t begin = 0;
    for (auto i = 1; i < s.size(); ++i) {
        if ((s[i] == '+' || s[i] == '-') && s[i - 1] != '^' && s[i-1] != 'e') {
            poly.push_back(s.substr(begin, i - begin));
            begin = i;
        }
    }
    poly.push_back(s.substr(begin, s.size()));

    for (const auto &i : poly) {
        stringstream ss1, ss2;
        int power;
        double coefficient;

        size_t j = 0;
        for (; j < i.size(); ++j) {
            if (i[j] == 'x')
                break;
            ss1 << i[j];
        }
        if (j == 0) // x
            ss1 << '1';
        if (j == 1 && ss1.str().size() == 1) {
            // 单个数字或 +x -x
            if (ss1.str()[0] == '+' || ss1.str()[0] == '-') // +x 或 -x
                ss1 << '1';
        }
        ss1 >> coefficient;
        if (!ss1.eof() || ss1.fail()) {
            std::cerr << "input polynomial error!" << std::endl;
            exit(-1);
        }

        if (j == i.size())
            power = 0;
        else if (j + 1 == i.size()) {
            power = 1;
        } else if (j + 2 < i.size()) {
            for (j += 2; j < i.size(); ++j) {
                ss2 << i[j];
            }
            ss2 >> power;
            if (!ss2.eof() || ss2.fail()) {
                std::cerr << "input polynomial error!" << std::endl;
                exit(-1);
            }
        } else {
            std::cerr << "input polynomial error!" << std::endl;
            exit(-1);
        }

        if (this->_polynome.find(power) == this->_polynome.end()) {
            this->_polynome.insert(std::make_pair(power, coefficient));
        } else {
            this->_polynome[power] += coefficient;
        }
    }
    // std::cout << to_string(this->_polynome) << std::endl;
}


string Polymerization::to_string() {
    string s;
    stringstream ss;
    for (auto i = this->_polynome.rbegin(); i != this->_polynome.rend(); ++i) {
        if ((*i).second > 0)
            ss << '+';
        ss << (*i).second << "x^" << (*i).first;
    }

    ss >> s;
    return s;
}

Polymerization Polymerization::operator+(const Polymerization &p) {
    poly_type res;
    for (auto i : this->_polynome) {
        if (p._polynome.find(i.first) == p._polynome.end()) {
            res.insert(i);
        } else {
            res.insert(std::make_pair(i.first, i.second + p._polynome.at(i.first)));
        }
    }
    for (auto i : p._polynome) {
        if (this->_polynome.find(i.first) == this->_polynome.end())
            res.insert(i);
    }

    return Polymerization(res);
}

Polymerization Polymerization::operator+(double n) {
    if (this->_polynome.find(0) == this->_polynome.end()) {
        this->_polynome.insert(make_pair(0, n));
    } else {
        this->_polynome[0] += n;
    }

    return Polymerization(this->_polynome);
}

Polymerization Polymerization::operator-(const Polymerization &p) {
    poly_type res;
    for (auto i : this->_polynome) {
        auto j = p._polynome.find(i.first);
        if (j == p._polynome.end()) {
            res.insert(i);
        } else {
            res.insert(std::make_pair(i.first, i.second - p._polynome.at(i.first)));
        }
    }
    for (auto i : p._polynome) {
        if (this->_polynome.find(i.first) == this->_polynome.end())
            res.insert(std::make_pair(i.first, -i.second));
    }

    return Polymerization(res);
}

Polymerization Polymerization::operator-(const double n) {
    if (this->_polynome.find(0) == this->_polynome.end()) {
        this->_polynome.insert(make_pair(0, -n));
    } else {
        this->_polynome[0] -= n;
    }

    return Polymerization(this->_polynome);
}

Polymerization Polymerization::operator*(const Polymerization &p) {
    poly_type res;
    for (auto i : this->_polynome) {
        for (auto j : p._polynome) {
            auto coefficient = i.second * j.second;
            auto power = i.first + j.first;
            if (!res.empty() && res.find(power) != res.end())
                res[power] += coefficient;
            else
                res.insert(std::make_pair(power, coefficient));
        }
    }

    return Polymerization(res);
}

Polymerization Polymerization::operator*(const double n) {
    poly_type res;
    for (auto i : this->_polynome) {
        res.insert(make_pair(i.first, i.second * n));
    }
    return Polymerization(res);
}

Polymerization &Polymerization::operator=(const Polymerization &p)
= default;

Polymerization &Polymerization::operator+=(const Polymerization &p) {
    for (auto i : p._polynome) {
        auto j = this->_polynome.find(i.first);
        if (j == this->_polynome.end()) {
            this->_polynome.insert(i);
        } else {
            this->_polynome[i.first] += p._polynome.at(i.first);
        }
    }

    return *this;
}

Polymerization &Polymerization::operator-=(const Polymerization &p) {
    for (auto i : p._polynome) {
        auto j = this->_polynome.find(i.first);
        if (j == this->_polynome.end()) {
            this->_polynome.insert(make_pair(i.first, -i.second));
        } else {
            this->_polynome[i.first] -= p._polynome.at(i.first);
        }
    }

    return *this;
}

Polymerization &Polymerization::operator*=(const Polymerization &p) {
    poly_type res;
    for (auto i : this->_polynome) {
        for (auto j : p._polynome) {
            auto coefficient = i.second * j.second;
            auto power = i.first + j.first;
            if (!res.empty() && res.find(power) != res.end())
                res[power] += coefficient;
            else
                res.insert(std::make_pair(power, coefficient));
        }
    }

    this->_polynome = res;
    return *this;
}

bool Polymerization::operator==(const Polymerization &p) {
    if (this->_polynome.size() == p._polynome.size()) {
        for (auto i : this->_polynome) {
            if (p._polynome.find(i.first) == p._polynome.end() || this->_polynome[i.first] != p._polynome.at(i.first)) {
                return false;
            }
        }
        return true;
    }
    return false;
}


Polymerization Polymerization::derivative() {
    poly_type res;
    for (auto i : this->_polynome) {
        if (i.first != 0) {
            auto coefficient = i.first * i.second;
            auto power = i.first - 1;
            res.insert(std::make_pair(power, coefficient));
        }
    }
    if (res.empty())
        res.insert(std::make_pair(0, 0));
    return Polymerization(res);
}

Polymerization Polymerization::second_derivative() {
    return this->derivative().derivative();
}

Polymerization Polymerization::third_derivative() {
    return this->second_derivative().derivative();
}

double Polymerization::value_at_point(double x) {
    double res = 0;
    for (auto i : this->_polynome) {
        res += i.second * pow(x, i.first);
    }
    return res;
}

Polymerization operator+(double n, Polymerization p) {
    return p + n;
}

Polymerization operator-(const double n, const Polymerization &p) {
    return Polymerization(std::to_string(n)) - p;
}

Polymerization operator*(const double n, Polymerization p) {
    return p * n;
}


