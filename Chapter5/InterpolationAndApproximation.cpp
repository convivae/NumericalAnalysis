//
// Created by convi on 2020/5/12.
//

#include "InterpolationAndApproximation.h"

void convivae::InterpolationAndApproximation::draw_dividing_line() {
    cout << "---------------------------------------------------------" << endl;
}

void convivae::InterpolationAndApproximation::read_data() {
    if (_initialized)
        return;

    ifstream in(_rdr);
    string s;

    getline(in, s);
    istringstream sx(s);
    for (f8 i; sx >> i;) {
        this->X.push_back(i);
    }

    getline(in, s);
    istringstream sy(s);
    for (f8 i; sy >> i;) {
        this->Y.push_back(i);
    }

    in.close();

    if (this->X.size() != this->Y.size()) {
        cerr << "data error" << endl;
        exit(-1);
    }

    this->_initialized = true;
    this->_n = this->X.size();
}

Polymerization convivae::InterpolationAndApproximation::Lagrange_polynomial(int n) {
    cout << "Lagrange 插值法" << endl;
    read_data();
    if (n > _n) {
        cerr << "error" << endl;
        exit(-1);
    }
    Polymerization res("0");
    string s;
    stringstream ss;
    for (auto k = 0; k <= n; ++k) {
        Polymerization tmp("1");
        for (auto j = 0; j <= n; ++j) {
            if (j == k)
                continue;
            double a = 1 / (X[k] - X[j]);
            double b = X[j] == 0 ? 0.0 : -X[j] / (X[k] - X[j]);
            ss << a << "x";
            b < 0.0 ? ss << b : ss << "+" << b;
            ss >> s;
            ss.clear();
            ss.str("");
            tmp *= Polymerization(s);
        }
        res += tmp * Y[k];
    }

    cout << res.to_string() << endl;
    draw_dividing_line();
    return res;
}

Polymerization convivae::InterpolationAndApproximation::Newton_polynomial(int n) {
    cout << "Newton 插值法" << endl;
    read_data();
    if (n > _n) {
        cerr << "error" << endl;
        exit(-1);
    }

    vector<vector<double>> table(n + 1, vector<double>(n + 2, 0));
    vector<vector<bool>> flag(n + 1, vector<bool>(n + 2, false));
    for (auto i = 0; i <= n; ++i) {
        table[i][0] = X[i];
        table[i][1] = Y[i];
        flag[i][0] = true;
        flag[i][1] = true;
    }

    for (auto i = 1; i <= n; ++i) {
        for (auto j = 2; j <= i + 1; ++j) {
            table[i][j] = differential_quotient(i - j + 1, i - 1, i, table, flag);
        }
    }

    cout << "table:" << endl;
    cout << "xi\tf[xi]\t";
    for (auto i = 0; i < n; ++i) {
        cout << "No." << i + 1 << "\t";
    }
    cout << endl;

    for (auto i = 0; i <= n; ++i) {
        for (auto j = 0; j < i + 2; ++j) {
            cout << table[i][j] << "\t";
        }
        cout << endl;
    }

    cout << "The Newton polynomial is:" << endl;
    Polymerization res(to_string(table[0][1]));
    Polymerization p("1");

    for (auto i = 1; i <= n; ++i) {
        auto s = to_string(-table[i - 1][0]);
        if (-table[i - 1][0] >= 0)
            s = string("x+").append(s);
        else
            s = string("x").append(s);
        p *= Polymerization(s);
        res += p * table[i][i + 1];

    }

    cout << res.to_string() << endl;
    draw_dividing_line();
    return res;
}

double
convivae::InterpolationAndApproximation::differential_quotient(int start_index, int before_end_index, int end_index,
                                                               vector<vector<double>> &table,
                                                               vector<vector<bool>> &flag) {
    if (start_index == end_index) { // f(xk) 零阶
        table[start_index][1] = Y[start_index];
        flag[start_index][1] = true;
        return Y[start_index];
    }
    if (before_end_index == start_index) { //f(x_start, x_end) 一阶
        double res = (Y[end_index] - Y[start_index]) / (X[end_index] - X[start_index]);
        if (end_index - before_end_index == 1) {
            table[end_index][2] = res;
            flag[end_index][2] = true;
        }
        return res;
    }

    // f(x0, x1, ... , x_(k-1), xk) k阶
    double res;
    int row = end_index, column = before_end_index - start_index + 2;
    if (end_index - before_end_index == 1 && flag[row][column]) {
        res = table[row][column];
    } else {
        double a = differential_quotient(start_index, before_end_index - 1, end_index, table, flag);
        double b = differential_quotient(start_index, before_end_index - 1, before_end_index, table, flag);
        res = (a - b) / (X[end_index] - X[before_end_index]);
        if (end_index - before_end_index == 1) {
            table[row][column] = res;
            flag[row][column] = true;
        }
    }

    return res;
}


