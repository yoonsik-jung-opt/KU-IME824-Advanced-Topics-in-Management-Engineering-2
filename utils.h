//
// Created by yoonsikjung on 2022/10/25.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <istream>
#include <numeric>
#include <algorithm>
#include <random>
#include <filesystem>
#include "LinAlg.h"

namespace fs = std::__fs::filesystem;
using namespace std;

template<typename T>
vector<double> readRow(istream &file, char delimiter);

template<typename T>
void readData(matx &ptr, string fname, char delimiter) {
    ifstream file(fname);
    vector<vector<double>> res;
    if (file.fail()) {
        cout << "File Not Exists" << endl;
    }
    while (file.good()) {
        vector<double> row = readRow<double>(file, delimiter);
        if (row.size() == 0)
            break;
        else
            res.push_back(row);
    }
    file.close();
    ptr = res;
}

template<typename T>
vector<double> readRow(istream &file, char delimiter) {
    vector<double> res;
    stringstream ss;
    string s;
    T value;
    bool inQuotes = false;
    while (file.good()) {
        char c = (char) file.get();
        if (c == '\r' || c == '\n' || c == '\xff') {
            break;
        } else if (c == delimiter) {
            continue;
        } else {
            while ((c != delimiter) & !(c == '\r' || c == '\n')) {
                s.push_back(c);
                c = (char) file.get();
            }

            if (std::is_same<T, int>() || std::is_same<T, bool>()) {
                res.push_back(stoi(s));
            } else if (std::is_same<T, float>() || std::is_same<T, double>()) {
                res.push_back(stod(s));
            }
            s.clear();
            file.unget();
        }
    }
    return res;
};