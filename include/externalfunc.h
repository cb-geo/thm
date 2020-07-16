#pragma once
#include <string>
#include <vector>

double interpolate1d(std::vector<std::vector<double>>& Data, double x,
                     bool extrapolate);
void get_parameter(std::string filename, std::vector<std::vector<double>>& data,
                   const int dim);
std::vector<double> linspace(double a, double b, int N);
