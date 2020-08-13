#pragma once
#include <string>
#include <vector>

double interpolate1d(std::vector<std::vector<double>>& Data, double x,
                     bool extrapolate);
void get_parameter(std::string filename, std::vector<std::vector<double>>& data,
                   const int dim);
int inOrNot(int poly_sides, double* poly_X, double* poly_Y, double x, double y);
std::vector<double> linspace(double a, double b, int N);

