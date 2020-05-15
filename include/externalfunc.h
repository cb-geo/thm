#pragma once
#include <vector>
#include <string>

double interpolate1d( std::vector<std::vector<double>>& Data, 
                      double x, bool extrapolate );
void get_parameter(std::string filename, 
                  std::vector<std::vector<double>>& data, const int dim);