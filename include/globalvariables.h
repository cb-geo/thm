#pragma once

#include <iostream>
#include <string>

#include "csv.h"

// std::vector<double> coord_perm(11,0);
// std::vector<double> perm_value(11,0);

// void get_parameter(std::string filename, std::vector<double>&
// coord,std::vector<double>& data);

// get_parameter("parameters.txt", coord_perm, perm_value);//1e-5/1000/9.8;

namespace EquationData {

double g_perm = 1e-5 / 1000 / 9.8;
const double g_c_T = 1.2e6;
const double g_lam = 1.2;
const double g_c_w = 1e6;
const double g_B_w = 1e5;

std::vector<std::vector<double>> g_perm_list = {{0, 0}};

// Temperature seetings
const int g_num_T_bnd_id = 4;
const int g_T_bnd_id[g_num_T_bnd_id] = {2, 3, 4, 5};
const double g_Tb_well = 273.15 + 5;
const double g_Tb_top = 273.15 + 5;
const double g_T_grad = 0.05;

// Pressure seetings
const int g_num_P_bnd_id = 1;
const int g_P_bnd_id[g_num_P_bnd_id] = {5};
const double g_Pb_top = 100000;
const double g_P_grad = 1000 * 9.8;
const double g_Qb_well = -0.0005;

}  // namespace EquationData
