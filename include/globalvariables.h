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

double g_perm = 1e-5 / 1000 / 9.8;  // permeability
const double g_c_T = 1.2e6;         // heat capacity of the mixture
const double g_lam = 1.2;           // heat conductivity
const double g_c_w = 1e6;           // heat capacity of water
const double g_B_w = 1e5;           // bulk modulus of pores

std::vector<std::vector<double>> g_perm_list = {
    {0, 0}};  // permeablity list used for interpolation

// Temperature seetings
const int g_num_T_bnd_id = 4;  // numbers of  temperature boudnary condition id
const int g_T_bnd_id[g_num_T_bnd_id] = {
    2, 3, 4, 5};                      //  temperature boudnary condition id
const double g_Tb_well = 273.15 + 5;  // wellbore temperature
const double g_Tb_top = 273.15 + 5;   // termperature at the top of model
const double g_T_grad = 0.05;  // temperature gradient in verital direction

// Pressure seetings
const int g_num_P_bnd_id = 1;  // numbers of  pressure boudnary condition id
const int g_P_bnd_id[g_num_P_bnd_id] = {5};  // pressure boundary condition id
const double g_Pb_top = 100000;              // pressure at the top of model
const double g_P_grad = 1000 * 9.8;  // pressure gradient in vertial direction
const double g_Qb_well = -0.0005;    // wellbore temperature

// solver settings
const double g_theta = 0.3;
const double g_period = 1 * 3600 * 24;     // simulation time
const double g_time_step = g_period / 20;  // time step

const unsigned int g_P_max_iteration_number = 1000;
const unsigned int g_T_max_iteration_number = 4000;
const double g_P_tol_residual = 1e-8;
const double g_T_tol_residual = 1e-8;

}  // namespace EquationData
