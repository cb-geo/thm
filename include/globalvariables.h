#pragma once

#include <iostream>
#include <string>

#include "csv.h"
#include <vector>

namespace EquationData {

double g_perm = 1e-6 / 1000 / 9.8;  // permeability
const double g_c_T = 1.2e6;         // heat capacity of the mixture
const double g_lam = 1.2;           // heat conductivity
const double g_c_w = 1e6;           // heat capacity of water
const double g_B_w = 1e5;           // bulk modulus of pores

// Pressure seetings

const int g_num_P_bnd_id = 8;  // numbers of  pressure boudnary condition id
const int g_P_bnd_id[g_num_P_bnd_id] = {
    1, 2, 4, 5, 7, 8, 9051, 9052};  // pressure boundary condition id

const double g_Pb_top = 0;           // pressure at the top of model
const double g_P_grad = 1000 * 9.8;  // pressure gradient in vertial direction
const double g_P_grad_x = 1000 * 9.8 * 10 / 2000;  // 2000m drops ten meters

// Velocity settings

const int g_num_QP_bnd_id = 4;  // numbers of velocity boudnary condition id
const int g_QP_bnd_id[g_num_QP_bnd_id] = {1, 4, 9051,
                                          9052};  // velocity  boundary

const double g_Qb_lateral = 0;  // wellbore temperature

// Temperature seetings
const int g_num_T_bnd_id = 11;
const int g_T_bnd_id[g_num_T_bnd_id] = {13, 9, 3, 1, 2, 4, 5, 7, 8, 9051, 9052};

const double g_Tb_well = 273.15 + 25;  // wellbore temperature
const double g_Tb_top = 273.15 + 15;   // termperature at the top of model
const double g_Tb_bottom = 273.15 + 18;
const double g_Tb_seabed_top = 273.15 + 15;  // termperature at the top of model
const double g_Tb_seabed_bottom = g_Tb_bottom;
const double g_T_grad = (g_Tb_bottom - g_Tb_top) /
                        60.96;  // temperature gradient in verital direction
const double g_T_seabed_grad = (g_Tb_seabed_bottom - g_Tb_seabed_top) / 60.96;

// Heat flow rate settins

const int g_num_QT_bnd_id = 2;  // numbers of velocity boudnary condition id
const int g_QT_bnd_id[g_num_QT_bnd_id] = {6, 10};  // velocity  boundary

const double g_QT_well = 0;  // wellbore temperature
const double g_QT_top = -g_lam * g_T_grad;
const double g_QT_bottom = g_lam * g_T_grad;

// solver settings

const bool is_linspace = false;
const double g_period = 86400 * 180;  // simulation time
const int g_n_time_step = 15;         // simulation time
std::vector<double> g_time_sequence = {0,  0.1, 0.2, 0.5, 1,   2,   5,  10,
                                       15, 30,  60,  90,  120, 150, 180};
const char time_unit = 'd';
const unsigned int n_g_P_max_iteration = 1000;
const unsigned int n_g_T_max_iteration = 4000;
const double g_P_tol_residual = 1e-8;
const double g_T_tol_residual = 1e-10;

// dimention of the input data file (parameters_for_interpolation.txt in
// inputfiles is used in the example)
const int dimension = 3;
// dimension in x, y and z directions
std::string file_name_interpolation =
    "inputfiles/parameters_for_interpolation.txt";
// number of sample point in x directin, y direction and z direction
const int dimension_x = 20;
const int dimension_y = 10;
const int dimension_z = 5;

}  // namespace EquationData
