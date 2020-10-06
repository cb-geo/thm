#pragma once

#include <iostream>
#include <string>

#include "csv.h"
#include <vector>

namespace EquationData {

double g_perm = 1e-6 / 1000 / 9.8;  // permeability
double g_c_T = 1.2e6;               // heat capacity of the mixture
double g_lam = 1.2;                 // heat conductivity
const double g_c_w = 1e6;           // heat capacity of water
const double g_B_w = 1e4;           // bulk modulus of pores

// Pressure seetings
const int g_num_P_bnd_id = 10;  // numbers of  pressure boudnary condition id
const int g_P_bnd_id[g_num_P_bnd_id] = {
    6, 3, 1, 2, 4, 5, 7, 8, 11311, 11312};  // pressure boundary condition id

const double g_Pb_top = 0;           // pressure at the top of model
const double g_P_grad = 1000 * 9.8;  // pressure gradient in vertial direction
const double g_P_grad_x = 1000 * 9.8 * 1 / 2000;  // 2000m drops ten meters

// Velocity settings
// const int g_num_QP_bnd_id = 4;  // numbers of velocity boudnary condition id
// const int g_QP_bnd_id[g_num_QP_bnd_id] = {2, 5, 7, 8};
const int g_num_QP_bnd_id = 0;  // numbers of velocity boudnary condition id
const int g_QP_bnd_id[1] = {-1};

// velocity  boundary
const double g_Qb_lateral = 0;

// Temperature seetings for injection
const int g_num_T_bnd_id = 3;
const int g_T_bnd_id[g_num_T_bnd_id] = {23, 6, 3};
// const int g_num_T_bnd_id = 11;
// const int g_T_bnd_id[g_num_T_bnd_id] = {23, 6, 3, 1, 2, 4, 5, 7, 8, 11311,
// 11312};

const double g_Tb_top = 273.15 + 14.3;  // termperature at the top of model
std::vector<std::vector<double>> g_Tb_top_seasonal = {
    {86400 * 30 * 0, 273.15 + 14.5},  {86400 * 30 * 1, 273.15 + 15.5},
    {86400 * 30 * 2, 273.15 + 15.5},  {86400 * 30 * 3, 273.15 + 16.5},
    {86400 * 30 * 4, 273.15 + 17},    {86400 * 30 * 5, 273.15 + 17},
    {86400 * 30 * 6, 273.15 + 14.5},  {86400 * 30 * 7, 273.15 + 11},
    {86400 * 30 * 8, 273.15 + 11},    {86400 * 30 * 9, 273.15 + 12.5},
    {86400 * 30 * 10, 273.15 + 13.5}, {86400 * 30 * 11, 273.15 + 13.5}};
const double g_T_grad = 0.03;  // temperature gradient in verital direction
const double g_T_seabed_grad = g_T_grad;
const double g_Tb_bottom = g_Tb_top * g_T_grad * 100;
const double g_Tb_seabed_top = g_Tb_top;  // termperature at the top of model
const double g_Tb_seabed_bottom = g_Tb_bottom;

const double g_Tb_well =
    (2 * g_Tb_top + g_T_grad * 60) / 2;  // average pile temperature

// Heat flow rate settins
// const int g_num_QT_bnd_id = 1;  // numbers of velocity boudnary condition id
// const int g_QT_bnd_id[g_num_QT_bnd_id] = {3};  // velocity  boundary
// const double g_QT_well = 0;  // wellbore temperature
// const double g_QT_top = -2.2 * g_T_grad;
// const double g_QT_bottom = 2.97 * g_T_grad;

// // Heat flow rate settins for injection
const int g_num_QT_bnd_id = 0;
const int g_QT_bnd_id[1] = {-1};  // velocity  boundary
// const double g_QT_well = 0;       // wellbore temperature
// const double g_QT_top = 0;
const double g_QT_bottom = 2.97 * g_T_grad;

// solver settings

const bool is_linspace = false;
const double g_total_time = 86400 * 360;  // simulation time
const double g_period = 86400 * 360;      // for periodic load
const int g_n_time_step = 25;             // simulation time
// std::vector<double> g_time_sequence = {
//     0,    90,   180,  270,  360,  450,  540,  630,  720,  810,  900,
//     990,  1080, 1170, 1260, 1350, 1440, 1530, 1620, 1710, 1800, 1890,
//     1980, 2070, 2160, 2250, 2340, 2430, 2520, 2610, 2700, 2790, 2880,
//     2970, 3060, 3150, 3240, 3330, 3420, 3510, 3600};
// std::vector<double> g_time_sequence = {0,   30,  60,  90,  120, 150, 180,
//                                        210, 240, 270, 300, 330, 360};
std::vector<double> g_time_sequence = {
    0,     0.5, 1,   2,   5,   7,   14,  30,  60,  90,  120, 150, 180,
    180.5, 181, 182, 185, 187, 194, 210, 240, 270, 300, 330, 360};

const char time_unit = 'd';
const unsigned int n_g_P_max_iteration = 1000;
const unsigned int n_g_T_max_iteration = 4000;
const double g_P_tol_residual = 1e-8;
const double g_T_tol_residual = 1e-10;

// dimention of the input data file (parameters_for_interpolation.txt in
// inputfiles is used in the example)
const int dimension = 3;
// dimension in x, y and z directions
// std::string mesh_file_name = "inputfiles/new_build.msh";
std::string mesh_file_name = "inputfiles/shortbeam_hexa_fracture_update.msh";

std::string perm_file_name_interpolation =
    "inputfiles/parameters_for_perm_interpolation.txt";
std::string therm_file_name_interpolation =
    "inputfiles/parameters_for_therm_interpolation.txt";
std::string capa_file_name_interpolation =
    "inputfiles/parameters_for_capa_interpolation.txt";
// number of sample point in x directin, y direction and z direction
const int dimension_x = 20;
const int dimension_y = 10;
const int dimension_z = 10;

}  // namespace EquationData
