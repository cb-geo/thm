#pragma once
#include <iostream>
#include "csv.h"
#include <string> 

// std::vector<double> coord_perm(11,0);
// std::vector<double> perm_value(11,0);

// void get_parameter(std::string filename, std::vector<double>& coord,std::vector<double>& data);

// get_parameter("parameters.txt", coord_perm, perm_value);//1e-5/1000/9.8;


namespace EquationData {

double perm = 1e-5/1000/9.8;
const double c_T = 1.2e6;
const double lam = 1.2;
const double c_w = 1e6;
const double B_w = 1e5;

std::vector<std::vector<double>> perm_list = {{0,0}};


// Temperature seetings
const int num_T_bnd_id = 4;
const int T_bnd_id[num_T_bnd_id] ={2,3,4,5};
const double Tb_well = 273.15+5;
const double Tb_top = 273.15+5;
const double T_grad = 0.05;

// Pressure seetings
const int num_P_bnd_id = 1;
const int P_bnd_id[num_P_bnd_id] ={5};
const double Pb_top = 100000;
const double P_grad = 1000*9.8;
const double Qb_well = -0.0005;


}  // namespace EquationData

