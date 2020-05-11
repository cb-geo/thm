#pragma once
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <utility>
#include <tuple>


namespace EquationData {
std::vector<double> coord_perm;
std::vector<double> perm;
auto [coord_perm, perm] = get_parameter("parameter.txt");//1e-5/1000/9.8;

const double c_T = 1.2e6;
// const double c_T = 1;
const double lam = 1.2;
// const double lam = 1;
const double c_w = 1e6;
// const double c_w = 6;
const double B_w = 1e5;

// Temperature seetings
const int num_T_bnd_id = 4;
int T_bnd_id[num_T_bnd_id] ={2,3,4,5};
const double Tb_well = 273.15+5;
const double Tb_top = 273.15+5;
const double T_grad = 0.05;

// Pressure seetings
const int num_P_bnd_id = 1;
int P_bnd_id[num_P_bnd_id] ={5};
const double Pb_top = 100000;
const double P_grad = 1000*9.8;
const double Qb_well = -0.0005;




std::tuple<std::vector<double>,std::vector<double>> get_parameter(std::string filename)
{
    std::vector<double> coord;
    std::vector<double> data;
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::file_source(filename));
    for (unsigned int line = 0; line < 83600; ++line)
    {
      try
        {
          double z, parameter;
          in >> z >> parameter;
          coord.push_back(z);
          data.push_back(parameter);
        }
      catch (...)
        {
          AssertThrow(false,
                      ExcMessage("Could not read all data points "));
        }
    }
  return make_tuple(coord, data);
}



}  // namespace EquationData
