#pragma once

#include <fstream>
#include <iostream>

// #include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

#include "globalvariables.h"

using namespace dealii;

namespace EquationData {

template <int dim>
class PressureInitialValues : public Function<dim> {
 public:
  PressureInitialValues() : Function<dim>(1) {}

  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;
  // virtual void vector_value(const Point<dim> &p,
  //                           Vector<double> &value) const;
};

template <int dim>
double PressureInitialValues<dim>::value(const Point<dim>& p,
                                         const unsigned int) const {
  const int poly_sides = 12;
  double poly_x[poly_sides] = {-2100.609, -3253.2102, -2031.5439, -211.4263, 763.4342, 
                575.0074, 279.3723, -423.9736, -564.616, -869.7101, -1054.3573, -1726.8217};

  double poly_y[poly_sides] = {4435.6898, 5455.6523, 6610.3567, 5117.4889, 3169.5471, 
                2851.097, 2913.4759, 2501.9842, 2612.1521, 3204.4684, 3462.1967, 4216.4399};

  int res = inOrNot(poly_sides, poly_x, poly_y, p[0], p[1]);

  double water_table = 0;
  if (res == 0) {
    water_table = ground_level_lower - 2;
  } else if (res == 1) {
    water_table = ground_level_upper - 2;
  }
  //return (water_table - p[2]) > 0 ? g_P_grad*p[2] : 0;
  return 0;
}

// template <int dim>
// void PressureInitialValues<dim>::vector_value(const Point<dim> &p,
//                                               Vector<double> &values) const
// {
//     for (unsigned int c = 0; c < this->n_components; ++c)
//         values(c) = PressureInitialValues<dim>::value(p, c);
// }

template <int dim>
class TemperatureInitialValues : public Function<dim> {
 public:
  TemperatureInitialValues() : Function<dim>(1) {}

  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;
  // virtual void vector_value(const Point<dim> &p,
  //                           Vector<double> &value) const;
};

template <int dim>
double TemperatureInitialValues<dim>::value(const Point<dim>& p,
                                            const unsigned int) const {
  const int poly_sides = 12;
  double poly_x[poly_sides] = {-2100.609, -3253.2102, -2031.5439, -211.4263, 763.4342, 
                575.0074, 279.3723, -423.9736, -564.616, -869.7101, -1054.3573, -1726.8217};

  double poly_y[poly_sides] = {4435.6898, 5455.6523, 6610.3567, 5117.4889, 3169.5471, 
                2851.097, 2913.4759, 2501.9842, 2612.1521, 3204.4684, 3462.1967, 4216.4399};

  int res = inOrNot(poly_sides, poly_x, poly_y, p[0], p[1]);

  double depth = 0;
  if (res == 0) {
    depth = ground_level_lower - p[2];
  } else if (res == 1) {
    depth = ground_level_upper - p[2];
  }
  if(depth < 1) return g_Ti_1m;
  else if(depth < 2) return g_Ti_2m;
  else if(depth < 3) return g_Ti_3m;
  else if(depth < 4) return g_Ti_4m;
  else return g_Ti_d;
}

// template <int dim>
// void TemperatureInitialValues<dim>::vector_value(const Point<dim> &p,
//                                                  Vector<double> &values)
//                                                  const
// {
//     for (unsigned int c = 0; c < this->n_components; ++c)
//         values(c) = TemperatureInitialValues<dim>::value(p, c);
// }
}  // namespace EquationData
