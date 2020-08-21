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
  return g_Pb_top + g_P_grad * (0. - p[2]) - g_P_grad_x * (0 - p[0]);
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
  const int poly_sides = 7;
  double poly_x[poly_sides] = {0,         0,         -289.56, -1371.6,
                               -1626.108, -1626.108, -1371.6};
  double poly_y[poly_sides] = {0,       710.184, 1000.658, 1000.658,
                               710.184, 279.806, 0};
  int res = inOrNot(poly_sides, poly_x, poly_y, p[0], p[1]);

  if (res == 0) {
    return g_Tb_seabed_top + g_T_seabed_grad * (0. - p[2]);
  } else if (res == 1) {
    return g_Tb_top + g_T_grad * (0. - p[2]);
  }
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
