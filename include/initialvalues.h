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
  return g_Pb_top + g_P_grad * (0. - p[2]);
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
  if (p[0] >= 0 ||
      (1000.658 - 710.184) / (-289.56) <= (p[1] - 710.184) / p[0] ||
      p[1] >= 1000.658 ||
      (1000.658 - 710.184) / (-1371.6 + 1626.108) <=
          (p[1] - 710.184) / (p[0] + 1626.108) ||
      p[0] <= -1626.108 ||
      279.806 / (-1626.108 + 1371.6) <= p[1] / (p[0] + 1371.6) || p[1] <= 0) {
    return g_Tb_seabed + g_T_grad * (0. - p[2]);
  } else {
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
