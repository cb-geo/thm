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
class PressureSourceTerm : public Function<dim> {
 public:
  PressureSourceTerm() : Function<dim>(), period_(0.2) {}

  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim>& p, Vector<double>& values) const;

 private:
  const double period_;
};

template <int dim>
double PressureSourceTerm<dim>::value(const Point<dim>& p,
                                      const unsigned int component) const {
  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());
  const double time = this->get_time();  // get time

  if (component == 0)
    return 0.;
  else if (component == 1)
    return 0.;
  else if (component == 2)
    // return -EquationData::B_w*EquationData::perm*EquationData::P_grad;
    return 0.;
}

template <int dim>
void PressureSourceTerm<dim>::vector_value(const Point<dim>& p,
                                           Vector<double>& values) const {
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = PressureSourceTerm<dim>::value(p, c);
}

template <int dim>
class TemperatureSourceTerm : public Function<dim> {
 public:
  TemperatureSourceTerm() : Function<dim>(), period_(0.2) {}

  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;

  // virtual void vector_value(const Point<dim> &p,
  //                           Vector<double> &value) const;

 private:
  const double period_;
};

template <int dim>
double TemperatureSourceTerm<dim>::value(const Point<dim>& p,
                                         const unsigned int) const {
  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());
  const double time = this->get_time();  // get time
  return 0.;
}

// template <int dim>
// void TemperatureSourceTerm<dim>::vector_value(const Point<dim> &p,
//                                                  Vector<double> &values)
//                                                  const
// {
//     for (unsigned int c = 0; c < this->n_components; ++c)
//         values(c) = TemperatureSourceTerm<dim>::value(p, c);
// }

}  // namespace EquationData
