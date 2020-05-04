
#pragma once
// #include <deal.II/numerics/vector_tools.h>
#include "globalvariables.h"
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iostream>

using namespace dealii;

namespace EquationData {
template <int dim>
class PressureBoundaryValues : public Function<dim> {
 public:
  PressureBoundaryValues() : Function<dim>(1), period(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary
  virtual void vector_value(const Point<dim>& p,  //之前的没有
                            Vector<double>& value) const;

 private:
  const double period;  // value
};

template <int dim>
double PressureBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());
  const double time = this->get_time();  // get time
  return p0;
}

template <int dim>
class TemperatureBoundaryValues : public Function<dim> {
 public:
  TemperatureBoundaryValues() : Function<dim>(1), period(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary
  virtual void vector_value(const Point<dim>& p,  //之前的没有
                            Vector<double>& value) const;

 private:
  const double period;  // value
};

template <int dim>
double TemperatureBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());

  const double time = this->get_time();
  return T0 + 10. * sin(time * 3.1415926);  // boundary value is set to zero in
                                            // this case
}

template <int dim>
void TemperatureBoundaryValues<dim>::vector_value(
    const Point<dim>& p, Vector<double>& values) const {
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = TemperatureBoundaryValues<dim>::value(p, c);
}
}  // namespace EquationData