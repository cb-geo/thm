#pragma once

#include <fstream>
#include <iostream>

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>

#include "globalvariables.h"

using namespace dealii;

namespace EquationData {
template <int dim>
class PressureDirichletBoundaryValues : public Function<dim> {
 private:
  const double period_;  // value
  int boundary_id_{-1};

 public:
  PressureDirichletBoundaryValues()
      : Function<dim>(), period_(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary

  // virtual void vector_value(const Point<dim>& p,  //放在这里没啥用
  //                           Vector<double>& value) const;
  virtual void set_boundary_id(int bnd_id) { boundary_id_ = bnd_id; }
};

template <int dim>
double PressureDirichletBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());
  const double time = this->get_time();  // get time

  // if (boundary_id_ == 4) {
  //   return g_Pb_top;
  // } else if (boundary_id_ == 3) {
  //   return g_Pb_top + g_P_grad * (0. - p[2]);
  // } else if (boundary_id_ == 5) {
  //   return g_Pb_top + g_P_grad * (0. - p[2]);
  // }
  if (boundary_id_ == 3 || boundary_id_ == 8 || boundary_id_ == 12 ||
      boundary_id_ == 13 || boundary_id_ == 14) {
    return g_Pb_top + g_P_grad * (0. - p[2]);
  }
}

template <int dim>
class TemperatureDirichletBoundaryValues : public Function<dim> {
 private:
  const double period_;  // value
  int boundary_id_{-1};

 public:
  TemperatureDirichletBoundaryValues()
      : Function<dim>(), period_(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary
  // virtual void vector_value(const Point<dim>& p,  //放在这里没啥用
  //                           Vector<double>& value) const;
  virtual void set_boundary_id(int bnd_id) { boundary_id_ = bnd_id; }
};

template <int dim>
double TemperatureDirichletBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());

  const double time = this->get_time();

  // if (boundary_id_ == 2) {
  //   return g_Tb_well; // boundary value is set to zero in this case
  // } else if (boundary_id_ == 4) {
  //   return g_Tb_top;
  // } else if (boundary_id_ == 3) {
  //   return g_Tb_top + g_T_grad * (0. - p[2]);
  // } else if (boundary_id_ == 5) {
  //   return g_Tb_top + g_T_grad * (0. - p[2]);
  // }
  if (boundary_id_ == 4) {
    return g_Tb_well;
  } else if (boundary_id_ == 3 || boundary_id_ == 8 || boundary_id_ == 12 ||
             boundary_id_ == 13 || boundary_id_ == 14) {
    return g_Tb_top + g_T_grad * (0. - p[2]);
  }
}

// template <int dim>
// void TemperatureDirichletBoundaryValues<dim>::vector_value(
//     const Point<dim>& p, Vector<double>& values) const {
//   for (unsigned int c = 0; c < this->n_components; ++c)
//     values(c) = TemperatureDirichletBoundaryValues<dim>::value(p, c);
// }

template <int dim>
class PressureNeumanBoundaryValues : public Function<dim> {
 private:
  const double period_;  // value
  int boundary_id_{-1};

 public:
  PressureNeumanBoundaryValues()
      : Function<dim>(), period_(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary

  // virtual void vector_value(const Point<dim>& p,  //放在这里没啥用
  //                           Vector<double>& value) const;
  virtual void set_boundary_id(int bnd_id) { boundary_id_ = bnd_id; }
};

template <int dim>
double PressureNeumanBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());
  double time = this->get_time();  // get time
  if (boundary_id_ == 4) {
    return g_Qb_well;
  }
}

template <int dim>
class TemperatureNeumanBoundaryValues : public Function<dim> {
 private:
  const double period;  // value
  int boundary_id_{-1};

 public:
  TemperatureNeumanBoundaryValues()
      : Function<dim>(), period(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary
  // virtual void vector_value(const Point<dim>& p,  //放在这里没啥用
  //                           Vector<double>& value) const;
  virtual void set_boundary_id(int bnd_id) { boundary_id_ = bnd_id; }
};

template <int dim>
double TemperatureNeumanBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());

  const double time = this->get_time();
  if (boundary_id_ == 4) {
    return g_QT_well;  // boundary value is set to zero in
                       // this case
  } else if (boundary_id_ == 3) {
    return g_QT_top;
  } else if (boundary_id_ == 8) {
    return g_QT_bottom;
  }
}

// template <int dim>
// void TemperatureBoundaryValues<dim>::vector_value(
//     const Point<dim>& p, Vector<double>& values) const {
//   for (unsigned int c = 0; c < this->n_components; ++c)
//     values(c) = TemperatureBoundaryValues<dim>::value(p, c);
// }

}  // namespace EquationData
