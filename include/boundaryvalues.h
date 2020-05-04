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
class PressureDirichletBoundaryValues : public Function<dim> {
 public:
  PressureDirichletBoundaryValues()
      : Function<dim>(), period(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary

  // virtual void vector_value(const Point<dim>& p,  //放在这里没啥用
  //                           Vector<double>& value) const;
  virtual void set_boundary_id(int bnd_id){
    boundary_id = bnd_id;
  }

 private:
  const double period;  // value
  int boundary_id;
  
};

template <int dim>
double PressureDirichletBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());
  const double time = this->get_time();  // get time
  // return Pb;
  if (boundary_id == 4)
  {
    return Pb_top; 
  } 
  else if (boundary_id == 3)
  {
    return Pb_top + P_grad*(0. - p[2]); 
  } 
  else if (boundary_id == 5)
  {
    return Pb_top + P_grad*(0. - p[2]);
  }

}

template <int dim>
class TemperatureDirichletBoundaryValues : public Function<dim> {
 public:
  TemperatureDirichletBoundaryValues()
      : Function<dim>(), period(0.2) {}  // 之前的没有

  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary
  // virtual void vector_value(const Point<dim>& p,  //放在这里没啥用
  //                           Vector<double>& value) const;
  virtual void set_boundary_id(int bnd_id){
    boundary_id = bnd_id;
  }

 private:
  const double period;  // value
  int boundary_id;
};

template <int dim>
double TemperatureDirichletBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());

  const double time = this->get_time();
  // return T0 + 10. * sin(time / (0.5*3600*24) * 3.1415926);  // boundary value is set to zero in
                                            // this case
  if (boundary_id == 2)
  {
    return Tb_well ;                                          
  } 
  else if (boundary_id == 4)
  {
    return Tb_top; 
  } 
  else if (boundary_id == 3)
  {
    return Tb_top + T_grad*(0. - p[2]);
  } 
  else if (boundary_id == 5)
  {
    return Tb_top + T_grad*(0. - p[2]);
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
 public:
  PressureNeumanBoundaryValues()
      : Function<dim>(), period(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary

  // virtual void vector_value(const Point<dim>& p,  //放在这里没啥用
  //                           Vector<double>& value) const;

 private:
  const double period;  // value
};

template <int dim>
double PressureNeumanBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());
  double time = this->get_time();  // get time
  return Qb_well;
}

template <int dim>
class TemperatureNeumanBoundaryValues : public Function<dim> {
 public:
  TemperatureNeumanBoundaryValues()
      : Function<dim>(), period(0.2) {}  // 之前的没有
  virtual double value(const Point<dim>& p,
                       const unsigned int component = 0) const;  // boundary
  // virtual void vector_value(const Point<dim>& p,  //放在这里没啥用
  //                           Vector<double>& value) const;

 private:
  const double period;  // value
};

template <int dim>
double TemperatureNeumanBoundaryValues<dim>::value(
    const Point<dim>& p, const unsigned int /*component*/) const {

  // (void)component;
  // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  // Assert(dim == 3, ExcNotImplemented());

  const double time = this->get_time();
  return 0.;  // boundary value is set to zero in
              // this case
}

// template <int dim>
// void TemperatureBoundaryValues<dim>::vector_value(
//     const Point<dim>& p, Vector<double>& values) const {
//   for (unsigned int c = 0; c < this->n_components; ++c)
//     values(c) = TemperatureBoundaryValues<dim>::value(p, c);
// }

}  // namespace EquationData
