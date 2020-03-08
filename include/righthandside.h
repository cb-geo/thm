
#pragma once
// #include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <iostream>
#include <fstream>
#include "globalvariables.h"

using namespace dealii;

namespace EquationData
{
template <int dim>
class PressureRightHandSide : public Function<dim>
{
public:
    PressureRightHandSide() : Function<dim>(1), period(0.2) {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const; // assign rhs
                                                                  // equal to zero
                                                                  // at each point

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &value) const;

private:
    const double period;
};

template <int dim>
double
PressureRightHandSide<dim>::value(const Point<dim> &p,
                                  const unsigned int component) const
{
    // (void)component;
    // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
    // Assert(dim == 3, ExcNotImplemented());
    const double time = this->get_time(); // get time
    return 0.;
}

template <int dim>
void PressureRightHandSide<dim>::vector_value(const Point<dim> &p,
                                              Vector<double> &values) const
{
    for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = PressureRightHandSide<dim>::value(p, c);
}

template <int dim>
class TemperatureRightHandSide : public Function<dim>
{
public:
    TemperatureRightHandSide() : Function<dim>(1), period(0.2) {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &value) const;

private:
    const double period;
};

template <int dim>
double
TemperatureRightHandSide<dim>::value(const Point<dim> &p,
                                     const unsigned int component) const
{
    // (void)component;
    // Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
    // Assert(dim == 3, ExcNotImplemented());
    const double time = this->get_time(); // get time
    return 0.;
}

template <int dim>
void TemperatureRightHandSide<dim>::vector_value(const Point<dim> &p,
                                                 Vector<double> &values) const
{
    for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = TemperatureRightHandSide<dim>::value(p, c);
}

} // namespace EquationData