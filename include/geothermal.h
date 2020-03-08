#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include <map>

namespace Geothermal
{
using namespace dealii;

namespace EquationData
{
const double perm = 1e-9;
const double eta = 1;
const double kappa = 1e-6;
const double beat = 10;
const double density = 1;
const double T0 = 278.15;
const double p0 = 100000.;

template <int dim>
class PressureInitialValues : public Function<dim>
{
public:
  PressureInitialValues() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;
  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const;
};

template <int dim>
double
PressureInitialValues<dim>::value(const Point<dim> &,
                                  const unsigned int) const
{
  return p0;
}

template <int dim>
void PressureInitialValues<dim>::vector_value(const Point<dim> &p,
                                              Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = PressureInitialValues<dim>::value(p, c);
}

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
  (void)component;
  Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  Assert(dim == 3, ExcNotImplemented());
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
class PressureBoundaryValues : public Function<dim>
{
public:
  PressureBoundaryValues() : Function<dim>(1), period(0.2) {} // 之前的没有
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const; // boundary
  virtual void vector_value(const Point<dim> &p,                //之前的没有
                            Vector<double> &value) const;

private:
  const double period; // value
};

template <int dim>
double PressureBoundaryValues<dim>::value(const Point<dim> &p,
                                          const unsigned int /*component*/) const
{

  (void)component;
  Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  Assert(dim == 3, ExcNotImplemented());
  const double time = this->get_time(); // get time
  return p0;
}

template <int dim>
class TemperatureInitialValues : public Function<dim>
{
public:
  TemperatureInitialValues() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;
  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const;
};

template <int dim>
double
TemperatureInitialValues<dim>::value(const Point<dim> &,
                                     const unsigned int) const
{
  return T0;
}

template <int dim>
void TemperatureInitialValues<dim>::vector_value(const Point<dim> &p,
                                                 Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = TemperatureInitialValues<dim>::value(p, c);
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
  (void)component;
  Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  Assert(dim == 3, ExcNotImplemented());
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

template <int dim>
class TemperatureBoundaryValues : public Function<dim>
{
public:
  TemperatureBoundaryValues() : Function<dim>(1) {} // 之前的没有
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const; // boundary
  virtual void vector_value(const Point<dim> &p,                //之前的没有
                            Vector<double> &value) const;

private:
  const double period; // value
};

template <int dim>
double TemperatureBoundaryValues<dim>::value(const Point<dim> &p,
                                          const unsigned int /*component*/) const
{

  (void)component;
  Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
  Assert(dim == 3, ExcNotImplemented());

  const double time = this->get_time();
  return T0 + 10. *
                  sin(time * 3.1415926); // boundary value is set to zero in this case
}

template <int dim>
void TemperatureBoundaryValues<dim>::vector_value(const Point<dim> &p,
                                                 Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = TemperatureBoundaryValues<dim>::value(p, c);
}

} // namespace EquationData

template <int dim>
class HeatEquation
{
public:
  HeatEquation();
  void run();

private:
  void grid_input();
  void setup_system();
  void solve_time_step();
  void output_results() const;

  Triangulation<dim> triangulation; // grid
  FE_Q<dim> fe;                     // element
  DoFHandler<dim> dof_handler;      // grid<->eleemnt

  ConstraintMatrix constraints; // hanging node

  SparsityPattern sparsity_pattern;    // sparsity
  SparseMatrix<double> mass_matrix;    // M
  SparseMatrix<double> laplace_matrix; // A
  SparseMatrix<double> system_matrix;  // M + k*theta*A

  Vector<double> solution;     // solution at n
  Vector<double> old_solution; // solution at n-1
  Vector<double> system_rhs;   // rhs

  double time;
  double time_step;
  unsigned int timestep_number;

  const double theta;
};

// template <int dim>
// class RightHandSide : public Function<dim>
// {
// public:
//   RightHandSide() : Function<dim>(), period(0.2) {}
//   virtual double value(const Point<dim> &p,
//                        const unsigned int component = 0) const; // assign rhs
//                                                                 // equal to zero
//                                                                 // at each point

// private:
//   const double period;
// };

// template <int dim>
// class BoundaryValues : public Function<dim>
// {
// public:
//   virtual double value(const Point<dim> &p,
//                        const unsigned int component = 0) const; // boundary
//                                                                 // value
// };

// template <int dim>
// double RightHandSide<dim>::value(const Point<dim> &p,
//                                  const unsigned int component) const
// {
//   (void)component;
//   Assert(component == 0, ExcIndexRange(component, 0, 1)); // for debug
//   Assert(dim == 3, ExcNotImplemented());

//   const double time = this->get_time(); // get time
//   return 0;
// }

// template <int dim>
// double BoundaryValues<dim>::value(const Point<dim> &p,
//                                   const unsigned int /*component*/) const
// {
//   // (void)component;
//   // Assert(component == 0, ExcIndexRange(component, 0, 1));
//   const double time = this->get_time();
//   return 10. *
//          sin(time * 3.1415926); // boundary value is set to zero in this case
// }

template <int dim>
HeatEquation<dim>::HeatEquation() // initialization
    : fe(1),
      dof_handler(triangulation),
      time(0.0),
      time_step(1. / 20), // a time step constant at 1/500 (remember that one
                          // period of the source on the right hand side was set
                          // to 0.2 above, so we resolve each period with 100
                          // time steps)
      timestep_number(0),
      theta(0.5)
{
}

template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;
  {
    std::map<unsigned int, unsigned int> boundary_count;
    typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();
    for (; cell != endc; ++cell)
    {
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
      {
        if (cell->face(face)->at_boundary())
          boundary_count[cell->face(face)->boundary_id()]++;
      }
    }

    std::cout << " boundary indicators: ";
    for (std::map<unsigned int, unsigned int>::iterator it =
             boundary_count.begin();
         it != boundary_count.end(); ++it)
    {
      std::cout << it->first << "(" << it->second << " times) ";
    }
    std::cout << std::endl;
  }

  std::ofstream out(filename.c_str());
  GridOut grid_out;
  grid_out.write_eps(triangulation, out);
  std::cout << " written to " << filename << std::endl
            << std::endl;
}

template <int dim>
void HeatEquation<dim>::grid_input()
{
  GridIn<dim> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f("mesh.msh");
  gridin.read_msh(f);

  print_mesh_info(triangulation, "grid-1.eps");
}

template <int dim>
void HeatEquation<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe); // distribute dofs to grid

  std::cout << std::endl
            << "===========================================" << std::endl
            << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl
            << std::endl;

  constraints.clear();
  DoFTools::make_hanging_node_constraints(
      dof_handler,
      constraints); // setting the hanging node
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs()); // sparsity
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints,
                                  /*keep_constrained_dofs = */ true);
  sparsity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparsity_pattern);    // initialize M using given sparsity
                                           // parttern
  laplace_matrix.reinit(sparsity_pattern); // initialize A using given sparsity
                                           // parttern
  system_matrix.reinit(sparsity_pattern);  // initialize M + k*theta*A using
                                           // given sparsity parttern

  MatrixCreator::create_mass_matrix(
      dof_handler, QGauss<dim>(fe.degree + 1),
      mass_matrix); // Assemble the mass matrix and a right hand side vector.
                    // If no coefficient is given (i.e., if the pointer to a
                    // function object is zero as it is by default), the
                    // coefficient is taken as being constant and equal to one.
  MatrixCreator::create_laplace_matrix(
      dof_handler, QGauss<dim>(fe.degree + 1),
      laplace_matrix); // Assemble the Laplace matrix.
                       // If no coefficient is given (i.e., if the pointer to a
                       // function object is zero as it is by default), the
                       // coefficient is taken as being constant and equal to
                       // one. In case you want to specify constraints and use
                       // the default argument for the coefficient you have to
                       // specify the (unused) coefficient argument as (const
                       // Function<spacedim> *const)nullptr.

  solution.reinit(dof_handler.n_dofs());
  old_solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim>
void HeatEquation<dim>::solve_time_step()
{
  SolverControl solver_control(1000,
                               1e-8 * system_rhs.l2_norm()); // setting for cg
  SolverCG<> cg(solver_control);                             // config cg

  PreconditionSSOR<> preconditioner;             // precond
  preconditioner.initialize(system_matrix, 1.0); // initialize precond

  cg.solve(system_matrix, solution, system_rhs,
           preconditioner); // solve eq

  constraints.distribute(solution); // make sure if the value is consistent at
                                    // the constraint point

  std::cout << "     " << solver_control.last_step() << " CG iterations."
            << std::endl;
}

// @sect4{<code>HeatEquation::output_results</code>}
//
// Neither is there anything new in generating graphical output:
template <int dim>
void HeatEquation<dim>::output_results() const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "U");

  data_out.build_patches();

  const std::string filename =
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
}

template <int dim>
void HeatEquation<dim>::run()
{
  grid_input();
  setup_system();

  Vector<double> tmp;           // this vector is for
  Vector<double> forcing_terms; // this vector is for forcing_terms

  tmp.reinit(solution.size());           // initialize tmp
  forcing_terms.reinit(solution.size()); // initialize forcing_terms

  VectorTools::interpolate(
      dof_handler, ZeroFunction<dim>(),
      old_solution); // interpolate the old solution based on dof_handler, here
                     // using interpolation because we refine the global
  solution =
      old_solution; // updating the solutin with sinterpolated old solution

  output_results(); // output

  while (time <= 1.)
  {
    time += time_step; // get new time point
    ++timestep_number; // get the No. of the new time step

    std::cout << "Time step " << timestep_number << " at t=" << time
              << std::endl;

    mass_matrix.vmult(system_rhs,
                      old_solution); // matrix multiplication system_rhs =
                                     // mass_matrix*old_solution

    laplace_matrix.vmult(tmp,
                         old_solution); // tmp = laplace_matrix*old_solution
    system_rhs.add(
        -(1 - theta) * time_step,
        tmp); // system_rhs = system_rhs -(1 - theta) * time_step*tmp
              // 注意，这里system_rhs是一个vector，所以add是将两个元素相乘了

    RightHandSide<dim> rhs_function;
    rhs_function.set_time(time);
    VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe.degree + 1),
                                        rhs_function, tmp);
    forcing_terms = tmp;
    forcing_terms *= time_step * theta;

    rhs_function.set_time(time - time_step);
    VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe.degree + 1),
                                        rhs_function, tmp);

    forcing_terms.add(
        time_step * (1 - theta),
        tmp); // 形成forcing term = f(x,n-1)*(1-theta) + f(x,n)*theta

    system_rhs +=
        forcing_terms; // system_rhs = system_rhs + forcing_term :
                       // sys_Old*U_Old + f(x,n-1)*(1-theta) + f(x,n)*theta

    system_matrix.copy_from(mass_matrix);
    system_matrix.add(theta * time_step, laplace_matrix); // sys = M +
                                                          // k*theta*A

    constraints.condense(system_matrix, system_rhs); // 压缩

    {
      EquationData::TemperatureBoundaryValues<dim> boundary_values_function; // creat boundary value
                                                    // object
      boundary_values_function.set_time(time);      // set the proper time

      std::map<types::global_dof_index, double> boundary_values;
      VectorTools::interpolate_boundary_values(
          dof_handler, // evaluate value by interpolation
          1, boundary_values_function, boundary_values);

      MatrixTools::apply_boundary_values(boundary_values, system_matrix,
                                         solution, system_rhs);
    }

    solve_time_step();

    output_results();

    old_solution = solution;
  }
}
} // namespace Geothermal
