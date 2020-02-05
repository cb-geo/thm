#ifndef THM_GEOTHERMAL_H_
#define THM_GEOTHERMAL_H_

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// C++ headers  files:
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

using namespace dealii;

// for setting parameters
namespace EquationData
{
const double eta = 1;
const double kappa = 1e-6;
const double beta = 10;
const double density = 1;

// iniital value for temparature --- value at point, forming vector for diff
// points
template <int dim>
class TemperatureInitialValues : public Function<dim>
{
public:
  TemperatureInitialValues() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim> &p, Vector<double> &value) const;
};

template <int dim>
double TemperatureInitialValues<dim>::value(const Point<dim> &,
                                            const unsigned int) const
{
  return 0;
}

template <int dim>
void TemperatureInitialValues<dim>::vector_value(const Point<dim> &p,
                                                 Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = TemperatureInitialValues<dim>::value(p, c);
}

// right hand side for temperature, value at point, forming vector for diff
// points
template <int dim>
class TemperatureRightHandSide : public Function<dim>
{
public:
  TemperatureRightHandSide() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim> &p, Vector<double> &value) const;
};

template <int dim>
double TemperatureRightHandSide<dim>::value(
    const Point<dim> &p, const unsigned int component) const
{
  return 0;
}

template <int dim>
void TemperatureRightHandSide<dim>::vector_value(const Point<dim> &p,
                                                 Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = TemperatureRightHandSide<dim>::value(p, c);
}

///////////////////////////////////////////
// boudnary value for temperature
template <int dim>
class TemperatureBoundaryValues : public Function<dim>
{
public:
  TemperatureBoundaryValues() : Function<dim>(1) {}
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;
};

template <int dim>
double TemperatureBoundaryValues<dim>::value(
    const Point<dim> &p, const unsigned int component) const
{
  return 0;
}

// right hand side for pressure, value at point, forming vector for diff points
template <int dim>
class PressureRightHandSide : public Function<dim>
{
public:
  PressureRightHandSide() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim> &p, Vector<double> &value) const;
};

template <int dim>
double PressureRightHandSide<dim>::value(const Point<dim> &p,
                                         const unsigned int component) const
{
  return 0;
}

template <int dim>
void PressureRightHandSide<dim>::vector_value(const Point<dim> &p,
                                              Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = PressureRightHandSide<dim>::value(p, c);
}

// boudnary value for Pressure
template <int dim>
class PressureBoundaryValues : public Function<dim>
{
public:
  PressureBoundaryValues() : Function<dim>(1) {}
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;
};

template <int dim>
double PressureBoundaryValues<dim>::value(const Point<dim> &p,
                                          const unsigned int component) const
{
  return 0;
}

// right hand side for Velocity, value at point, forming vector for diff points
template <int dim>
class VelocityRightHandSide : public Function<dim>
{
public:
  VelocityRightHandSide() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim> &p, Vector<double> &value) const;
};

template <int dim>
double VelocityRightHandSide<dim>::value(const Point<dim> &p,
                                         const unsigned int component) const
{
  return 0;
}

template <int dim>
void VelocityRightHandSide<dim>::vector_value(const Point<dim> &p,
                                              Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = VelocityRightHandSide<dim>::value(p, c);
}

// boudnary value for Velocity
template <int dim>
class VelocityBoundaryValues : public Function<dim>
{
public:
  VelocityBoundaryValues() : Function<dim>(1) {}
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim> &p, Vector<double> &value) const;
};

template <int dim>
double VelocityBoundaryValues<dim>::value(const Point<dim> &p,
                                          const unsigned int component) const
{
  return 0;
}

template <int dim>
void VelocityBoundaryValues<dim>::vector_value(const Point<dim> &p,
                                               Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = VelocityRightHandSide<dim>::value(p, c);
}

//////////////////////////////////////////////

} // namespace EquationData

namespace LinearSolvers
{

// class for calcuating inverse matrix
template <class MatrixType, class PreconditionerType>
class InverseMatrix : public Subscriptor
{
public:
  // inverse matrix
  InverseMatrix(const MatrixType &m, const PreconditionerType &preconditioner);

  // vector multiply
  template <typename VectorType>
  void vmult(VectorType &dst, const VectorType &src) const;

private:
  const SmartPointer<const MatrixType> matrix;
  const PreconditionerType &preconditioner;
};

template <class MatrixType, class PreconditionerType>
InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix(
    const MatrixType &m, const PreconditionerType &preconditioner)
    : matrix(&m), preconditioner(preconditioner) {}

// using cg to calculate the inverse of the matrix
template <class MatrixType, class PreconditionerType>
template <typename VectorType>
void InverseMatrix<MatrixType, PreconditionerType>::vmult(
    VectorType &dst, const VectorType &src) const
{
  SolverControl solver_control(src.size(), 1e-7 * src.l2_norm());
  SolverCG<VectorType> cg(solver_control);

  dst = 0;

  try
  {
    cg.solve(*matrix, dst, src, preconditioner);
  }
  catch (std::exception &e)
  {
    Assert(false, ExcMessage(e.what()));
  }
}

// Schur complement preconditioner
template <class PreconditionerType>
class SchurComplement : public Subscriptor
{
public:
  SchurComplement(
      const BlockSparseMatrix<double> &system_matrix,
      const InverseMatrix<SparseMatrix<double>, PreconditionerType> &M_inverse);

  void vmult(Vector<double> &dst, const Vector<double> &src) const;

private:
  const SmartPointer<const BlockSparseMatrix<double>> system_matrix;
  const SmartPointer<const InverseMatrix<SparseMatrix<double>, PreconditionerType>> M_inverse;

  mutable Vector<double> tmp1, tmp2;
};

template <class PreconditionerType>
SchurComplement<PreconditionerType>::SchurComplement(
    const BlockSparseMatrix<double> &system_matrix,
    const InverseMatrix<SparseMatrix<double>, PreconditionerType> &M_inverse)
    : system_matrix(&system_matrix),
      M_inverse(&M_inverse),
      tmp1(system_matrix.block(0, 0).m()),
      tmp2(system_matrix.block(0, 0).m()) {}

template <class PreconditionerType>
void SchurComplement<PreconditionerType>::vmult(
    Vector<double> &dst, const Vector<double> &src) const
{
  system_matrix->block(0, 1).vmult(tmp1, src);
  M_inverse->vmult(tmp2, tmp1);
  system_matrix->block(1, 0).vmult(dst, tmp2);
}

} // namespace LinearSolvers

// Geothermal
template <int dim>
class Geothermal
{
public:
  Geothermal();
  void run(){};

private:
  void grid_input();
  void setup_dofs();
  void assemble_preconditioner();
  void build_preconditioner();
  void assemble_flow_system();
  void assemble_temperature_system(const double maximal_velocity);
  void assemble_temperature_matrix();
  double get_maximal_velocity() const;
  std::pair<double, double> get_extrapolated_temperature_range() const;
  void solve();
  void output_results() const;
  void refine_mesh(const unsigned int max_grid_level);

  double compute_viscosity(
      const std::vector<double> &old_temperature,
      const std::vector<double> &old_old_temperature,
      const std::vector<Tensor<1, dim>> &old_temperature_grads,
      const std::vector<Tensor<1, dim>> &old_old_temperature_grads,
      const std::vector<double> &old_temperature_laplacians,
      const std::vector<double> &old_old_temperature_laplacians,
      const std::vector<Tensor<1, dim>> &old_velocity_values,
      const std::vector<Tensor<1, dim>> &old_old_velocity_values,
      const std::vector<double> &gamma_values, const double global_u_infty,
      const double global_T_variation, const double cell_diameter) const;

  Triangulation<dim> triangulation;
  double global_Omega_diameter;

  const unsigned int flow_degree;
  FESystem<dim> flow_fe;
  DoFHandler<dim> flow_dof_handler;

  BlockSparseMatrix<double> flow_matrix;
  BlockSparseMatrix<double> flow_preconditioner_matrix;

  BlockVector<double> flow_solution;
  BlockVector<double> old_flow_solution;
  BlockVector<double> flow_rhs;

  const unsigned int temperature_degree;
  FE_Q<dim> temperature_fe;
  DoFHandler<dim> temperature_dof_handler;

  SparseMatrix<double> temperature_mass_matrix;
  SparseMatrix<double> temperature_stiffness_matrix;
  SparseMatrix<double> temperature_matrix;

  Vector<double> temperature_solution;
  Vector<double> old_temperature_solution;
  Vector<double> old_old_temperature_solution;
  Vector<double> temperature_rhs;

  double time_step;
  double old_time_step;
  unsigned int timestep_number;

  bool rebuild_flow_matrix;
  bool rebuild_temperature_matrices;
  bool rebuild_flow_preconditioner;
};

// Geothermal constructor

template <int dim>
Geothermal<dim>::Geothermal()
    : triangulation(Triangulation<dim>::maximum_smoothing),
      global_Omega_diameter(std::numeric_limits<double>::quiet_NaN()),
      flow_degree(1),
      flow_fe(FE_Q<dim>(flow_degree + 1), dim, FE_Q<dim>(flow_degree), 1),
      flow_dof_handler(triangulation),

      temperature_degree(2),
      temperature_fe(temperature_degree),
      temperature_dof_handler(triangulation),

      time_step(0),
      old_time_step(0),
      timestep_number(0),
      rebuild_flow_matrix(true),
      rebuild_temperature_matrices(true),
      rebuild_flow_preconditioner(true)
{
}

template <int dim>
void Geothermal<dim>::grid_input()
{
  GridIn<dim> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f("mesh.msh");
  gridin.read_msh();
  print_mesh_inf(triangulation, "grid-1.eps");
}

template <int dim>
void Geothermal<dim>::setup_dofs()
{

  //flow matrix
  std::vector<unsigned int> flow_sub_block(dim + 1, 0);
  flow_sub_block[dim] = 1;

  flow_dof_handler.distribute_dofs(flow_fe);
  DoFRenumbering::component_wise(flow_dof_handler, flow_sub_block);

  temperature_dof_handler.distribute_dofs(temperature_fe);
  temperature_constraints.clear();

  std::vector<types::global_dof_index> flow_dofs_per_block(2);
  DoFTools::count_dofs_per_block(flow_dof_handler, flow_dofs_per_block, flow_sub_block);

  const unsigned int n_u = flow_dofs_per_block[0],
                     n_p = flow_dofs_per_block[1],
                     n_T = temperature_dof_handler.n_dofs();
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << " (on " << triangulation.n_levels() << " levels)" << std::endl
            << "Number of degrees of freedom: " << n_u + n_p + n_T << " ("
            << n_u << '+' << n_p << '+' << n_T << ')' << std::endl
            << std::endl;
  BlockDynamicSparsityPattern dsp(2, 2);
  dsp.block(0, 0).reinit(n_u, n_u);
  dsp.block(1, 0).reinit(n_p, n_u);
  dsp.block(0, 1).reinit(n_u, n_p);
  dsp.block(1, 1).reinit(n_p, n_p);
  dsp.collect_sizes();

  Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
  for (unsigned int c = 0; c < dim + 1; ++c)
    for (unsigned int d = 0; d < dim + 1; ++d)
      if (!((c == dim) && (d == dim)))
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::none;
  DoFTools::make_sparsity_pattern(
      flow_dof_handler, coupling, dsp, constraints, false);
  sparsity_pattern.copy_from(dsp);
  flow_matrix.reinit(sparsity_pattern);

  //flow preconidtioner matrix
  BlockDynamicSparsityPattern preconditioner_dsp(2, 2);
  preconditioner_dsp.block(0, 0).reinit(n_u, n_u);
  preconditioner_dsp.block(1, 0).reinit(n_p, n_u);
  preconditioner_dsp.block(0, 1).reinit(n_u, n_p);
  preconditioner_dsp.block(1, 1).reinit(n_p, n_p);
  preconditioner_dsp.collect_sizes();
  Table<2, DoFTools::Coupling> preconditioner_coupling(dim + 1, dim + 1);
  for (unsigned int c = 0; c < dim + 1; ++c)
    for (unsigned int d = 0; d < dim + 1; ++d)
      if (((c == dim) && (d == dim)))
        preconditioner_coupling[c][d] = DoFTools::always;
      else
        preconditioner_coupling[c][d] = DoFTools::none;
  DoFTools::make_sparsity_pattern(flow_dof_handler,
                                  preconditioner_coupling,
                                  preconditioner_dsp,
                                  constraints,
                                  false);
  preconditioner_sparsity_pattern.copy_from(preconditioner_dsp);
  flow_preconditioner_matrix.reinit(preconditioner_sparsity_pattern);

  // temperature matrix
  DynamicSparsityPattern dsp(n_T, n_T);
  DoFTools::make_sparsity_pattern(temperature_dof_handler,
                                  dsp,
                                  temperature_constraints,
                                  false);
  temperature_matrix.reinit(dsp);
  temperature_mass_matrix.reinit(temperature_matrix);
  temperature_stiffness_matrix.reinit(temperature_matrix);

  //initialization of solution and rhs
  solution.reinit(2);
  solution.block(0).reinit(n_u);
  solution.block(1).reinit(n_p);
  solution.collect_sizes();
  system_rhs.reinit(2);
  system_rhs.block(0).reinit(n_u);
  system_rhs.block(1).reinit(n_p);
  system_rhs.collect_sizes();
}

template <int dim>
double Geothermal<dim>::get_maximal_velocity() const
{
  const QIterated<dim> quadrature_formula(QTrapez<1>(), flow_degree + 1);
  const unsigned int n_q_points = quadrature_formula.size();
  FEValues<dim> fe_values(flow_fe, quadrature_formula, update_values);
  std::vector<Tensor<1, dim>> velocity_values(n_q_points);
  double max_velocity = 0;
  const FEValuesExtractors::Vector velocities(0);
  for (const auto &cell : flow_dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    fe_values[velocities].get_function_values(stokes_solution,
                                              velocity_values);
    for (unsigned int q = 0; q < n_q_points; ++q)
      max_velocity = std::max(max_velocity, velocity_values[q].norm());
  }
  return max_velocity;
}

#endif

template <int dim>
std::pair<double, double>
Geothermal<dim>::get_extrapolated_temperature_range() const
{
  const QIterated<dim> quadrature_formula(QTrapez<1>(), temperature_degree);
  const unsigned int n_q_points = quadrature_formula.size();
  FEValues<dim> fe_values(temperature_fe, quadrature_formula, update_values);
  std::vector<double> old_temperature_values(n_q_points);
  std::vector<double> old_old_temperature_values(n_q_points);
  if (timestep_number != 0)
  {
    double min_temperature = std::numeric_limits<double>::max(),
           max_temperature = -std::numeric_limits<double>::max();
    for (const auto &cell : temperature_dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(old_temperature_solution,
                                    old_temperature_values);
      fe_values.get_function_values(old_old_temperature_solution,
                                    old_old_temperature_values);
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const double temperature =
            (1. + time_step / old_time_step) * old_temperature_values[q] -
            time_step / old_time_step * old_old_temperature_values[q];
        min_temperature = std::min(min_temperature, temperature);
        max_temperature = std::max(max_temperature, temperature);
      }
    }
    return std::make_pair(min_temperature, max_temperature);
  }
  else
  {
    double min_temperature = std::numeric_limits<double>::max(),
           max_temperature = -std::numeric_limits<double>::max();
    for (const auto &cell : temperature_dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(old_temperature_solution,
                                    old_temperature_values);
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const double temperature = old_temperature_values[q];
        min_temperature = std::min(min_temperature, temperature);
        max_temperature = std::max(max_temperature, temperature);
      }
    }
    return std::make_pair(min_temperature, max_temperature);
  }
}

template <int dim>
double Geothermal<dim>::compute_viscosity(
    const std::vector<double> &old_temperature,
    const std::vector<double> &old_old_temperature,
    const std::vector<Tensor<1, dim>> &old_temperature_grads,
    const std::vector<Tensor<1, dim>> &old_old_temperature_grads,
    const std::vector<double> &old_temperature_laplacians,
    const std::vector<double> &old_old_temperature_laplacians,
    const std::vector<Tensor<1, dim>> &old_velocity_values,
    const std::vector<Tensor<1, dim>> &old_old_velocity_values,
    const std::vector<double> &gamma_values,
    const double global_u_infty,
    const double global_T_variation,
    const double cell_diameter) const
{
  constexpr double beta = 0.017 * dim;
  constexpr double alpha = 1.0;
  if (global_u_infty == 0)
    return 5e-3 * cell_diameter;
  const unsigned int n_q_points = old_temperature.size();
  double max_residual = 0;
  double max_velocity = 0;
  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    const Tensor<1, dim> u =
        (old_velocity_values[q] + old_old_velocity_values[q]) / 2;
    const double dT_dt =
        (old_temperature[q] - old_old_temperature[q]) / old_time_step;
    const double u_grad_T =
        u * (old_temperature_grads[q] + old_old_temperature_grads[q]) / 2;
    const double kappa_Delta_T =
        EquationData::kappa *
        (old_temperature_laplacians[q] + old_old_temperature_laplacians[q]) /
        2;
    const double residual =
        std::abs((dT_dt + u_grad_T - kappa_Delta_T - gamma_values[q]) *
                 std::pow((old_temperature[q] + old_old_temperature[q]) / 2,
                          alpha - 1.));
    max_residual = std::max(residual, max_residual);
    max_velocity = std::max(std::sqrt(u * u), max_velocity);
  }
  const double c_R = std::pow(2., (4. - 2 * alpha) / dim);
  const double global_scaling = c_R * global_u_infty * global_T_variation *
                                std::pow(global_Omega_diameter, alpha - 2.);
  return (
      beta * max_velocity *
      std::min(cell_diameter,
               std::pow(cell_diameter, alpha) * max_residual / global_scaling));
}
