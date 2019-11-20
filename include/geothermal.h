#ifndef THM_GEOTHERMAL_H_
#define THM_GEOTHERMAL_H_

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

// C++ headers  files:
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

using namespace dealii;

// for setting parameters
namespace EquationData
{
const double eta = 1;
const double kappa = 1e-6;
const double beta = 10;
const double density = 1;

// iniital value for temparature --- value at point, forming vector for diff points
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
  return 0;
}

template <int dim>
void TemperatureInitialValues<dim>::vector_value(const Point<dim> &p,
                                                 Vector<double> &values) const
{
  for (unsigned int c = 0; c < this->n_components; ++c)
    values(c) = TemperatureInitialValues<dim>::value(p, c);
}

// right hand side for temperature, value at point, forming vector for diff points
template <int dim>
class TemperatureRightHandSide : public Function<dim>
{
public:
  TemperatureRightHandSide() : Function<dim>(1) {}

  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const;
};

template <int dim>
double TemperatureRightHandSide<dim>::value(const Point<dim> &p,
                                            const unsigned int component) const
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
double TemperatureBoundaryValues<dim>::value(const Point<dim> &p,
                                             const unsigned int component) const
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

  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const;
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

//boudnary value for Pressure
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

  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const;
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

//boudnary value for Velocity
template <int dim>
class VelocityBoundaryValues : public Function<dim>
{
public:
  VelocityBoundaryValues() : Function<dim>(1) {}
  virtual double value(const Point<dim> &p,
                       const unsigned int component = 0) const;

  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const;
};

template <int dim>
double VelocityBoundaryValues<dim>::value(const Point<dim> &p,
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
  InverseMatrix(const MatrixType &m,
                const PreconditionerType &preconditioner);

  // vector multiply
  template <typename VectorType>
  void vmult(VectorType &dst,
             const VectorType &src) const;

private:
  const SmartPointer<const MatrixType> matrix;
  const PreconditionerType &preconditioner;
};

template <class MatrixType, class PreconditionerType>
InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix(const MatrixType &m,
                                                             const PreconditionerType &preconditioner)
    : matrix(&m),
      preconditioner(preconditioner)
{
}

// using cg to calculate the inverse of the matrix
template <class MatrixType, class PreconditionerType>
template <typename VectorType>
void InverseMatrix<MatrixType, PreconditionerType>::vmult(VectorType &dst,
                                                          const VectorType &src) const
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
  SchurComplement(const BlockSparseMatrix<double> &system_matrix,
                  const InverseMatrix<SparseMatrix<double>, PreconditionerType> &M_inverse);

  void vmult(Vector<double> &dst,
             const Vector<double> &src) const;

private:
  const SmartPointer<const BlockSparseMatrix<double>> system_matrix;
  const SmartPointer<const InverseMatrix<SparseMatrix<double>, PreconditionerType>> M_inverse;

  mutable Vector<double> tmp1, tmp2;
};

template <class PreconditionerType>
SchurComplement<PreconditionerType>::SchurComplement(const BlockSparseMatrix<double> &system_matrix,
                                                     const InverseMatrix<SparseMatrix<double>, PreconditionerType> &M_inverse)
    : system_matrix(&system_matrix),
      M_inverse(&M_inverse),
      tmp1(system_matrix.block(0, 0).m()),
      tmp2(system_matrix.block(0, 0).m())
{
}

template <class PreconditionerType>
void SchurComplement<PreconditionerType>::vmult(Vector<double> &dst,
                                                const Vector<double> &src) const
{
  system_matrix->block(0, 1).vmult(tmp1, src);
  M_inverse->vmult(tmp2, tmp1);
  system_matrix->block(1, 0).vmult(dst, tmp2);
}

} // namespace LinearSolvers

// Geothermal
template <int dim>
class Goethermal
{
public:
  Goethermal();
  void run();

private:
  void make_grid();
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

  double
  compute_viscosity(const std::vector<double> &old_temperature,
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
                    const double cell_diameter) const;

  Triangulation<dim> triangulation;
  double global_Omega_diameter;

  const unsigned int flow_degree;
  FESystem<dim> flow_fe;
  DoFHandler<dim> flow_dof_handler;
  ConstraintMatrix flow_constraints;

  BlockSparseMatrix<double> flow_matrix;
  BlockSparseMatrix<double> flow_preconditioner_matrix;

  BlockVector<double> flow_solution;
  BlockVector<double> old_flow_solution;
  BlockVector<double> flow_rhs;

  const unsigned int temperature_degree;
  FE_Q<dim> temperature_fe;
  DoFHandler<dim> temperature_dof_handler;
  ConstraintMatrix temperature_constraints;

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
Goethermal<dim>::Goethermal()
    : triangulation(Triangulation<dim>::maximum_smoothing),
      global_Omega_diameter(std::numeric_limits<double>::quiet_NaN()),
      flow_degree(1),
      flow_fe(FE_Q<dim>(flow_degree + 1), dim,
              FE_Q<dim>(flow_degree), 1),
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

#endif 
