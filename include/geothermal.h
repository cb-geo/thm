#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <vector>

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
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

using namespace dealii;

template <int dim>
class CoupledTH {
 public:
  CoupledTH(const unsigned int degree);
  void run();

 private:
  void make_grid_and_dofs();
  void assemble_T_system();
  void assemble_P_system();
  void linear_solve_P();
  void linear_solve_T();
  void output_results(Vector<double>&, std::string) const;

  Triangulation<dim> triangulation;  // grid
  const unsigned int T_degree;       // element degree
  const unsigned int P_degree;

  FE_Q<dim> P_fe;               // element type
  FE_Q<dim> T_fe;               // element type
  DoFHandler<dim> dof_handler;  // grid<->eleemnt

  // ConstraintMatrix constraints;  // hanging node

  SparsityPattern sparsity_pattern;      // sparsity
  SparseMatrix<double> P_system_matrix;  // M_P + K_P
  SparseMatrix<double> T_system_matrix;  // M_T + K_T

  Vector<double> P_solution;      // P solution at n
  Vector<double> T_solution;      // T solution at n
  Vector<double> old_P_solution;  // P solution at n-1
  Vector<double> old_T_solution;  // T solution at n-1
  Vector<double> P_system_rhs;    // right hand side of P system
  Vector<double> T_system_rhs;    // right hand side of T system

  double time;                   // t
  unsigned int timestep_number;  // n_t
  std::vector<double> time_sequence;

  double period;
  int n_time_step;
  double time_step;

  unsigned int P_iteration_namber;
  unsigned int T_iteration_namber;
  unsigned int n_P_max_iteration = EquationData::n_g_P_max_iteration;
  unsigned int n_T_max_iteration = EquationData::n_g_T_max_iteration;
  double P_tol_residual = EquationData::g_P_tol_residual;
  double T_tol_residual = EquationData::g_T_tol_residual;
};

template <int dim>
CoupledTH<dim>::CoupledTH(const unsigned int degree)  // initialization
    : T_degree(degree),
      P_degree(degree),
      P_fe(P_degree),
      T_fe(T_degree),
      dof_handler(triangulation),

      time(0.0),
      timestep_number(0),
      T_iteration_namber(0),
      P_iteration_namber(0) {
  if (EquationData::is_linspace) {
    period = EquationData::g_period;
    n_time_step = EquationData::g_n_time_step;
    time_sequence = linspace(0.0, period, n_time_step);
    time_step = time_sequence[1] - time_sequence[0];
  } else {
    time_sequence = EquationData::g_time_sequence;
    n_time_step = time_sequence.size();
    period = time_sequence[n_time_step - 1];
    time_step = time_sequence[1] - time_sequence[0];
  }
}

// CHECK THE GRID INPUT
template <int dim>
void print_mesh_info(const Triangulation<dim>& triangulation,
                     const std::string& filename) {
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;
  {
    std::map<unsigned int, unsigned int> boundary_count;
    typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();
    for (; cell != endc; ++cell) {
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face) {
        if (cell->face(face)->at_boundary())
          boundary_count[cell->face(face)->boundary_id()]++;
      }
    }

    std::cout << " boundary indicators: ";
    for (std::map<unsigned int, unsigned int>::iterator it =
             boundary_count.begin();
         it != boundary_count.end(); ++it) {
      std::cout << it->first << "(" << it->second << " times) ";
    }
    std::cout << std::endl;
  }

  std::ofstream out(filename.c_str());
  GridOut grid_out;
  grid_out.write_eps(triangulation, out);
  std::cout << " written to " << filename << std::endl << std::endl;
}

template <int dim>
void CoupledTH<dim>::make_grid_and_dofs() {

  cbgeo::Clock timer;
  timer.tick();
  GridIn<dim> gridin;  // instantiate a gridinput
  gridin.attach_triangulation(triangulation);
  std::ifstream f("inputfiles/mesh.msh");
  gridin.read_msh(f);

  print_mesh_info(triangulation, "outputfiles/grid-1.eps");

  dof_handler.distribute_dofs(P_fe);  // distribute dofs to grid
  dof_handler.distribute_dofs(T_fe);  // distribute dofs to grid

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << " (on " << triangulation.n_levels() << " levels)" << std::endl
            << "Number of degrees of freedom: " << 2 * dof_handler.n_dofs()
            << " (" << dof_handler.n_dofs() << '+' << dof_handler.n_dofs()
            << ')' << std::endl
            << std::endl;

  // sparsity pattern for T
  DynamicSparsityPattern T_dsp(dof_handler.n_dofs());  // sparsity
  DoFTools::make_sparsity_pattern(
      dof_handler, T_dsp /*constraints, keep_constrained_dofs = ,true*/);
  sparsity_pattern.copy_from(T_dsp);

  // forming system matrixes and initialize these matrixes
  T_system_matrix.reinit(sparsity_pattern);
  T_system_rhs.reinit(dof_handler.n_dofs());
  T_solution.reinit(dof_handler.n_dofs());
  old_T_solution.reinit(dof_handler.n_dofs());

  // sparsity pattern for P
  DynamicSparsityPattern P_dsp(dof_handler.n_dofs());  // sparsity
  DoFTools::make_sparsity_pattern(
      dof_handler, P_dsp /*constraints, keep_constrained_dofs = ,true*/);
  sparsity_pattern.copy_from(P_dsp);

  // forming system matrixes and initialize these matrixes
  P_system_matrix.reinit(sparsity_pattern);
  P_system_rhs.reinit(dof_handler.n_dofs());
  P_solution.reinit(dof_handler.n_dofs());
  old_P_solution.reinit(dof_handler.n_dofs());

  // constraints.clear();
  // DoFTools::make_hanging_node_constraints(
  //     dof_handler,
  //     constraints); // setting the hanging node
  // constraints.close();

  timer.tock("Grid_and_dofs");
}

template <int dim>
void CoupledTH<dim>::assemble_P_system() {
  cbgeo::Clock timer;
  timer.tick();

  // reset matreix to zero
  P_system_matrix = 0;
  P_system_rhs = 0;
  P_solution = 0;

  // Getting T values
  QGauss<dim> T_quadrature_formula(T_degree + 1);

  QGauss<dim - 1> T_face_quadrature_formula(T_degree + 1);

  FEValues<dim> T_fe_values(T_fe, T_quadrature_formula,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

  FEFaceValues<dim> T_fe_face_values(T_fe, T_face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                         update_quadrature_points |
                                         update_JxW_values);

  // Getting P values
  QGauss<dim> P_quadrature_formula(P_degree + 1);

  QGauss<dim - 1> P_face_quadrature_formula(P_degree + 1);

  FEValues<dim> P_fe_values(P_fe, P_quadrature_formula,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

  FEFaceValues<dim> P_fe_face_values(P_fe, P_face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                         update_quadrature_points |
                                         update_JxW_values);

  // define loop number for P
  const unsigned int P_dofs_per_cell = P_fe.dofs_per_cell;
  const unsigned int P_n_q_points = P_quadrature_formula.size();
  const unsigned int P_n_face_q_points = P_face_quadrature_formula.size();

  // define loop number for T
  const unsigned int T_dofs_per_cell = T_fe.dofs_per_cell;
  const unsigned int T_n_q_points = T_quadrature_formula.size();
  const unsigned int T_n_face_q_points = T_face_quadrature_formula.size();

  // store the value at previous step at q_point for T
  std::vector<double> old_T_sol_values(T_n_q_points);
  std::vector<Tensor<1, dim>> old_T_sol_grads(T_n_q_points);

  // store the value at previous step at q_point for P
  std::vector<double> old_P_sol_values(P_n_q_points);
  std::vector<Tensor<1, dim>> old_P_sol_grads(P_n_q_points);

  // store the rhs and bd and old solution value at q_point of element for P
  std::vector<double> P_source_values(P_n_q_points);
  std::vector<double> QP_bd_values(P_n_face_q_points);

  // store the coordinate of gauss point
  // const Point<dim> P_quadrature_coord;

  //  local element matrix
  FullMatrix<double> P_local_mass_matrix(P_dofs_per_cell, P_dofs_per_cell);
  FullMatrix<double> P_local_stiffness_matrix(P_dofs_per_cell, P_dofs_per_cell);
  Vector<double> P_local_rhs(P_dofs_per_cell);
  std::vector<types::global_dof_index> P_local_dof_indices(P_dofs_per_cell);

  // boudnary condition and source term
  EquationData::PressureSourceTerm<dim> P_source_term;
  EquationData::PressureNeumanBoundaryValues<dim> QP_boundary;
  EquationData::PressureDirichletBoundaryValues<dim> P_boundary;

  // loop for cell
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell) {
    // initialization
    P_local_mass_matrix = 0;
    P_local_stiffness_matrix = 0;
    P_local_rhs = 0;
    T_fe_values.reinit(cell);
    P_fe_values.reinit(cell);

    // get the values at gauss point old_T_sol_values from the system
    // old_T_solution
    T_fe_values.get_function_values(old_T_solution, old_T_sol_values);
    T_fe_values.get_function_gradients(old_T_solution, old_T_sol_grads);
    P_fe_values.get_function_values(old_P_solution, old_P_sol_values);
    P_fe_values.get_function_gradients(old_P_solution, old_P_sol_grads);

    // get source term value at the gauss point
    P_source_term.set_time(time);
    P_source_term.value_list(P_fe_values.get_quadrature_points(),
                             P_source_values);  // 一列q个

    // loop for q_point ASSMBLING CELL METRIX (weak form equation writing)
    for (unsigned int q = 0; q < P_n_q_points; ++q) {

      const auto P_quadrature_coord = P_fe_values.quadrature_point(q);
      EquationData::g_perm = interpolate1d(
          EquationData::g_perm_list, P_quadrature_coord[2], false);  // step-5

      for (unsigned int i = 0; i < P_dofs_per_cell; ++i) {
        const Tensor<1, dim> grad_phi_i_P = P_fe_values.shape_grad(i, q);
        const double phi_i_P = P_fe_values.shape_value(i, q);
        for (unsigned int j = 0; j < P_dofs_per_cell; ++j) {
          const Tensor<1, dim> grad_phi_j_P = P_fe_values.shape_grad(j, q);
          const double phi_j_P = P_fe_values.shape_value(j, q);
          P_local_mass_matrix(i, j) += (phi_i_P * phi_j_P * P_fe_values.JxW(q));
          P_local_stiffness_matrix(i, j) +=
              time_step * (EquationData::g_perm * EquationData::g_B_w *
                           grad_phi_i_P * grad_phi_j_P * P_fe_values.JxW(q));
        }
        P_local_rhs(i) += (time_step * phi_i_P * P_source_values[q] +
                           time_step * grad_phi_i_P * (Point<dim>(0, 0, 1)) *
                               (-EquationData::g_B_w * EquationData::g_perm *
                                EquationData::g_P_grad) +
                           phi_i_P * old_P_sol_values[q]) *
                          P_fe_values.JxW(q);
      }
    }

    // APPLIED NEWMAN BOUNDARY CONDITION
    for (int bd_i = 0; bd_i < EquationData::g_num_QP_bnd_id; bd_i++) {
      for (unsigned int face_no = 0;
           face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        if (cell->at_boundary(face_no) && cell->face(face_no)->boundary_id() ==
                                              EquationData::g_QP_bnd_id[bd_i]) {
          P_fe_face_values.reinit(cell, face_no);

          // get boundary condition
          QP_boundary.set_time(time);
          QP_boundary.set_boundary_id(*(EquationData::g_QP_bnd_id + bd_i));
          QP_boundary.value_list(P_fe_face_values.get_quadrature_points(),
                                 QP_bd_values);

          for (unsigned int q = 0; q < P_n_face_q_points; ++q) {

            const auto P_face_quadrature_coord =
                P_fe_face_values.quadrature_point(q);
            EquationData::g_perm =
                interpolate1d(EquationData::g_perm_list,
                              P_face_quadrature_coord[2], false);  // step-5

            for (unsigned int i = 0; i < P_dofs_per_cell; ++i) {
              P_local_rhs(i) += -time_step * EquationData::g_B_w *
                                (P_fe_face_values.shape_value(i, q) *
                                 QP_bd_values[q] * P_fe_face_values.JxW(q));
            }
          }
        }
      }
    }

    // local ->globe
    cell->get_dof_indices(P_local_dof_indices);

    for (unsigned int i = 0; i < P_dofs_per_cell; ++i) {
      for (unsigned int j = 0; j < P_dofs_per_cell; ++j) {
        P_system_matrix.add(P_local_dof_indices[i], P_local_dof_indices[j],
                            P_local_mass_matrix(i, j));
        P_system_matrix.add(P_local_dof_indices[i], P_local_dof_indices[j],
                            P_local_stiffness_matrix(i, j));
      }
      P_system_rhs(P_local_dof_indices[i]) += P_local_rhs(i);
    }
  }

  // // ADD DIRICLET BOUNDARY
  {

    for (int bd_i = 0; bd_i < EquationData::g_num_P_bnd_id; bd_i++) {

      P_boundary.set_time(time);
      P_boundary.set_boundary_id(*(EquationData::g_P_bnd_id + bd_i));
      std::map<types::global_dof_index, double> P_bd_values;
      VectorTools::interpolate_boundary_values(
          dof_handler, *(EquationData::g_P_bnd_id + bd_i), P_boundary,
          P_bd_values);  // i表示边界的index
      MatrixTools::apply_boundary_values(P_bd_values, P_system_matrix,
                                         P_solution, P_system_rhs);
    }
  }
  timer.tock("assemble_P_system");
}

template <int dim>
void CoupledTH<dim>::assemble_T_system() {
  cbgeo::Clock timer;
  timer.tick();
  // reset matreix to zero
  T_system_matrix = 0;
  T_system_rhs = 0;
  T_solution = 0;

  // Getting T values
  QGauss<dim> T_quadrature_formula(T_degree + 1);

  QGauss<dim - 1> T_face_quadrature_formula(T_degree + 1);

  FEValues<dim> T_fe_values(T_fe, T_quadrature_formula,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

  FEFaceValues<dim> T_fe_face_values(T_fe, T_face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                         update_quadrature_points |
                                         update_JxW_values);

  // Getting P values
  QGauss<dim> P_quadrature_formula(P_degree + 1);

  QGauss<dim - 1> P_face_quadrature_formula(P_degree + 1);

  FEValues<dim> P_fe_values(P_fe, P_quadrature_formula,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

  FEFaceValues<dim> P_fe_face_values(P_fe, P_face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                         update_quadrature_points |
                                         update_JxW_values);

  // define loop number
  const unsigned int P_dofs_per_cell = P_fe.dofs_per_cell;
  const unsigned int P_n_q_points = P_quadrature_formula.size();
  const unsigned int P_n_face_q_points = P_face_quadrature_formula.size();
  const unsigned int T_dofs_per_cell = T_fe.dofs_per_cell;
  const unsigned int T_n_q_points = T_quadrature_formula.size();
  const unsigned int T_n_face_q_points = T_face_quadrature_formula.size();

  // store the value at previous step at q_point for T
  std::vector<double> old_T_sol_values(T_n_q_points);
  std::vector<Tensor<1, dim>> old_T_sol_grads(T_n_q_points);

  // store the value at previous step at q_point for P
  std::vector<double> old_P_sol_values(P_n_q_points);
  std::vector<Tensor<1, dim>> old_P_sol_grads(P_n_q_points);

  // store the source and bd value at q_point
  std::vector<double> T_source_values(T_n_q_points);
  std::vector<double> QT_bd_values(T_n_face_q_points);

  //  local element matrix
  FullMatrix<double> T_local_mass_matrix(T_dofs_per_cell, T_dofs_per_cell);
  FullMatrix<double> T_local_stiffness_matrix(T_dofs_per_cell, T_dofs_per_cell);
  FullMatrix<double> T_local_convection_matrix(T_dofs_per_cell,
                                               T_dofs_per_cell);
  Vector<double> T_local_rhs(T_dofs_per_cell);
  std::vector<types::global_dof_index> T_local_dof_indices(T_dofs_per_cell);

  // boudnary condition
  EquationData::TemperatureSourceTerm<dim> T_source_term;
  EquationData::TemperatureNeumanBoundaryValues<dim> QT_boundary;
  EquationData::TemperatureDirichletBoundaryValues<dim> T_boundary;

  // loop for cell
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell) {
    // initialization
    T_local_mass_matrix = 0;
    T_local_stiffness_matrix = 0;
    T_local_convection_matrix = 0;
    T_local_rhs = 0;
    T_fe_values.reinit(cell);
    P_fe_values.reinit(cell);  //// may combine in the same cell
    // get the values at gauss point
    T_fe_values.get_function_values(old_T_solution, old_T_sol_values);
    T_fe_values.get_function_gradients(old_T_solution, old_T_sol_grads);
    P_fe_values.get_function_values(old_P_solution, old_P_sol_values);
    P_fe_values.get_function_gradients(old_P_solution, old_P_sol_grads);
    // get right hand side
    T_source_term.set_time(time);
    T_source_term.value_list(T_fe_values.get_quadrature_points(),
                             T_source_values);

    // loop for q_point ASSMBLING CELL METRIX (weak form equation writing)
    for (unsigned int q = 0; q < T_n_q_points; ++q) {
      const auto T_quadrature_coord = T_fe_values.quadrature_point(q);
      EquationData::g_perm = interpolate1d(
          EquationData::g_perm_list, T_quadrature_coord[2], false);  // step-5

      for (unsigned int i = 0; i < T_dofs_per_cell; ++i) {
        const Tensor<1, dim> grad_phi_i_T = T_fe_values.shape_grad(i, q);
        const double phi_i_T = T_fe_values.shape_value(i, q);
        for (unsigned int j = 0; j < T_dofs_per_cell; ++j) {
          const Tensor<1, dim> grad_phi_j_T = T_fe_values.shape_grad(j, q);
          const double phi_j_T = T_fe_values.shape_value(j, q);
          T_local_mass_matrix(i, j) += (phi_i_T * phi_j_T * T_fe_values.JxW(q));
          T_local_stiffness_matrix(i, j) +=
              time_step * (EquationData::g_lam / EquationData::g_c_T *
                           grad_phi_i_T * grad_phi_j_T * T_fe_values.JxW(q));
          T_local_convection_matrix(i, j) +=
              time_step * EquationData::g_c_w / EquationData::g_c_T * phi_i_T *
              (-EquationData::g_perm *
               (old_P_sol_grads[q] +
                (Point<dim>(0, 0, 1)) * EquationData::g_P_grad) *
               grad_phi_j_T * T_fe_values.JxW(q));
        }

        T_local_rhs(i) +=
            (time_step * T_source_values[q] + old_T_sol_values[q]) * phi_i_T *
            T_fe_values.JxW(q);
      }
    }

    // APPLIED NEUMAN BOUNDARY CONDITION
    for (int bd_i = 0; bd_i < EquationData::g_num_QT_bnd_id; bd_i++) {
      for (unsigned int face_no = 0;
           face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        if (cell->at_boundary(face_no) && cell->face(face_no)->boundary_id() ==
                                              EquationData::g_QT_bnd_id[bd_i]) {
          T_fe_face_values.reinit(cell, face_no);

          // get boundary condition
          QT_boundary.set_time(time);
          QT_boundary.set_boundary_id(*(EquationData::g_QT_bnd_id + bd_i));
          QT_boundary.value_list(T_fe_face_values.get_quadrature_points(),
                                 QT_bd_values);

          for (unsigned int q = 0; q < T_n_face_q_points; ++q) {

            const auto T_face_quadrature_coord =
                T_fe_face_values.quadrature_point(q);
            EquationData::g_perm =
                interpolate1d(EquationData::g_perm_list,
                              T_face_quadrature_coord[2], false);  // step-5

            for (unsigned int i = 0; i < T_dofs_per_cell; ++i) {
              T_local_rhs(i) += -time_step / EquationData::g_c_T *
                                T_fe_face_values.shape_value(i, q) *
                                QT_bd_values[q] * T_fe_face_values.JxW(q);
            }
          }
        }
      }
    }
    // local ->globe
    cell->get_dof_indices(T_local_dof_indices);
    for (unsigned int i = 0; i < T_dofs_per_cell; ++i) {
      for (unsigned int j = 0; j < T_dofs_per_cell; ++j) {
        T_system_matrix.add(T_local_dof_indices[i], T_local_dof_indices[j],
                            T_local_mass_matrix(i, j));  // sys_mass_matrix
        T_system_matrix.add(
            T_local_dof_indices[i], T_local_dof_indices[j],
            T_local_stiffness_matrix(i, j));  // sys_stiff_matrix
        T_system_matrix.add(
            T_local_dof_indices[i], T_local_dof_indices[j],
            T_local_convection_matrix(i, j));  // sys_convect_matrix
      }
      T_system_rhs(T_local_dof_indices[i]) += T_local_rhs(i);
    }
  }

  // ADD DIRICHLET BOUNDARY
  {

    for (int bd_i = 0; bd_i < EquationData::g_num_T_bnd_id; bd_i++) {

      T_boundary.set_time(time);
      T_boundary.set_boundary_id(*(EquationData::g_T_bnd_id + bd_i));
      std::map<types::global_dof_index, double> T_bd_values;
      VectorTools::interpolate_boundary_values(
          dof_handler, *(EquationData::g_T_bnd_id + bd_i), T_boundary,
          T_bd_values);  // i表示边界的index
      MatrixTools::apply_boundary_values(T_bd_values, T_system_matrix,
                                         T_solution, T_system_rhs);
    }
  }
  timer.tock("assemble_T_system");
}

template <int dim>
void CoupledTH<dim>::linear_solve_P() {
  cbgeo::Clock timer;
  timer.tick();
  SolverControl solver_control(
      n_P_max_iteration,
      P_tol_residual * P_system_rhs.l2_norm());     // setting for cg
  SolverCG<> cg(solver_control);                    // config cg
  PreconditionSSOR<> preconditioner;                // precond
  preconditioner.initialize(P_system_matrix, 1.0);  // initialize precond
  cg.solve(P_system_matrix, P_solution, P_system_rhs,
           preconditioner);  // solve eq

  P_iteration_namber = solver_control.last_step();
  // constraints.distribute(solution);  // make sure if the value is
  // consistent at
  // the constraint point

  std::cout << "\nCG iterations: " << P_iteration_namber << std::endl;

  timer.tock("linear_solve_P");
}

template <int dim>
void CoupledTH<dim>::linear_solve_T() {
  cbgeo::Clock timer;
  timer.tick();
  SolverControl solver_control(
      std::max<std::size_t>(n_T_max_iteration, T_system_rhs.size() / 10),
      T_tol_residual * T_system_rhs.l2_norm());     // setting for solver
  SolverGMRES<> solver(solver_control);             // config solver
  PreconditionJacobi<> preconditioner;              // precond
  preconditioner.initialize(T_system_matrix, 1.0);  // initialize precond
  solver.solve(T_system_matrix, T_solution, T_system_rhs,
               preconditioner);  // solve eq

  Vector<double> residual(dof_handler.n_dofs());
  T_system_matrix.vmult(residual, T_solution);
  residual -= T_system_rhs;
  T_iteration_namber = solver_control.last_step();

  std::cout << "  Iterations required for convergence:    "
            << T_iteration_namber << "\n"
            << "  Max norm of residual:      " << residual.linfty_norm()
            << "\n";

  // constraints.distribute(solution);  // make sure if the value is
  // consistent at the constraint point

  timer.tock("linear_solve_T");
}

// @sect4{<code>CoupledTH::output_results</code>}
//
// Neither is there anything new in generating graphical output:
template <int dim>
void CoupledTH<dim>::output_results(Vector<double>& solution,
                                    std::string var_name) const {
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, var_name);

  data_out.build_patches();

  const std::string filename = "outputfiles/" + var_name + "-solution-" +
                               Utilities::int_to_string(timestep_number, 3) +
                               ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
}

template <int dim>
void CoupledTH<dim>::run() {

  unsigned int binary_search_number;
  double initial_time_step;
  double theta;

  make_grid_and_dofs();

  VectorTools::interpolate(dof_handler,
                           EquationData::TemperatureInitialValues<dim>(),
                           old_T_solution);
  VectorTools::interpolate(
      dof_handler, EquationData::PressureInitialValues<dim>(), old_P_solution);

  output_results(old_T_solution, "T");
  output_results(old_P_solution, "P");

  do {
    std::cout << "\nTimestep " << timestep_number;

    binary_search_number = 1;
    initial_time_step =
        time_sequence[timestep_number + 1] - time_sequence[timestep_number];
    time_step = initial_time_step / 2;
    theta = 0;

    do {

      assemble_P_system();

      linear_solve_P();

      old_P_solution = P_solution;

      assemble_T_system();

      linear_solve_T();

      old_T_solution = T_solution;

      time += time_step;

      theta += pow(0.5, binary_search_number);

      if (P_iteration_namber > n_P_max_iteration / 2 ||
          T_iteration_namber > n_T_max_iteration / 2) {
        time_step = time_step / 2;
        ++binary_search_number;
      }

      std::cout << "\nt=" << time << ", dt=" << time_step << '.' << std::endl;

    } while ((1 - theta) > 0.00001);

    timestep_number += 1;
    output_results(T_solution, "T");
    output_results(P_solution, "P");

    std::cout << "\n" << std::endl << std::endl;

    // MatrixOut matrix_out;
    // std::ofstream out ("2rhs_matrix_at_"+std::to_string(time));
    // // matrix_out.build_patches (system_matrix, "system_matrix");
    // // matrix_out.write_gnuplot (out);
    // // system_matrix.print_formatted(out);
    // system_rhs.print(out);

  } while (time < period);
}
