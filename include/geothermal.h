#pragma once
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_sparsity_pattern.h>
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

#include "globalvariables.h"

using namespace dealii;

template <int dim>
class CoupledTH {
 public:
  CoupledTH();
  void run();

 private:
  void make_grid_and_dofs();
  void assemble_T_system();
  void assemble_p_system();
  void solve_p_time_step();
  void solve_T_time_step();
  void output_results() const;

  Triangulation<dim> triangulation;  // grid
  const unsigned int T_degree;
  const unsigned int p_degree;

  FE_Q<dim> p_fe;                 // element
  FE_Q<dim> T_fe;                 // element
  DoFHandler<dim> p_dof_handler;  // grid<->eleemnt
  DoFHandler<dim> T_dof_handler;  // grid<->eleemnt

  // ConstraintMatrix constraints;  // hanging node

  SparsityPattern sparsity_pattern;         // sparsity
  SparseMatrix<double> p_mass_matrix;       // M
  SparseMatrix<double> T_mass_matrix;       // M
  SparseMatrix<double> p_stiffness_matrix;  // A
  SparseMatrix<double> T_stiffness_matrix;  // A
  SparseMatrix<double> p_system_matrix;            // M + k*theta*A
  SparseMatrix<double> T_system_matrix;            // M + k*theta*A

  Vector<double> p_solution;      // solution at n
  Vector<double> T_solution;      // solution at n
  Vector<double> old_p_solution;  // solution at n-1
  Vector<double> old_T_solution;  // solution at n-1
  Vector<double> p_system_rhs;    // rhs
  Vector<double> T_system_rhs;    // rhs

  double time;
  double time_step;
  unsigned int timestep_number;

  const double theta;
};

template <int dim>
CoupledTH<dim>::CoupledTH()  // initialization
    : p_fe(1),
      T_fe(1),
      T_degree(T_degree),
      p_degree(p_degree),
      p_dof_handler(triangulation),
      T_dof_handler(triangulation),

      time(0.0),
      time_step(1. / 20),  // a time step constant at 1/500 (remember that one
                           // period of the source on the right hand side was
                           // set to 0.2 above, so we resolve each period with
                           // 100 time steps)
      timestep_number(0),
      theta(0.5) {}

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
  GridIn<dim> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f("mesh.msh");
  gridin.read_msh(f);

  print_mesh_info(triangulation, "grid-1.eps");

  p_dof_handler.distribute_dofs(p_fe);  // distribute dofs to grid
  T_dof_handler.distribute_dofs(T_fe);  // distribute dofs to grid

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << " (on " << triangulation.n_levels() << " levels)" << std::endl
            << "Number of degrees of freedom: "
            << p_dof_handler.n_dofs() + T_dof_handler.n_dofs() << " (" < < < <
      p_dof_handler.n_dofs()
          << '+' << T_dof_handler.n_dofs() << ')' << std::endl
          << std::endl;
  DynamicSparsityPattern dsp(T_dof_handler.n_dofs());  // sparsity
  DynamicSparsityPattern dsp(p_dof_handler.n_dofs());  // sparsity
  DoFTools::make_sparsity_pattern(
      T_dof_handler, dsp /*constraints, /*keep_constrained_dofs = ,true*/);
  DoFTools::make_sparsity_pattern(
      p_dof_handler, dsp /*constraints, /*keep_constrained_dofs = ,true*/);
  sparsity_pattern.copy_from(dsp);

// initialize MATRIX using given sparsity
  T_system_matrix.reinit(sparsity_pattern);  
  p_system_matrix.reinit(sparsity_pattern); 

  T_mass_matrix.reinit(sparsity_pattern);
  p_mass_matrix.reinit(sparsity_pattern);

  T_stiffness_matrix.reinit(sparsity_pattern);
  p_stiffness_matrix.reinit(sparsity_pattern);

  T_solution.reinit(T_dof_handler.n_dofs());
  p_solution.reinit(p_dof_handler.n_dofs());

  old_T_solution.reinit(T_dof_handler.n_dofs());
  old_p_solution.reinit(p_dof_handler.n_dofs());

  T_system_rhs.reinit(T_dof_handler.n_dofs());
  p_system_rhs.reinit(p_dof_handler.n_dofs());

  // constraints.clear();
  // DoFTools::make_hanging_node_constraints(
  //     dof_handler,
  //     constraints); // setting the hanging node
  // constraints.close();
}

template <int dim>
void CoupledTH<dim>::assemble_p_system() {

  // reset matreix to zero
  // p_mass_matrix = 0;
  // p_stiffness_matrix = 0;

  QGauss<dim> T_quadrature_formula(T_degree + 1);
  QGauss<dim> p_quadrature_formula(p_degree + 1);

  QGauss<dim - 1> T_face_quadrature_formula(T_degree + 1);
  QGauss<dim - 1> p_face_quadrature_formula(p_degree + 1);

  FEValues<dim> T_fe_values(T_fe, T_quadrature_formula,
      update_values | update_gradients | update_JxW_values);
  FEFaceValues<dim> T_fe_face_values(T_fe, T_face_quadrature_formula,
      update_values | update_normal_vectors |
      update_quadrature_points | update_JxW_values);

  FEValues<dim> p_fe_values(
      p_fe, p_quadrature_formula,
      update_values | update_gradients | update_JxW_values);
  FEFaceValues<dim> p_fe_face_values(p_fe, p_face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                         update_quadrature_points |
                                         update_JxW_values);

  const unsigned int T_dofs_per_cell = T_fe.dofs_per_cell;
  const unsigned int p_dofs_per_cell = p_fe.dofs_per_cell;

  const unsigned int T_n_q_points = T_quadrature_formula.size();
  const unsigned int p_n_q_points = p_quadrature_formula.size();

  const unsigned int T_n_face_q_points = T_face_quadrature_formula.size();
  const unsigned int p_n_face_q_points = p_face_quadrature_formula.size();

  FullMatrix<double> p_local_mass_matrix(p_dofs_per_cell, p_dofs_per_cell);
  FullMatrix<double> p_local_stiffness_matrix(p_dofs_per_cell, p_dofs_per_cell);
  Vector<double> p_local_rhs(p_dofs_per_cell);
  std::vector<types::global_dof_index> p_local_dof_indices(p_dofs_per_cell);

  // boudnary condition
  const EquationData::PressureRightHandSide<dim> p_right_hand_side;
  const EquationData::PressureBoundaryValues<dim> p_boundary;

  // store the rhs and bd value at q_point
  std::vector<double> p_rhs_values(p_n_q_points);
  std::vector<double> p_bd_values(p_n_face_q_points);

  // store the value at previous step at q_point
  std::vector<double> old_T_sol_values(T_n_q_points);
  std::vector<double> old_p_sol_values(p_n_q_points);
  std::vector<Tensor<1, dim>> old_T_sol_grads(T_n_q_points);
  std::vector<Tensor<1, dim>> old_p_sol_grads(p_n_q_points);

  // loop for cell
  cell = p_dof_handler.begin_active(), endc = p_dof_handler.end();
  for (; cell != endc; ++cell) {
    // initialization
    p_local_mass_matrix = 0;
    p_local_stiffness_matrix = 0;
    p_local_rhs = 0;
    p_fe_values.reinit(cell);

    // get teh values at gauss point
    T_fe_values.get_function_values(old_T_solution, old_T_sol_values);
    p_fe_values.get_function_values(old_p_solution, old_p_sol_values);
    p_right_hand_side.value_list(p_fe_values.get_quadrature_points(),
                                 p_rhs_values);

    // loop for q_point ASSMBLING CELL METRIX
    for (unsigned int q = 0; q < p_n_q_points; ++q) {
      for (unsigned int i = 0; i < p_dofs_per_cell; ++i) {
        const double grad_phi_i_p = p_fe_values.shape_grad(i, q);
        const double phi_i_p = p_fe_values.shape_value(i, q);
        for (unsigned int j = 0; j < p_dofs_per_cell; ++j) {
          const double grad_phi_j_p = p_fe_values.shape_grad(j, q);
          const double phi_j_p = p_fe_values.shape_value(j, q);
          p_local_mass_matrix(i, j) += (phi_i_p * phi_j_p * p_fe_values.JxW(q));
          p_local_stiffness_matrix(i, j) +=
              (EquationData::kappa * grad_phi_i_p * grad_phi_j_p *
               p_fe_values.JxW(q));
        }
        p_local_rhs(i) +=
            (-phi_i_p * pressure_rhs_values[q]) * fe_values.JxW(q);
      }
    }

    // APPLIED BOUNDARY CONDITION
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
         ++face_no)
      if (cell->at_boundary(face_no)) {
        p_fe_face_values.reinit(cell, face_no);
        p_boundary.value_list(p_fe_face_values.get_quadrature_points(),
                              p_bd_values);

        for (unsigned int q = 0; q < p_n_face_q_points; ++q) {
          for (unsigned int i = 0; i < p_dofs_per_cell; ++i) {
            const double phi_i_p = p_fe_face_values.value(i, q);
            p_local_rhs(i) += -(phi_i_p * p_fe_face_values.normal_vector(q) *
                                p_bd_values[q] * p_fe_face_values.JxW(q))
          }
        }
      }

    // local ->globe
    cell->get_dof_indices(p_local_dof_indices);

    for (unsigned int i = 0; i < p_dofs_per_cell; ++i)
      for (unsigned int j = 0; j < p_dofs_per_cell; ++j)
        p_mass_matrix.add(p_local_dof_indices[i], p_local_dof_indices[j],
                          p_local_mass_matrix(i, j));
    p_stiffness_matrix.add(p_local_dof_indices[i], p_local_dof_indices[j],
                           p_local_stiffness_matrix(i, j));

    for (unsigned int i = 0; i < p_dofs_per_cell; ++i)
      p_system_rhs(local_dof_indices[i]) += p_local_rhs(i);

    //
  }
}

template <int dim>
void CoupledTH<dim>::assemble_T_system() {

  T_mass_matrix = 0;
  T_stiffness_matrix = 0;

  QGauss<dim> T_quadrature_formula(T_degree + 1);
  QGauss<dim> p_quadrature_formula(p_degree + 1);

  QGauss<dim - 1> T_face_quadrature_formula(T_degree + 1);
  QGauss<dim - 1> p_face_quadrature_formula(p_degree + 1);

  FEValues<dim> T_fe_values(
      T_fe, T_quadrature_formula,
      update_values | update_gradients | update_JxW_values);
  FEFaceValues<dim> T_fe_face_values(T_fe, T_face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                         update_quadrature_points |
                                         update_JxW_values);
  FEValues<dim> p_fe_values(
      p_fe, p_quadrature_formula,
      update_values | update_gradients | update_JxW_values);
  FEFaceValues<dim> p_fe_face_values(p_fe, p_face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                         update_quadrature_points |
                                         update_JxW_values);

  const unsigned int T_dofs_per_cell = T_fe.dofs_per_cell;
  const unsigned int p_dofs_per_cell = p_fe.dofs_per_cell;

  const unsigned int T_n_q_points = T_quadrature_formula.size();
  const unsigned int p_n_q_points = p_quadrature_formula.size();

  const unsigned int T_n_face_q_points = T_face_quadrature_formula.size();
  const unsigned int p_n_face_q_points = p_face_quadrature_formula.size();

  FullMatrix<double> T_local_mass_matrix(T_dofs_per_cell, T_dofs_per_cell);
  FullMatrix<double> T_local_stiffness_matrix(T_dofs_per_cell, T_dofs_per_cell);
  Vector<double> T_local_rhs(T_dofs_per_cell);
  std::vector<types::global_dof_index> T_local_dof_indices(T_dofs_per_cell);

  // boudnary condition
  const EquationData::TemperatureRightHandSide<dim> T_right_hand_side;
  const EquationData::TemperatureBoundaryValues<dim> T_boundary;

  // store the rhs and bd value at q_point
  std::vector<double> T_rhs_values(T_n_q_points);
  std::vector<double> T_bd_values(T_n_face_q_points);

  // store the value at previous step at q_point
  std::vector<double> old_T_sol_values(T_n_q_points);
  std::vector<double> old_p_sol_values(p_n_q_points);
  std::vector<Tensor<1, dim>> old_T_sol_grads(T_n_q_points);
  std::vector<Tensor<1, dim>> old_p_sol_grads(p_n_q_points);

  // loop for cell
  cell = T_dof_handler.begin_active(), endc = T_dof_handler.end();
  for (; cell != endc; ++cell) {
    // initialization
    T_local_mass_matrix = 0;
    T_local_stiffness_matrix = 0;
    T_local_rhs = 0;
    T_fe_values.reinit(cell);

    // get teh values at gauss point
    T_fe_values.get_function_values(old_T_solution, old_T_sol_values);
    p_fe_values.get_function_values(old_p_solution, old_p_sol_values);
    T_right_hand_side.value_list(T_fe_values.get_quadrature_points(),
                                 T_rhs_values);

    // loop for q_point ASSMBLING CELL METRIX
    for (unsigned int q = 0; q < T_n_q_points; ++q) {
      for (unsigned int i = 0; i < T_dofs_per_cell; ++i) {
        const double grad_phi_i_T = T_fe_values.shape_grad(i, q);
        const double phi_i_T = T_fe_values.shape_value(i, q);
        for (unsigned int j = 0; j < T_dofs_per_cell; ++j) {
          const double grad_phi_j_T = T_fe_values.shape_grad(j, q);
          const double phi_j_T = T_fe_values.shape_value(j, q);
          T_local_mass_matrix(i, j) +=
              (phi_i_T * grad_phi_j_T * T_fe_values.JxW(q));
          T_local_stiffness_matrix(i, j) +=
              (EquationData::kappa * grad_phi_i_T * grad_phi_j_T *
               T_fe_values.JxW(q));
        }
        T_local_rhs(i) +=
            (-phi_T[i] * temperature_rhs_values[q]) * T_fe_values.JxW(q);
      }
    }

    // APPLIED BOUNDARY CONDITION
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
         ++face_no)
      if (cell->at_boundary(face_no)) {
        T_fe_face_values.reinit(cell, face_no);
        T_boundary.value_list(T_fe_face_values.get_quadrature_points(),
                              T_bd_values);

        for (unsigned int q = 0; q < T_n_face_q_points; ++q)
          for (unsigned int i = 0; i < T_dofs_per_cell; ++i) {
            const double phi_i_T = T_fe_face_values.value(i, q);
            local_rhs(i) += -(phi_i_T * T_fe_face_values.normal_vector(q) *
                              T_bd_values[q] * T_fe_face_values.JxW(q))
          }
      }

    // local ->globe
    cell->get_dof_indices(T_local_dof_indices);

    for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
      for (unsigned int j = 0; j < T_dofs_per_cell; ++j)
        T_mass_matrix.add(T_local_dof_indices[i], T_local_dof_indices[j],
                          T_local_mass_matrix(i, j));
    T_stiffness_matrix.add(T_local_dof_indices[i], T_local_dof_indices[j],
                           T_local_stiffness_matrix(i, j));

    for (unsigned int i = 0; i < T_dofs_per_cell; ++i)
      T_system_rhs(local_dof_indices[i]) += T_local_rhs(i);

    //
  }
}

template <int dim>
void CoupledTH<dim>::solve_p_time_step() {
  SolverControl solver_control(1000,
                               1e-8 * p_system_rhs.l2_norm());  // setting for cg
  SolverCG<> cg(solver_control);                              // config cg

  PreconditionSSOR<> preconditioner;              // precond
  preconditioner.initialize(p_system_matrix, 1.0);  // initialize precond

  cg.solve(p_system_matrix, p_solution, p_system_rhs,
           preconditioner);  // solve eq

  // constraints.distribute(solution);  // make sure if the value is
  // consistent at
  // the constraint point

  std::cout << "     " << solver_control.last_step() << " CG iterations."
            << std::endl;
}

// @sect4{<code>CoupledTH::output_results</code>}
//
// Neither is there anything new in generating graphical output:
template <int dim>
void CoupledTH<dim>::output_results() const {
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
void CoupledTH<dim>::run() {
  make_grid_and_dofs();
  assemble_pressure_system();

  Vector<double> tmp;            // this vector is for
  Vector<double> forcing_terms;  // this vector is for forcing_terms

  tmp.reinit(solution.size());            // initialize tmp
  forcing_terms.reinit(solution.size());  // initialize forcing_terms

  VectorTools::interpolate(
      dof_handler, ZeroFunction<dim>(),
      old_solution);  // interpolate the old solution based on dof_handler,
                      // here using interpolation because we refine the
                      // global
  solution =
      old_solution;  // updating the solutin with sinterpolated old solution

  output_results();  // output

  while (time <= 1.) {
    time += time_step;  // get new time point
    ++timestep_number;  // get the No. of the new time step

    std::cout << "Time step " << timestep_number << " at t=" << time
              << std::endl;

    mass_matrix.vmult(system_rhs,
                      old_solution);  // matrix multiplication system_rhs =
                                      // mass_matrix*old_solution

    laplace_matrix.vmult(tmp,
                         old_solution);  // tmp = laplace_matrix*old_solution
    system_rhs.add(
        -(1 - theta) * time_step,
        tmp);  // system_rhs = system_rhs -(1 - theta) * time_step*tmp
               // 注意，这里system_rhs是一个vector，所以add是将两个元素相乘了

    EquationData::TemperatureRightHandSide<dim> rhs_function;
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
        tmp);  // 形成forcing term = f(x,n-1)*(1-theta) + f(x,n)*theta

    system_rhs +=
        forcing_terms;  // system_rhs = system_rhs + forcing_term :
                        // sys_Old*U_Old + f(x,n-1)*(1-theta) + f(x,n)*theta

    system_matrix.copy_from(mass_matrix);
    system_matrix.add(theta * time_step, laplace_matrix);  // sys = M +
                                                           // k*theta*A

    // constraints.condense(system_matrix, system_rhs);  // 压缩

    {
      EquationData::TemperatureBoundaryValues<dim>
          boundary_values_function;             // creat boundary value
                                                // object
      boundary_values_function.set_time(time);  // set the proper time

      std::map<types::global_dof_index, double> boundary_values;
      VectorTools::interpolate_boundary_values(
          dof_handler,  // evaluate value by interpolation
          1, boundary_values_function, boundary_values);

      MatrixTools::apply_boundary_values(boundary_values, system_matrix,
                                         solution, system_rhs);
    }

    solve_time_step();

    output_results();

    old_solution = solution;
  }
}
