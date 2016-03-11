#ifndef _STOKESPROBLEM_HH_
#define _STOKESPROBLEM_HH_

#include <functional>
#include <vector>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/synchronous_iterator.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nothing.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include "CopyData.hh"
#include "ExactSolutionEx1.hh"
#include "ExactSolutionEx2.hh"
#include "ExactSolutionEx3.hh"
#include "ExactSolutionEx4.hh"
#include "Parameters.hh"
#include "RightHandSideEx1.hh"
#include "RightHandSideEx2.hh"
#include "RightHandSideEx3.hh"
#include "RightHandSideEx4.hh"
#include "ScratchData.hh"

using namespace dealii;

enum PRIMALDUAL{primal,dual};

template <int dim>
class StokesProblem
{
public:
  typedef typename hp::DoFHandler<dim>::active_cell_iterator DoFHandler_active_cell_iterator;
  typedef typename hp::DoFHandler<dim>::cell_iterator DoFHandler_cell_iterator;
  typedef typename Triangulation<dim>::active_cell_iterator Triangulation_active_cell_iterator;
  typedef typename Triangulation<dim>::cell_iterator Triangulation_cell_iterator;

  StokesProblem (Parameters const &parameters);

  void generate_mesh();

  void setup_system(PRIMALDUAL primal_dual);

  void assemble_system(PRIMALDUAL primal_dual);

  void solve(PRIMALDUAL primal_dual);

  void compute_error();

  void compute_error_estimator();

  std::vector<std::pair<double, DoFHandler_active_cell_iterator>> 
    compute_goal_oriented_error_estimator();

  void mark_cells_goal_oriented(const unsigned int cycle, const double theta,
      std::vector<std::pair<double, DoFHandler_active_cell_iterator>> 
      const &go_error_estimator_square);

  void mark_cells(const unsigned int cycle, const double theta);

  void output_results (const unsigned int cycle);

  void refine_in_h_p();

  unsigned int n_active_cells();

  types::global_dof_index n_dofs();

  double error_l2_norm();

  double error_estimate_l2_norm();

  double pressure_mean_value();

  // The following functions are only for testing purposes
  Triangulation<dim>& get_triangulation();

private:
  void set_active_fe_indices (hp::DoFHandler<dim> &local_dof_handler,
      std::map<Triangulation_active_cell_iterator,DoFHandler_active_cell_iterator>
      &patch_to_global_tria_map);

  void patch_output (unsigned int patch_number, const unsigned int cycle,
                     hp::DoFHandler<dim> &local_dof_handler, BlockVector<double> &local_solu);

  void p_refinement(hp::DoFHandler<dim> &local_dof_handler,
      std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator>
      &patch_to_global_tria_map, unsigned int level_p_refine, BlockVector<double> &local_solution);

  void h_refinement(Triangulation<dim> &local_triangulation,
                    hp::DoFHandler<dim> &local_dof_handler,
                    unsigned int level_h_refine,
                    BlockVector<double> &local_solution);


  void patch_assemble_system(hp::DoFHandler<dim> const &local_dof_handler,
      ConstraintMatrix const &constraints_patch, BlockVector<double> const &local_solu,
      BlockSparseMatrix<double> &patch_system, BlockVector<double> &patch_rhs,
      PRIMALDUAL primal_dual);

  void patch_solve(hp::DoFHandler<dim> &local_dof_handler,
                   unsigned int patch_number, unsigned int cycle, BlockVector<double> &local_solu,
                   double &conv_est, double &workload_num);

  void patch_convergence_estimator(const unsigned int cycle,
                                   SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
                                   std::vector<unsigned int>::iterator>> const &synch_iterator,
                                   ScratchData &scratch_data, CopyData<dim> &copy_data);

  std::vector<DoFHandler_cell_iterator> get_cells_at_coarsest_common_level (
    const std::vector<DoFHandler_active_cell_iterator> &patch);


  bool sort_decreasing_order (const std::pair<double, DoFHandler_active_cell_iterator> &i,
                              const std::pair<double, DoFHandler_active_cell_iterator > &j);

  void copy_to_refinement_maps(CopyData<dim> const &copy_data);

  std::vector<DoFHandler_active_cell_iterator> get_patch_around_cell(
    const DoFHandler_active_cell_iterator &cell);

  void build_triangulation_from_patch (
    const std::vector<DoFHandler_active_cell_iterator>  &patch,
    Triangulation<dim> &local_triangulation, unsigned int &level_h_refine,
    unsigned int &level_p_refine,
    std::map<Triangulation_active_cell_iterator,
    DoFHandler_active_cell_iterator> &patch_to_global_tria_map);

  void compute_local_dual_residual(
      hp::DoFHandler<dim> &local_dof_handler,
      unsigned int level_p_refine,
      std::vector<unsigned int> const &block_component_patch,
      std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator>
        &patch_to_global_tria_map,
      hp::DoFHandler<dim> &local_dual_dof_handler,
      BlockVector<double> &local_dual_solu,
      BlockVector<double> &dual_residual_solu);

  double compute_local_go_error_estimator_square(
      hp::DoFHandler<dim> const &dual_local_dof_handler,
      BlockVector<double> const &local_dual_solution,
      BlockVector<double> const &dual_residual_solution);

  std::unique_ptr<Function<dim>> exact_solution;
  std::unique_ptr<Function<dim>> rhs_function;
  std::unique_ptr<FunctionParser<dim>> dual_source;

  hp::FECollection<dim> fe_collection;
  hp::FECollection<dim> dual_fe_collection;

  Triangulation<dim> triangulation;
  hp::DoFHandler<dim> dof_handler;

  hp::QCollection<dim> quadrature_collection;
  hp::QCollection<dim-1> face_quadrature_collection;

  hp::QCollection<dim> quadrature_collection_error;

  ConstraintMatrix constraints;
  BlockSparsityPattern sparsity_pattern;
  BlockSparseMatrix<double> system_matrix;

  BlockVector<double> solution;
  BlockVector<double> system_rhs;
  BlockVector<double> dual_solution;

  bool verbose;
  EXAMPLE example;
  REFINEMENT refinement;
  const unsigned int max_degree;

  Vector<float> marked_cells;
  Vector<double> h_Conv_Est;
  Vector<double> p_Conv_Est;
  Vector<double> hp_Conv_Est;
  Vector<double> convergence_est_per_cell;
  Vector<double> error_per_cell;
  Vector<double> est_per_cell;
  Vector<double> Vect_Pressure_Err;
  Vector<double> Vect_grad_Velocity_Err;
  Vector<double> Vect_Velocity_Err;
  std::vector<std::pair<double, DoFHandler_active_cell_iterator>> primal_indicator_cell;
  std::vector<DoFHandler_active_cell_iterator> candidate_cell_set;
  std::map<DoFHandler_active_cell_iterator, bool> p_ref_map;
};


template <int dim>
inline
unsigned int StokesProblem<dim>::n_active_cells()
{
  return triangulation.n_active_cells();
}


template <int dim>
inline
types::global_dof_index StokesProblem<dim>::n_dofs()
{
  return dof_handler.n_dofs();
}


template <int dim>
inline
double StokesProblem<dim>::error_l2_norm()
{
  return error_per_cell.l2_norm();
}


template <int dim>
inline
double StokesProblem<dim>::error_estimate_l2_norm()
{
  return est_per_cell.l2_norm();
}


template <int dim>
inline
Triangulation<dim>& StokesProblem<dim>::get_triangulation()
{
  return triangulation;
}

#endif
