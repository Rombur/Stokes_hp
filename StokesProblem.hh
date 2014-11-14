#ifndef _STOKESPROBLEM_HH_
#define _STOKESPROBLEM_HH_

#include <vector>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nothing.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/base/std_cxx1x/bind.h>

#include "ExactSolution.hh"
#include "RightHandSide.hh"

using namespace dealii;

template <int dim>
class StokesProblem
{
  public:
    StokesProblem ();
    ~StokesProblem();
    void run();

    std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> get_patch_around_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell);

    void build_triangulation_from_patch (const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch, Triangulation<dim> &tria_patch, unsigned int &level_h_refine, unsigned int &level_p_refine);

  private:

    const RightHandSide<dim> rhs_function;
    const ExactSolution<dim> exact_solution;

    void generate_mesh ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    double pressure_mean_value () const;
    void compute_error ();
    void estimate (Vector<double> &est_per_cell);
    void set_active_fe_indices (hp::DoFHandler<dim> &dof_handler_patch);

    void h_patch_conv_load_no (double &h_convergence_est_per_cell, unsigned int &h_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell);
    void p_patch_conv_load_no (double &p_convergence_est_per_cell, unsigned int &p_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell );

    std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> get_cells_at_coarsest_common_level ( const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch);

    bool decreasing (const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator> &i, const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator > &j);
   
    
    void marking_cells (const unsigned int cycle,  Vector<float> & marked_cells, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set, 
    std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map);
    void output_results (const unsigned int cycle,Vector<float> & marked_cells);
    void refine_in_h_p (const unsigned int cycle, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set, 
    std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map );

    Triangulation<dim> triangulation;
    hp::DoFHandler<dim> dof_handler;

    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim> quadrature_collection;
    hp::QCollection<dim-1> face_quadrature_collection;

    ConstraintMatrix constraints;
    BlockSparsityPattern sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;

    const unsigned int max_degree;
    const double Tolerance;
};

#endif
