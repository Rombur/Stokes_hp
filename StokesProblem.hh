#ifndef _STOKESPROBLEM_HH_
#define _STOKESPROBLEM_HH_

#include <vector>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/synchronous_iterator.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/precondition_block.h>

#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

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
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/base/std_cxx1x/bind.h>

#include "ExactSolutionEx1.hh"
#include "ExactSolutionEx2.hh"
#include "ExactSolutionEx3.hh"
#include "ExactSolutionEx4.hh"

#include "CopyData.hh"
#include "Parameters.hh"
#include "RightHandSide.hh"
#include "ScratchData.hh"


using namespace dealii;

template <int dim>
class StokesProblem
{
  public:
    typedef typename hp::DoFHandler<dim>::active_cell_iterator DoFHandler_active_cell_iterator;
    typedef typename hp::DoFHandler<dim>::cell_iterator DoFHandler_cell_iterator;
    typedef typename Triangulation<dim>::active_cell_iterator Triangulation_active_cell_iterator;
    typedef typename Triangulation<dim>::cell_iterator Triangulation_cell_iterator;

    StokesProblem (bool verbose, EXAMPLE example, QUADRATURE quadrature, 
        unsigned int max_degree, unsigned int n_max_cycles, double theta,
        double tolerance);
    
    ~StokesProblem();
    
    void run();

    // The following functions are put as public to be able to test it.
    std::vector<DoFHandler_active_cell_iterator> get_patch_around_cell(
        const DoFHandler_active_cell_iterator &cell);

    void build_triangulation_from_patch (
        const std::vector<DoFHandler_active_cell_iterator>  &patch, 
        Triangulation<dim> &local_triangulation, unsigned int &level_h_refine, 
        unsigned int &level_p_refine, 
        std::map<Triangulation_active_cell_iterator, 
        DoFHandler_active_cell_iterator> &patch_to_global_tria_map);

  private:
    void generate_mesh ();

    void set_global_active_fe_indices (hp::DoFHandler<dim> &dof_handler);
    
    void setup_system ();
    
    void assemble_system ();
    
    void solve();
    
    double pressure_mean_value () const;
    
    double exact_pressure_mean_value () const;

    void compute_error (Vector<double> &error_per_cell, 
        Vector<double> &Vect_Pressure_Err, Vector<double> &Vect_grad_Velocity_Err, 
        Vector<double> &Vect_Velocity_Err);

    void estimate (Vector<double> &est_per_cell);

    void set_active_fe_indices (hp::DoFHandler<dim> &local_dof_handler, 
        std::map<Triangulation_active_cell_iterator, 
        DoFHandler_active_cell_iterator> &patch_to_global_tria_map);

    void patch_output (unsigned int patch_number, const unsigned int cycle, 
        hp::DoFHandler<dim> &local_dof_handler, BlockVector<double> &local_solu);

    void p_refinement(hp::DoFHandler<dim> &local_dof_handler, 
        std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator>
        &patch_to_global_tria_map, unsigned int level_p_refine, BlockVector<double> &local_solution);

    void h_refinement(Triangulation<dim> &local_triangulation,
        hp::DoFHandler<dim> &local_dof_handler, unsigned int level_h_refine, 
        BlockVector<double> &local_solution);

    void patch_assemble_system(hp::DoFHandler<dim> const &local_dof_handler,
        ConstraintMatrix const &constraints_patch, BlockVector<double> const &local_solu, 
        BlockSparseMatrix<double> &patch_system, BlockVector<double> &patch_solution, 
        BlockVector<double> &patch_rhs);

    void patch_solve(hp::DoFHandler<dim> &local_dof_handler, 
        unsigned int patch_number, unsigned int cycle, BlockVector<double> &local_solu,
        double &conv_est, double &workload_num);

    void patch_conv_load_no(const unsigned int cycle,
        SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
        std::vector<unsigned int>::iterator>> const &synch_iterator,
        ScratchData &scratch_data, CopyData<dim> &copy_data);

    std::vector<DoFHandler_cell_iterator> get_cells_at_coarsest_common_level (
        const std::vector<DoFHandler_active_cell_iterator> &patch);

    bool decreasing (const std::pair<double, DoFHandler_active_cell_iterator> &i, 
        const std::pair<double, DoFHandler_active_cell_iterator > &j);

    void copy_to_refinement_maps(Vector<double> const &est_per_cell, 
        CopyData<dim> const &copy_data);

    void marking_cells (const unsigned int cycle,  Vector<float> & marked_cells, 
        std::vector<DoFHandler_active_cell_iterator> &candidate_cell_set);

    void output_results (const unsigned int cycle, Vector<float> & marked_cells, 
        Vector<double> &est_per_cell, Vector<double> &error_per_cell, 
        Vector<double> &Vect_Pressure_Err, Vector<double> &Vect_grad_Velocity_Err);

    void refine_in_h_p (std::vector<DoFHandler_active_cell_iterator> &candidate_cell_set);

    Function<dim>* exact_solution;

    Triangulation<dim> triangulation;
    hp::DoFHandler<dim> dof_handler;

    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim> quadrature_collection;
    hp::QCollection<dim-1> face_quadrature_collection;

    hp::QCollection<dim> quadrature_collection_Err;
    hp::QCollection<dim-1> face_quadrature_collection_Err;

    ConstraintMatrix constraints;
    BlockSparsityPattern sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;

    bool verbose;
    EXAMPLE example;
    const unsigned int max_degree;
    const unsigned int max_n_cycles;
    double theta;
    const double tolerance;

    Vector<double> h_Conv_Est;
    Vector<double> p_Conv_Est;
    Vector<double> hp_Conv_Est;
    Vector<double> convergence_est_per_cell;
    std::vector<std::pair<double, DoFHandler_active_cell_iterator> > to_be_sorted;
    std::map<DoFHandler_active_cell_iterator, bool > p_ref_map;
};

#endif
