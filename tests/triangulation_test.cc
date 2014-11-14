#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <set>

#include "../StokesProblem.hh"
#include "deal.II/grid/grid_generator.h"

using namespace dealii;

TEST_CASE("Build triangulation patch","[patch]")
{
  StokesProblem<2> stokes_problem;

  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);
  hp::DoFHandler<2> dof_handler(triangulation);

  typename hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
  for (unsigned int i=0; i<6; ++i)
    ++cell;

  Triangulation<2> tria_patch;
  std::vector<typename hp::DoFHandler<2>::active_cell_iterator> patch = 
    stokes_problem.get_patch_around_cell(cell);		

  std::set<int> coord;
  coord.insert(-2);
  coord.insert(0);
  coord.insert(2);

  for (unsigned int k=0; k<patch.size(); ++k)
  {
    int i((patch[k]->center()(0)-cell->center()(0))/0.125);
    int j((patch[k]->center()(1)-cell->center()(1))/0.125);

    REQUIRE(coord.count(i));
    REQUIRE(coord.count(j));
  }


  unsigned int level_h_refine;
  unsigned int level_p_refine;
  stokes_problem.build_triangulation_from_patch (patch, tria_patch, level_h_refine, level_p_refine);

  hp::DoFHandler<2> dof_handler_patch(tria_patch);

  typename hp::DoFHandler<2>::active_cell_iterator patch_cell;
  for (patch_cell = dof_handler_patch.begin(); patch_cell!=dof_handler_patch.end(); ++patch_cell)
  {
    int i((patch_cell->center()(0)-cell->center()(0))/0.125);
    int j((patch_cell->center()(1)-cell->center()(1))/0.125);

    REQUIRE(coord.count(i));
    REQUIRE(coord.count(j));
  }
}
