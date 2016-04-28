#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "Parameters.hh"

TEST_CASE("parameters", "[parameters]")
{
  Parameters param("test_project.inp");
  REQUIRE(param.do_goal_oriented() == 0);
  REQUIRE(param.get_verbose() == 0);
  REQUIRE(param.get_example() == 0);
  REQUIRE(param.get_quadrature() == 1);
  REQUIRE(param.get_refinement() == 0);
  REQUIRE(param.get_dim() == 2);
  REQUIRE(param.get_max_degree() == 10);
  REQUIRE(param.get_max_n_cycles() == 15);
  REQUIRE(param.get_theta() == 1.0);
  REQUIRE(param.get_tolerance() == 1e-15);
  REQUIRE(param.get_dual_source().compare("if(x>0.5,y>0.5,10.);if(x>0.5,y>0.5,10.);0.") == 0);
}
