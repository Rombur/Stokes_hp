#ifndef _PARAMETERS_HH_
#define _PARAMETERS_HH_

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

// Example 1 is the first example of "hp-Adaptive Discontinuous Galerkin Finite
// Element Methods for the Stokes Problem."
// Example 2 is the second example of "An adaptive Uzawa FEM for the Stokes
// Problem: Convergence without the Inf-Sup Condition"
// Example 3 is the first example of "An adaptive Uzawa FEM for the Stokes
// Problem: Convergence without the Inf-Sup Condition"
// Example 4 is the third example of "An adaptive Uzawa FEM for the Stokes
// Problem: Convergence without the Inf-Sup Condition"
// Example 5 is the same as example 3 but in 3D. This example is for debugging
// purpose
enum EXAMPLE {example_1,example_2,example_3,example_4,example_5};
enum QUADRATURE {gauss_legendre,gauss_lobatto};
enum REFINEMENT {h_refine,p_refine,hp_refine};

class Parameters
{
  public :
    Parameters(const std::string &input_filename);

    bool do_goal_oriented() const;

    bool get_verbose() const;

    EXAMPLE get_example() const;

    QUADRATURE get_quadrature() const;

    REFINEMENT get_refinement() const;

    unsigned int get_dim() const;

    unsigned int get_max_degree() const;

    unsigned int get_max_n_cycles() const;

    double get_theta() const;

    double get_tolerance() const;

    std::string get_dual_source() const;

    // The following functions are only for testing purposes
    void set_example(EXAMPLE ex);

  private :
    void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);

    ParameterHandler prm;

    bool goal_oriented;
    bool verbose;
    EXAMPLE example;
    QUADRATURE quadrature;
    REFINEMENT refinement;
    unsigned int dim;
    unsigned int max_degree;
    unsigned int max_n_cycles;
    double theta;
    double tolerance;
    std::string dual_source;
};


inline
bool Parameters::do_goal_oriented() const
{
  return goal_oriented;
}


inline
bool Parameters::get_verbose() const
{
  return verbose;
}


inline
EXAMPLE Parameters::get_example() const
{
  return example;
}


inline
QUADRATURE Parameters::get_quadrature() const
{
  return quadrature;
}


inline
REFINEMENT Parameters::get_refinement() const
{
  return refinement;
}


inline
unsigned int Parameters::get_dim() const
{
  return dim;
}


inline
unsigned int Parameters::get_max_degree() const
{
  return max_degree;
}


inline
unsigned int Parameters::get_max_n_cycles() const
{
  return max_n_cycles;
}


inline
double Parameters::get_theta() const
{
  return theta;
}


inline
double Parameters::get_tolerance() const
{
  return tolerance;
}


inline
std::string Parameters::get_dual_source() const
{
  return dual_source;
}


inline
void Parameters::set_example(EXAMPLE ex) 
{
  example = ex;
}

#endif
