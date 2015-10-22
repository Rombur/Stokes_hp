#ifndef _EXACTSOLUTIONEX2_HH_
#define _EXACTSOLUTIONEX2_HH_

#include <vector>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

using namespace dealii;

template <int dim>
class ExactSolutionEx2 : public Function<dim>
{
public:
  ExactSolutionEx2 () : Function<dim>(dim+1) {}
  virtual void vector_value (const Point<dim> &p,
                             Vector<double>   &value) const;
  virtual void vector_gradient (const Point<dim> &p,
                                std::vector< Tensor< 1, dim > > &gradients) const;
};

#endif
