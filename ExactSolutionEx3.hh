#ifndef _EXACTSOLUTIONEX3_HH_
#define _EXACTSOLUTIONEX3_HH_

#include <vector>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h> 

using namespace dealii;

template <int dim>
class ExactSolutionEx3 : public Function<dim>
{
public:
	ExactSolutionEx3 () : Function<dim>(dim+1) {};

	void vector_value (const Point<dim> &p, Vector<double> &value) const;
	
  void vector_gradient (const Point<dim> &p, std::vector<Tensor<1,dim>> &gradients) const;
};

#endif
