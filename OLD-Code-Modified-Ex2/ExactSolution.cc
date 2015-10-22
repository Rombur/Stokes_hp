#include "ExactSolution.hh"
#include <cmath>

#include <deal.II/lac/vector.h>

template <int dim>
void ExactSolution <dim>::vector_value (const Point<dim> &p,
                                        Vector<double>   &values) const
{
  /*
  //square domain in [-1,1]/ zero boundary for the first component of velocity is expected
  values(0)=(std::pow(p(0),2)-1)*(std::pow(p(1),2)-1);
  values(1)=(-2*p(0))*(1/3* std::pow(p(1),3)-p(1));
  values(2)=0.;
  */


  /*
  //L-shape domain zero velocity boundary for the first component

   values(0)=(p(0)-std::pow(p(0),3))*(p(1)-std::pow(p(1),3));
   values(1)=(3*std::pow(p(0),2)-1)*(1/2*std::pow(p(1),2)-1/4*std::pow(p(1),4));
   values(2)=0.; //p(0)+p(1);
  */



//L-shape/ actual problem / solution from Houston, et.al's paper

  values(0) = -1.0 * std::exp (p(0)) * (p(1) * std::cos (p(1)) + std::sin (p(1)));
  values(1) =    std::exp(p(0)) * p(1) * std::sin (p(1));
  values(2) = 2.0 * std::exp (p(0)) * std::sin (p(1)) - (2.0 * (1.0 - std::exp(1.0)) * (std::cos (1.0) - 1.0)) / 3.0;

}


template <int dim>
void ExactSolution <dim>::vector_gradient  (const Point<dim> &p,
                                            std::vector< Tensor< 1, dim > > &gradients) const
{
  /*
  gradients[0][0]=2*p(0)*(std::pow(p(1),2)-1);
  gradients[0][1]= 2*p(1) * (std::pow(p(0),2)-1);
  gradients[1][0]= -2*(1/3 * std::pow(p(1), 3)-p(1));
  gradients[1][1]= (-2*p(0))*(std::pow(p(1),2)-1);
  gradients[2][0]=0.;
  gradients[2][1]=0.;
  */

  /*
  //L-shape domain zero velocity boundary for the first component
  gradients[0][0]= (p(1)-std::pow(p(1),3))* (1-3*std::pow(p(0),2));
  gradients[0][1]= (p(0)-std::pow(p(0),3))* (1-3*std::pow(p(1),2));
  gradients[1][0]=  (1/2*std::pow(p(1),2)-1/4*std::pow(p(1),4))*(6*p(0));
  gradients[1][1]=  (p(1)-std::pow(p(1),3))*(3*std::pow(p(0),2)-1);
  gradients[2][0]=  0.;//1.;
  gradients[2][1]=  0.;//1.;
  */


//L-shape/ actual problem / solution from Houston, et.al's paper

  gradients[0][0]= -1.0 * std::exp (p(0)) * (p(1) * std::cos (p(1)) + std::sin (p(1)));
  gradients[0][1]= -1.0 * std::exp (p(0)) * ( 2.0 * std::cos (p(1)) - p(1) * std::sin (p(1)));
  gradients[1][0]=  std::exp(p(0)) * p(1) * std::sin (p(1));
  gradients[1][1]=  std::exp(p(0)) * (std::sin (p(1)) + p(1) * std::cos(p(1)));
  gradients[2][0]=  2.0 * std::exp (p(0)) * std::sin (p(1));
  gradients[2][1]=  2.0 * std::exp (p(0)) * std::cos (p(1));

}

//Explicit initialization

template class ExactSolution<2>;
