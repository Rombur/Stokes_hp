#include "ExactSolutionEx3.hh"
#include <cmath>


template <>
void ExactSolutionEx3<2>::vector_value(const Point<2> &p, Vector<double> &values) const
{
  values(0) = 2*p(1)*std::cos(std::pow(p(0),2)+std::pow(p(1),2));
  values(1) =  -2*p(0)*std::cos(std::pow(p(0),2)+std::pow(p(1),2));
  values(2) = std::exp(-10*(std::pow(p(0),2)+std::pow(p(1),2)));	
}


template <>
void ExactSolutionEx3<2>::vector_gradient(const Point<2> &p,
    std::vector<Tensor<1,2>> &gradients) const
{	
  gradients[0][0]= -4*p(0)*p(1)*std::sin(std::pow(p(0),2)+std::pow(p(1),2));		
  gradients[0][1]= 2*std::cos(std::pow(p(0),2)+std::pow(p(1),2))-4*std::pow(p(1),2)*
    std::sin(std::pow(p(0),2)+std::pow(p(1),2));
  gradients[1][0]= -2*std::cos(std::pow(p(0),2)+std::pow(p(1),2) ) + 
    4*std::pow(p(0),2)*std::sin(std::pow(p(0),2)+std::pow(p(1),2));
  gradients[1][1]=  4*p(0)*p(1)*std::sin(std::pow(p(0),2)+std::pow(p(1),2)) ;
  gradients[2][0]= std::exp(-10*(std::pow(p(0),2)+std::pow(p(1),2))) * -20*p(0);
  gradients[2][1]= std::exp(-10*(std::pow(p(0),2)+std::pow(p(1),2))) * -20*p(1);
}



template <>
void ExactSolutionEx3<3>::vector_value(const Point<3> &p, Vector<double> &values) const
{
  // Silence warnings
  (void) p;
  (void) values;
}


template <>
void ExactSolutionEx3<3>::vector_gradient(const Point<3> &p,
    std::vector<Tensor<1,3>> &gradients) const
{
  // Silence warnings
  (void) p;
  (void) gradients;
}
