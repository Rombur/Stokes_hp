#include "ExactSolutionEx3.hh"
#include <cmath>

using std::cos;
using std::sin;
using std::pow;
using std::exp;

template <>
void ExactSolutionEx3<2>::vector_value(const Point<2> &p, Vector<double> &values) const
{
  const double x(p[0]);
  const double y(p[1]);

  values(0) = 2.*y*cos(pow(x,2)+pow(y,2));
  values(1) = -2.*x*cos(pow(x,2)+pow(y,2));
  values(2) = exp(-10.*(pow(x,2)+pow(y,2)))- 0.3142;
// 0.3142 : is the pressure mean value
//(1/10)* ( pow( std::erf(std::sqrt(10)*std::sqrt(log(exp(1)))),2) * pi / log(exp(1)) );	
}


template <>
void ExactSolutionEx3<2>::vector_gradient(const Point<2> &p,
    std::vector<Tensor<1,2>> &gradients) const
{	
  const double x(p[0]);
  const double y(p[1]);

  gradients[0][0]= -4.*x*y*sin(pow(x,2)+pow(y,2));		
  gradients[0][1]= 2.*cos(pow(x,2)+pow(y,2)) -4*pow(y,2)*
    sin(pow(x,2)+pow(y,2));
  gradients[1][0]= -2.*cos(pow(x,2)+pow(y,2)) +4*pow(x,2)*
    sin(pow(x,2)+pow(y,2));
  gradients[1][1]= 4*x*y*sin(pow(x,2)+pow(y,2));
  gradients[2][0]= (-20.)*x*exp(-10.*(pow(x,2)+pow(y,2)));
  gradients[2][1]= (-20.)*y*exp(-10.*(pow(x,2)+pow(y,2)));
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
