#include "ExactSolutionEx5.hh"
#include <cmath>

using std::exp;
using std::pow;
using std::sin;
using std::cos;

template <>
void ExactSolutionEx5<2>::vector_value(const Point<2> &p, Vector<double> &values) const
{
  (void) p;
  (void) values;
}


template <>
void ExactSolutionEx5<2>::vector_gradient(const Point<2> &p,
                                          std::vector<Tensor<1,2>> &gradients) const
{
  (void) p;
  (void) gradients;
}



template <>
void ExactSolutionEx5<3>::vector_value(const Point<3> &p, Vector<double> &values) const
{
  const double x(p[0]);
  const double y(p[1]);

  values(0) = 2.*y*cos(pow(x,2)+pow(y,2));
  values(1) = -2.*x*cos(pow(x,2)+pow(y,2));
  values(2) = 0.;
  values(3) = exp(-10.*(pow(x,2)+pow(y,2))) - 
    2.*((1./40.)*M_PI*std::pow(std::erf(std::sqrt(10.)),2));

}


template <>
void ExactSolutionEx5<3>::vector_gradient(const Point<3> &p,
                                          std::vector<Tensor<1,3>> &gradients) const
{
  const double x(p[0]);
  const double y(p[1]);

  gradients[0][0]= -4.*x*y*sin(pow(x,2)+pow(y,2));
  gradients[0][1]= 2.*cos(pow(x,2)+pow(y,2)) -4*pow(y,2)*
                   sin(pow(x,2)+pow(y,2));
  gradients[0][2]= 0.;
  gradients[1][0]= -2.*cos(pow(x,2)+pow(y,2)) +4*pow(x,2)*
                   sin(pow(x,2)+pow(y,2));
  gradients[1][1]= 4*x*y*sin(pow(x,2)+pow(y,2));
  gradients[1][2]= 0.;
  gradients[2][0]= 0.;
  gradients[2][1]= 0.;
  gradients[2][2]= 0.;
  gradients[3][0]= (-20.)*x*exp(-10.*(pow(x,2)+pow(y,2)));
  gradients[3][1]= (-20.)*y*exp(-10.*(pow(x,2)+pow(y,2)));
  gradients[3][2]= 0;
}
