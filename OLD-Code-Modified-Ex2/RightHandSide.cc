#include "RightHandSide.hh"
#include <cmath>

template <int dim>
double RightHandSide <dim>::value (const Point<dim> &p, const unsigned int  component) const
{
  if (component==0)
    // example-1 and example-2 , paper Houston et.al
    return 0.;
  /*
    // example 3.1, paper Morin-Nocheto- Uzawa
    return 16*p(1)*std::sin(std::pow(p(0),2)+std::pow(p(1),2))+ 8*std::pow(p(0),2)*p(1)*std::cos( std::pow(p(0),2)+std::pow(p(1),2))+8*std::pow(p(1),3)*std::cos( std::pow(p(0),2)
    +std::pow(p(1),2))-20*p(0)*std::exp(-10*( std::pow(p(0),2)+std::pow(p(1),2)));
  */

  //for assembling system as step-22 with 1/2 as coefficient in front
  //return std::exp (p(0)) * std::sin (p(1));

  //Manufacturred solution, L-shape domain without 1/2 in front of laplacian.
  //return 12*p(0)*p(1)-6*std::pow(p(0),3)*p(1)-6*std::pow(p(1),3)*p(0);

  //for square D, [-1,1]
  //return 2-std::pow(p(0),2)-std::pow(p(1),2);

  if (component==1)
    // using gradi_phi_u and rhs==0 without 1/2 in front of laplacian. ( // example-1 and example-2 , paper Houston et.al )
    return 0.;
  /*
    // example 3.1, paper Morin-Nocheto- Uzawa
    return  -16*p(0)*std::sin(  std::pow(p(0),2)+std::pow(p(1),2)) - 8*p(0)*std::pow(p(1),2)*std::cos(  std::pow(p(0),2)+std::pow(p(1),2))-8*std::pow(p(0),3)*std::cos( std::pow(p(0),2)+std::pow(p(1),2))
    -20*p(1)*std::exp(-10*( std::pow(p(0),2)+std::pow(p(1),2)));
  */


  //for assembling system as step-22 with 1/2 as coefficient in front
  //return std::exp (p(0)) * std::cos (p(1));


  //Manufacturred solution, L-shape domain without 1/2 in front of laplacian.
  //return 1+9*std::pow(p(0),2)*std::pow(p(1),2)+3/2*std::pow(p(1),4)-6*std::pow(p(1),2)-3*std::pow(p(0),2);


  //for square D, [-1,1]
  //return 2*p(0);

}


template <int dim>
void
RightHandSide <dim>::vector_value (const Point<dim> &p, Vector<double> &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = RightHandSide::value (p, c);
}

//Explicit initialization

template class RightHandSide<2>;
