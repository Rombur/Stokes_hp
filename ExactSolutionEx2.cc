#include "ExactSolutionEx2.hh"
#include <cmath>

#include <deal.II/lac/vector.h> 


template <int dim>
void ExactSolutionEx2 <dim>::vector_value (const Point<dim> &p,
		Vector<double>   &values) const
		{
	const double PI = std::atan(1.0)*4;
	const double Lambda= 0.54448373678246;
	const double Omega=3/2*(PI);
	double r= std::sqrt(std::pow(p(0),2)+std::pow(p(1),2));

	double Psi= (std::sin((1+Lambda)*(std::atan(p(1)/p(0))))*std::cos(Lambda*Omega))/(1+Lambda)-std::cos((1+Lambda)*std::atan(p(1)/p(0)))
	-(std::sin((1-Lambda)*std::atan(p(1)/p(0))) * std::cos(Lambda*Omega))/(1-Lambda) +std::cos((1-Lambda)*std::atan(p(1)/p(0)));

	double Deriv_Psi= (std::cos(Lambda*Omega) * std::cos((1+Lambda)*std::atan(p(1)/p(0))))+ (1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0)))-
			std::cos(Lambda*Omega)*std::cos(((1-Lambda)*std::atan(p(1)/p(0))))-(1-Lambda)*std::sin((1-Lambda)*std::atan(p(1)/p(0)));

	double third_Deriv_Psi= - (std::pow((1+Lambda),2)) * std::cos(Lambda*Omega)*std::cos((1+Lambda)*std::atan(p(1)/p(0)))- std::pow((1+Lambda),3)*std::sin((1+Lambda)*std::atan(p(1)/p(0)))+ std::pow((1-Lambda),2)*std::cos(Lambda*Omega)*std::cos((1-Lambda)*std::atan(p(1)/p(0)))+ std::pow((1-Lambda),3)*std::sin((1-Lambda)* std::atan(p(1)/p(0)));

	//L-shape/ actual problem / solution from Houston, et.al's paper

	values(0) = std::pow(r, Lambda) *((1+Lambda)*std::sin(std::atan(p(1)/p(0))) * Psi + std::cos(std::atan(p(1)/p(0))) *
			Deriv_Psi ) ;
	values(1) =  std::pow(r, Lambda) * (std::sin(std::atan(p(1)/p(0)))* Deriv_Psi - (1+Lambda) * std::cos(std::atan(p(1)/p(0))) * Psi);
	values(2) =  -std::pow(r, (Lambda-1)) * (std::pow((1+Lambda),2)*Deriv_Psi + third_Deriv_Psi )/ (1-Lambda);
		}

template <int dim>
void ExactSolutionEx2 <dim>::vector_gradient  (const Point<dim> &p,
		std::vector< Tensor< 1, dim > > &gradients) const

		{
const double PI = std::atan(1.0)*4;
	const double Lambda= 0.54448373678246;
	const double Omega=3/2*(PI);
	double r= std::sqrt(std::pow(p(0),2)+std::pow(p(1),2));

	//L-shape/ actual problem / solution from Houston, et.al's paper


	gradients[0][0]= 1/(std::pow(r,2)) *(std::pow(r, Lambda)*Lambda*p(0)*( (1+Lambda)*p(1)* ( (std::sin((1+Lambda)*atan(p(1)/p(0)))* std::cos(Lambda*Omega)) / (1+Lambda) - std::cos((1+Lambda)* std::atan(p(1)/p(0)) ) - ( (std::sin((1-Lambda)* std::atan(p(1)/p(0)) )*std::cos(Lambda*Omega)) / (1-Lambda) )+ ( std::cos((1-Lambda)*std::atan(p(1)/p(0)))) )  /  (p(0)* std::sqrt(1+ std::pow( p(1)/p(0),2)))+
			( std::cos(Lambda*Omega)* std::cos((1+Lambda)*std::atan(p(1)/p(0))) + (1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0)))-  std::cos(Lambda*Omega)*std::cos((1-Lambda)* std::atan(p(1)/p(0))) - (1-Lambda)*std::sin((1-Lambda)*std::atan(p(1)/p(0))) )  /(std::sqrt(1+ std::pow (p(1)/p(0),2))) ))+
			std::pow(r, Lambda)*( -( (1+Lambda)*p(1)* ( (std::sin((1+Lambda)*std::atan(p(1)/p(0)))* std::cos(Lambda*Omega) ) /(1+Lambda) -   std::cos((1+Lambda)*std::atan(p(1)/p(0))) -  ((std::sin((1-Lambda)*std::atan(p(1)/p(0))))*std::cos(Lambda*Omega))/(1-Lambda) + std::cos((1-Lambda)*std::atan(p(1)/p(0))) ) )/(std::pow(p(0),2) * std::sqrt (1+ std::pow(p(1)/p(0),2))) +
  ( (1+Lambda)*std::pow(p(1),3)* ((std::sin((1+Lambda)*std::atan(p(1)/p(0)))* std::cos(Lambda*Omega))  /(1+Lambda)   -  std::cos((1+Lambda)*std::atan(p(1)/p(0)))   -  (std::sin((1-Lambda)*std::atan(p(1)/p(0))) * std::cos(Lambda*Omega))/(1-Lambda)   +  std::cos((1-Lambda)*std::atan(p(1)/p(0))) ))/( std::pow(p(0),4)* std::pow( (1+std::pow(p(1)/p(0) ,2)) ,(3/2) )) +
( (1+Lambda)*p(1)* ( -( std::cos((1+Lambda)*std::atan(p(1)/p(0)))*p(1)* std::cos(Lambda*Omega) )/(std::pow(p(0),2)* (1+ std::pow((p(1)/p(0)) ,2))  )  -
 (std::sin((1+Lambda)*std::atan(p(1)/p(0)) )* (1+Lambda)*p(1) )/(std::pow(p(0),2)* (1+ std::pow((p(1)/p(0)) ,2))  ) +
(std::cos((1-Lambda)* std::atan(p(1)/p(0)))* p(1)*std::cos(Lambda*Omega) )/(std::pow(p(0),2)* (1+ std::pow((p(1)/p(0)) ,2))  ) +
(std::sin((1-Lambda)*std::atan(p(1)/p(0))) * (1-Lambda)* p(1))/(std::pow(p(0),2)* (1+ std::pow((p(1)/p(0)) ,2))  ) ) )/( p(0)*std::sqrt(1+std::pow(p(1)/p(0) , 2)) ) +
( std::pow(p(1),2)* ( std::cos(Lambda*Omega)* std::cos((1+Lambda)*std::atan(p(1)/p(0))) +  (1+Lambda)* std::sin((1+Lambda)*std::atan(p(1)/p(0))) - std::cos(Lambda*Omega)*std::cos((1-Lambda)*std::atan(p(1)/p(0))) - (1-Lambda)*std::sin((1-Lambda)*std::atan(p(1)/p(0))) ))/(std::pow(p(0),3)* std::pow( 1+ std::pow(p(1)/p(0) ,2) , (3/2)) ) +

( ( std::cos(Lambda*Omega)*(1+Lambda)*p(1)*std::sin((1+Lambda)* std::atan(p(1)/p(0))) )/(std::pow(p(0),2)*(1+std::pow(p(1)/p(0),2)) ) -
( std::pow((1+Lambda),2)* p(1)* std::cos((1+Lambda)* std::atan(p(1)/p(0))) )/(std::pow(p(0),2)*(1+std::pow(p(1)/p(0),2)) )-
((1-Lambda)*p(1)*std::cos(Lambda*Omega)*std::sin((1-Lambda)*(std::atan(p(1)/p(0))) ))/(std::pow(p(0),2)*(1+std::pow(p(1)/p(0),2)) )+
(std::pow((1-Lambda),2)*p(1)* std::cos((1-Lambda)* std::atan(p(1)/p(0))))/(std::pow(p(0),2)*(1+std::pow(p(1)/p(0),2)) )     )/(std::sqrt(1+std::pow(p(1)/p(0),2))) );

 


gradients[0][1]= 1/(std::pow(r,2))* (  (std::pow(r,Lambda))*Lambda*p(1)* ( ((1+Lambda)*p(1)* ((std::cos(Lambda*Omega)* std::sin((1+Lambda)*std::atan(p(1))/p(0)))/((1+Lambda))-  std::cos((1+Lambda)*std::atan(p(1)/p(0))) -  ( std::cos(Lambda*Omega)* std::sin((1-Lambda)*std::atan(p(1)/p(0)) ))/(1-Lambda) + std::cos((1-Lambda)* std::atan(p(1)/p(0))) ) )/( p(0)*std::sqrt(1+std::pow(p(1)/p(0) ,2)))  +
( std::cos(Lambda*Omega)*std::cos((1+Lambda)*std::atan(p(1)/p(0))) + (1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0)) ) -
					std::cos(Lambda*Omega)*std::cos((1-Lambda)* std::atan(p(1)/p(0))) -
					(1-Lambda)*std::sin((1-Lambda)*std::atan(p(1)/p(0))) )/(std::sqrt(1+std::pow(p(1)/p(0) ,2)))  ) )+ 
std::pow(r,Lambda)* (( (1+Lambda)* ( (std::sin((1+Lambda)*std::atan(p(1)/p(0)) )*std::cos(Lambda*Omega)  )/(1+Lambda) -  std::cos((1+Lambda)*std::atan(p(1)/p(0))) -(std::cos(Lambda*Omega)* std::sin((1-Lambda)* std::atan(p(1)/p(0)) ))/(1-Lambda) + std::cos((1-Lambda)*std::atan(p(1)/p(0))) ) )/(p(0)* std::sqrt(1+std::pow(p(1)/p(0) ,2)))
-((1+Lambda)*std::pow(p(1),2)*((std::cos(Lambda*Omega)* std::sin((1+Lambda)*std::atan(p(1)/p(0))))/(1+Lambda) -  std::cos((1+Lambda)*std::atan(p(1)/p(0))) -  (std::cos(Lambda*Omega)* std::sin((1-Lambda)*std::atan(p(1)/p(0))))/(1-Lambda) + std::cos((1-Lambda)*std::atan(p(1)/p(0)))  ))/( std::pow(p(0),3)* std::pow((1+std::pow(p(1)/p(0),2)), 3/2)  )
+((1+Lambda)*p(1)* ( ( std::cos(Lambda*Omega)* std::cos((1+Lambda)*std::atan(p(1)/p(0))) )/(p(0)* (1+std::pow(p(1)/p(0),2))  ) + ( (1+Lambda)* std::sin((1+Lambda)* std::atan(p(1)/p(0)) ) )/(p(0)* (1+std::pow(p(1)/p(0),2))  ) -
( (1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0)) )  )/(p(0)* (1+std::pow(p(1)/p(0),2))  )  -
( std::cos(Lambda*Omega)*std::cos((1-Lambda)*std::atan(p(1)/p(0)) ))/(p(0)* (1+std::pow(p(1)/p(0),2))  ) ) )/( p(0)* std::sqrt(1+ std::pow(p(1)/p(0),2))  )
-((p(1)*( std::cos(Lambda*Omega)*std::cos((1+Lambda)*std::atan(p(1)/p(0))) + (1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0))) - std::cos(Lambda*Omega)*std::cos((1-Lambda)*std::atan(p(1)/p(0)))  -  (1-Lambda)*std::sin((1-Lambda)* std::atan(p(1)/p(0))))))/(std::pow(p(0),2)*std::pow( (1+p(1)/p(0)) , (3/2) )  )
+( -( std::cos(Lambda*Omega)* (1+Lambda)* std::sin((1+Lambda)*std::atan(p(1)/p(0))) )/(p(0)*(1+std::pow(p(1)/p(0),2)) )  + 
( std::pow((1+Lambda),2)*std::cos((1+Lambda)*std::atan(p(1)/p(0)))  )/(p(0)*(1+std::pow(p(1)/p(0),2)) ) 
+((1-Lambda)* std::cos(Lambda*Omega)* std::sin( (1-Lambda)* std::atan(p(1)/p(0)))   )/(p(0)*(1+std::pow(p(1)/p(0),2)) )  -(( std::pow((1-Lambda),2)* std::cos((1-Lambda)* std::atan(p(1)/p(0)))))/(p(0)*(1+std::pow(p(1)/p(0),2)) ) )/( std::sqrt(1+ std::pow(p(1)/p(0),2)) ));






	gradients[1][0]= 1/(std::pow(r,2))* ( std::pow(r,Lambda)*Lambda*p(0)*( (p(1)* (  std::cos(Lambda*Omega)*std::cos((1+Lambda)*std::atan(p(1)/p(0))) +
			(1+Lambda)*std::sin( (1+Lambda)*std::atan(p(1)/p(0)) ) -
			std::cos(Lambda*Omega)*std::cos( (1-Lambda)*std::atan(p(1)/p(0)) )  - (1-Lambda)*
			std::sin( (1-Lambda)*std::atan(p(1)/p(0)) )  ))/(p(0)*std::sqrt(1+ std::pow(p(1)/p(0) ,2)))-
			( (1+Lambda)* ( (std::cos(Lambda*Omega)*std::sin((1+Lambda)*std::atan(p(1)/p(0))))/(1+Lambda) -
					std::cos( (1+Lambda)*std::atan(p(1)/p(0)) )  -( std::cos(Lambda*Omega)* std::sin((1-Lambda)*std::atan(p(1)/p(0))))/(1-Lambda) +
					std::cos( (1-Lambda)*std::atan(p(1)/p(0)) ) )   )/(p(0)*std::sqrt(1+ std::pow(p(1)/p(0) ,2))) ))+
					std::pow(r,Lambda) *  ( -(p(1)*(std::cos(Lambda*Omega)*std::cos((1+Lambda)*std::atan(p(1)/p(0)))
	+(1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0)))- std::cos(Lambda*Omega)*std::cos((1-Lambda)*std::atan(p(1)/p(0))) -
	(1-Lambda)* std::sin((1-Lambda)*std::atan(p(1)/p(0))) ))/(std::pow(p(0),2)* std::sqrt(1+std::pow(p(1)/p(0) ,2))) +
							(std::pow(p(1),3)*( std::cos(Lambda*Omega)*std::cos((1+Lambda)* std::atan(p(1)/p(0)))+ (1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0))) -
									std::cos(Lambda*Omega)*std::cos((1-Lambda)*std::atan(p(1)/p(0)))-(1-Lambda)*std::sin((1-Lambda)*std::atan(p(1)/p(0)))))/( std::pow(p(0),2)*std::pow( (1+p(1)/p(0)) , 3/2)) +
									(p(1)*( ((1+Lambda)*p(1)* std::cos(Lambda*Omega)* std::sin((1+Lambda)* std::atan(p(1)/p(0))) )/(std::pow(p(0),2)* (1+ std::pow( p(1)/p(0) ,2))) -
											( p(1)* std::pow((1+Lambda) ,2)* std::cos( (1+Lambda)*std::atan(p(1)/p(0)) ) )/(std::pow(p(0),2)* (1+ std::pow( p(1)/p(0) ,2)))
											-( (1-Lambda)*p(1)*std::cos(Lambda*Omega)* std::sin((1-Lambda)*std::atan(p(1)/p(0))) )/(std::pow(p(0),2)* (1+ std::pow( p(1)/p(0) ,2)))
											+(p(1)*std::pow( (1-Lambda) ,2)*std::cos((1-Lambda)*std::atan(p(1)/p(0))))/(std::pow(p(0),2)* (1+ std::pow( p(1)/p(0) ,2))) ))/(p(0)* std::sqrt(1+ std::pow(p(1)/p(0), 2) )) -
											( (1+Lambda)*std::pow(p(1),2)* ( (std::cos(Lambda*Omega)* std::sin( (1+Lambda)*std::atan(p(1)/p(0))) )/(1+Lambda) -  std::cos((1+Lambda)*std::atan(p(1)/p(0))) - (std::cos(Lambda*Omega)*std::sin((1-Lambda)*std::atan(p(1)/p(0))))/(1-Lambda) + std::cos((1-Lambda)*std::atan(p(1)/p(0)))  ))/( std::pow(p(1),3)*  std::pow(  1+std::pow( p(1)/p(0) ,2) ,3/2)) -
											((1+Lambda)* ( -(p(1)*std::cos(Lambda*Omega)*std::cos((1+Lambda)*std::atan(p(1)/p(0))))/(std::pow(p(0),2)* (1+ std::pow( p(1)/p(0) ,2))) -
													((1+Lambda)*p(1)* std::sin((1+Lambda)*std::atan(p(1)/p(0))) )/(std::pow(p(0),2)* (1+ std::pow( p(1)/p(0) ,2)))
													+(p(1)*std::cos(Lambda*Omega)* std::cos((1-Lambda)*std::atan(p(1)/p(0))))/(std::pow(p(0),2)* (1+ std::pow( p(1)/p(0) ,2)))
													+((1-Lambda)*p(1)* std::sin((1-Lambda)*std::atan(p(1)/p(0))))/(std::pow(p(0),2)* (1+ std::pow( p(1)/p(0) ,2))) ) )/(std::sqrt(1+std::pow(p(1)/p(0),2))) );



	gradients[1][1]= 1/(std::pow(r,2))* ( std::pow(r,Lambda) *Lambda* p(1)* ( (p(1)*( std::cos(Lambda*Omega)* std::cos((1+Lambda)*std::atan(p(1)/p(0)))+(1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0)))-std::cos(Lambda*Omega)*std::cos((1-Lambda)*std::atan(p(1)/p(0)))-(1-Lambda)*std::sin((1-Lambda)*std::atan(p(1)/p(0)))  ) )/(p(0)*std::sqrt(1+std::pow(p(1)/p(0),2)))
			- ( (1+Lambda)* ((std::sin((1+Lambda)*std::atan(p(1)/p(0)))*  std::cos(Lambda*Omega))/(1+Lambda) - std::cos((1+Lambda)*std::atan(p(1)/p(0)))  -(std::cos(Lambda*Omega)*std::sin((1-Lambda)*std::atan(p(1)/p(0)))  )/(1-Lambda)    + std::cos((1-Lambda)*std::atan(p(1)/p(0)))  ) )/(std::sqrt(1+std::pow(p(1)/p(0),2)))  ))

			+std::pow(r,Lambda)* ( (std::cos(Lambda*Omega)* std::cos((1+Lambda)*std::atan(p(1)/p(0))) + (1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0))) -std::cos(Lambda*Omega)* std::cos((1-Lambda)*std::atan(p(1)/p(0))) - (1-Lambda)*std::sin((1-Lambda)*std::atan(p(1)/p(0))) )/(p(0)*std::sqrt(1+ std::pow(p(1)/p(0) ,2))  )
					- ( std::pow(p(1),2)* (  std::cos(Lambda*Omega)* std::cos((1+Lambda)*std::atan(p(1)/p(0)))   + (1+Lambda)*std::sin((1+Lambda) *std::atan(p(1)/p(0))) - std::cos(Lambda*Omega)*std::cos((1-Lambda)*std::atan(p(1)/p(0)))  - (1-Lambda)*std::sin((1-Lambda)*std::atan(p(1)/p(0)))  ))/( std::pow(p(0),3)* std::pow( (1+ std::pow(p(1)/p(0),2) ), (3/2) ))
					+ (p(1)* (-( (1+Lambda)*std::cos(Lambda*Omega)*std::sin((1+Lambda)*std::atan(p(1)/p(0))))/(p(0)*(std::sqrt(1+std::pow(p(1)/p(0),2))))
							+ (std::pow((1+Lambda),2)*std::cos( (1+Lambda)* std::atan(p(1)/p(0)))  )/(p(0)*(std::sqrt(1+std::pow(p(1)/p(0),2))))
							+((1-Lambda)*std::cos(Lambda*Omega)*std::sin((1-Lambda)*std::atan(p(1)/p(0)))  )/(p(0)*(std::sqrt(1+std::pow(p(1)/p(0),2))))
							- ( std::pow((1-Lambda),2) * std::cos((1-Lambda)*std::atan(p(1)/p(0)))  )/(p(0)*(std::sqrt(1+std::pow(p(1)/p(0),2)))) ) )/(p(0)* std::sqrt(1+ std::pow(p(1)/p(0),2)))
							+( (1+Lambda)*p(1)  * ((std::cos(Lambda*Omega)* std::sin((1+Lambda)*std::atan(p(1)/p(0))) )/(1+Lambda) - std::cos((1+Lambda)*std::atan(p(1)/p(0)))  - ( std::cos(Lambda*Omega)* std::sin((1-Lambda)*std::atan(p(1)/p(0))) )/(1-Lambda)  + std::cos((1-Lambda)*std::atan(p(1)/p(0)))))/(std::pow(p(0),2)* std::pow( (1+std::pow(p(1)/p(0),2)) , (3/2)) )
							-(1+Lambda)*( (std::cos(Lambda*Omega)* std::cos((1+Lambda)*std::atan(p(1)/p(0))) )/( p(0)*(1+std::pow(p(1)/p(0),2)) ) +
									((1+Lambda)*std::sin((1+Lambda)*std::atan(p(1)/p(0))))/(p(0)*(1+std::pow(p(1)/p(0),2))) -
									(std::cos(Lambda*Omega)* std::cos((1-Lambda)*std::atan(p(1)/p(0))) )/(p(0)*(1+std::pow(p(1)/p(0),2))) -
									( (1-Lambda)* std::sin((1-Lambda)*std::atan(p(1)/p(0))) )/(p(0)*(1+std::pow(p(1)/p(0),2)))   )/(std::sqrt(1+std::pow(p(1)/p(0) ,2))));


	gradients[2][0]= 0 ;
	gradients[2][1]=  0 ;

		}



//Explicit initialization

template class ExactSolutionEx2<2>;
