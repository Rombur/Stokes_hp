#include "StokesProblem.hh"
#include "Parameters.hh"

using namespace dealii;

int main(int argc, char* argv[])
{
	try
	{
		deallog.depth_console (0);
    Parameters parameters(argv[1]);

    if (parameters.get_dim()==2)
    {
      StokesProblem<2> stokesproblem(parameters.get_verbose(),
          parameters.get_exemple(), parameters.get_quadrature(), 
          parameters.get_max_degree(), parameters.get_max_n_cycles(),
          parameters.get_theta(), parameters.get_tolerance());
      stokesproblem.run();
    }
    else
      AssertThrow(false,ExcMessage("3D not implemented yet."));
	  
	}
	catch (std::exception &exc)
	{
		std::cerr << std::endl << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
			<< exc.what() << std::endl
			<< "Aborting!" << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		return 1;
	}
	catch (...)
	{
		std::cerr << std::endl << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
			<< "Aborting!" << std::endl
			<< "----------------------------------------------------"
			<< std::endl;
		return 1;
	}
	return 0;
}
