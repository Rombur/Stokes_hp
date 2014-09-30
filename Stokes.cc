// This code works up to function estimate ()


#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <set>

#include <vector>
#include <algorithm>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nothing.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#define PRINT(expr) {std::cout << ""#expr" = " << std::setprecision(15) << (expr) << std::endl;}
/*----------------------------------------------------------------------------------------*/

namespace hp_Stokes
{
	using namespace dealii;
	using namespace std;

	 template <int dim>
	class ExactSolution : public Function<dim>
	{
	public:
		ExactSolution () : Function<dim>(dim+1) {}
		virtual void vector_value (const Point<dim> &p,
			Vector<double>   &value) const;
		virtual void vector_gradient (const Point<dim> &p,
			vector< Tensor< 1, dim > > &gradients) const;
	};
	  template <int dim>
	void ExactSolution <dim>::vector_value (const Point<dim> &p,
		Vector<double>   &values) const
	{
		values(0) = -1.0 * exp (p(0)) * (p(1) * cos (p(1)) + sin (p(1)));
		values(1) = exp(p(0)) * p(1) * sin (p(1));
		values(2) = 2.0 * exp (p(0)) * sin (p(1)) - (2.0 * (1.0 - exp(1.0)) * (cos (1.0) - 1.0)) / 3.0;
	}
	  template <int dim>
	void ExactSolution <dim>::vector_gradient  (const Point<dim> &p,
		vector< Tensor< 1, dim > > &gradients) const
	{
		gradients[0][0]= -1.0 * exp (p(0)) * (p(1) * cos (p(1)) + sin (p(1)));
		gradients[0][1]= -1.0 * exp (p(0)) * ( 2.0 * cos (p(1)) - p(1) * sin (p(1)));
		gradients[1][0]=  exp(p(0)) * p(1) * sin (p(1));
		gradients[1][1]=  exp(p(0)) * (sin (p(1)) + p(1) * cos(p(1)));
		gradients[2][0]=  2.0 * exp (p(0)) * sin (p(1));
		gradients[2][1]=  2.0 * exp (p(0)) * cos (p(1));
	}
	/*......................................................................................*/
	// Calculate Right Hand Side
	template <int dim>
	class RightHandSide : public Function<dim>
	{
	public:
		RightHandSide () : Function<dim>(dim+1) {}
		virtual double value (const Point<dim> &p, const unsigned int component = 0) const;
		virtual void vector_value (const Point<dim> &p, Vector<double> &value) const;
	};
	template <int dim>
	double RightHandSide <dim>::value (const Point<dim> & p, const unsigned int  component) const
	{
		return 0;
	}
	template <int dim>
	void
		RightHandSide <dim>::vector_value (const Point<dim> &p, Vector<double> &values) const
	{
		for (unsigned int c=0; c<this->n_components; ++c)
			values(c) = RightHandSide::value (p, c);
	}
	/*......................................................................................*/
	template <int dim>
	class StokesProblem
	{
	public:
		StokesProblem ();
		~StokesProblem();
		void run();

	private:
		//enum
		//{
		//	patch_on,
		//	patch_off
		//};

		//static bool cell_is_on_patch (const typename hp::DoFHandler<2>::cell_iterator &cell);
		//static bool cell_is_NOT_on_patch (const typename hp::DoFHandler<2>::cell_iterator &cell);

		const RightHandSide<dim> rhs_function;
		const ExactSolution<dim> exact_solution;

		void generate_mesh ();
		void setup_system ();
		void assemble_system ();
		void solve ();
		double pressure_mean_value () const;
		void compute_error ();
		void estimate (Vector<double> &est_per_cell);

		//void h_patch_conv_load_no (double &h_convergence_est_per_cell, unsigned int &h_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell) const;
		//void p_patch_conv_load_no (double &p_convergence_est_per_cell, unsigned int &p_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell ) const;

		//vector<hp::DoFHandler<dim>::active_cell_iterator> get_patch_around_cell( typename hp::DoFHandler<dim>::active_cell_iterator &cell);
		//vector<hp::DoFHandler<dim>::active_cell_iterator> get_cells_at_coarsest_common_level ( vector<hp::DoFHandler<dim>::active_cell_iterator>  &patch);
		//void build_triangulation_from_patch (vector<hp::DoFHandler<dim>::active_cell_iterator>  &patch, Triangulation<dim> &tria_patch, unsigned int &level_h_refine, unsigned int &level_p_refine);

		//bool decreasing (pair<double,hp::DoFHandler<dim> > &i, pair<double,hp::DoFHandler<dim> > &j);

		//void postprocess (const unsigned int cycle);

		//const double Tolerance;
		

		Triangulation<dim> triangulation;
		hp::DoFHandler<dim> dof_handler;

		//FESystem<2>     stokes_fe;
		//FESystem<2>     fe_no;
		//Triangulation<2> tria_patch;
		//hp::DoFHandler<2> dof_handler_patch;

		hp::FECollection<dim> fe_collection;
		hp::QCollection<dim> quadrature_collection;
		hp::QCollection<dim-1> face_quadrature_collection;

		ConstraintMatrix constraints;
		BlockSparsityPattern sparsity_pattern;
		BlockSparseMatrix<double> system_matrix;

		BlockVector<double> solution;
		BlockVector<double> system_rhs;

                const unsigned int max_degree;

		double L2_norm_est;
		double L1_norm_est;
		Vector<double> est_per_cell;

	};
	/*......................................................................................*/

	class SchurComplement : public Subscriptor
	{
	public:
		SchurComplement (const BlockSparseMatrix<double> &system_matrix,
			const SparseDirectUMFPACK &A_inverse);
		void vmult (Vector<double>       &dst,
			const Vector<double> &src) const;
	private:
		const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
		const SmartPointer<const SparseDirectUMFPACK> A_inverse;
		mutable Vector<double> tmp1, tmp2;
	};
	SchurComplement::
		SchurComplement (const BlockSparseMatrix<double> &system_matrix,
		const SparseDirectUMFPACK &A_inverse)
		:
	system_matrix (&system_matrix),
		A_inverse (&A_inverse),
		tmp1 (system_matrix.block(0,0).m()),
		tmp2 (system_matrix.block(0,0).m())
	{}
	void SchurComplement::vmult (Vector<double>       &dst,
		const Vector<double> &src) const
	{
		system_matrix->block(0,1).vmult (tmp1, src);
		A_inverse->vmult (tmp2, tmp1);
		system_matrix->block(1,0).vmult (dst, tmp2);
	}
	/*......................................................................................*/
	// constructor and destructor
	//Tolerance (0.0001),

	template <int dim>
	StokesProblem<dim>::StokesProblem(): dof_handler(triangulation),  max_degree (5)


		//stokes_fe (FE_Q<2>(QGaussLobatto<1> (unsigned int degree + 2)), 2, FE_Q<2> (QGaussLobatto<1> (unsigned int degree + 1)), 1),
		//fe_no (FE_Nothing<2>(), 3)
	{
		for (unsigned int degree=1; degree<=max_degree; ++degree)
		{

			fe_collection.push_back(FESystem<dim>(FE_Q<dim> (QGaussLobatto<1> (degree + 2)), dim, FE_Q<dim> (QGaussLobatto<1> (degree + 1)), 1));//QGaussLobatto<1>?
			quadrature_collection.push_back(QGauss<dim> (degree+2));
			face_quadrature_collection.push_back (QGauss<dim-1> (degree+1));
		}
		//fe_collection.push_back(FESystem<dim> (FE_Nothing<dim>(), dim+1) );
		//quadrature_collection.push_back(QGauss<dim>(0));
		//face_quadrature_collection.push_back (QGauss<dim-1>(0));
	}
	/*.....................................................................................*/
	template <int dim>
	StokesProblem <dim>::~StokesProblem() {
		dof_handler.clear();
	}
	/*......................................................................................*/
	/*
	template <int dim>
	bool StokesProblem <dim>::decreasing (pair<double,hp::DoFHandler<dim> > &i, pair<double,hp::DoFHandler<dim> > &j)
	{
		return ((i.first) > (j.first));
	}
	*/
	/*......................................................................................*/
	/*
	template <int dim>
	bool StokesProblem <dim>::cell_is_on_patch (const typename hp::DoFHandler<dim>::cell_iterator &cell)
	{
		return (cell->material_id() == patch_on);
	}
	*/
	/*......................................................................................*/
	/*
	template <int dim>
	bool StokesProblem <dim>::cell_is_NOT_on_patch (const typename hp::DoFHandler<dim>::cell_iterator &cell)
	{
		return (cell->material_id() == patch_off);
	}
	*/
	/*......................................................................................*/
	// Generate mesh
	template <int dim>
	void StokesProblem <dim>::generate_mesh(){

		vector<Point<dim>> vertices (8);

		vertices [0]=Point<dim> (-1,-1);
		vertices [1]=Point<dim> (0,-1);
		vertices [2]=Point<dim> (-1,0);
		vertices [3]=Point<dim> (0,0);
		vertices [4]=Point<dim> (1,0);
		vertices [5]=Point<dim> (-1,1);
		vertices [6]=Point<dim> (0,1);
		vertices [7]=Point<dim> (1,1);

		const unsigned int n_cells=3;
		vector<CellData<dim> > cell(n_cells);
		cell[0].vertices[0]=0;
		cell[0].vertices[1]=1;
		cell[0].vertices[2]=2;
		cell[0].vertices[3]=3;

		cell[1].vertices[0]=2;
		cell[1].vertices[1]=3;
		cell[1].vertices[2]=5;
		cell[1].vertices[3]=6;

		cell[2].vertices[0]=3;
		cell[2].vertices[1]=4;
		cell[2].vertices[2]=6;
		cell[2].vertices[3]=7;

		triangulation.create_triangulation(vertices,cell,SubCellData());
		triangulation.refine_global (1);
		ofstream out ("grid-L-Shape.eps");
		GridOut grid_out;
		grid_out.write_eps (triangulation, out);
		cout<<"Number of active cells: "<< triangulation.n_active_cells() << endl;
		cout<<"Total number of cells: " << triangulation.n_cells() << endl ;
	}
	/*......................................................................................*/
	// setup system()
	template <int dim>
	void StokesProblem <dim>::setup_system(){
		dof_handler.distribute_dofs (fe_collection);

		vector<unsigned int> block_component (dim+1, 0);
		block_component[dim]=1;
		DoFRenumbering::component_wise(dof_handler, block_component);

		{
			constraints.clear ();
			FEValuesExtractors::Vector velocities(0);
			DoFTools::make_hanging_node_constraints (dof_handler, constraints);
			VectorTools::interpolate_boundary_values (dof_handler,0,exact_solution,constraints, fe_collection.component_mask(velocities));
		}
		constraints.close();

		vector<types::global_dof_index> dofs_per_block (2);
		DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
		const unsigned int n_u=dofs_per_block[0], n_p=dofs_per_block[1];

		cout<< "Number of degrees of freedom: " << dof_handler. n_dofs()<<
			"(" << n_u << "+" << n_p << ")" << endl;
		{
			BlockCompressedSetSparsityPattern csp (dofs_per_block,dofs_per_block);

			DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
			sparsity_pattern.copy_from(csp);
		}

		system_matrix.reinit (sparsity_pattern);
		solution.reinit (dofs_per_block);
		system_rhs.reinit (dofs_per_block);

	}// End of function setup_system()
	/*......................................................................................*/
	// assemble system
	template <int dim>
	void StokesProblem <dim>::assemble_system () {
		hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients);


		FullMatrix<double> local_matrix;
		Vector<double> local_rhs;
		vector<types::global_dof_index> local_dof_indices;

		vector<Vector<double> >  rhs_values;

		const FEValuesExtractors::Vector velocities (0);
		const FEValuesExtractors::Scalar pressure (dim);

		vector<Tensor<2,dim> > grad_phi_u;
		vector<double> div_phi_u;
		vector<Tensor<1,dim> > phi_u;
		vector<double> phi_p;

		typename hp::DoFHandler<dim>::active_cell_iterator
			cell = dof_handler.begin_active(),
			endc = dof_handler.end();
		for (; cell!=endc; ++cell)
		{
			const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
			local_matrix.reinit (dofs_per_cell, dofs_per_cell);
			local_rhs.reinit (dofs_per_cell);

			hp_fe_values.reinit (cell);
			const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
			const vector<double>& JxW_values = fe_values.get_JxW_values ();
			const unsigned int n_q_points = fe_values.n_quadrature_points;

			rhs_values.resize(n_q_points, Vector<double>(dim+1));
			rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

			grad_phi_u.resize(dofs_per_cell);
			div_phi_u.resize(dofs_per_cell);
			phi_u.resize (dofs_per_cell);
			phi_p.resize(dofs_per_cell);

			for (unsigned int q=0; q<n_q_points; ++q)
			{
				for (unsigned int k=0; k<dofs_per_cell; ++k)
				{
					grad_phi_u[k] = fe_values[velocities].gradient (k, q);
					div_phi_u[k] = fe_values[velocities].divergence (k, q);
					phi_u[k] = fe_values[velocities].value (k, q);
					phi_p[k] = fe_values[pressure].value (k, q);
				}
				for (unsigned int i=0; i<dofs_per_cell; ++i)
				{
					for (unsigned int j=0; j<=i; ++j)
					{
						local_matrix(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
							- div_phi_u[i] * phi_p[j]
						- phi_p[i] * div_phi_u[j]
						+ phi_p[i] * phi_p[j])
							* JxW_values[q];
					} // end of loop over 'j'
					local_rhs(i) += (phi_u[i][0] * rhs_values[q](0) + phi_u[i][1] * rhs_values [q](1)) * JxW_values[q];//?
				} // end of loop 'i'
			} // end of loop 'q'
			// to create the upper triangle of local_matrix

			for (unsigned int i=0; i<dofs_per_cell; ++i)
				for (unsigned int j=i+1; j<dofs_per_cell; ++j)
					local_matrix(i,j) = local_matrix(j,i);
			//local system to global system
			local_dof_indices.resize (dofs_per_cell);
			cell->get_dof_indices (local_dof_indices);
			constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
		} // end of iteration for cells
	}
	/*......................................................................................*/
	template <int dim>
	void StokesProblem <dim>::solve ()
	{
		SparseDirectUMFPACK A_inverse;
		A_inverse.initialize (system_matrix.block(0,0),
			SparseDirectUMFPACK::AdditionalData());
		Vector<double> tmp (solution.block(0).size());
		{
			Vector<double> schur_rhs (solution.block(1).size());
			A_inverse.vmult (tmp, system_rhs.block(0));
			system_matrix.block(1,0).vmult (schur_rhs, tmp);
			schur_rhs -= system_rhs.block(1);

			SchurComplement schur_complement (system_matrix, A_inverse);
			SolverControl solver_control (solution.block(1).size(),
					1e-4);

			SolverCG<>    cg (solver_control);

			SparseDirectUMFPACK preconditioner;
			preconditioner.initialize (system_matrix.block(1,1),
					SparseDirectUMFPACK::AdditionalData());

			cg.solve (schur_complement, solution.block(1), schur_rhs,
				preconditioner);
			//cout<<" residuals of each step " << solver_control.enable_history_data() << endl;
			constraints.distribute (solution);
			cout << "  "
				<< solver_control.last_step()
				<< " outer CG Schur complement iterations for pressure"
				<< endl;
		}
		system_matrix.block(0,1).vmult (tmp, solution.block(1));
		tmp *= -1.0;
		tmp += system_rhs.block(0);
		A_inverse.vmult (solution.block(0), tmp);

		constraints.distribute (solution);
		solution.block (1).add (-1.0 * pressure_mean_value ());
		constraints.distribute (solution);
	}
	/*......................................................................................*/
	template <int dim>
	double StokesProblem <dim>::pressure_mean_value () const
	{
		// get pressure such that satisfies mean value property:
		hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_JxW_values);
		const FEValuesExtractors::Scalar pressure (dim);

		vector<double> values;
		double domain_mean_val_p=0;
		// double measure_domain=0;
		typename hp::DoFHandler<dim>::active_cell_iterator
			cell = dof_handler.begin_active(),
			endc = dof_handler.end();
		for (; cell!=endc; ++cell)
		{
			hp_fe_values.reinit (cell);

			const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
			const vector<double>& JxW_values = fe_values.get_JxW_values ();
			const unsigned int n_q_points = fe_values.n_quadrature_points;
			values.resize(n_q_points);
			fe_values[pressure].get_function_values(solution, values);
			for (unsigned int q=0; q<n_q_points; ++q)
			{
				domain_mean_val_p += values[q]*JxW_values[q];
				// measure_domain += JxW_values[q];
			}//q
		}//cell
		// 3 here is the area corresponding to our 3 cells.
		return domain_mean_val_p/3.0;
		// return domain_mean_val_p / measure_domain;
	}
	/*.................................................................................................*/
	template <int dim>
	void StokesProblem <dim>::compute_error ()
	{
		hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients);
		const FEValuesExtractors::Vector velocities (0);
		const FEValuesExtractors::Scalar pressure (dim);

		vector<double> values;
		vector<Tensor<2,dim>> gradients;
		vector<vector<Tensor<1,dim>>> exact_solution_gradients;
		vector<Vector<double>> exact_solution_values;

		Vector<double> error_per_cell(triangulation.n_active_cells());

		typename hp::DoFHandler<dim>::active_cell_iterator
			cell = dof_handler.begin_active(),
			endc = dof_handler.end();
		unsigned int cell_index=0;
		for (; cell!=endc; ++cell,++cell_index)
		{
			hp_fe_values.reinit (cell);
			const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
			const vector<double>& JxW_values = fe_values.get_JxW_values ();
			const vector<Point<dim>>& quadrature_points = fe_values.get_quadrature_points();
			const unsigned int n_q_points = fe_values.n_quadrature_points;

			gradients.resize(n_q_points);
			values.resize(n_q_points);
			exact_solution_gradients.resize(n_q_points , vector<Tensor<1,dim>> (dim+1));//?
			exact_solution_values.resize(n_q_points, Vector<double> (dim+1));//?

			fe_values[velocities].get_function_gradients(solution, gradients);
			fe_values[pressure].get_function_values(solution, values);

			exact_solution.vector_gradient_list(quadrature_points, exact_solution_gradients);
			exact_solution.vector_value_list(quadrature_points, exact_solution_values);

			double subtract_p=0;
			double grad_u_vals=0;
			for (unsigned int q=0; q<n_q_points; ++q)
			{
				values[q] -= exact_solution_values[q](dim);//?
				subtract_p +=values[q]*values[q]* JxW_values[q];
				for (unsigned int i=0; i<dim; ++i)
					gradients[q][i]-=exact_solution_gradients[q][i];
				grad_u_vals +=double_contract(gradients[q],gradients[q])* JxW_values[q];
			} // q
			error_per_cell(cell_index) =(sqrt(subtract_p) + sqrt(grad_u_vals));

		}// cell
		cout<< "Vector of Compute Error per Cell: " << error_per_cell<< endl ;
		double L1_norm=error_per_cell.l1_norm();
		cout<< "L1_norm of ERROR is: "<< L1_norm << endl;
		double L2_norm=error_per_cell.l2_norm();
		cout<< "L2_norm of ERROR is: "<< L2_norm << endl;
	}
	/*......................................................................................*/
	
	template <int dim>
	void StokesProblem <dim>::estimate (Vector<double> &est_per_cell)  {
		hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);
		hp::FEFaceValues<dim> hp_fe_face_values(fe_collection, face_quadrature_collection, update_JxW_values|update_gradients|update_normal_vectors);
		hp::FEFaceValues<dim> hp_neighbor_face_values(fe_collection, face_quadrature_collection, update_gradients);
		hp::FESubfaceValues<dim> hp_subface_values(fe_collection, face_quadrature_collection, update_JxW_values|update_gradients|update_normal_vectors);
		hp::FESubfaceValues<dim> hp_neighbor_subface_values(fe_collection, face_quadrature_collection, update_gradients);

		vector<Tensor<1,dim>> gradients_p;
		vector<double> divergences;
		vector<Tensor<1,dim>> laplacians;

		vector<Tensor<2,dim>> gradients;
		vector<Tensor<2,dim>> neighbor_gradients;

		const FEValuesExtractors::Vector velocities (0);
		const FEValuesExtractors::Scalar pressure (dim);

		const RightHandSide<dim> rhs_function;
		vector<Vector<double> >  rhs_values;

		//est_per_cell.reinit (triangulation.n_active_cells());
		Vector<double> res_est_per_cell(triangulation.n_active_cells());
		Vector<double> Jump_est_per_cell(triangulation.n_active_cells());

		typename hp::DoFHandler<dim>::active_cell_iterator
			cell = dof_handler.begin_active(),
			endc = dof_handler.end();
		unsigned int cell_index=0;
		for (; cell!=endc; ++cell,++cell_index)
		{
			hp_fe_values.reinit (cell);
			const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
			const vector<double>& JxW_values = fe_values.get_JxW_values ();
			const unsigned int n_q_points = fe_values.n_quadrature_points;

			rhs_values.resize(n_q_points, Vector<double>(dim+1));
			rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);


			divergences.resize(n_q_points);
			gradients_p.resize(n_q_points) ;
			laplacians.resize(n_q_points);

			fe_values[pressure].get_function_gradients(solution, gradients_p);
			fe_values[velocities].get_function_divergences(solution, divergences);
			fe_values[velocities].get_function_laplacians(solution, laplacians);
			//PRINT(gradients_p[0]);

			double term2=0;//divergence term in estimator definition
			double term1=0;// for the residual term in estimator definition
			for (unsigned int q=0; q<n_q_points; ++q){
				term2 += (divergences[q])*(divergences[q])*JxW_values[q];

				for (unsigned int i=0; i<2; ++i)
					gradients_p[q][i]-= (rhs_values[q](i)+laplacians[q][i]);

				term1+= contract(gradients_p[q],gradients_p[q])*JxW_values[q];
			}// q
			res_est_per_cell(cell_index)= pow((cell->diameter())/(cell->get_fe().degree), 2.0 ) * (term1) + term2;

			// ........................................... compute jump_est_per_cell..............................................................
			double term3=0;//jumpped part of the estimator
			for (unsigned int face_number=0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)

				if ((cell->face(face_number)->at_boundary()==false)	&& (cell->face(face_number)->has_children() == false) && (cell->face(face_number)->level() == cell->level()))
				{
					const unsigned int q_index = std::max (cell->active_fe_index(),
						cell->neighbor(face_number)->active_fe_index());

					hp_fe_face_values.reinit       (cell,                        face_number,                             q_index);
					hp_neighbor_face_values.reinit (cell->neighbor(face_number), cell->neighbor_of_neighbor(face_number), q_index);

					const FEFaceValues<2> &neighbor_face_values =hp_neighbor_face_values.get_present_fe_values ();//?
					const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values ();//?

					const vector<double>& JxW_values = fe_face_values.get_JxW_values ();

					const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

					gradients.resize(n_face_q_points);
					neighbor_gradients.resize(n_face_q_points);

					neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
					fe_face_values[velocities].get_function_gradients(solution, gradients);

					vector<Tensor<1,dim>> jump_per_face;//?
					jump_per_face.resize(n_face_q_points);
					double jump_val=0;
					for (unsigned int q=0; q<n_face_q_points; ++q)
					{
						for (unsigned int i=0; i<2; ++i){
							for (unsigned int j=0; j<2; ++j){
							jump_per_face[q][i] = (gradients[q][i][j]-neighbor_gradients[q][i][j]) *(fe_face_values.normal_vector(q)[j]);
							}
						}
						jump_val += contract(jump_per_face[q],jump_per_face[q])*JxW_values[q];
					}//q_per_face
					term3 +=(cell->face(face_number)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;
				} //if level

				// if the neighbor has children

				else if ( (cell->face(face_number)->at_boundary()==false) && (cell->face(face_number)->has_children() == true))
				{
					for (unsigned int subface=0;
						subface< cell->face(face_number)->n_children(); ++subface)
					{
						const unsigned int q_index = std::max(cell->neighbor_child_on_subface (face_number, subface)->active_fe_index(), cell->active_fe_index());

						hp_neighbor_face_values.reinit (cell->neighbor_child_on_subface (face_number, subface), cell->neighbor_of_neighbor(face_number), q_index);
						hp_subface_values.reinit (cell,face_number, subface, q_index);

						const FEFaceValues<2> &neighbor_face_values  = hp_neighbor_face_values.get_present_fe_values ();
						const FESubfaceValues<2> &fe_subface_values = hp_subface_values.get_present_fe_values ();


						const vector<double>& JxW_values = fe_subface_values.get_JxW_values ();

						const unsigned int n_subface_q_points = fe_subface_values.n_quadrature_points;//?

						gradients.resize(n_subface_q_points);
						neighbor_gradients.resize(n_subface_q_points);

						neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
						fe_subface_values[velocities].get_function_gradients(solution, gradients);

						vector<Tensor<1,dim>> jump_per_subface;
						jump_per_subface.resize(n_subface_q_points);//?

						double jump_val=0;
						for (unsigned int q=0; q<n_subface_q_points; ++q)
						{
							for (unsigned int i=0; i<2; ++i){
								for (unsigned int j=0; j<2; ++j){
								jump_per_subface[q][j] += (gradients[q][i][j]- neighbor_gradients[q][i][j])*(fe_subface_values.normal_vector(q)[j]);
								}//j
							}// i
							jump_val += contract(jump_per_subface[q],jump_per_subface[q])*(JxW_values[q]);
						}//q_per_subface
						term3 +=(cell->face(face_number)->child(subface)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;
					}// subface
				}//else if

				// if the neighbor is coarser

				else if ( (cell->face(face_number)->at_boundary()==false) && (cell->neighbor_is_coarser(face_number)))

				{
					const unsigned int q_index = std::max(cell->active_fe_index(),cell->neighbor(face_number)->active_fe_index());
					hp_fe_face_values.reinit(cell, face_number,q_index);
					hp_neighbor_subface_values.reinit(cell->neighbor(face_number),cell->neighbor_of_coarser_neighbor(face_number).first, cell->neighbor_of_coarser_neighbor(face_number).second,q_index);

					const FEFaceValues<dim> &fe_face_values  = hp_fe_face_values.get_present_fe_values ();
					const FESubfaceValues<dim> &neighbor_subface_values = hp_neighbor_subface_values.get_present_fe_values ();


					const vector<double>& JxW_values = fe_face_values.get_JxW_values ();

					const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

					gradients.resize(n_face_q_points);
					neighbor_gradients.resize(n_face_q_points);

					neighbor_subface_values[velocities].get_function_gradients(solution, neighbor_gradients);
					fe_face_values[velocities].get_function_gradients(solution, gradients);

					vector<Tensor<1,dim>> jump_per_face;
					jump_per_face.resize(n_face_q_points);
					double jump_val=0;
					for (unsigned int q=0;
						q<n_face_q_points; ++q)
					{
						for (unsigned int i=0; i<2; ++i){
							for (unsigned int j=0; j<2; ++j){
							jump_per_face[q][i] += (gradients[q][i][j]- neighbor_gradients[q][i][j])*(fe_face_values.normal_vector(q)[j]);
							}//j
						}// i
						jump_val += contract(jump_per_face[q],jump_per_face[q])*JxW_values[q];
					}//q_per_face
					term3 +=(cell->face(face_number)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;

				} // else if coarse neighbor

				Jump_est_per_cell(cell_index) = term3;
				est_per_cell(cell_index)=sqrt(Jump_est_per_cell(cell_index)+res_est_per_cell(cell_index));
		}//cell
		
		cout<< "Vector of Error Estimate: "<< est_per_cell << endl;
		L1_norm_est= est_per_cell.l1_norm();
		cout<< "L1_norm of ERROR Estimate is: "<< L1_norm_est << endl;
		L2_norm_est= est_per_cell.l2_norm();	
		cout<< "L2_norm of ERROR Estimate is: "<< L2_norm_est << endl;
		

	}// func.estimate ()
	
	/*......................................................................................*/
	/*
	template <int dim>
	vector<hp::DoFHandler<dim>::active_cell_iterator> StokesProblem <dim>::get_patch_around_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell)
	{
		vector<hp::DoFHandler<dim>::active_cell_iterator> patch;
		patch.push_back (cell);
		for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
			if (cell->face(face_number)->at_boundary()==false)
			{
				if (cell->face(face_number)->has_children() == false)
					patch.push_back (cell->neighbor(face_number));
				else
					for (unsigned int subface=0; subface< cell->face(face_number)->n_children(); ++subface)
						patch.push_back (cell->neighbor_child_on_subface (face_number, subface));
			}
			return patch;
	}
	*/
	/*..............................................................................................................................................*/
	/*
	template <int dim>
	vector<hp::DoFHandler<dim>::active_cell_iterator> StokesProblem <dim>::get_cells_at_coarsest_common_level (const vector<hp::DoFHandler<dim>::active_cell_iterator>  &patch)
	{
		Assert (patch.size() > 0, ExcMessage("vector containing patch cells should not be an empty vector!"));//?
		unsigned int min_level = patch[0]->level();
		unsigned int max_level = patch[0]->level();
		for (unsigned int i=0; i<patch.size();++i)
		{
			min_level = std::min (min_level, patch[i]->level());
			max_level = std::max (max_level, patch[i]->level());
		}
		if (min_level == max_level)
			return patch;
		else
		{
			set<hp::DoFHandler<dim>::active_cell_iterator>  uniform_cells;
			typename hp::DoFHandler<dim>::active_cell_iterator  patch_c=patch.begin();
			const typename hp::DoFHandler<dim>::active_cell_iterator end_c=patch.end ();
			for (; patch_c!=end_c ; ++patch_c){
				if (patch_c->level() == min_level)
					uniform_cells.push_back (patch_c);
				else
				{
					hp::DoFHandler<dim>::active_cell_iterator parent = patch_c;

					while (parent->level() > min_level)
						parent = parent->parent();
					uniform_cells.insert (parent);
				}
			}

			return vector<hp::DoFHandler<dim>::active_cell_iterator> (uniform_cells.begin(),
				uniform_cells.end());
		}// else
	}	//get_cells_at_coarsest_common_level
	*/
	/*.....................................................................................................................................*/
/*
	template <int dim>
	void StokesProblem <dim>::build_triangulation_from_patch (const vector<hp::DoFHandler<dim>::active_cell_iterator>  &patch,
		Triangulation<dim> &tria_patch, unsigned int &level_h_refine, unsigned int &level_p_refine )
	{
		const vector<hp::DoFHandler<dim>::active_cell_iterator> uniform_cells = get_cells_at_coarsest_common_level (patch);

		level_h_refine=patch[0]->level();  //...............................................................................................?
		level_p_refine=patch[0]->active_fe_index();//............................................................................................?

		tria_patch.clear();
		vector<Point<dim>> vertices;
		const unsigned int n_uniform_cells=uniform_cells.size();
		vector<CellData<dim> > cells(n_uniform_cells);
		unsigned int k=0;// for enumerating cells
		unsigned int i=0;// for enumerating vertices
		typename hp::DoFHandler<dim>::active_cell_iterator  uniform_c=uniform_cells.begin();
		const typename hp::DoFHandler<dim>::active_cell_iterator end_uniform_c=uniform_cells.end();

		for (; uniform_c!=end_uniform_c; ++uniform_c)
		{
			bool repeat_vertex;
			for (unsigned int j=0;  j<GeometryInfo<2>::vertices_per_cell; ++j)
			{
				Point<dim>& position=uniform_c->vertex (j);//.......................?
				repeat_vertex=false;
				for (unsigned int m=0; m<i; ++m)
				{
					//double difference= std::sqrt(std::pow((vertices[m][0]-position[0]),2)+std::pow((vertices[m][1]-position[1]),2));
					//if (difference < 1e-10)
					//{
						 if (*position == vertices[m]){ //...............?
						repeat_vertex=true;
						cells[k].vertices[j]=m ;
						break;//.............?
					}//if
				}//for  m

				if (repeat_vertex==false)
				{
					vertices.push_back(position);
					cells[k].vertices[j]=i;
					i=i+1;
				}
			}//for vertices_per_cell
			k=k+1;
		}//uniform_c
		tria_patch.create_triangulation(vertices,cells,SubCellData());

		Assert (tria_patch.n_active_cells() == uniform_cells.size(), ExcInternalError());

		std::map<typename Triangulation<dim>::cell_iterator, typename DoFHandler<dim>::cell_iterator> patch_to_global_tria_map;
		unsigned int index=0;
		for (typename Triangulation<dim>::cell_iterator cell = tria_patch.begin(); cell != tria_patch.end(); ++cell, ++index)
			patch_to_global_tria_map.insert (std::make_pair(cell, uniform_cells[index]));

		bool refinement_necessary;
		do
		{
			refinement_necessary = false;

			for (typename Triangulation<dim>::active_cell_iterator cell = tria_patch.begin_active(); cell != tria_patch.end(); ++cell)
				if (patch_to_global_tria_map[cell]->has_children())
				{
					cell->set_refine_flag();
					refinement_necessary = true;
				}

			if (refinement_necessary)
			{
				tria_patch.execute_coarsening_and_refinement ();

				for (typename Triangulation<dim>::cell_iterator cell = tria_patch.begin(); cell != tria_patch.end(); ++cell)
					if (cell->has_children())
						// these children may not yet be in the map
						for (unsigned int c=0; c<cell->n_children(); ++c)
							if (patch_to_global_tria_map.find(cell->child(c)) == patch_to_global_tria_map.end())
								patch_to_global_tria_map.insert (std::make_pair(cell->child(c),
											patch_to_global_tria_map[cell]->child(c)));
			}
		}
		while (refinement_necessary);

		// mark the cells in 'tria' that exist in the patch: go
		// through all cells in 'tria' and see whether the
		// corresponding cell in the global triangulation is part of
		// the 'patch' list of cells
		for (typename Triangulation<dim>::cell_iterator cell = tria_patch.begin(); cell != tria_patch.end(); ++cell)
		{
			typename DoFHandler<dim>::cell_iterator global_cell = patch_to_global_tria_map[cell];
			bool global_cell_is_in_patch = false;
			for (vector<hp::DoFHandler<dim>::active_cell_iterator>::const_iterator p=patch.begin();
					p != patch.end(); ++p)
				if (*p == global_cell)
				{
					global_cell_is_in_patch = true;
					break;
				}

			if (global_cell_is_in_patch)
				cell->set_material_id (patch_on);
			else
				cell->set_material_id (patch_off);
		}

	} // build_triangulation
*/
	/*..........................    Compute h_convergence_estimator   &   h_workload_number  for each patch around cell   ..........................*/

/*	template <int dim>
	void StokesProblem <dim>::h_patch_conv_load_no (double &h_convergence_est_per_cell, unsigned int &h_workload_num, typename hp::DoFHandler<dim>::active_cell_iterator const &cell) const
	{
		Triangulation<dim> tria_patch;//
		hp::DoFHandler<dim> dof_handler_patch(tria_patch);//
		h_convergence_est_per_cell=0.;
		double h_solu_norm_per_patch=0.;
		vector<hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell(cell);
		//PRINT(tria_patch.n_cells());
		build_triangulation_from_patch (patch, tria_patch, level_h_refine, level_p_refine);
		//PRINT(tria_patch.n_cells());
		//GridOut::write_gnuplot (tria_patch, out);
		//................................................set_fe_nothing........................................................................//
		//set_active_fe_indices...
		typename hp::DoFHandler<dim>::active_cell_iterator
			cell = dof_handler_patch.begin_active(),//.....................................................tria.begin_active().....................................???
			endc = dof_handler_patch.end();
		for (; cell!=endc; ++cell){

			if (cell_is_on_patch(cell))
				cell->set_active_fe_index (0);//............?
			else if (cell_is_NOT_on_patch(cell))
				cell->set_active_fe_index (max_degree);//.............?
			else
				Assert (false, ExcNotImplemented());
		}//cell
		//...........................................................  h_refinement of patch cells ..................................................................//
		bool need_to_refine = false;
		do
		{
			need_to_refine = false;
			for (typename hp::DoFHandler<dim>::active_cell_iterator
				cell = dof_handler_patch.begin_active();
				cell != dof_handler_patch.end(); ++cell)
				for (; cell!=endc; ++cell){
					if (cell_is_on_patch(cell)){  //.....................................................?
						if ( (cell_is_on_patch(cell)) && (cell->level() <  level_h_refine+1 ))  {
							need_to_refine = true;
							cell->set_refine_flag();
						}
					}
					}//cell
					if (need_to_refine == true)
						tria_patch.execute_coarsening_and_refinement ();
				}//do
				while (need_to_refine == true);
				//..................................................  setup_h_patch_system and  patch_rhs  .......................... //
				dof_handler_patch.distribute_dofs (fe_collection);// fe_collection
				unsigned int h_workload_num = dof_handler_patch. n_dofs();
				unsigned int local_system_size = dof_handler_patch. n_dofs();

				vector<unsigned int> block_component (dim+1, 0);
				block_component[dim]=1;
				DoFRenumbering::component_wise(dof_handler_patch, block_component);
				{
					constraints.clear ();
					FEValuesExtractors::Vector velocities(0);
					DoFTools::make_hanging_node_constraints <hp::DoFHandler<dim>> (dof_handler_patch, constraints);

				}
				//......................... Zero_Bdry_Condition_on_Patch .......................................//
				{
					std::vector<types::global_dof_index> local_face_dof_indices (fe_collection.dofs_per_face);//stokes_fe
					for (typename hp::DoFHandler<dim>::active_cell_iterator
						cell = dof_handler_patch.begin_active();
						cell != dof_handler_patch.end(); ++cell)
						if (cell_is_on_patch (cell)){
							for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
								if (!cell->at_boundary(f))
								{
									bool face_is_on_patch_Bdry = false;
									if ((cell->neighbor(f)->has_children() == false)
										&&
										(cell_is_NOT_on_patch (cell->neighbor(f))))
										face_is_on_patch_Bdry = true;
									else if (cell->neighbor(f)->has_children() == true)
									{
										for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf)
											if (cell_is_NOT_on_patch (cell->neighbor_child_on_subface
												(f, sf)))
											{
												face_is_on_patch_Bdry = true;
												break;
											}
									}
									if (face_is_on_patch_Bdry)
									{

										cell->face(f)->get_dof_indices (local_face_dof_indices, 0);
										for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
											if (fe_collection.face_system_to_component_index(i).first < dim)//stokes_fe
												constraints.add_line (local_face_dof_indices[i]);
									}
								}
						}// if (cell_is_on_patch (cell))
				}
				//

				constraints.close();

				{
					BlockCompressedSetSparsityPattern csp (local_system_size,local_system_size);

					DoFTools::make_sparsity_pattern (dof_handler_patch, csp, constraints, false);
					sparsity_pattern.copy_from(csp);
				}

				BlockSparseMatrix<double> patch_system (sparsity_pattern);
				BlockVector<double> patch_solution (local_system_size);
				BlockVector<double> patch_rhs (local_system_size);

				// ..................................................  assemble  patch_system  and patch_rhs .............................. //

				hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);

				FullMatrix<double> local_matrix;
				Vector<double> local_rhs;
				Vector<double> local_rhs1;
				Vector<double> local_rhs2;
				vector<types::global_dof_index> local_dof_indices;

				vector<Vector<double> >  rhs_values;

				const FEValuesExtractors::Vector velocities (0);
				const FEValuesExtractors::Scalar pressure (dim);

				vector<Tensor<2,dim> > grad_phi_u;
				vector<double> div_phi_u;
				vector<Tensor<1,dim> > phi_u;
				vector<double> phi_p;

				vector<Tensor<1,dim>> gradients_p;
				vector<double> divergences;
				vector<Tensor<1,dim>> laplacians;

				vector<double> values;
				vector<Tensor<2,dim>> gradients;


				typename hp::DoFHandler<dim>::active_cell_iterator
					cell = dof_handler_patch.begin_active(),//...................................................which
					endc = dof_handler_patch.end();
				for (; cell!=endc; ++cell){

					const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
					local_matrix.reinit (dofs_per_cell, dofs_per_cell);
					local_rhs.reinit (dofs_per_cell);
					local_rhs1.reinit (dofs_per_cell);
					local_rhs2.reinit (dofs_per_cell);

					hp_fe_values.reinit (cell);
					const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
					const vector<double>& JxW_values = fe_values.get_JxW_values ();
					const unsigned int n_q_points = fe_values.n_quadrature_points;

					rhs_values.resize(n_q_points, Vector<double>(dim+1));
					rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

					grad_phi_u.resize(dofs_per_cell);
					div_phi_u.resize(dofs_per_cell);
					phi_u.resize (dofs_per_cell);
					phi_p.resize(dofs_per_cell);

					divergences.resize(n_q_points);
					gradients_p.resize(n_q_points) ;
					laplacians.resize(n_q_points);

					fe_values[pressure].get_function_gradients(solution, gradients_p);
					fe_values[velocities].get_function_divergences(solution, divergences);
					fe_values[velocities].get_function_laplacians(solution, laplacians);

					gradients.resize(n_q_points);
					values.resize(n_q_points);


					for (unsigned int q=0; q<n_q_points; ++q)
					{
						for (unsigned int k=0; k<dofs_per_cell; ++k)
						{
							grad_phi_u[k] = fe_values[velocities].gradient (k, q);
							div_phi_u[k] = fe_values[velocities].divergence (k, q);
							phi_u[k] = fe_values[velocities].value (k, q);
							phi_p[k] = fe_values[pressure].value (k, q);
						}
						for (unsigned int i=0; i<dofs_per_cell; ++i)
						{
							for (unsigned int j=0; j<dofs_per_cell; ++j)
							{
								local_matrix(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
									- div_phi_u[i] * phi_p[j]
								- phi_p[i] * div_phi_u[j]
								+ phi_p[i] * phi_p[j])
									* JxW_values[q];
							} // end of loop over 'j'

							local_rhs1(i)+= ((rhs_values[q](0)+laplacians[q][0]-gradients_p[q][0])*(phi_u[i][0])+ (rhs_values[q](1)+laplacians[q][1]-gradients_p[q][1])*(phi_u[i][1]))* JxW_values[q];
							local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];
							local_rhs(i)= local_rhs1(i)-local_rhs2(i);
						}// i
					}//q
					local_dof_indices.resize (dofs_per_cell);
					cell->get_dof_indices (local_dof_indices);
					constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, patch_system, patch_rhs);
				}// for cell

				// .....................................solve patch_system and patch_rhs ............get  patch_solution ...................... //

					SparseDirectUMFPACK A_inverse;
					A_inverse.initialize (patch_system.block(0,0), SparseDirectUMFPACK::AdditionalData());
					Vector<double> tmp (patch_solution.block(0).size());
					{
						Vector<double> schur_rhs (patch_solution.block(1).size());
						A_inverse.vmult (tmp, patch_rhs.block(0));
						patch_system.block(1,0).vmult (schur_rhs, tmp);
						schur_rhs -= patch_rhs.block(1);

						SchurComplement schur_complement (patch_system, A_inverse);
						SolverControl solver_control (patch_solution.block(1).size(), 1e-6);
						SolverCG<>    cg (solver_control);
						SparseDirectUMFPACK preconditioner;
						preconditioner.initialize (patch_system.block(1,1),
							SparseDirectUMFPACK::AdditionalData());

						cg.solve (schur_complement, patch_solution.block(1), schur_rhs,
							preconditioner);
						constraints.distribute (patch_solution);
						cout << "  "
							<< solver_control.last_step()
							<< " outer CG Schur complement iterations for pressure"
							<< endl;
					}

					{
						patch_system.block(0,1).vmult (tmp, patch_solution.block(1));
						tmp *= -1.0;
						tmp += patch_rhs.block(0);
						A_inverse.vmult (patch_solution.block(0), tmp);
						constraints.distribute (patch_solution);
					}

					//.......................................  get the L2 norm of the gradient of velocity solution and pressure value  .....................//
					fe_values[velocities].get_function_gradients(patch_solution, gradients);
					fe_values[pressure].get_function_values(patch_solution, values);

					double pressure_val=0;
					double grad_u_val=0;
					typename hp::DoFHandler<dim>::active_cell_iterator
						cell = dof_handler_patch.begin_active(),
						endc = dof_handler_patch.end();
					for (; cell!=endc; ++cell)
					{
					hp_fe_values.reinit (cell);
					const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
					const vector<double>& JxW_values = fe_values.get_JxW_values ();
					const unsigned int n_q_points = fe_values.n_quadrature_points;
						for (unsigned int q=0; q<n_q_points; ++q)
						{
							pressure_val +=values[q]*values[q]* JxW_values[q];

							for (unsigned int i=0; i<2; ++i)

								grad_u_val += contract(gradients[q][i],gradients[q][i])* JxW_values[q];
						} // q
						h_solu_norm_per_patch +=(sqrt(pressure_val) + sqrt(grad_u_val));//?
					}// cells on patch


				h_convergence_est_per_cell = h_solu_norm_per_patch ;

		}// function h_patch_con_loadnum
*/
		/*..........................    Compute p_convergence_estimator   &   p_workload_number  for each patch around cell   ..........................*/
/*
	template <int dim>
		void StokesProblem <dim>::p_patch_conv_load_no (double &p_convergence_est_per_cell, unsigned int &p_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell ) const
		{
			Triangulation<dim> tria_patch;//
			hp::DoFHandler<dim> dof_handler_patch;//

			p_convergence_est_per_cell=0.;
			p_solu_norm_per_patch=0.;

			vector<hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell (cell);
			build_triangulation_from_patch (patch, tria_patch,level_h_refine, level_p_refine);
			dof_handler_patch(tria_patch);//
			//................................................set_fe_nothing........................................................................//
			typename DoFHandler<dim>::active_cell_iterator
				cell = dof_handler_patch.begin_active(),
				endc = dof_handler_patch.end();
			for (; cell!=endc; ++cell){

				if (cell_is_on_patch(cell))
					cell->set_active_fe_index (0);
				else if (cell_is_NOT_on_patch(cell))
					cell->set_active_fe_index (max_degree);
				else
					Assert (false, ExcNotImplemented());
			}//cell
			//......................................................................................................................................//
			typename DoFHandler<dim>::active_cell_iterator
				cell = dof_handler_patch.begin_active(),
				endc = dof_handler_patch.end();
			for (; cell!=endc; ++cell){
				if (cell_is_on_patch(cell)){
					if (cell->active_fe_index()+1 < fe_collection.size() && cell->active_fe_index() < level_p_refine+1)
						cell->set_active_fe_index (cell->active_fe_index() + 1);
				}
			}//cell
			//..................................................  setup   p_patch_system and  patch_rhs  .......................... //
			dof_handler_patch.distribute_dofs (fe_collection);
			unsigned int h_workload_num = dof_handler_patch. n_dofs();//
			unsigned int local_system_size = dof_handler_patch . n_dofs();//

			vector<unsigned int> block_component (dim+1, 0);
			block_component[dim]=1;
			DoFRenumbering::component_wise(dof_handler_patch, block_component);

			{
				constraints.clear ();
				FEValuesExtractors::Vector velocities(0);
				DoFTools::make_hanging_node_constraints <hp::DoFHandler<dim>> (dof_handler_patch, constraints);

			}

			//...................... Zero_Bdry_Condition_on_Patch ................................................//
			{
				std::vector<types::global_dof_index> local_face_dof_indices (stokes_fe.dofs_per_face);
				for (typename hp::DoFHandler<dim>::active_cell_iterator
					cell = dof_handler_patch.begin_active();
					cell != dof_handler_patch.end(); ++cell)
					if (cell_is_on_patch (cell)){
						for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
							if (!cell->at_boundary(f))
							{
								bool face_is_on_patch_Bdry = false;
								if ((cell->neighbor(f)->has_children() == false)
									&&
									(cell_is_NOT_on_patch (cell->neighbor(f))))
									face_is_on_patch_Bdry = true;
								else if (cell->neighbor(f)->has_children() == true)
								{
									for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf)
										if (cell_is_NOT_on_patch (cell->neighbor_child_on_subface
											(f, sf)))
										{
											face_is_on_patch_Bdry = true;
											break;
										}
								}
								if (face_is_on_patch_Bdry)
								{
									cell->face(f)->get_dof_indices (local_face_dof_indices, 0);
									for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
										if (fe_collection.face_system_to_component_index(i).first < 2)
											constraints.add_line (local_face_dof_indices[i]);
								}
							}
					}//  if (cell_is_on_patch (cell))
			}
			//

			constraints.close();

			{
				BlockCompressedSetSparsityPattern csp (local_system_size,local_system_size);

				DoFTools::make_sparsity_pattern (dof_handler_patch, csp, constraints, false);
				sparsity_pattern.copy_from(csp);
			}

			BlockSparseMatrix<double> patch_system (sparsity_pattern);
			BlockVector<double> patch_solution (local_system_size);
			BlockVector<double> patch_rhs (local_system_size);

			// ..................................................  assemble  patch_system  and patch_rhs .............................. //
			hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);

			FullMatrix<double> local_matrix;
			Vector<double> local_rhs;
			Vector<double> local_rhs1;
			Vector<double> local_rhs2;
			vector<types::global_dof_index> local_dof_indices;

			vector<Vector<double> >  rhs_values;

			const FEValuesExtractors::Vector velocities (0);
			const FEValuesExtractors::Scalar pressure (dim);

			vector<Tensor<2,dim> > grad_phi_u;
			vector<double> div_phi_u;
			vector<Tensor<1,dim> > phi_u;
			vector<double> phi_p;

			vector<Tensor<1,dim>> gradients_p;
			vector<double> divergences;
			vector<Tensor<1,dim>> laplacians;

			vector<double> values;
			vector<Tensor<2,dim>> gradients;

			typename hp::DoFHandler<dim>::active_cell_iterator
				cell = dof_handler_patch.begin_active(),
				endc = dof_handler_patch.end();
			for (; cell!=endc; ++cell)
			{

				const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
				local_matrix.reinit (dofs_per_cell, dofs_per_cell);
				local_rhs.reinit (dofs_per_cell);
				local_rhs1.reinit (dofs_per_cell);
				local_rhs2.reinit (dofs_per_cell);

				hp_fe_values.reinit (cell);
				const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
				const vector<double>& JxW_values = fe_values.get_JxW_values ();
				const unsigned int n_q_points = fe_values.n_quadrature_points;

				rhs_values.resize(n_q_points, Vector<double>(dim+1));
				rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

				grad_phi_u.resize(dofs_per_cell);
				div_phi_u.resize(dofs_per_cell);
				phi_u.resize (dofs_per_cell);
				phi_p.resize(dofs_per_cell);

				divergences.resize(n_q_points);
				gradients_p.resize(n_q_points) ;
				laplacians.resize(n_q_points);

				fe_values[pressure].get_function_gradients(solution, gradients_p);
				fe_values[velocities].get_function_divergences(solution, divergences);
				fe_values[velocities].get_function_laplacians(solution, laplacians);

				gradients.resize(n_q_points);
				values.resize(n_q_points);


				for (unsigned int q=0; q<n_q_points; ++q)
				{
					for (unsigned int k=0; k<dofs_per_cell; ++k)
					{
						grad_phi_u[k] = fe_values[velocities].gradient (k, q);
						div_phi_u[k] = fe_values[velocities].divergence (k, q);
						phi_u[k] = fe_values[velocities].value (k, q);
						phi_p[k] = fe_values[pressure].value (k, q);
					}
					for (unsigned int i=0; i<dofs_per_cell; ++i)
					{
						for (unsigned int j=0; j<dofs_per_cell; ++j)
						{
							local_matrix(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
								- div_phi_u[i] * phi_p[j]
							- phi_p[i] * div_phi_u[j]
							+ phi_p[i] * phi_p[j])
								* JxW_values[q];
						} // end of loop over 'j'

						local_rhs1(i)+= ((rhs_values[q](0)+laplacians[q][0]-gradients_p[q][0])*(phi_u[i][0])+ (rhs_values[q](1)+laplacians[q][1]-gradients_p[q][1])*(phi_u[i][1]))* JxW_values[q];//?
						local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];//?
						local_rhs(i)= local_rhs1(i)-local_rhs2(i);//?
					}// i
				}//q


				local_dof_indices.resize (dofs_per_cell);
				cell->get_dof_indices (local_dof_indices);
				constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, patch_system, patch_rhs);
			}// for patch_c

			// .....................................solve patch_system and patch_rhs ............get  patch_solution ...................... //
			{
				SparseDirectUMFPACK A_inverse;
				A_inverse.initialize (patch_system.block(0,0), SparseDirectUMFPACK::AdditionalData());
				Vector<double> tmp (patch_solution.block(0).size());
				{
					Vector<double> schur_rhs (patch_solution.block(1).size());
					A_inverse.vmult (tmp, patch_rhs.block(0));
					patch_system.block(1,0).vmult (schur_rhs, tmp);
					schur_rhs -= patch_rhs.block(1);

					SchurComplement schur_complement (patch_system, A_inverse);
					SolverControl solver_control (patch_solution.block(1).size(), 1e-6);
					SolverCG<>    cg (solver_control);
					SparseDirectUMFPACK preconditioner;
					preconditioner.initialize (patch_system.block(1,1),
						SparseDirectUMFPACK::AdditionalData());

					cg.solve (schur_complement, patch_solution.block(1), schur_rhs,
						preconditioner);
					constraints.distribute (patch_solution);
					cout << "  "
						<< solver_control.last_step()
						<< " outer CG Schur complement iterations for pressure"
						<< endl;
				}

				{
					patch_system.block(0,1).vmult (tmp, patch_solution.block(1));
					tmp *= -1.0;
					tmp += patch_rhs.block(0);
					A_inverse.vmult (patch_solution.block(0), tmp);
					constraints.distribute (patch_solution);
				}

				//.........................  get the L2 norm of the gradient of velocity solution and pressure value  .....................//

				fe_values[velocities].get_function_gradients(patch_solution, gradients);
				fe_values[pressure].get_function_values(patch_solution, values);

				double pressure_val=0;
				double grad_u_val=0;

				typename hp::DoFHandler<2>::active_cell_iterator
					cell = dof_handler_patch.begin_active(),
					endc = dof_handler_patch.end();
				for (; cell!=endc; ++cell)
				{
					hp_fe_values.reinit (cell);
					const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
					const vector<double>& JxW_values = fe_values.get_JxW_values ();
					const unsigned int n_q_points = fe_values.n_quadrature_points;

					for (unsigned int q=0; q<n_q_points; ++q)
					{
						pressure_val +=values[q]*values[q]* JxW_values[q];

						for (unsigned int i=0; i<2; ++i)

							grad_u_val += contract(gradients[q][i],gradients[q][i])* JxW_values[q];
					} // q

				}// cells on patch
				p_solu_norm_per_patch +=(sqrt(pressure_val) + sqrt(grad_u_val));  //?
			}// solve//......................................................................?
			p_convergence_est_per_cell = p_solu_norm_per_patch ;

		} // function p_patch_conv_load_no
		*/
		/*..................................................................POSTPROCESS()....................................................*/
	/*
	template <int dim>
	void StokesProblem <dim>:: postprocess (const unsigned int cycle){

			const double theta= 0.5;

			vector<pair<double, hp::DoFHandler<dim>::active_cell_iterator> > to_be_sorted;

			Vector<double> est_per_cell (triangulation.n_active_cells());
			estimate(est_per_cell);

//............................................................................................................................
{
			Vector<float> fe_degrees (triangulation.n_active_cells());
  	    {
  			      typename hp::DoFHandler<dim>::active_cell_iterator
  			      cell = dof_handler.begin_active(),
   				     endc = dof_handler.end();
				        for (unsigned int index=0; cell!=endc; ++cell, ++index)
  				        fe_degrees(index)= fe_collection[cell->active_fe_index()].degree;
            }


      DataOut<dim,hp::DoFHandler<dim> > data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.add_data_vector (est_per_cell, "error");
      data_out.add_data_vector (fe_degrees, "fe_degree");
      data_out.build_patches ();
      const std::string filename = "solution-" +
                                   Utilities::int_to_string (cycle, 2) +
                                   ".vtk";
      std::ofstream output (filename.c_str());
      data_out.write_vtk (output);
    }

//..............................................................................................................................
			Vector<double> convergence_est_per_cell (triangulation.n_active_cells());
			bool need_to_p_refine = false;
			bool need_to_h_refine = false;
			unsigned int cell_index=0;
			typename hp::DoFHandler<dim>::active_cell_iterator
				cell = dof_handler.begin_active(),
				endc = dof_handler.end();
			for (; cell!=endc; ++cell,++cell_index)
			{
				bool need_to_p_refine = false;
				bool need_to_h_refine = false;
				double indicator_per_cell =0.;

				double h_convergence_est_per_cell;
				double h_workload_num;
				h_patch_conv_load_no (h_convergence_est_per_cell,h_workload_num, cell);
				h_convergence_est_per_cell = h_convergence_est_per_cell /est_per_cell(cell_index);

				double p_convergence_est_per_cell;
				double p_workload_num;
				p_patch_conv_load_no (p_convergence_est_per_cell,p_workload_num, cell);
				p_convergence_est_per_cell = p_convergence_est_per_cell /est_per_cell(cell_index);

				double h_ratio= h_convergence_est_per_cell /  h_workload_num ;
				double p_ratio= p_convergence_est_per_cell /  p_workload_num ;

				if (h_ratio < p_ratio) {
					convergence_est_per_cell(cell_index)=p_convergence_est_per_cell;
					indicator_per_cell=convergence_est_per_cell(cell_index)*est_per_cell(cell_index);
					need_to_p_refine = true;

				}
				else{
					convergence_est_per_cell(cell_index)=h_convergence_est_per_cell;
					indicator_per_cell=convergence_est_per_cell(cell_index)*est_per_cell(cell_index);
					need_to_h_refine = true;

				}

				to_be_sorted.push_back(pair<double, hp::DoFHandler<dim>::active_cell_iterator> > (indicator_per_cell,cell));//.............?
			}// cell
			sort (to_be_sorted.first.begin(), to_be_sorted.first.end(), decreasing);


			vector<hp::DoFHandler<dim>::active_cell_iterator> candidate_cell_set;
			double L2_norm=est_per_cell.l2_norm();
			typename hp::DoFHandler<dim>::active_cell_iterator  cell_sort=to_be_sorted.second.begin(), endc_sort=to_be_sorted.second.end(); //...............?
			double sum=0;
			unsigned int index_cell=0;
			for (; cell_sort!=endc_sort; ++cell_sort,++index_cell)
			{
				sum +=to_be_sorted.first(index_cell)*to_be_sorted.first(index_cell);
				if (sum < (theta*(L2_norm))*(theta*(L2_norm)))
					candidate_cell_set.push_back (cell_sort);
			}// cell_sort

			//......................................  Output results ....................................................//
				{
			std::vector<std::string> solution_names (dim, "velocity");
			solution_names.push_back ("pressure");

			std::vector<DataComponentInterpretation::DataComponentInterpretation>
			data_component_interpretation
			(dim, DataComponentInterpretation::component_is_part_of_vector);
			data_component_interpretation
			.push_back (DataComponentInterpretation::component_is_scalar);

			DataOut<2,hp::DoFHandler<2> > data_out;
			data_out.attach_dof_handler (dof_handler);
			data_out.add_data_vector (solution, solution_names,
			DataOut<2,hp::DoFHandler<2> >::type_dof_data,
			data_component_interpretation);

			data_out.add_data_vector (est_per_cell, "error_estimator");
			data_out.add_data_vector (est_per_cell, "error_estimator");

			data_out.build_patches ();

			std::ostringstream filename;
			filename << "solution-"
			<< Utilities::int_to_string (cycle, 2)
			<< ".vtk";

			std::ofstream output (filename.str().c_str());
			data_out.write_vtk (output);
			}

			//............................................................................................................//


			typename hp::DoFHandler<dim>::active_cell_iterator  cell_candidate=candidate_cell_set.begin(), endc_candidate=candidate_cell_set.end();
			for (; cell_candidate!=endc_candidate; ++ cell_candidate)
			{
				Triangulation<dim> tria_patch;//
				hp::DoFHandler<dim> dof_handler_patch;//
				dof_handler_patch(tria_patch);//

				if (need_to_h_refine == true){

					vector<hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell (cell_candidate);
					build_triangulation_from_patch (patch, tria_patch,level_h_refine, level_p_refine);
					//................................................set_fe_nothing...........................................//
					typename DoFHandler<dim>::active_cell_iterator

						cell_p = dof_handler_patch.begin_active(),
						endc_p = dof_handler_patch.end();
					for (; cell_p!=endc_p; ++cell_p){

						if (cell_is_on_patch(cell_p))
							cell_p->set_active_fe_index (0);
						else if (cell_is_NOT_on_patch(cell_p))
							cell_p->set_active_fe_index (max_degree);
						else
							Assert (false, ExcNotImplemented());
					}//cell_p
					//.........................................................................................................//
					bool need_to_refine;
					do
					{
						need_to_refine = false;
						typename hp::DoFHandler<2>::active_cell_iterator
							cell_p = dof_handler_patch.begin_active(),
							endc_p = dof_handler_patch.end();
						for (; cell_p!=endc_p; ++cell_p){

							if ( (cell_is_on_patch(cell_p)) &&  (cell_p->level() <  level_h_refine+1) )  {
								need_to_refine = true;
								cell_p->set_refine_flag(); ///..........................?
							
						}
					}//cell_p
					if (need_to_refine == true)
						tria_patch.execute_coarsening_and_refinement ();
				}//do
				while (need_to_refine == true);

			}//if    need_to_h_refine
			//......................................................................................................................................................//
			if(need_to_p_refine == true){

				vector<hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell (cell_candidate);
				build_triangulation_from_patch (patch, tria_patch,level_h_refine, level_p_refine);

				//................................................set_fe_nothing...........................................//
				typename DoFHandler<dim>::active_cell_iterator
					cell_p = dof_handler_patch.begin_active(),
					endc_p = dof_handler_patch.end();
				for (; cell_p!=endc_p; ++cell_p){

					if (cell_is_on_patch(cell_p))
						cell_p->set_active_fe_index (0);
					else if (cell_is_NOT_on_patch(cell_p))
						cell_p->set_active_fe_index (max_degree);
					else
						Assert (false, ExcNotImplemented());
				}//cell
				//.........................................................................................................//
				typename DoFHandler<dim>::active_cell_iterator
					cell_p = dof_handler_patch.begin_active(),
					endc_p = dof_handler_patch.end();
				for (; cell_p!=endc_p; ++cell_p){
					if (cell_is_on_patch(cell_p)){
						if (cell_p->active_fe_index()+1 < fe_collection.size() && cell_p->active_fe_index() < level_p_refine+1)
							cell_P->set_active_fe_index (cell_p->active_fe_index() + 1);
					}//if
				}//cell_p
			}// if  need_to_p_refine
			//......................................................................................................................................................//
		}// for... candidate_cells





	}// postprocess
*/
	/*......................................................................................................................................................*/
template <int dim>
	void StokesProblem <dim>::run(){

		//for (unsigned int cycle=0; cycle<6; ++cycle)
		//{
			//cout << "Cycle " << cycle << ':' << endl;
			//if (cycle == 0)

				generate_mesh();

			setup_system ();
			assemble_system();
			solve ();
			compute_error ();
                        Vector<double> est_per_cell (triangulation.n_active_cells());
			estimate(est_per_cell);
			//L2_norm_est= est_per_cell.l2_norm();
			//cout<< "L2_norm Estimate is: "<< L2_norm_est << endl;

			//if (L2_norm_est < Tolerance) break;

			//postprocess(cycle);
		//}
	}//run
}//  namespace hp-Stokes
/*......................................................................................*/
int main ()
{

	try
	{
		using namespace dealii;
		using namespace hp_Stokes;
		deallog.depth_console (0);
		StokesProblem<2> stokesproblem;
		stokesproblem.run();
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

