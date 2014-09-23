#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <set>

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

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
/*----------------------------------------------------------------------------------------*/
using namespace std;
using namespace dealii;

class ExactSolution : public Function<2>
{
public:
	ExactSolution () : Function<2>(2+1) {}
	virtual void vector_value (const Point<2> &p,
		Vector<double>   &value) const;
	virtual void vector_gradient (const Point<2> &p,
		vector< Tensor< 1, 2 > > &gradients) const;
};

void ExactSolution::vector_value (const Point<2> &p,
	Vector<double>   &values) const
{
	values(0) = -1.0 * exp (p(0)) * (p(1) * cos (p(1)) + sin (p(1)));
	values(1) = exp(p(0)) * p(1) * sin (p(1));
	values(2) = 2.0 * exp (p(0)) * sin (p(1)) - (2.0 * (1.0 - exp(1.0)) * (cos (1.0) - 1.0)) / 3.0;
}

void ExactSolution::vector_gradient  (const Point<2> &p,
	vector< Tensor< 1, 2 > > &gradients) const
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
class RightHandSide : public Function<2>
{
public:
	RightHandSide () : Function<2>(2+1) {}
	virtual double value (const Point<2> &p, const unsigned int component = 0) const;
	virtual void vector_value (const Point<2> &p, Vector<double> &value) const;
};
double RightHandSide::value (const Point<2> & p, const unsigned int  component) const
{
	return 0;
}
void
	RightHandSide::vector_value (const Point<2> &p, Vector<double> &values) const
{
	for (unsigned int c=0; c<this->n_components; ++c)
		values(c) = RightHandSide::value (p, c);
}
/*......................................................................................*/
class StokesProblem
{
public:  
	StokesProblem ();
	~StokesProblem(); 
	void run();

private:
	RightHandSide rhs_function;
	ExactSolution exact_solution;
	double pressure_mean_value () const;
	void generate_mesh ();
	void setup_system ();
	void assemble_system ();
	void solve ();
	void compute_error ();
	void estimate (Vector<double> &est_per_cell) const;

	void convergence_estimate ();
	void postprocess (const unsigned int cycle);

	const double Tolerance;
	const unsigned int max_degree;

	Triangulation<2> triangulation;
	hp::DoFHandler<2> dof_handler;
	hp::FECollection<2> fe_collection; 
	hp::QCollection<2> quadrature_collection;
	hp::QCollection<2-1> face_quadrature_collection;

	ConstraintMatrix constraints;
	BlockSparsityPattern sparsity_pattern;
	BlockSparseMatrix<double> system_matrix;

	BlockVector<double> solution;
	BlockVector<double> system_rhs;
	double L1_norm_est;
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
StokesProblem::StokesProblem(): max_degree (5), dof_handler(triangulation), Tolerance (0.00001){
	for (unsigned int degree=1; degree<=max_degree; ++degree)
	{
		fe_collection.push_back(FESystem<2> (FE_Q<2> (QGaussLobatto<1> (degree + 2)), 2, FE_Q<2> (QGaussLobatto<1> (degree + 1)), 1));
		quadrature_collection.push_back(QGauss<2> (degree+2));
		face_quadrature_collection.push_back (QGauss<2-1> (degree+1));
	}
} 
/*.....................................................................................*/
StokesProblem::~StokesProblem() {
	dof_handler.clear();
}
/*......................................................................................*/
// Generate mesh
void StokesProblem::generate_mesh(){
	vector<Point<2> > vertices (8);

	vertices [0]=Point<2> (-1,-1);
	vertices [1]=Point<2> (0,-1);
	vertices [2]=Point<2> (-1,0);
	vertices [3]=Point<2> (0,0);
	vertices [4]=Point<2> (1,0);
	vertices [5]=Point<2> (-1,1);
	vertices [6]=Point<2>(0,1);
	vertices [7]=Point<2>(1,1);

	const unsigned int n_cells=3;
	vector<CellData<2> > cell(n_cells);
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
	triangulation.refine_global (2);
	ofstream out ("grid-L-Shape.eps");
	GridOut grid_out;
	grid_out.write_eps (triangulation, out);
	cout<<"Number of active cells: "<< triangulation.n_active_cells() << endl;
	cout<<"Total number of cells: " << triangulation.n_cells() << endl ;
} 
/*......................................................................................*/
// setup system()
void StokesProblem::setup_system(){
	dof_handler.distribute_dofs (fe_collection);

	vector<unsigned int> block_component (2+1, 0);
	block_component[2]=1;
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
void StokesProblem::assemble_system () {
	hp::FEValues<2> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients);


	FullMatrix<double> local_matrix;
	Vector<double> local_rhs;
	vector<types::global_dof_index> local_dof_indices;

	vector<Vector<double> >  rhs_values;

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (2);

	vector<Tensor<2,2> > grad_phi_u;
	vector<double> div_phi_u;
	vector<Tensor<1,2> > phi_u;
	vector<double> phi_p;

	typename hp::DoFHandler<2>::active_cell_iterator
		cell = dof_handler.begin_active(),
		endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
		local_matrix.reinit (dofs_per_cell, dofs_per_cell);
		local_rhs.reinit (dofs_per_cell);

		hp_fe_values.reinit (cell);
		const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
		const vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;

		rhs_values.resize(n_q_points, Vector<double>(3));
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
				local_rhs(i) += (phi_u[i][0] * rhs_values[q](0) + phi_u[i][1] * rhs_values [q](1)) * JxW_values[q];
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
void StokesProblem::solve ()
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
			1e-6);
		SolverCG<>    cg (solver_control);
		SparseDirectUMFPACK preconditioner;
		preconditioner.initialize (system_matrix.block(1,1),
			SparseDirectUMFPACK::AdditionalData());

		cg.solve (schur_complement, solution.block(1), schur_rhs,
			preconditioner);
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
double StokesProblem::pressure_mean_value () const
{
	// get pressure such that satisfies mean value property:
	hp::FEValues<2> hp_fe_values (fe_collection, quadrature_collection, update_values|update_JxW_values);
	const FEValuesExtractors::Scalar pressure (2);

	vector<double> values;
	double domain_mean_val_p=0;
	// double measure_domain=0;
	typename hp::DoFHandler<2>::active_cell_iterator
		cell = dof_handler.begin_active(),
		endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		hp_fe_values.reinit (cell);

		const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
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
	return domain_mean_val_p/3;
	// return domain_mean_val_p / measure_domain;
}
/*.................................................................................................*/
void StokesProblem::compute_error () {
	hp::FEValues<2> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients);
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (2);

	vector<double> values;
	vector<Tensor<2,2>> gradients;
	vector<vector<Tensor<1,2>>> exact_solution_gradients;
	vector<Vector<double>> exact_solution_values;

	Vector<double> error_per_cell(triangulation.n_active_cells());

	typename hp::DoFHandler<2>::active_cell_iterator
		cell = dof_handler.begin_active(),
		endc = dof_handler.end();
	unsigned cell_index=0;
	for (; cell!=endc; ++cell,++cell_index)
	{
		hp_fe_values.reinit (cell);
		const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
		const vector<double>& JxW_values = fe_values.get_JxW_values ();
		const vector<Point<2>>& quadrature_points = fe_values.get_quadrature_points();
		const unsigned int n_q_points = fe_values.n_quadrature_points;

		gradients.resize(n_q_points);
		values.resize(n_q_points);
		exact_solution_gradients.resize(n_q_points , vector<Tensor<1,2>> (3));
		exact_solution_values.resize(n_q_points, Vector<double> (3));

		fe_values[velocities].get_function_gradients(solution, gradients);
		fe_values[pressure].get_function_values(solution, values);

		exact_solution.vector_gradient_list(quadrature_points, exact_solution_gradients);
		exact_solution.vector_value_list(quadrature_points, exact_solution_values);

		double subtract_p=0;
		double grad_u_vals=0;
		for (unsigned int q=0; q<n_q_points; ++q)
		{
			values[q] -= exact_solution_values[q](2);
			subtract_p +=values[q]*values[q]* JxW_values[q];
			for (unsigned int i=0; i<2; ++i)
				gradients[q][i]-=exact_solution_gradients[q][i];	
			grad_u_vals +=double_contract(gradients[q],gradients[q])* JxW_values[q];
		} // q
		error_per_cell(cell_index) =(sqrt(subtract_p) + sqrt(grad_u_vals));

	}// cell
	double L1_norm=error_per_cell.l1_norm();
	cout<< "L1_norm is: "<< L1_norm << endl;
}
/*......................................................................................*/
void StokesProblem::estimate (Vector<double> &est_per_cell) const {
	hp::FEValues<2> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);
	hp::FEFaceValues<2> hp_fe_face_values(fe_collection, face_quadrature_collection, update_JxW_values|update_gradients|update_normal_vectors);
	hp::FEFaceValues<2> hp_neighbor_face_values(fe_collection, face_quadrature_collection, update_gradients);
	hp::FESubfaceValues<2> hp_subface_values(fe_collection, face_quadrature_collection, update_JxW_values|update_gradients|update_normal_vectors);
	hp::FESubfaceValues<2> hp_neighbor_subface_values(fe_collection, face_quadrature_collection, update_gradients);

	vector<Tensor<1,2>> gradients_p;
	vector<double> divergences;
	vector<Tensor<1,2>> laplacians;

	vector<Tensor<2,2>> gradients;
	vector<Tensor<2,2>> neighbor_gradients;

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (2);

	const RightHandSide rhs_function;
	vector<Vector<double> >  rhs_values;

	Vector<double> res_est_per_cell(triangulation.n_active_cells());
	Vector<double> Jump_est_per_cell(triangulation.n_active_cells());


	typename hp::DoFHandler<2>::active_cell_iterator
		cell = dof_handler.begin_active(),
		endc = dof_handler.end();
	unsigned cell_index=0;
	for (; cell!=endc; ++cell,++cell_index)
	{
		hp_fe_values.reinit (cell);
		const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
		const vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;

		rhs_values.resize(n_q_points, Vector<double>(3));
		rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

		divergences.resize(n_q_points);
		gradients_p.resize(n_q_points) ;
		laplacians.resize(n_q_points);

		fe_values[pressure].get_function_gradients(solution, gradients_p);
		fe_values[velocities].get_function_divergences(solution, divergences);
		fe_values[velocities].get_function_laplacians(solution, laplacians);

		double term2=0;
		double term1=0;
		for (unsigned int q=0; q<n_q_points; ++q){
			term2 += (divergences[q])*(divergences[q])*JxW_values[q];

			for (unsigned int i=0; i<2; ++i)
				gradients_p[q][i]-= (rhs_values[q](i)+laplacians[q][i]);

			term1+= contract(gradients_p[q],gradients_p[q])*JxW_values[q];
		}// q
		res_est_per_cell(cell_index)=sqrt((pow((cell->diameter())/(cell->get_fe().degree), 2.0 )*(term1)) + term2);

		// ******************************************** compute jump_est_per_cell**********************************************************************
		double term3=0;
		for (unsigned int face_number=0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)

			if ((cell->face(face_number)->at_boundary()==false)	&& (cell->face(face_number)->has_children() == false) && (cell->face(face_number)->level() == cell->level()))
			{
				hp_fe_face_values.reinit (cell, face_number,cell->active_fe_index());
				hp_neighbor_face_values.reinit (cell->neighbor(face_number), cell->neighbor_of_neighbor(face_number)); //cell->active_fe_index()

				const FEFaceValues<2> &neighbor_face_values =hp_neighbor_face_values.get_present_fe_values ();
				const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values ();

				const vector<double>& JxW_values = fe_face_values.get_JxW_values ();

				const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

				gradients.resize(n_face_q_points);
				neighbor_gradients.resize(n_face_q_points);

				neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
				fe_face_values[velocities].get_function_gradients(solution, gradients);

				vector<Tensor<1,2>> jump_per_face;
				jump_per_face.resize(n_face_q_points);
				double jump_val=0;
				for (unsigned int q=0; q<n_face_q_points; ++q)
				{
					for (unsigned int i=0; i<2; ++i){
						jump_per_face[q] += (gradients[q][i]- neighbor_gradients[q][i])*(fe_face_values.normal_vector(q));		
					}// i
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

					hp_neighbor_face_values.reinit (cell->neighbor_child_on_subface (face_number, subface), cell->neighbor_of_neighbor(face_number));
					hp_subface_values.reinit (cell,face_number, subface, cell->active_fe_index());

					const FEFaceValues<2> &neighbor_face_values  = hp_neighbor_face_values.get_present_fe_values ();
					const FESubfaceValues<2> &fe_subface_values = hp_subface_values.get_present_fe_values ();


					const vector<double>& JxW_values = fe_subface_values.get_JxW_values ();

					const unsigned int n_subface_q_points = fe_subface_values.n_quadrature_points;//?

					gradients.resize(n_subface_q_points);
					neighbor_gradients.resize(n_subface_q_points);

					neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
					fe_subface_values[velocities].get_function_gradients(solution, gradients);

					vector<Tensor<1,2>> jump_per_subface;
					jump_per_subface.resize(n_subface_q_points);//?
					double jump_val=0;
					for (unsigned int q=0; q<n_subface_q_points; ++q)
					{
						for (unsigned int i=0; i<2; ++i){
							jump_per_subface[q] += (gradients[q][i]- neighbor_gradients[q][i])*(fe_subface_values.normal_vector(q));		
						}// i
						jump_val += contract(jump_per_subface[q],jump_per_subface[q])*(JxW_values[q]);
					}//q_per_subface
					term3 +=(face(face_number)->child(subface)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;
				}// subface
			}//else if

			// if the neighbor is coarser

			else if ( (cell->face(face_number)->at_boundary()==false) && (cell->neighbor_is_coarser(face_number)))
			{
				hp_fe_face_values.reinit(cell, face_number,cell->active_fe_index());
				hp_neighbor_subface_values.reinit(cell->neighbor(face_number),cell->neighbor_of_coarser_neighbor(face_number).first, cell->neighbor_of_coarser_neighbor(face_number).second);

				const FEFaceValues<2> &fe_face_values  = hp_fe_face_values.get_present_fe_values ();
				const FESubfaceValues<2> &neighbor_subface_values = hp_neighbor_subface_values.get_present_fe_values ();


				const vector<double>& JxW_values = fe_face_values.get_JxW_values ();

				const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

				gradients.resize(n_face_q_points);
				neighbor_gradients.resize(n_face_q_points);

				neighbor_subface_values[velocities].get_function_gradients(solution, neighbor_gradients);
				fe_face_values[velocities].get_function_gradients(solution, gradients);

				vector<Tensor<1,2>> jump_per_face;
				jump_per_face.resize(n_face_q_points);
				double jump_val=0;
				for (unsigned int q=0; q<n_face_q_points; ++q)
				{
					for (unsigned int i=0; i<2; ++i){
						jump_per_face[q] += (gradients[q][i]- neighbor_gradients[q][i])*(fe_face_values.normal_vector(q));		
					}// i
					jump_val += contract(jump_per_face[q],jump_per_face[q])*JxW_values[q];
				}//q_per_face
				term3 +=(cell->face(face_number)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;

			} // else if coarse neighbor

			Jump_est_per_cell(cell_index) = sqrt(term3);
			est_per_cell(cell_index)=Jump_est_per_cell(cell_index)+res_est_per_cell(cell_index);
	}//cell
	L1_norm_est= est_per_cell.l1_norm();
	cout<< "L1_norm is: "<< L1_norm_est << endl;
}// func.estimate ()
/*......................................................................................*/
vector<hp::DoFHandler<2>::active_cell_iterator> void StokesProblem::get_patch_around_cell(const typename hp::DoFHandler<2>::active_cell_iterator &cell)   
		{
			vector<hp::DoFHandler<2>::active_cell_iterator> patch;  
			patch.push_back (cell);
			for (unsigned int face_number=0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)
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
/*......................................................................................*/
unsigned int StokesProblem::count_dofs_on_patch(const vector<hp::DoFHandler<2>::active_cell_iterator> &patch)
		{
			set<types::global_dof_index> dofs_on_patch;
			vector<types::global_dof_index> local_dof_indices;
			typename hp::DoFHandler<2>::active_cell_iterator  cell=patch.begin();
			const typename hp::DoFHandler<2>::active_cell_iterator endc=patch.end ();
			for(; cell!=endc ; ++cell)
			{
				local_dof_indices.resize (cell->get_fe().dofs_per_cell);
				cell->get_dof_indices (local_dof_indices);
				dofs_on_patch.insert (local_dof_indices.begin(),local_dof_indices.end());
			}
			return dofs_on_patch.size();  
		}
/*......................................................................................*/
map<types::global_dof_index,unsigned int> StokesProblem::map_global_dofs_to_patch_indices (const vector<hp::DoFHandler<2>::active_cell_iterator> &patch)
		{
			map<types::global_dof_index,unsigned int> dofs_mapping;
			vector<types::global_dof_index> local_dof_indices;

			unsigned int next_unused_patch_dof_index = 0;
			typename hp::DoFHandler<2>::active_cell_iterator  cell=patch.begin();
			const typename hp::DoFHandler<2>::active_cell_iterator endc=patch.end ();
			for(; cell!=endc ; ++cell)
			{
				local_dof_indices.resize (cell->get_fe().dofs_per_cell);
				cell->get_dof_indices(local_dof_indices);

				for (unsigned int i=0; i < cell->get_fe().dofs_per_cell ; ++i)
					if (dofs_mapping.find(local_dof_indices[i]) != dofs_mapping.end())  
					{
						dofs_mapping.insert(make_pair(local_dof_indices[i], next_unused_patch_dof_index));
						++next_unused_patch_dof_index;
					}
			}
			return dof_mapping;
		}
/*......................................................................................*/


/*......................................................................................*/





/*......................................................................................*/
void StokesProblem::convergence_estimate () 
{
	Vector<double> est_per_cell(triangulation.n_active_cells());
	estimate(est_per_cell);

	//Vector<double> convergence_est_per_cell(triangulation.n_active_cells());  //? reinit in postprocess
	//Vector<double> solu_norm_per_patch(triangulation.n_active_cells());
	double convergence_est_per_cell;
	double solu_norm_per_patch;


	unsigned cell_index=0;
	typename hp::DoFHandler<2>::active_cell_iterator
		cell = dof_handler.begin_active(),
		endc = dof_handler.end();
	unsigned cell_index=0;
	for (; cell!=endc; ++cell,++cell_index)
	{
		

		


		

		/**************************************************************************************************/
		solve_on_each_patch ()
		{ 	
			vector<hp::DoFHandler<2>::active_cell_iterator> patch = get_patch_around_cell <hp::DoFHandler<2>> (cell);

			unsigned int local_system_size = count_dofs_on_patch <hp::DoFHandler<2>> (patch);

			// setup_patch_system and  patch_rhs                /////////////////////////////////////////////////

			dof_handler.distribute_dofs (fe_collection);
			vector<unsigned int> block_component (2+1, 0);
			block_component[2]=1;
			DoFRenumbering::component_wise(dof_handler, block_component);

			{
				constraints.clear ();
				make_hanging_node_constraints <hp::DoFHandler<2>> (dof_handler, constraints);
			}
			constraints.close();

			{
				BlockCompressedSetSparsityPattern csp (local_system_size,local_system_size);

				DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
				sparsity_pattern.copy_from(csp);
			}

			BlockSparseMatrix<double> patch_system (sparsity_pattern);
			BlockVector<double> patch_solution (local_system_size);
			BlockVector<double> patch_rhs (local_system_size);

			map<types::global_dof_index,unsigned int> global_to_patch_index_map;
			global_to_patch_index_map = map_global_dofs_to_patch_indices <hp::DoFHandler<2>> (patch); 

			// assemble  patch_system  and patch_rhs 

			hp::FEValues<2> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);

			FullMatrix<double> local_matrix;
			Vector<double> local_rhs;
			Vector<double> local_rhs1;
			Vector<double> local_rhs2;
			vector<types::global_dof_index> local_dof_indices;


			vector<Vector<double> >  rhs_values;

			const FEValuesExtractors::Vector velocities (0);
			const FEValuesExtractors::Scalar pressure (2);

			vector<Tensor<2,2> > grad_phi_u;
			vector<double> div_phi_u;
			vector<Tensor<1,2> > phi_u;
			vector<double> phi_p;

			vector<Tensor<1,2>> gradients_p;
			vector<double> divergences;
			vector<Tensor<1,2>> laplacians;

			vector<double> values;
			vector<Tensor<2,2>> gradients;



			typename hp::DoFHandler<2>::active_cell_iterator  cell=patch.begin();
			const typename hp::DoFHandler<2>::active_cell_iterator endc=patch.end ();
			for(; cell!=endc ; ++cell)
			{
				const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
				local_matrix.reinit (dofs_per_cell, dofs_per_cell);
				local_rhs.reinit (dofs_per_cell);
				local_rhs1.reinit (dofs_per_cell);
				local_rhs2.reinit (dofs_per_cell);


				hp_fe_values.reinit (cell);
				const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
				const vector<double>& JxW_values = fe_values.get_JxW_values ();
				const unsigned int n_q_points = fe_values.n_quadrature_points;

				rhs_values.resize(n_q_points, Vector<double>(3));
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

				fe_values[velocities].get_function_gradients(solution, gradients);
				fe_values[pressure].get_function_values(solution, values);

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
						for (unsigned int l=0; l<2; ++l)
							local_rhs1(i)+= contract(rhs_values[q](l)+laplacians[q][l]-gradients_p[q][l],phi_u[i][l])* JxW_values[q];
						local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];
						local_rhs(i)= local_rhs1(i)-local_rhs2(i);
					}// i
				}//q

				vector<types::global_dof_index> patch_dof_indices;/////////////////////////////////////////////////
				patch_dof_indices.resize (cell->get_fe().dofs_per_cell);

				for (unsigned int i=0; i<dofs_per_cell; ++i)
					patch_dof_indices[i] =global_to_patch_index_map[local_dof_indices[i]];
				constraints.distribute_local_to_global (local_matrix, local_rhs, patch_dof_indices, patch_system, patch_rhs);  /////////////////////////////////////////////////

				/*					for (unsigned int i=0; i<dofs_per_cell; ++i)
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				patch_system.add (global_to_patch_index_map[local_dof_indices[i]],
				global_to_patch_index_map[local_dof_indices[j]],
				local_matrix(i,j));
				patch_rhs(global_to_patch_index_map[local_dof_indices[i]]) += 
				local_rhs(i);
				}
				*/
				// solve patch_system and patch_rhs ............get  patch_solution		
				SparseDirectUMFPACK A_inverse;
				A_inverse.initialize (patch_system.block(0,0),
					SparseDirectUMFPACK::AdditionalData());
				Vector<double> tmp (patch_solution.block(0).size());
				{
					Vector<double> schur_rhs (patch_solution.block(1).size());
					A_inverse.vmult (tmp, patch_rhs.block(0));
					patch_system.block(1,0).vmult (schur_rhs, tmp);
					schur_rhs -= patch_rhs.block(1);

					SchurComplement schur_complement (patch_system, A_inverse);
					SolverControl solver_control (patch_solution.block(1).size(),
						1e-6);
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

				patch_system.block(0,1).vmult (tmp, patch_solution.block(1));
				tmp *= -1.0;
				tmp += patch_rhs.block(0);
				A_inverse.vmult (patch_solution.block(0), tmp);
				constraints.distribute (patch_solution);

				// get the L2 norm of the gradient of velocity solution and pressure value

				double pressure_val=0;
				double grad_u_val=0;
				for (unsigned int q=0; q<n_q_points; ++q)
				{
					pressure_val +=values[q]*values[q]* JxW_values[q];

					for (unsigned int i=0; i<2; ++i)

						grad_u_val +=double_contract(gradients[q],gradients[q])* JxW_values[q];
				} // q

				solu_norm_per_patch(cell_index)+=(sqrt(pressure_val) + sqrt(grad_u_val));

			}// cells on patch

		}// solve
		convergence_est_per_cell(cell_index)=solu_norm_per_patch(cell_index)/est_per_cell(cell_index);
	}// cell

}// convergence estimator

/*......................................................................................*/
/*
void StokesProblem:: postprocess (const unsigned int cycle){

Vector<double> est_per_cell (triangulation.n_active_cells());
estimate(est_per_cell);

Vector<double> convergence_est_per_cell (triangulation.n_active_cells());
convergence_estimate(convergence_est_per_cell);


}
*/
/*......................................................................................*/

void StokesProblem::run(){

	for (unsigned int cycle=0; cycle<6; ++cycle)
	{
		cout << "Cycle " << cycle << ':' << endl;
		if (cycle == 0)

			generate_mesh();

		setup_system ();
		assemble_system();
		solve ();
		if (L1_norm_est < Tolerance) break;
		//compute_error();
		//estimate ();
		//convergence_estimate ();
		postprocess(cycle);
	}
}
/*......................................................................................*/
int main()
{
	StokesProblem stokesproblem;
	stokesproblem.run();
	return 0;
}
