#include "StokesProblem.hh"
#include<cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include "SchurComplement.hh"

/*......................................................................................*/
// constructor 

template <int dim>
StokesProblem<dim>::StokesProblem(): dof_handler(triangulation),  max_degree (6), Tolerance (0.001)

{
  for (unsigned int degree=1; degree<=max_degree; ++degree)
  {
   //fe_collection.push_back (FESystem<dim>(FE_Q<dim> (QGaussLobatto<1> (degree + 2)), dim,
   //	  FE_Q<dim> (QGaussLobatto<1> (degree + 1)), 1));
    fe_collection.push_back (FESystem<dim>(FE_Q<dim> (degree + 2), dim,
	  FE_Q<dim> (degree + 1), 1));
    quadrature_collection.push_back(QGauss<dim> (degree+3));
    face_quadrature_collection.push_back (QGauss<dim-1> (degree+2));
  }

  fe_collection.push_back (FESystem<dim>(FE_Nothing<dim>(), dim, 
          FE_Nothing<dim>(), 1));
  quadrature_collection.push_back(QGauss<dim>(1));
  face_quadrature_collection.push_back (QGauss<dim-1>(1));
}

/*.....................................................................................*/
template <int dim>
StokesProblem <dim>::~StokesProblem() {
  dof_handler.clear();
}

/*......................................................................................*/
  template <int dim>
bool StokesProblem <dim>::decreasing (const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator > &i, const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator > &j)
{
  return ((i.first) > (j.first));
}

/*......................................................................................*/
// Generate mesh


template <int dim>
void StokesProblem <dim>::generate_mesh(){

  std::vector<Point<dim> > vertices (8);

  vertices [0]=Point<dim> (-1,-1);
  vertices [1]=Point<dim> (0,-1);
  vertices [2]=Point<dim> (-1,0);
  vertices [3]=Point<dim> (0,0);
  vertices [4]=Point<dim> (1,0);
  vertices [5]=Point<dim> (-1,1);
  vertices [6]=Point<dim> (0,1);
  vertices [7]=Point<dim> (1,1);

  const unsigned int n_cells=3;
  std::vector<CellData<dim> > cell(n_cells);
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

  std::ofstream out ("grid-L-Shape.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);
 // std::cout<<"Number of active cells: "<< triangulation.n_active_cells() << std::endl;
 // std::cout<<"Total number of cells: " << triangulation.n_cells() << std::endl ;
}
//.....................................................................................
template <int dim>
void StokesProblem <dim>::set_global_active_fe_indices (hp::DoFHandler<dim> &dof_handler)
{
  typename hp::DoFHandler<dim>::active_cell_iterator cell= dof_handler.begin_active(), end_cell = dof_handler.end();
  for (; cell!=end_cell; ++cell)
// if (static_cast<unsigned int> (cell->active_fe_index())  < (fe_collection.size()-1) )
//cell->set_active_fe_index (static_cast<unsigned int> (cell->active_fe_index()));

     cell->set_active_fe_index (0);
}
/*......................................................................................*/
// setup system()

template <int dim>
void StokesProblem <dim>::setup_system(){

  system_matrix.clear();

  set_global_active_fe_indices(dof_handler);
  dof_handler.distribute_dofs (fe_collection);
  DoFRenumbering::Cuthill_McKee (dof_handler);

  std::vector<unsigned int> block_component (dim+1, 0);
  block_component[dim]=1;
  DoFRenumbering::component_wise(dof_handler, block_component);

  {
    constraints.clear ();
    FEValuesExtractors::Vector velocities(0);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler,0,exact_solution,constraints, fe_collection.component_mask(velocities));
  }
  constraints.close();

  std::vector<types::global_dof_index> dofs_per_block (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
  const unsigned int n_u=dofs_per_block[0], n_p=dofs_per_block[1];

  std::cout<< "Number of degrees of freedom: " << dof_handler. n_dofs()<<
    "(" << n_u << "+" << n_p << ")" << std::endl;


  {
    BlockCompressedSetSparsityPattern csp (dofs_per_block,dofs_per_block);

    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    sparsity_pattern.copy_from(csp);
  }

  system_matrix.reinit (sparsity_pattern);
  solution.reinit (dofs_per_block);
  system_rhs.reinit (dofs_per_block);
}
/*......................................................................................*/
// assemble system

template <int dim>
void StokesProblem <dim>::assemble_system () {
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients);

  FullMatrix<double> local_matrix;
  Vector<double> local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  //const RightHandSide<dim> rhs_function;
  std::vector<Vector<double> >  rhs_values;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);


  std::vector<SymmetricTensor<2,dim> > symgrad_phi_u; 
 // std::vector<Tensor<2,dim> > grad_phi_u;
  std::vector<double> div_phi_u;
  std::vector<Tensor<1,dim> > phi_u;
  std::vector<double> phi_p;

  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
	 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
    local_matrix.reinit (dofs_per_cell, dofs_per_cell);
    local_rhs.reinit (dofs_per_cell);
    local_matrix=0;
    local_rhs=0;

    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
    const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
    const unsigned int n_q_points = fe_values.n_quadrature_points;

    rhs_values.resize(n_q_points, Vector<double>(dim+1));
    rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);
    
    symgrad_phi_u.resize(dofs_per_cell);
    //grad_phi_u.resize(dofs_per_cell);
    div_phi_u.resize(dofs_per_cell);
    phi_u.resize (dofs_per_cell);
    phi_p.resize(dofs_per_cell);

    for (unsigned int q=0; q<n_q_points; ++q)
    {
      for (unsigned int k=0; k<dofs_per_cell; ++k)
      {
        symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
	//grad_phi_u[k] = fe_values[velocities].gradient (k, q);
	div_phi_u[k] = fe_values[velocities].divergence (k, q);
	phi_u[k] = fe_values[velocities].value (k, q);
	phi_p[k] = fe_values[pressure].value (k, q);
      }
      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
                   for (unsigned int j=0; j<=i; ++j)
					{

                                          local_matrix(i,j) += (symgrad_phi_u[i] * symgrad_phi_u[j]
                                          - div_phi_u[i] * phi_p[j]
                                          - phi_p[i] * div_phi_u[j]
                                          + phi_p[i] * phi_p[j])
                                          * fe_values.JxW(q);
                                           

/*
						local_matrix(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
							- div_phi_u[i] * phi_p[j]
						- phi_p[i] * div_phi_u[j]
						+ phi_p[i] * phi_p[j])
  							* JxW_values[q];
  */ 
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
for (unsigned int l=0; l< dof_handler.n_dofs() ; ++l)
  for (unsigned int m=0; m< dof_handler.n_dofs() ; ++m)
  //std::cout<<  "[" << l<< " , "  << m << "]= "  << patch_system.el (l,m) << "............." <<std::endl;
   
 std::cout << "[" << l << " , "  << m << "]= " << system_matrix.el (l,m) << std::endl;
}

/*......................................................................................*/
// Solve

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
	1e-6*schur_rhs.l2_norm());

    SolverCG<>    cg (solver_control);

    SparseDirectUMFPACK preconditioner;
    preconditioner.initialize (system_matrix.block(1,1),
	SparseDirectUMFPACK::AdditionalData());

    cg.solve (schur_complement, solution.block(1), schur_rhs,
	preconditioner);
    //cout<<" residuals of each step " << solver_control.enable_history_data() << endl;
    constraints.distribute (solution);
    //  std::cout << "  "
    //<< solver_control.last_step()
    //<< " outer CG Schur complement iterations for pressure"
    //<< std::endl;
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

  std::vector<double> values;
  double domain_mean_val_p=0;
  double measure_domain=0;
  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
	 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    hp_fe_values.reinit (cell);

    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
    const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
    const unsigned int n_q_points = fe_values.n_quadrature_points;
    values.resize(n_q_points);
    fe_values[pressure].get_function_values(solution, values);
    for (unsigned int q=0; q<n_q_points; ++q)
    {
      domain_mean_val_p += values[q]*JxW_values[q];
      measure_domain += JxW_values[q];
    }//q
  }//cell
  // 3 here is the area corresponding to our 3 cells.
 // return domain_mean_val_p/3.0;
  return domain_mean_val_p / measure_domain;
}
/*......................................................................................*/
// compute_error

  template <int dim>
void StokesProblem <dim>::compute_error (Vector<double> &error_per_cell)
{
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients);
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<double> values;
  std::vector<Tensor<2,dim> > gradients;
  std::vector<std::vector<Tensor<1,dim> > > exact_solution_gradients;
  std::vector<Vector<double> > exact_solution_values;

  //error_per_cell(triangulation.n_active_cells());

  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
	 endc = dof_handler.end();
  unsigned int cell_index=0;
  for (; cell!=endc; ++cell,++cell_index)
  {
    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
    const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
    const std::vector<Point<dim> >& quadrature_points = fe_values.get_quadrature_points();
    const unsigned int n_q_points = fe_values.n_quadrature_points;

    gradients.resize(n_q_points);
    values.resize(n_q_points);
    exact_solution_gradients.resize(n_q_points , std::vector<Tensor<1,dim> > (dim+1));//?
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
  //  std::cout<< "Vector of Compute Error per Cell: " << error_per_cell<< std::endl ;
  double L1_norm=error_per_cell.l1_norm();
  std::cout<< "L1_norm of ERROR is: "<< L1_norm << std::endl;
  double L2_norm=error_per_cell.l2_norm();
  std::cout<< "L2_norm of ERROR is: "<< L2_norm << std::endl;
}


/*......................................................................................*/
// compute_estimator

template <int dim>
void StokesProblem <dim>::estimate (Vector<double> &est_per_cell)  {
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);
  hp::FEFaceValues<dim> hp_fe_face_values(fe_collection, face_quadrature_collection, update_JxW_values|update_gradients|update_normal_vectors);
  hp::FEFaceValues<dim> hp_neighbor_face_values(fe_collection, face_quadrature_collection, update_gradients);
  hp::FESubfaceValues<dim> hp_subface_values(fe_collection, face_quadrature_collection, update_JxW_values|update_gradients|update_normal_vectors);
  hp::FESubfaceValues<dim> hp_neighbor_subface_values(fe_collection, face_quadrature_collection, update_gradients);

  std::vector<Tensor<1,dim> > gradients_p;
  std::vector<double> divergences;
  std::vector<Tensor<1,dim> > laplacians;

  std::vector<Tensor<2,dim> > gradients;
  std::vector<Tensor<2,dim> > neighbor_gradients;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  const RightHandSide<dim> rhs_function;
  std::vector<Vector<double> >  rhs_values;

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
    const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
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

	const std::vector<double>& JxW_values = fe_face_values.get_JxW_values ();

	const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

	gradients.resize(n_face_q_points);
	neighbor_gradients.resize(n_face_q_points);

	neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
	fe_face_values[velocities].get_function_gradients(solution, gradients);

	std::vector<Tensor<1,dim> > jump_per_face;//?
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


	  const std::vector<double>& JxW_values = fe_subface_values.get_JxW_values ();

	  const unsigned int n_subface_q_points = fe_subface_values.n_quadrature_points;//?

	  gradients.resize(n_subface_q_points);
	  neighbor_gradients.resize(n_subface_q_points);

	  neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
	  fe_subface_values[velocities].get_function_gradients(solution, gradients);

	  std::vector<Tensor<1,dim> > jump_per_subface;
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


	const std::vector<double>& JxW_values = fe_face_values.get_JxW_values ();

	const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

	gradients.resize(n_face_q_points);
	neighbor_gradients.resize(n_face_q_points);

	neighbor_subface_values[velocities].get_function_gradients(solution, neighbor_gradients);
	fe_face_values[velocities].get_function_gradients(solution, gradients);

	std::vector<Tensor<1,dim> > jump_per_face;
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
}// func.estimate ()


/*......................................................................................*/
// get_layers_of_patch_around_cell

  template <int dim>
std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> StokesProblem <dim>::get_patch_around_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell)
{

  std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch;
  std::set<typename hp::DoFHandler<dim>::active_cell_iterator> cells_done;

  patch.push_back (cell);
  cells_done.insert(cell);
  //  i counter for the number of patch layers ... n_layers
  for (unsigned int i=0; i<1; ++i)
  {
    const unsigned int patch_size = patch.size();
    for (unsigned int j=0; j<patch_size; ++j)
    { 
      for (unsigned int face_number=0; face_number< GeometryInfo<dim>::faces_per_cell; ++face_number)
      {
	if (patch[j]->face(face_number)->at_boundary()==false)
	{
	  if (patch[j]->face(face_number)->has_children() == false)
	  {
	    typename hp::DoFHandler<dim>::active_cell_iterator celll = patch[j]->neighbor(face_number);
	    if (cells_done.count(celll)==0)
	    {
	      patch.push_back(celll);
	      cells_done.insert(celll);
	    }
	  }
	  else
	    for (unsigned int subface=0; subface< patch[j]->face(face_number)->n_children(); ++subface)
	    {                
	      typename hp::DoFHandler<dim>::active_cell_iterator child_cell = patch[j]->neighbor_child_on_subface (face_number, subface);
	      if (cells_done.count(child_cell)==0)
	      {
		patch.push_back(child_cell);
		cells_done.insert(child_cell);
	      }
	    }// for subface             
	} // if at_boundary()==false
      }// face_number
    } // for j
  } // for i
  return patch;
}

/*......................................................................................*/
//get_cells_at_coarsest_common_level

template <int dim>
std::vector<typename hp::DoFHandler<dim>::cell_iterator> StokesProblem <dim>::get_cells_at_coarsest_common_level (const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch)
{
  Assert (patch.size() > 0, ExcMessage("vector containing patch cells should not be an empty vector!"));
  //Assert (patch.size() > 0, ExcInternalError());
  unsigned int min_level = static_cast<unsigned int> (patch[0]->level());
  unsigned int max_level = static_cast<unsigned int> (patch[0]->level());
  for (unsigned int i=0; i<patch.size();++i)
  {
    min_level = std::min (min_level, static_cast<unsigned int> (patch[i]->level()) );
    max_level = std::max (max_level, static_cast<unsigned int> (patch[i]->level()) );
  }

    std::set<typename hp::DoFHandler<dim>::cell_iterator>  uniform_cells;
 
    typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::const_iterator  patch_c;

    for (patch_c=patch.begin(); patch_c!=patch.end () ; ++patch_c){
      if (static_cast<unsigned int>((*patch_c)->level()) == min_level)
	uniform_cells.insert (*patch_c);
      else
      {
	typename hp::DoFHandler<dim>::cell_iterator parent = *patch_c;

	while (static_cast<unsigned int> (parent->level()) > min_level)
	  parent = parent-> parent();
	uniform_cells.insert (parent);
      }
    }

    return std::vector<typename hp::DoFHandler<dim>::cell_iterator> (uniform_cells.begin(), uniform_cells.end());

}	//get_cells_at_coarsest_common_level	
    
/*......................................................................................*/
//build_triangulation_from_patch
template <int dim>
void StokesProblem <dim>::build_triangulation_from_patch (const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch,
    Triangulation<dim> &local_triangulation, unsigned int &level_h_refine, unsigned int &level_p_refine, std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> & patch_to_global_tria_map )
{
  std::vector<typename hp::DoFHandler<dim>::cell_iterator> uniform_cells = get_cells_at_coarsest_common_level (patch); // uniform_cells as const vector?

  level_h_refine=static_cast<unsigned int> (patch[0]->level());  
  level_p_refine=static_cast<unsigned int> (patch[0]->active_fe_index());

  local_triangulation.clear();
  std::vector<Point<dim> > vertices;
  const unsigned int n_uniform_cells=uniform_cells.size();
  std::vector<CellData<dim> > cells(n_uniform_cells);
  unsigned int k=0;// for enumerating cells
  unsigned int i=0;// for enumerating vertices

  typename std::vector<typename hp::DoFHandler<dim>::cell_iterator>::const_iterator uniform_c;

  for (uniform_c=uniform_cells.begin(); uniform_c!=uniform_cells.end(); ++uniform_c)
  {
    bool repeat_vertex;
    for (unsigned int j=0;  j< GeometryInfo<dim>::vertices_per_cell; ++j)
    {
      Point<dim> position=(*uniform_c)->vertex (j);
      repeat_vertex=false;

      for (unsigned int m=0; m<i; ++m)
      {

	if (position == vertices[m]){ 
	  repeat_vertex=true;
	  cells[k].vertices[j]=m ;
	  break;
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
  local_triangulation.create_triangulation(vertices,cells,SubCellData());

  Assert (local_triangulation.n_active_cells() == uniform_cells.size(), ExcInternalError());

  local_triangulation.clear_user_flags ();
  unsigned int index=0;
  std::map<typename Triangulation<dim>::cell_iterator, typename hp::DoFHandler<dim>::cell_iterator> patch_to_global_tria_map_tmp;
  for (typename Triangulation<dim>::cell_iterator coarse_cell_t = local_triangulation.begin(); coarse_cell_t != local_triangulation.end(); ++coarse_cell_t, ++index)
{
      patch_to_global_tria_map_tmp.insert (std::make_pair(coarse_cell_t, uniform_cells[index]));
    AssertThrow (( std::fabs(coarse_cell_t->center()(0) - uniform_cells[index]->center()(0))<1e-12 && std::fabs(coarse_cell_t->center()(1) - uniform_cells[index]->center()(1)) <1e-12)  , ExcInternalError());
}

  bool refinement_necessary;
  do
  {
    refinement_necessary = false;
    for (typename Triangulation<dim>::active_cell_iterator cell_tt = local_triangulation.begin_active(); cell_tt != local_triangulation.end(); ++cell_tt)
   {
     // if (patch_to_global_tria_map.count(cell_tt)==0)
     if (patch_to_global_tria_map_tmp[cell_tt]->has_children())
      {
	cell_tt -> set_refine_flag();
	refinement_necessary = true;
        
      }
      else for (unsigned int i=0; i<patch.size(); ++i){
	if (patch_to_global_tria_map_tmp[cell_tt]==patch[i]){
// this flag shows that this cell is on the patch
	  cell_tt->set_user_flag();
        break;}
      }
    }//for tria_patch.begin...

    if (refinement_necessary)
    {
      local_triangulation.execute_coarsening_and_refinement ();

      for (typename Triangulation<dim>::cell_iterator cell_ttt = local_triangulation.begin(); cell_ttt != local_triangulation.end(); ++cell_ttt)  
	if (cell_ttt-> has_children()){
	     // Since the cell got children, then it should not be in the map anymore...children may be added into map
              

             // these children may not yet be in the map
	  for (unsigned int c=0; c< cell_ttt ->n_children(); ++c){
	    if (patch_to_global_tria_map_tmp.find(cell_ttt->child(c)) == patch_to_global_tria_map_tmp.end())
                {
	          patch_to_global_tria_map_tmp.insert (std::make_pair(cell_ttt ->child(c), patch_to_global_tria_map_tmp[cell_ttt]->child(c)));
                   AssertThrow( (std::fabs (cell_ttt ->child(c)->center()(0)- patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(0)) < 1e-12 && std::fabs (cell_ttt ->child(c)->center()(1)-           patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(1)) < 1e-12) , ExcInternalError());
                }
            }
          patch_to_global_tria_map_tmp.erase(cell_ttt);  
      }

// how do we know that here they pair the cells with the same center?
    }
  }
  while (refinement_necessary);

//for (typename Triangulation<dim>::cell_iterator tria_cell = local_triangulation.begin(); tria_cell != local_triangulation.end(); ++tria_cell)
//Assert ((tria_cell.center()(0) == patch_to_global_tria_map_tmp[tria_cell].center()(0) && tria_cell.center()(1) == patch_to_global_tria_map_tmp[tria_cell].center()(1))  , ExcInternalError());

typename std::map<typename Triangulation<dim>::cell_iterator, typename hp::DoFHandler<dim>::cell_iterator>::iterator map_tmp_it = patch_to_global_tria_map_tmp.begin(),
map_tmp_end = patch_to_global_tria_map_tmp.end();
 
for (; map_tmp_it!=map_tmp_end; ++map_tmp_it)
patch_to_global_tria_map[map_tmp_it->first] = map_tmp_it->second;

} // build_triangulation
/*......................................................................................................................................*/
//set_active_fe_indices for cells on and out of each patch


// mark the cells in 'tria' that exist in the patch: go
// through all cells in 'tria' and see whether the
// corresponding cell in the global triangulation is part of
// the 'patch' list of cells

  template <int dim>
void StokesProblem <dim>::set_active_fe_indices (hp::DoFHandler<dim> &local_dof_handler)
{
  typename hp::DoFHandler<dim>::active_cell_iterator patch_c= local_dof_handler.begin_active(), end_patch_c = local_dof_handler.end();
  for (; patch_c!=end_patch_c; ++patch_c){

    if (patch_c->user_flag_set()==true)
    {
// if (static_cast<unsigned int> (patch_c->get_active_fe_indices()) +1 < fe_collection.size()-1 )
// patch_c->set_active_fe_index (static_cast<unsigned int> (patch_c->get_active_fe_indices()) + 1);

      patch_c->set_active_fe_index (0);
    }
    else if (patch_c->user_flag_set()==false)
    {	
      patch_c->set_active_fe_index (max_degree);
    }
    else
      Assert (false, ExcNotImplemented());
  }

}

//.................................................................................................................................
template <int dim>
void StokesProblem <dim>:: patch_output (unsigned int patch_number, const unsigned int cycle, hp::DoFHandler<dim> &local_dof_handler, BlockVector<double> &local_solu)
{

  std::vector<std::string> solution_names (dim, "patch_velocity");
  solution_names.push_back ("patch_pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

  DataOut<dim,hp::DoFHandler<dim> > patch_data_out;

  patch_data_out.attach_dof_handler (local_dof_handler);


  patch_data_out.add_data_vector (local_solu, solution_names, DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, data_component_interpretation);
  patch_data_out.build_patches ();

  std::string filename = "patch_solution-" +
  Utilities::int_to_string (cycle, 2) +
  +"-"+Utilities::int_to_string (patch_number, 2) +".vtu";
  std::ofstream output (filename.c_str());
  patch_data_out.write_vtu (output);
}
/*......................................................................................................................................*/
//Compute h_convergence_estimator   &   h_workload_number  for each patch around cell   

  template <int dim>
void StokesProblem <dim>:: h_patch_conv_load_no ( const unsigned int cycle , double &h_convergence_est_per_cell, unsigned int &h_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell, unsigned int & patch_number)
{
  
  std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell(cell);
  Triangulation<dim> local_triangulation;		
  unsigned int level_h_refine;
  unsigned int level_p_refine;
  std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> patch_to_global_tria_map;

  build_triangulation_from_patch (patch, local_triangulation, level_h_refine, level_p_refine, patch_to_global_tria_map);
  hp::DoFHandler<dim> local_dof_handler(local_triangulation);

  std::cout<< "patch_number is :" << " " << patch_number << std:: endl ;
  
  set_active_fe_indices (local_dof_handler);
  local_dof_handler.distribute_dofs (fe_collection);

  std::cout<< "number of dofs on the current patch " << local_dof_handler.n_dofs() << std::endl;

  h_convergence_est_per_cell=0.;
  double h_solu_norm_per_patch=0.;		

  ConstraintMatrix constraints_patch;		
  BlockSparsityPattern sparsity_pattern_patch;		

  std::vector<unsigned int> block_component_patch (dim+1, 0);
  block_component_patch[dim]=1;
  DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

  std::vector<types::global_dof_index> dofs_per_block_patch (2);
  DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);


  //......................................................................................................
BlockVector<double> local_solu (dofs_per_block_patch);


typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
    for (; patch_cl !=end_patch_cl; ++patch_cl)
{

const unsigned int   dofs_per_cl = patch_cl->get_fe().dofs_per_cell;
if (dofs_per_cl!=0){
   
Vector<double> local_solution_values(dofs_per_cl); 
   typename hp::DoFHandler<dim>::active_cell_iterator global_cell = patch_to_global_tria_map[patch_cl];
//std::cout << "dofs_per_cell :" << dofs_per_cl<< std::endl;
//std::cout << "global_cell :" << global_cell->get_fe().dofs_per_cell << std::endl;
//std::cout << " solutin size:" << solution.size()<< std::endl;
//std::cout << " global dofs number:" << dof_handler.n_dofs()<< std::endl;

    global_cell-> get_dof_values (solution,local_solution_values);
    patch_cl->set_dof_values (local_solution_values, local_solu);
}
}
  
  //......................  Solution Transfer .......................  h_refinement of patch cells ................................. // 

  bool need_to_refine = false;

  typename hp::DoFHandler<dim>::active_cell_iterator patch_cc= local_dof_handler.begin_active(), end_patch_cc = local_dof_handler.end();
  for (; patch_cc!=end_patch_cc; ++patch_cc){
    if (static_cast<unsigned int> (patch_cc->level()) <  (level_h_refine+1) )  {
      need_to_refine = true;
      patch_cc->set_refine_flag();
    }
  }//patch_cc

  if (need_to_refine == true)
  {
    local_triangulation.prepare_coarsening_and_refinement ();
    SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim>> solution_transfer(local_dof_handler);
    solution_transfer.prepare_for_pure_refinement();

    local_triangulation.execute_coarsening_and_refinement (); 

    local_dof_handler.distribute_dofs (fe_collection);

   // take a copy of the tmp vector
   // BlockVector<double>   local_solu(temp);  
    
    DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
    // resize the vector temp to the correct size
    BlockVector<double> temp (dofs_per_block_patch);
    solution_transfer.refine_interpolate(local_solu , temp);
    local_solu = temp;
  }

  //......................  setup_h_patch_system and  patch_rhs .............. //


  unsigned int local_system_size = local_dof_handler. n_dofs();
  h_workload_num = local_dof_handler. n_dofs(); 

  //patch_output (patch_number , cycle, local_dof_handler, local_solu);

  {
  constraints_patch.clear ();
  DoFTools::make_hanging_node_constraints(local_dof_handler, constraints_patch);  
  //DoFTools::make_hanging_node_constraints < hp::DoFHandler<dim> > (local_dof_handler, constraints_patch);  

  
  //......................... Zero_Bdry_Condition_on_Patch .......................................//
{	
    typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
    for (; patch_cl !=end_patch_cl; ++patch_cl){
      std::vector<types::global_dof_index> local_face_dof_indices ((patch_cl->get_fe()).dofs_per_face);
     if (patch_cl->user_flag_set() == true)
         {
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                bool face_is_on_patch_Bdry = false;
                if ((patch_cl->face(f)->at_boundary()) ||  ((patch_cl->neighbor(f)->has_children() == false)
	         &&
	        (patch_cl->neighbor(f)->user_flag_set() == false)))
                  face_is_on_patch_Bdry = true;
                else if ((patch_cl->face(f)->at_boundary()) || (patch_cl->neighbor(f)->has_children() == true))
                   {
                     for (unsigned int sf=0; sf< patch_cl->face(f)->n_children(); ++sf)
                         if (patch_cl->neighbor_child_on_subface (f, sf) -> user_flag_set() == false)
                          {
		            face_is_on_patch_Bdry = true;
		               break;
	                  }
                   }// else if
              
                     if (face_is_on_patch_Bdry)
	              {        
	              patch_cl->face(f)->get_dof_indices (local_face_dof_indices, patch_cl->active_fe_index());
	              for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
	              if ((patch_cl->get_fe()).face_system_to_component_index(i).first < dim)
		      constraints_patch.add_line (local_face_dof_indices[i]);
	              }    
             }// face f
         }// if user_flag_true

     }// patch_cl
 }

}
 constraints_patch.close();

//std::cout<< "dofs_per_block_patch:   " << dofs_per_block_patch[0] <<std::endl;
//std::cout<< "dofs_per_block_patch:   " << dofs_per_block_patch[1] <<std::endl;
//std::cout << " local_system_size:  " << local_system_size <<std::endl;

  //std::vector<types::global_dof_index> dofs_per_block_patch (2);
  DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
  const unsigned int n_u=dofs_per_block_patch[0], n_p=dofs_per_block_patch[1];

  {
    BlockCompressedSetSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);

    DoFTools::make_sparsity_pattern (local_dof_handler, csp, constraints_patch, false);
    sparsity_pattern_patch.copy_from(csp);
  }
  BlockSparseMatrix<double> patch_system (sparsity_pattern_patch);
  BlockVector<double> patch_solution (dofs_per_block_patch);
  BlockVector<double> patch_rhs (dofs_per_block_patch);
  // .........................................  assemble  patch_system  and patch_rhs .............................. //

  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);

  FullMatrix<double> local_matrix_patch;
  Vector<double> local_rhs_patch;
  Vector<double> local_rhs1;
  Vector<double> local_rhs2;
  std::vector<types::global_dof_index> local_dof_indices;

  std::vector<Vector<double> >  rhs_values;
  const RightHandSide<dim> rhs_function;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<Tensor<2,dim> > grad_phi_u;
  std::vector<double> div_phi_u;
  std::vector<Tensor<1,dim> > phi_u;
  std::vector<double> phi_p;

  std::vector<Tensor<1,dim> > gradients_p;
  std::vector<double> divergences;
  std::vector<Tensor<1,dim> > laplacians;

  std::vector<double> values;
  std::vector<Tensor<2,dim> > gradients;


  typename hp::DoFHandler<dim>::active_cell_iterator patch_cll= local_dof_handler.begin_active(), end_patch_cll = local_dof_handler.end();
  for (; patch_cll!=end_patch_cll; ++patch_cll){

    const unsigned int   dofs_per_cell = patch_cll->get_fe().dofs_per_cell;
    if (dofs_per_cell!=0) {
    local_matrix_patch.reinit (dofs_per_cell, dofs_per_cell);
    local_rhs_patch.reinit (dofs_per_cell);
    local_rhs1.reinit (dofs_per_cell);
    local_rhs2.reinit (dofs_per_cell);

    hp_fe_values.reinit (patch_cll);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
    const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
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

    fe_values[pressure].get_function_gradients(local_solu, gradients_p);
    fe_values[velocities].get_function_divergences(local_solu, divergences);
    fe_values[velocities].get_function_laplacians(local_solu, laplacians);


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
	  local_matrix_patch(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
	      - div_phi_u[i] * phi_p[j]
	      - phi_p[i] * div_phi_u[j]
	      + phi_p[i] * phi_p[j])
	    * JxW_values[q];
	} // end of loop over 'j'

	local_rhs1(i)+= ((rhs_values[q](0)+laplacians[q][0]-gradients_p[q][0])*(phi_u[i][0])+ (rhs_values[q](1)+laplacians[q][1]-gradients_p[q][1])*(phi_u[i][1]))* JxW_values[q];
	local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];
	local_rhs_patch(i)= local_rhs1(i)-local_rhs2(i);
      }// i
    }//q
 
  for (unsigned int i=0; i<dofs_per_cell; ++i)
       for (unsigned int j=i+1; j<dofs_per_cell; ++j)
	local_matrix_patch(i,j) = local_matrix_patch(j,i);    

    local_dof_indices.resize (dofs_per_cell);
    patch_cll->get_dof_indices (local_dof_indices);
    constraints_patch.distribute_local_to_global (local_matrix_patch, local_rhs_patch, local_dof_indices, patch_system, patch_rhs);
}
  }// for patch_cll

//  ................................................  symmetricity check  ...........................................................
{
BlockVector<double> temp1 (dofs_per_block_patch);
BlockVector<double> temp2 (dofs_per_block_patch);
BlockVector<double> unit_vector (dofs_per_block_patch);
for (unsigned int i=0; i<local_dof_handler.n_dofs(); ++i)
{	
  for (unsigned int j=0; j<local_dof_handler.n_dofs(); ++j)
  {
    if (i==j)
  	unit_vector[j]=1;
    else
      unit_vector[j]=0;
  } // for j

  patch_system.vmult (temp1,unit_vector );
  patch_system.Tvmult (temp2, unit_vector);

  for (unsigned int k=0; k<local_dof_handler.n_dofs(); ++k)
  {
    if(std::fabs(temp1[k]-temp2[k]) > 1e-8)
    //if (temp1[k]!=temp2[k])
     {
       std::cout<< "symmetricity check  "<< std::endl;
       std::cout<< "comparison between row and column number: " << i << std::endl;
       std::cout<< "vector temp1 [ "<< k << " ]=" << temp1[k] << std::endl;
       std::cout<< "vector temp2 [ "<< k << " ]=" << temp2[k] << std::endl; 
     }
  }// for k 

}// for i
}
  // .....................................solve patch_system and patch_rhs ............get  patch_solution ...................... //

SparseDirectUMFPACK  A_direct;
A_direct.initialize(patch_system);
A_direct.vmult (patch_solution, patch_rhs);



// pressure mass matrix to precondition the Schur complement
//  Later, when solving, we then precondition the Schur complement with M^{−1}_p by doing a few CG iterations on the well-conditioned pressure mass matrix Mp stored in the (1,1) block.

// For the Stokes equation we consider here, the Schur complement is BA^{−1}B^{T} where the matrix A is related to the Laplace operator 
// Thus, solving with A is a lot more complicated: the matrix is badly conditioned and we know that we need many iterations unless we have a very good preconditioner

//in 2d, we use the ultimate preconditioner, namely a direct sparse LU decomposition of the matrix. This is implemented using the SparseDirectUMFPACK class that uses the UMFPACK direct solver to compute the decomposition.

/*
 SparseDirectUMFPACK A_inverse;
  A_inverse.initialize (patch_system.block(0,0), SparseDirectUMFPACK::AdditionalData());
  Vector<double> tmp (patch_solution.block(0).size());
  {
    Vector<double> schur_rhs (patch_solution.block(1).size());
    A_inverse.vmult (tmp, patch_rhs.block(0));
    patch_system.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= patch_rhs.block(1);

    SchurComplement schur_complement (patch_system, A_inverse);
    SolverControl solver_control (patch_solution.block(1).size(), 1e-18);
//1e-6*schur_rhs.l2_norm()
    SolverCG<>    cg (solver_control);
 
//  Now to the preconditioner to the Schur complement. As explained in the introduction, the preconditioning is done by a mass matrix in the pressure variable.

//   the inversion of a mass matrix is a rather cheap and straight-forward operation (compared to, e.g., a Laplace matrix). The CG method with ILU preconditioning converges in 5-10 steps, independently on the mesh size. This is precisely what we do here: We choose another ILU preconditioner and take it along to the InverseMatrix object via the corresponding template parameter. A CG solver is then called within the vmult operation of the inverse matrix.


    SparseDirectUMFPACK preconditioner;
    preconditioner.initialize (patch_system.block(1,1),
    SparseDirectUMFPACK::AdditionalData());


    cg.solve (schur_complement, patch_solution.block(1), schur_rhs,
	preconditioner);
    constraints_patch.distribute (patch_solution);
    //    std::cout << "  "
    //<< solver_control.last_step()
    //<< " outer CG Schur complement iterations for pressure"
    //<< std::endl;
  }

  {
    patch_system.block(0,1).vmult (tmp, patch_solution.block(1));
    tmp *= -1.0;
    tmp += patch_rhs.block(0);
    A_inverse.vmult (patch_solution.block(0), tmp);
    constraints_patch.distribute (patch_solution);

  }
*/
  //.......................................  get the L2 norm of the gradient of velocity solution and pressure value  .....................//

  double pressure_val=0;
  double grad_u_val=0;
  typename hp::DoFHandler<dim>::active_cell_iterator patch_cel= local_dof_handler.begin_active(), end_patch_cel = local_dof_handler.end();
  for (; patch_cel!=end_patch_cel; ++patch_cel){
    const unsigned int   dofs_per_cel = patch_cel->get_fe().dofs_per_cell;
    if (dofs_per_cel!=0) {
    hp_fe_values.reinit (patch_cel);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
    const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
    const unsigned int n_q_points = fe_values.n_quadrature_points;
    gradients.resize(n_q_points);
    values.resize(n_q_points);

    fe_values[velocities].get_function_gradients(patch_solution, gradients);
    fe_values[pressure].get_function_values(patch_solution, values);

    for (unsigned int q=0; q<n_q_points; ++q)
    {
      pressure_val +=values[q]*values[q]* JxW_values[q];

      for (unsigned int i=0; i<dim; ++i)

	grad_u_val += contract(gradients[q][i],gradients[q][i])* JxW_values[q];//.... double contract?
    } // q
    h_solu_norm_per_patch +=(sqrt(pressure_val) + sqrt(grad_u_val));// ?
}
  }// cells on patch

  h_convergence_est_per_cell = h_solu_norm_per_patch ;

}  // function h_patch_con_load_num


/*......................................................................................*/

//Compute p_convergence_estimator   &   p_workload_number  for each patch around cell 



//
//	template <int dim>
//		void StokesProblem <dim>::p_patch_conv_load_no (double &p_convergence_est_per_cell, unsigned int &p_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell ) 
//		{
//
//		Triangulation<dim> local_triangulation;
//    std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell(cell);		
//		unsigned int level_h_refine;
//		unsigned int level_p_refine;
//		build_triangulation_from_patch (patch, local_triangulation, level_h_refine, level_p_refine);
//		hp::DoFHandler<dim> local_dof_handler(local_triangulation);//
//
//
//		p_convergence_est_per_cell=0.;
//		double p_solu_norm_per_patch=0.;	
//	
//		ConstraintMatrix constraints_patch;		
//		BlockSparsityPattern sparsity_pattern_patch;

//			//............................  p_ refinement on each patch ...................................//
//			typename hp::DoFHandler<dim>::active_cell_iterator patch_c= local_dof_handler.begin_active(), end_patch_c = local_dof_handler.end();
//						for (; patch_c!=end_patch_c; ++patch_c){
//
//				
//					if (static_cast<unsigned int> (patch_c->active_fe_index()) +1 < fe_collection.size() && static_cast<unsigned int> (patch_c->active_fe_index()) < (level_p_refine+1))
//						patch_c->set_active_fe_index (static_cast<unsigned int> (patch_c->active_fe_index()) + 1);
//				
//			}//for patch_c
//			//......................... setup   p_patch_system and  patch_rhs  .......................... //
//			local_dof_handler.distribute_dofs (fe_collection);
//		 	p_workload_num = local_dof_handler. n_dofs();//
//			unsigned int local_system_size = local_dof_handler.n_dofs();//
//
//      std::vector<unsigned int> block_component_patch(dim+1, 0);
//			block_component_patch[dim]=1;
//			DoFRenumbering::component_wise(local_dof_handler, block_component_patch);
//
//			{
//					constraints_patch.clear ();
//					FEValuesExtractors::Vector velocities(0);
//					DoFTools::make_hanging_node_constraints <hp::DoFHandler<dim> > (local_dof_handler, constraints_patch);
//			}
//			//...................... Zero_Bdry_Condition_on_Patch ........................................//
//			{
//				typename hp::DoFHandler<dim>::active_cell_iterator patch_c= local_dof_handler.begin_active(), end_patch_c = local_dof_handler.end();
//						for (; patch_c!=end_patch_c; ++patch_c){
//					std::vector<types::global_dof_index> local_face_dof_indices (patch_c->get_fe().dofs_per_face);  
//					
//						for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
//							if (! (patch_c->at_boundary(f)))
//							{
//								bool face_is_on_patch_Bdry = false;
//								if ((patch_c->neighbor(f)->has_children() == false)
//									&&
//									(patch_c->neighbor(f)->user_flag_set() == false))
//									face_is_on_patch_Bdry = true;
//								else if (patch_c->neighbor(f)->has_children() == true)
//								{
//									for (unsigned int sf=0; sf< patch_c->face(f)->n_children(); ++sf)
//										if (patch_c->neighbor_child_on_subface (f, sf)->user_flag_set() == false)
//										{
//											face_is_on_patch_Bdry = true;
//											break;
//										}
//								}
//								if (face_is_on_patch_Bdry)
//								{
//									patch_c->face(f)->get_dof_indices (local_face_dof_indices, 0);
//									for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
//										if (patch_c->get_fe().face_system_to_component_index(i).first < 2)
//											constraints_patch.add_line (local_face_dof_indices[i]);
//								}
//							}
//					
//				}//for patch_c
//			}
//
//			constraints_patch.close();
//      std::vector<types::global_dof_index> dofs_per_block_patch (2);
//				DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
//				const unsigned int n_u=dofs_per_block_patch[0], n_p=dofs_per_block_patch[1];
//
//				{
//				  BlockCompressedSetSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);
//
//					DoFTools::make_sparsity_pattern (local_dof_handler, csp, constraints_patch, false);
//					sparsity_pattern_patch.copy_from(csp);
//				}
//				BlockSparseMatrix<double> patch_system (sparsity_pattern_patch);
//				BlockVector<double> patch_solution (dofs_per_block_patch);
//				BlockVector<double> patch_rhs (dofs_per_block_patch);
//
//			// ..................................................  assemble  patch_system  and patch_rhs .............................. //
//			hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);
//
//			FullMatrix<double> local_matrix_patch;
//			Vector<double> local_rhs_patch;
//			Vector<double> local_rhs1;
//			Vector<double> local_rhs2;
//      std::vector<types::global_dof_index> local_dof_indices;
//
//      std::vector<Vector<double> >  rhs_values;
//
//			const FEValuesExtractors::Vector velocities (0);
//			const FEValuesExtractors::Scalar pressure (dim);
//
//      std::vector<Tensor<2,dim> > grad_phi_u;
//			std::vector<double> div_phi_u;
//			std::vector<Tensor<1,dim> > phi_u;
//			std::vector<double> phi_p;
//
//      std::vector<Tensor<1,dim> > gradients_p;
//			std::vector<double> divergences;
//			std::vector<Tensor<1,dim> > laplacians;
//
//      std::vector<double> values;
//			std::vector<Tensor<2,dim> > gradients;
//
//			typename hp::DoFHandler<dim>::active_cell_iterator patch_cc= local_dof_handler.begin_active(), end_patch_cc = local_dof_handler.end();
//						for (; patch_cc!=end_patch_cc; ++patch_cc){
//
//
//				const unsigned int   dofs_per_cell = patch_cc->get_fe().dofs_per_cell;
//				local_matrix_patch.reinit (dofs_per_cell, dofs_per_cell);
//				local_rhs_patch.reinit (dofs_per_cell);
//				local_rhs1.reinit (dofs_per_cell);
//				local_rhs2.reinit (dofs_per_cell);
//
//				hp_fe_values.reinit (patch_cc);
//				const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
//				const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
//				const unsigned int n_q_points = fe_values.n_quadrature_points;
//
//				rhs_values.resize(n_q_points, Vector<double>(dim+1));
//				rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);
//
//				grad_phi_u.resize(dofs_per_cell);
//				div_phi_u.resize(dofs_per_cell);
//				phi_u.resize (dofs_per_cell);
//				phi_p.resize(dofs_per_cell);
//
//				divergences.resize(n_q_points);
//				gradients_p.resize(n_q_points) ;
//				laplacians.resize(n_q_points);
//
//				fe_values[pressure].get_function_gradients(patch_solution, gradients_p);
//				fe_values[velocities].get_function_divergences(patch_solution, divergences);
//				fe_values[velocities].get_function_laplacians(patch_solution, laplacians);
//
//				gradients.resize(n_q_points);
//				values.resize(n_q_points);
//
//
//				for (unsigned int q=0; q<n_q_points; ++q)
//				{
//					for (unsigned int k=0; k<dofs_per_cell; ++k)
//					{
//						grad_phi_u[k] = fe_values[velocities].gradient (k, q);
//						div_phi_u[k] = fe_values[velocities].divergence (k, q);
//						phi_u[k] = fe_values[velocities].value (k, q);
//						phi_p[k] = fe_values[pressure].value (k, q);
//					}
//					for (unsigned int i=0; i<dofs_per_cell; ++i)
//					{
//						for (unsigned int j=0; j<dofs_per_cell; ++j)
//						{
//							local_matrix_patch(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
//								- div_phi_u[i] * phi_p[j]
//							- phi_p[i] * div_phi_u[j]
//							+ phi_p[i] * phi_p[j])
//								* JxW_values[q];
//						} // end of loop over 'j'
//
//						local_rhs1(i)+= ((rhs_values[q](0)+laplacians[q][0]-gradients_p[q][0])*(phi_u[i][0])+ (rhs_values[q](1)+laplacians[q][1]-gradients_p[q][1])*(phi_u[i][1]))* JxW_values[q];//?
//						local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];//?
//						local_rhs_patch(i)= local_rhs1(i)-local_rhs2(i);//?
//					}// i
//				}//q
//
//
//				local_dof_indices.resize (dofs_per_cell);
//				patch_cc->get_dof_indices (local_dof_indices);
//				constraints_patch.distribute_local_to_global (local_matrix_patch, local_rhs_patch, local_dof_indices, patch_system, patch_rhs);
//			}// for patch_cc
//
//			// .....................................solve patch_system and patch_rhs ............get  patch_solution ...................... //
//			
//				SparseDirectUMFPACK A_inverse;
//				A_inverse.initialize (patch_system.block(0,0), SparseDirectUMFPACK::AdditionalData());
//				Vector<double> tmp (patch_solution.block(0).size());
//				{
//					Vector<double> schur_rhs (patch_solution.block(1).size());
//					A_inverse.vmult (tmp, patch_rhs.block(0));
//					patch_system.block(1,0).vmult (schur_rhs, tmp);
//					schur_rhs -= patch_rhs.block(1);
//
//					SchurComplement schur_complement (patch_system, A_inverse);
//					SolverControl solver_control (patch_solution.block(1).size(), 1e-6);
//					SolverCG<>    cg (solver_control);
//					SparseDirectUMFPACK preconditioner;
//					preconditioner.initialize (patch_system.block(1,1),
//						SparseDirectUMFPACK::AdditionalData());
//
//					cg.solve (schur_complement, patch_solution.block(1), schur_rhs,
//						preconditioner);
//					constraints_patch.distribute (patch_solution);
//        //  std::cout << "  "
//						//<< solver_control.last_step()
//						//<< " outer CG Schur complement iterations for pressure"
//						//<< std::endl;
//				}
//
//				{
//					patch_system.block(0,1).vmult (tmp, patch_solution.block(1));
//					tmp *= -1.0;
//					tmp += patch_rhs.block(0);
//					A_inverse.vmult (patch_solution.block(0), tmp);
//					constraints_patch.distribute (patch_solution);
//				}
//
//				//.........................  get the L2 norm of the gradient of velocity solution and pressure value  .....................//
//
//				
//				double pressure_val=0;
//				double grad_u_val=0;
//
//				typename hp::DoFHandler<dim>::active_cell_iterator patch_ccc= local_dof_handler.begin_active(), end_patch_ccc = local_dof_handler.end();
//						for (; patch_ccc!=end_patch_ccc; ++patch_ccc){
//					hp_fe_values.reinit (patch_ccc);
//					const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values ();
//					const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
//					const unsigned int n_q_points = fe_values.n_quadrature_points;
//					gradients.resize(n_q_points);
//					values.resize(n_q_points);
//
//					fe_values[velocities].get_function_gradients(patch_solution, gradients);
//					fe_values[pressure].get_function_values(patch_solution, values);
//
//					for (unsigned int q=0; q<n_q_points; ++q)
//					{
//						pressure_val +=values[q]*values[q]* JxW_values[q];
//
//						for (unsigned int i=0; i<2; ++i)
//
//							grad_u_val += contract(gradients[q][i],gradients[q][i])* JxW_values[q];// ......double_contract?
//					} // q
//
//				}// cells on patch
//				p_solu_norm_per_patch +=(sqrt(pressure_val) + sqrt(grad_u_val));  //?
//			
//			p_convergence_est_per_cell = p_solu_norm_per_patch ;
//
//		} // function p_patch_conv_load_no
//


/*..............................................................................................................*/
//marking_cells()

  template <int dim>
void StokesProblem <dim>:: marking_cells (const unsigned int cycle, Vector<float> & marked_cells, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set, 
    std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map)
{

  const double theta= 0.3;
// based on Dorfler's paper, this \theta should be chosen as a number close to minimum of convergence estimator!

  std::vector<std::pair<double, typename hp::DoFHandler<dim>::active_cell_iterator> > to_be_sorted;

  Vector<double> est_per_cell (triangulation.n_active_cells());
  estimate(est_per_cell);

  Vector<double> convergence_est_per_cell (triangulation.n_active_cells());

  //std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > p_ref_map;

  unsigned int cell_index=0;
  unsigned int patch_number=0; 

  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
	 endc = dof_handler.end();
  for (; cell!=endc; ++cell , ++cell_index, ++patch_number)
  {
   
    double indicator_per_cell =0.0;

    double h_convergence_est_per_cell;
    unsigned int h_workload_num;
    h_patch_conv_load_no (cycle ,h_convergence_est_per_cell,h_workload_num, cell, patch_number);

/*
convergence_est_per_cell(cell_index)=h_convergence_est_per_cell;
std::cout << "convergence_est_per_cell: "<< convergence_est_per_cell(cell_index) << std::endl;
*/

    h_convergence_est_per_cell = h_convergence_est_per_cell /est_per_cell(cell_index);
    std::cout << "error_est_per_cell: "<< est_per_cell(cell_index) << std::endl;
    
    //	double p_convergence_est_per_cell;
    //	unsigned int p_workload_num;
    //	p_patch_conv_load_no (p_convergence_est_per_cell,p_workload_num, cell);
    //	p_convergence_est_per_cell = p_convergence_est_per_cell /est_per_cell(cell_index);

    double h_ratio= h_convergence_est_per_cell /  h_workload_num ;
    //	double p_ratio= p_convergence_est_per_cell /  p_workload_num ;

    //	if (h_ratio < p_ratio) {
    //		convergence_est_per_cell(cell_index)=p_convergence_est_per_cell;
    //		indicator_per_cell=convergence_est_per_cell(cell_index)*est_per_cell(cell_index);
    //		//p_ref_map[cell] = true;
    //		p_ref_map.insert (std::make_pair(cell, true));
    //	}
    //	else{



    convergence_est_per_cell(cell_index)=h_convergence_est_per_cell;
    std::cout << "convergence_est_per_cell: "<< convergence_est_per_cell(cell_index) << std::endl;
    indicator_per_cell=convergence_est_per_cell(cell_index)*est_per_cell(cell_index);
    //p_ref_map[cell] = false;
    p_ref_map.insert (std::make_pair(cell, false));
    //	}

    //to_be_sorted.push_back(pair<double, typename hp::DoFHandler<dim>::active_cell_iterator> (indicator_per_cell,cell));........//?
    to_be_sorted.push_back(std::make_pair(indicator_per_cell,cell));
  }// cell
  std::sort (to_be_sorted.begin(), to_be_sorted.end(), std_cxx1x::bind(&StokesProblem<dim>::decreasing,this,std_cxx1x::_1,std_cxx1x::_2));

  double L2_norm=est_per_cell.l2_norm();


  double sum=0;
  for (unsigned int i=0; i< to_be_sorted. size(); ++i) {
    typename hp::DoFHandler<dim>::active_cell_iterator  cell_sort=to_be_sorted[i].second;
    sum+= (to_be_sorted[i].first)*(to_be_sorted[i].first);
    candidate_cell_set.push_back (cell_sort);
    if (sum >= (theta*(L2_norm))*(theta*(L2_norm)))
    break;
  }
  //..............................................................................................................................

  marked_cells =0.;
  unsigned int index=0;
  typename hp::DoFHandler<dim>::active_cell_iterator  celll = dof_handler.begin_active(),
	   endcl = dof_handler.end();
  for (; celll!=endcl; ++celll,++index)
  {
    typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::iterator  mark_candidate;
    for (mark_candidate=candidate_cell_set.begin(); mark_candidate!=candidate_cell_set.end(); ++ mark_candidate)
    {
      if (celll == (*mark_candidate))
	marked_cells(index)=1;

    }
  }
}// marking_cells ()

//.................................................................................................................................
//Output_result
  template <int dim>
void StokesProblem <dim>:: output_results (const unsigned int cycle, Vector<float> &marked_cells, Vector<double> &est_per_cell , Vector<double> &error_per_cell)
{

  Vector<float> fe_degrees (triangulation.n_active_cells());
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
	   endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      fe_degrees(index)= fe_collection[cell->active_fe_index()].degree;
  }



  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);


  DataOut<dim,hp::DoFHandler<dim> > data_out;

  data_out.attach_dof_handler (dof_handler);


  data_out.add_data_vector (solution, solution_names,DataOut<dim,hp::DoFHandler<dim> >::type_dof_data,data_component_interpretation);

 
  data_out.add_data_vector (marked_cells, "marked_cells");
  data_out.add_data_vector (fe_degrees, "fe_degree");
  data_out.add_data_vector (est_per_cell, "error_estimator");
  data_out.add_data_vector (error_per_cell, "error");

 

  data_out.build_patches ();
  std::string filename = "solution-" +
   Utilities::int_to_string (cycle, 2) +".vtu";
  std::ofstream output (filename.c_str());
  data_out.write_vtu (output);
}

/*..............................................................................................................*/
//refine_in_h_p()

  template <int dim>
void StokesProblem <dim>:: refine_in_h_p (const unsigned int cycle, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set, 
    std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map )

{
  bool need_to_h_refine=false;
  typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::iterator  cell_candidate;
  for (cell_candidate=candidate_cell_set.begin(); cell_candidate!=candidate_cell_set.end(); ++ cell_candidate)
  {
    std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell (*cell_candidate);
    //	Triangulation<dim> tria_patch;//
    //hp::DoFHandler<dim> local_dof_handler;//
    unsigned int level_h_refine = static_cast<unsigned int>((*cell_candidate)->level());
    unsigned int level_p_refine = static_cast<unsigned int>((*cell_candidate)->active_fe_index());
    //


    //build_triangulation_from_patch (patch, local_triangulation, level_h_refine, level_p_refine, patch_to_global_tria_map);
    //hp::DoFHandler<dim> local_dof_handler(local_triangulation);

    if (p_ref_map[*cell_candidate]==false)
    {	

       for (unsigned int i=0; i<patch.size(); ++i){
          if (static_cast<unsigned int> (patch[i]->level()) <  (level_h_refine+1) )  {
           need_to_h_refine = true;
           patch[i]->set_refine_flag();
           }// if
       }//patch[i]

if ( need_to_h_refine==true)
    triangulation.execute_coarsening_and_refinement();
    }//if 
    //......................................................................................................................................................//
    /*if(p_ref_map[*cell_candidate]==true){

      for (unsigned int i=0; i<patch.size(); ++i){

      if (static_cast<unsigned int> (patch[i]->active_fe_index())+ 1 < fe_collection.size() && static_cast<unsigned int> (patch[i]->active_fe_index()) < level_p_refine+1)
      patch[i]->set_active_fe_index (static_cast<unsigned int> (patch[i]->active_fe_index()) + 1);

      }//for i
      }// if 
     */

  }// for... candidate_cells
 
}// refine_in_h_p
/*......................................................................................................................................................*/
template <int dim>
void StokesProblem <dim>::run(){

  for (unsigned int cycle=0; cycle<1; ++cycle)
  {
    std::cout << "Cycle " << cycle << ':' << std::endl;
    if (cycle == 0)
      generate_mesh();
/*
    {
      std::vector<unsigned int> subdivisions (dim, 1);
      subdivisions[0] = 1;
      const Point<dim> bottom_left = (dim == 2 ?
                                      Point<dim>(-2,-1) :
                                      Point<dim>(-2,0,-1));
      const Point<dim> top_right   = (dim == 2 ?
                                      Point<dim>(2,0) :
                                      Point<dim>(2,1,0));
      GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                 subdivisions,
                                                 bottom_left,
                                                 top_right);
    }

    triangulation.refine_global (1);
*/
     std::cout<<"Number of active cells: "<< triangulation.n_active_cells() << std::endl;
    
     std::cout<<"Total number of cells: " << triangulation.n_cells() << std::endl ;
    setup_system ();
    assemble_system();

    solve ();
/*
    std::cout<<" *****************************SOLVE global System "<< std::endl;
    Vector<double> error_per_cell (triangulation.n_active_cells());
    compute_error (error_per_cell);
    Vector<double> est_per_cell (triangulation.n_active_cells());
    estimate(est_per_cell);

    //  std::cout<< "Vector of Error Estimate: "<< est_per_cell << std::endl;
    double L1_norm_est= est_per_cell.l1_norm();
    std::cout<< "L1_norm of ERROR Estimate is: "<< L1_norm_est << std::endl;
    double L2_norm_est= est_per_cell.l2_norm();	
    std::cout<< "L2_norm of ERROR Estimate is: "<< L2_norm_est << std::endl;	


    Vector<float> marked_cells(triangulation.n_active_cells());
    std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> candidate_cell_set; 
    std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > p_ref_map;
    marking_cells(cycle,  marked_cells, candidate_cell_set, p_ref_map);
    output_results(cycle, marked_cells, est_per_cell, error_per_cell);
    refine_in_h_p(cycle,  candidate_cell_set, p_ref_map);

    if (L2_norm_est < Tolerance) 
      break;
*/
  }


}//run


//Explicit initialization

template class StokesProblem<2>;
