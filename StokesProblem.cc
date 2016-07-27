#include "StokesProblem.hh"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <tuple>


template <int dim>
StokesProblem<dim>::StokesProblem(Parameters const &parameters):
  exact_solution(nullptr),
  rhs_function(nullptr),
  dof_handler(triangulation),
  verbose(parameters.get_verbose()),
  example(parameters.get_example()),
  refinement(parameters.get_refinement()),
  max_degree(parameters.get_max_degree())
{
  // Set the exact solution used in the Dirichlet boundary condition and to
  // compute the error. Also set the rhs.
  switch (example)
  {
    case (example_1) :
      {
        exact_solution.reset(new ExactSolutionEx1<dim>());
        rhs_function.reset(new RightHandSideEx1<dim>());

        break;
      }
    case (example_2) :
      {
        exact_solution.reset(new ExactSolutionEx2<dim>());
        rhs_function.reset(new RightHandSideEx2<dim>());

        break;
      }
    case (example_3) :
      {
        exact_solution.reset(new ExactSolutionEx3<dim>());
        rhs_function.reset(new RightHandSideEx3<dim>());

        break;
      }
    case (example_4) :
      {
        exact_solution.reset(new ExactSolutionEx4<dim>());
        rhs_function.reset(new RightHandSideEx4<dim>());

        break;
      }
    case (example_5) :
      {
        exact_solution.reset(new ExactSolutionEx5<dim>());
        rhs_function.reset(new RightHandSideEx5<dim>());

        break;
      }
    default :
      {
        AssertThrow(false, ExcMessage("Unknow Example"));
      }
  }

  // Set the quadrature and the fe_collection
  if (parameters.get_quadrature()==gauss_lobatto)
  {
    // One using GaussLobatto, we are seeting the number of points not the
    // degree of the polynomial -> we need to increase the "degree" by one
    for (unsigned int degree=1; degree<=max_degree; ++degree)
      fe_collection.push_back(FESystem<dim>(FE_Q<dim> (QGaussLobatto<1> (degree+3)), dim,
            FE_Q<dim> (QGaussLobatto<1> (degree+2)), 1));

    unsigned int n_quad_points = 4;
    for (unsigned int degree=1; degree<=max_degree; ++degree)
    {
      quadrature_collection.push_back(QGauss<dim> (n_quad_points));
      face_quadrature_collection.push_back(QGauss<dim-1> (n_quad_points));
      quadrature_collection_error.push_back(QGauss<dim> (n_quad_points+1));
      ++n_quad_points;
    }
  }
  else
  {
    for (unsigned int degree=1; degree<=max_degree; ++degree)
      fe_collection.push_back(FESystem<dim>(FE_Q<dim> (degree+2), dim,
            FE_Q<dim> (degree+1), 1));

    unsigned int n_quad_points(4);
    for (unsigned int degree=1; degree<=max_degree; ++degree)
    {
      quadrature_collection.push_back(QGauss<dim> (n_quad_points));
      face_quadrature_collection.push_back(QGauss<dim-1> (n_quad_points));
      quadrature_collection_error.push_back(QGauss<dim> (n_quad_points+1));
      ++n_quad_points;
    }
  }

  fe_collection.push_back(FESystem<dim>(FE_Nothing<dim>(), dim,
        FE_Nothing<dim>(), 1));
  quadrature_collection.push_back(QGauss<dim>(1));
  face_quadrature_collection.push_back (QGauss<dim-1>(1));


  // Set the dual source if necessary
  if (parameters.do_goal_oriented()==true)
  {
    dual_source.reset(new FunctionParser<dim>(dim+1));
    dual_source->initialize(FunctionParser<dim>::default_variable_names(),
        parameters.get_dual_source(), typename FunctionParser<dim>::ConstMap());
  }
}


template <int dim>
bool StokesProblem<dim>::sort_decreasing_order (
  const std::pair<double,DoFHandler_active_cell_iterator > &i,
  const std::pair<double,DoFHandler_active_cell_iterator > &j)
{
  return ((i.first) > (j.first));
}


template <int dim>
void StokesProblem<dim>::generate_mesh()
{
  // Generate a square/cube mesh or a L-shape mesh.
  if ((example==example_3) || (example==example_4))
    {
      GridGenerator::hyper_cube(triangulation, -1, 1);
      triangulation.refine_global(1);
    }
  else
    {
      std::vector<Point<dim>> vertices(8);

      vertices[0] = Point<dim>(-1,-1);
      vertices[1] = Point<dim>(0,-1);
      vertices[2] = Point<dim>(-1,0);
      vertices[3] = Point<dim>(0,0);
      vertices[4] = Point<dim>(1,0);
      vertices[5] = Point<dim>(-1,1);
      vertices[6] = Point<dim>(0,1);
      vertices[7] = Point<dim>(1,1);

      const unsigned int n_cells(3);
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
      triangulation.refine_global(1);
    }

  // Set the degree of the FE to the lowest possible on every cell.
  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  end_cell = dof_handler.end();
  for (; cell!=end_cell; ++cell)
    cell->set_active_fe_index(0);
}


template <int dim>
void StokesProblem<dim>::setup_system(PRIMALDUAL primal_dual)
{
  system_matrix.clear();

  dof_handler.distribute_dofs(fe_collection);

  std::vector<unsigned int> block_component (dim+1, 0);
  block_component[dim]=1;
  DoFRenumbering::component_wise(dof_handler, block_component);

  constraints.clear();
  FEValuesExtractors::Vector velocities(0);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  //TODO Here there is a problem in choosing the boundary condition for the dual
  //problem. For now, we do the easiest and choose homogeneous Dirichlet
  //boundary.
  if (primal_dual==primal)
  {
    // Project the boundary values. This allows to get the same results even
    // when changing the FE type.
    typename FunctionMap <dim>::type boundary_functions;
    boundary_functions[0]= exact_solution.get();
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::project_boundary_values (dof_handler, boundary_functions,
        face_quadrature_collection, boundary_values);
    IndexSet boundary_dofs;
    DoFTools::extract_boundary_dofs(dof_handler,
        fe_collection.component_mask(velocities), boundary_dofs);
    for (auto & p : boundary_values)
      if (boundary_dofs.is_element(p.first))
      {
        constraints.add_line(p.first);
        constraints.set_inhomogeneity(p.first, p.second);
      }
  }
  else
    // Since we impose homogenous boundary conditions, we can simply interpolate
    // the boundary values
    VectorTools::interpolate_boundary_values(dof_handler, 0, ZeroFunction<dim>(dim+1),
        constraints, fe_collection.component_mask(velocities));

  // Since with Dirichlet velocity Bdry, pressure will be defined up to a constant,
  // in order to make the solution to be unique, we need to add an additional
  // constraint for Pressure. We choose for example the first cell of triangulation
  // and do as follow:
  DoFHandler_active_cell_iterator first_cell = dof_handler.begin_active();
  std::vector<types::global_dof_index> local_dof_indices (first_cell->get_fe().dofs_per_cell);
  first_cell->get_dof_indices(local_dof_indices);
  types::global_dof_index first_pressure_dof =
    local_dof_indices[first_cell->get_fe().component_to_system_index(dim,0)];

  // component_to_system_index: "Compute the shape function for the given vector
  // component and index."
  constraints.add_line (first_pressure_dof);
  constraints.set_inhomogeneity (first_pressure_dof,0.);

  constraints.close();

  std::vector<types::global_dof_index> dofs_per_block (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
  const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];

  BlockDynamicSparsityPattern csp(2,2);
  csp.block(0,0).reinit(n_u,n_u);
  csp.block(1,0).reinit(n_p,n_u);
  csp.block(0,1).reinit(n_u,n_p);
  csp.block(1,1).reinit(n_p,n_p);
  csp.collect_sizes();
  DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
  sparsity_pattern.copy_from (csp);

  system_matrix.reinit(sparsity_pattern);
  solution.reinit(2);
  solution.block(0).reinit(n_u);
  solution.block(1).reinit(n_p);
  solution.collect_sizes();
  system_rhs.reinit(2);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.collect_sizes ();
}


template <int dim>
void StokesProblem<dim>::assemble_system(PRIMALDUAL primal_dual)
{
  system_matrix = 0.;
  system_rhs = 0.;

  hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
      update_values|update_quadrature_points|update_JxW_values|update_gradients);

  FullMatrix<double> local_matrix;
  Vector<double> local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  std::vector<Vector<double>> rhs_values;

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<Tensor<2,dim>> grad_phi_u;
  std::vector<double> div_phi_u;
  std::vector<Tensor<1,dim>> phi_u;
  std::vector<double> phi_p;

  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
  {
    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
    local_matrix.reinit(dofs_per_cell, dofs_per_cell);
    local_rhs.reinit(dofs_per_cell);
    local_matrix = 0;
    local_rhs = 0;

    hp_fe_values.reinit(cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    const std::vector<double> &JxW_values = fe_values.get_JxW_values();
    const unsigned int n_q_points = fe_values.n_quadrature_points;

    rhs_values.resize(n_q_points, Vector<double>(dim+1));
    if (primal_dual==primal)
      rhs_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);
    else
      dual_source->vector_value_list(fe_values.get_quadrature_points(), rhs_values);

    grad_phi_u.resize(dofs_per_cell);
    div_phi_u.resize(dofs_per_cell);
    phi_u.resize (dofs_per_cell);
    phi_p.resize(dofs_per_cell);

    for (unsigned int q=0; q<n_q_points; ++q)
    {
      for (unsigned int k=0; k<dofs_per_cell; ++k)
      {
        grad_phi_u[k] = fe_values[velocities].gradient(k,q);
        div_phi_u[k] = fe_values[velocities].divergence(k,q);
        phi_u[k] = fe_values[velocities].value(k,q);
        phi_p[k] = fe_values[pressure].value(k,q);
      }

      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        if (primal_dual==primal)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            local_matrix(i,j) += (scalar_product(grad_phi_u[i], grad_phi_u[j])
                - div_phi_u[i]*phi_p[j] - phi_p[i]*div_phi_u[j]) * JxW_values[q];
        else
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            local_matrix(i,j) += (scalar_product(grad_phi_u[i], grad_phi_u[j])
                - div_phi_u[i]*phi_p[j] + phi_p[i]*div_phi_u[j]) * JxW_values[q];

        double fe_rhs(0.);
        for (unsigned int d=0; d<dim; ++d)
          fe_rhs += phi_u[i][d]*rhs_values[q][d];
        
        local_rhs[i] += fe_rhs * JxW_values[q];
      }
    }

    //local system to global system
    local_dof_indices.resize (dofs_per_cell);
    cell->get_dof_indices (local_dof_indices);
    constraints.distribute_local_to_global(local_matrix, local_rhs,
        local_dof_indices, system_matrix, system_rhs);
  }
}


template <int dim>
void StokesProblem<dim>::solve(PRIMALDUAL primal_dual)
{
  SparseDirectUMFPACK A_inverse;
  A_inverse.initialize(system_matrix, SparseDirectUMFPACK::AdditionalData());
  if (primal_dual==primal)
  {
    A_inverse.vmult(solution, system_rhs);
    constraints.distribute(solution);
    solution.block(1).add(-pressure_mean_value());
  }
  else
  {
    A_inverse.vmult(dual_solution, system_rhs);
    constraints.distribute(dual_solution);
  }
}


// compute_error(): Computes the Energy norm of error which is the square
// root of (L-2 norm of error in velocity gradients + L-2 norm of pressure)
template <int dim>
void StokesProblem<dim>::compute_error()
{
  error_per_cell.reinit(triangulation.n_active_cells());
  Vect_Pressure_Err.reinit(triangulation.n_active_cells());
  Vect_grad_Velocity_Err.reinit(triangulation.n_active_cells());
  Vect_Velocity_Err.reinit(triangulation.n_active_cells());

  hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection_error,
      update_values|update_quadrature_points|update_JxW_values|update_gradients|
      update_hessians);
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<double> values;
  std::vector<Tensor<2,dim>> gradients;
  std::vector<Tensor<1,dim>> velocity_values;

  std::vector<std::vector<Tensor<1,dim>>> exact_solution_gradients;
  std::vector<Vector<double>> exact_solution_values;

  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  endc = dof_handler.end();

  for (unsigned int cell_index=0; cell!=endc; ++cell,++cell_index)
    {
      double subtract_p = 0.;
      double grad_u_vals = 0.;
      double u_vals = 0.;
      hp_fe_values.reinit(cell);
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
      const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
      const std::vector<Point<dim> > &quadrature_points = fe_values.get_quadrature_points();
      const unsigned int n_q_points = fe_values.n_quadrature_points;

      velocity_values.resize(n_q_points);
      gradients.resize(n_q_points);
      values.resize(n_q_points);
      exact_solution_gradients.resize(n_q_points, std::vector<Tensor<1,dim>>(dim+1));
      exact_solution_values.resize(n_q_points, Vector<double> (dim+1));

      fe_values[velocities].get_function_values(solution, velocity_values);
      fe_values[velocities].get_function_gradients(solution, gradients);
      fe_values[pressure].get_function_values(solution, values);

      exact_solution->vector_gradient_list(quadrature_points, exact_solution_gradients);
      exact_solution->vector_value_list(quadrature_points, exact_solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          values[q] -= exact_solution_values[q](dim);
          subtract_p += values[q]*values[q]* JxW_values[q];

          for (unsigned int i=0; i<dim; ++i)
            {
              velocity_values[q][i] -= exact_solution_values[q](i);
              gradients[q][i] -= exact_solution_gradients[q][i];
            }

          grad_u_vals += gradients[q].norm_square() * JxW_values[q];
          u_vals += velocity_values[q].norm_square() * JxW_values[q];
        }

      error_per_cell[cell_index] = std::sqrt(subtract_p + grad_u_vals);

      Vect_Pressure_Err[cell_index] = std::sqrt(subtract_p);
      Vect_grad_Velocity_Err[cell_index] = std::sqrt(grad_u_vals);
      Vect_Velocity_Err[cell_index] = std::sqrt(u_vals);
    }

  if (verbose)
    {
      std::cout<< std::endl;
      double L2_norm_grad_velocity_Err = Vect_grad_Velocity_Err.l2_norm();
      std::cout<< "L2_norm_grad_velocity_Err: "<< L2_norm_grad_velocity_Err << std::endl;
      std::cout<< std::endl;
      double L2_norm_velocity_Err = Vect_Velocity_Err.l2_norm();
      std::cout<< "L2_norm_velocity_Err: "<< L2_norm_velocity_Err << std::endl;
      std::cout<< std::endl;
      std::cout<< std::endl;
      double L2_norm_pressure_Err = Vect_Pressure_Err.l2_norm();
      std::cout<< "L2_norm_pressure_Err: "<< L2_norm_pressure_Err << std::endl;
      std::cout<< std::endl;
      std::cout<< std::endl;
      double L2_norm_total_Err = sqrt(std::pow(L2_norm_grad_velocity_Err,2) +
                                      std::pow (L2_norm_pressure_Err,2));
      std::cout<< "L2_norm of Total_ERROR is: "<< L2_norm_total_Err << std::endl;
      std::cout<< std::endl;
    }
}

// compute_error_estimator(): Computes the residual based error estimator in
// hp- adaptive FEM approach
template <int dim>
void StokesProblem<dim>::compute_error_estimator()
{
  est_per_cell.reinit(triangulation.n_active_cells());
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection,
      update_values|update_quadrature_points|update_JxW_values|update_gradients|
      update_hessians);
  hp::FEFaceValues<dim> hp_fe_face_values(fe_collection, face_quadrature_collection,
      update_JxW_values|update_gradients|update_normal_vectors);
  hp::FEFaceValues<dim> hp_neighbor_face_values(fe_collection,
      face_quadrature_collection, update_gradients);
  hp::FESubfaceValues<dim> hp_subface_values(fe_collection, face_quadrature_collection,
      update_JxW_values|update_gradients|update_normal_vectors);
  hp::FESubfaceValues<dim> hp_neighbor_subface_values(fe_collection,
      face_quadrature_collection, update_gradients);

  std::vector<Tensor<1,dim> > gradients_p;
  std::vector<double> divergences;
  std::vector<Tensor<1,dim> > laplacians;

  std::vector<Tensor<2,dim> > gradients;
  std::vector<Tensor<2,dim> > neighbor_gradients;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<Vector<double>> rhs_values;

  Vector<double> res_est_per_cell(triangulation.n_active_cells());
  Vector<double> Jump_est_per_cell(triangulation.n_active_cells());

  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  endc = dof_handler.end();
  unsigned int cell_index=0;
  for (; cell!=endc; ++cell,++cell_index)
  {
    hp_fe_values.reinit (cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
    const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
    const unsigned int n_q_points = fe_values.n_quadrature_points;

    rhs_values.resize(n_q_points, Vector<double>(dim+1));
    rhs_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);

    divergences.resize(n_q_points);
    gradients_p.resize(n_q_points);
    laplacians.resize(n_q_points);

    fe_values[pressure].get_function_gradients(solution, gradients_p);
    fe_values[velocities].get_function_divergences(solution, divergences);
    fe_values[velocities].get_function_laplacians(solution, laplacians);


    double term2=0;//divergence term in estimator definition
    double term1=0;// for the residual term in estimator definition
    for (unsigned int q=0; q<n_q_points; ++q)
    {
      term2 += (divergences[q])*(divergences[q])*JxW_values[q];

      for (unsigned int i=0; i<dim; ++i)
        gradients_p[q][i] -= (rhs_values[q](i)+laplacians[q][i]);

      term1 += scalar_product(gradients_p[q],gradients_p[q])*JxW_values[q];
    }
    res_est_per_cell[cell_index] = pow((cell->diameter())/(cell->get_fe().degree), 2.0)*(term1) + term2;

    //  compute jump_est_per_cell
    double term3=0;//jump part of the estimator
    for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
    {
      // Check that the face is not at the boundary
      if (cell->face(face_number)->at_boundary()==false)
      {
        // If the neighbor cell is active and on the same level
        if ((cell->neighbor(face_number)->active()) &&
            (cell->neighbor(face_number)->level() == cell->level()))
        {
          const unsigned int q_index = std::max (cell->active_fe_index(),
              cell->neighbor(face_number)->active_fe_index());

          hp_fe_face_values.reinit(cell, face_number, q_index);
          hp_neighbor_face_values.reinit (cell->neighbor(face_number),
              cell->neighbor_of_neighbor(face_number), q_index);

          const FEFaceValues<dim> &neighbor_face_values = hp_neighbor_face_values.get_present_fe_values();
          const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values ();
          const std::vector<double> &JxW_values = fe_face_values.get_JxW_values ();
          const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

          gradients.resize(n_face_q_points);
          neighbor_gradients.resize(n_face_q_points);

          neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
          fe_face_values[velocities].get_function_gradients(solution, gradients);

          std::vector<Tensor<1,dim> > jump_per_face;
          jump_per_face.resize(n_face_q_points);
          double jump_val=0;

          for (unsigned int q=0; q<n_face_q_points; ++q)
          {
            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int j=0; j<dim; ++j)
                jump_per_face[q][i] += (gradients[q][i][j]-neighbor_gradients[q][i][j]) *
                  (fe_face_values.normal_vector(q)[j]);
            jump_val += scalar_product(jump_per_face[q],jump_per_face[q])*JxW_values[q];
          }

          unsigned int min_fe_degree = std::min(cell->get_fe().degree,
              cell->neighbor(face_number)->get_fe().degree);


          term3 += (cell->face(face_number)->diameter())/(2.0*min_fe_degree)*jump_val;
        }
        // If the neighbor has children
        else if (cell->neighbor(face_number)->has_children() == true)
        {
          for (unsigned int subface=0; subface< cell->face(face_number)->n_children(); ++subface)
          {
            const unsigned int q_index = std::max(cell->neighbor_child_on_subface(
                  face_number, subface)->active_fe_index(), cell->active_fe_index());


            hp_neighbor_face_values.reinit(cell->neighbor_child_on_subface(
                  face_number, subface), cell->neighbor_of_neighbor(face_number), q_index);
            hp_subface_values.reinit (cell,face_number, subface, q_index);

            const FEFaceValues<dim> &neighbor_face_values = hp_neighbor_face_values.get_present_fe_values();
            const FESubfaceValues<dim> &fe_subface_values = hp_subface_values.get_present_fe_values();
            const std::vector<double> &JxW_values = fe_subface_values.get_JxW_values ();
            const unsigned int n_subface_q_points = fe_subface_values.n_quadrature_points;

            gradients.resize(n_subface_q_points);
            neighbor_gradients.resize(n_subface_q_points);

            neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
            fe_subface_values[velocities].get_function_gradients(solution, gradients);

            std::vector<Tensor<1,dim>> jump_per_subface;
            jump_per_subface.resize(n_subface_q_points);

            double jump_val=0;
            for (unsigned int q=0; q<n_subface_q_points; ++q)
            {
              for (unsigned int i=0; i<dim; ++i)
                for (unsigned int j=0; j<dim; ++j)
                  jump_per_subface[q][j] += (gradients[q][i][j] - neighbor_gradients[q][i][j]) *
                    (fe_subface_values.normal_vector(q)[j]);
              jump_val += scalar_product(jump_per_subface[q],jump_per_subface[q])*(JxW_values[q]);
            }

            unsigned int min_fe_degree = std::min(cell->get_fe().degree,
                cell->neighbor_child_on_subface(face_number,subface)->get_fe().degree);

            term3 +=(cell->face(face_number)->child(subface)->diameter()) / (2.0 *
                min_fe_degree)*jump_val;
          }
        }
        // The neighbor is coarser
        else
        {
          Assert(cell->neighbor_is_coarser(face_number),
              ExcMessage("Problem while computing the error estimator."));
          const unsigned int q_index = std::max(cell->active_fe_index(),
              cell->neighbor(face_number)->active_fe_index());
          hp_fe_face_values.reinit(cell, face_number,q_index);
          hp_neighbor_subface_values.reinit(cell->neighbor(face_number),
              cell->neighbor_of_coarser_neighbor(face_number).first,
              cell->neighbor_of_coarser_neighbor(face_number).second,q_index);

          const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values();
          const FESubfaceValues<dim> &neighbor_subface_values =
            hp_neighbor_subface_values.get_present_fe_values();

          const std::vector<double> &JxW_values = fe_face_values.get_JxW_values();

          const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

          gradients.resize(n_face_q_points);
          neighbor_gradients.resize(n_face_q_points);

          neighbor_subface_values[velocities].get_function_gradients(solution, neighbor_gradients);
          fe_face_values[velocities].get_function_gradients(solution, gradients);

          std::vector<Tensor<1,dim> > jump_per_face;
          jump_per_face.resize(n_face_q_points);
          double jump_val=0;
          for (unsigned int q=0; q<n_face_q_points; ++q)
          {
            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int j=0; j<dim; ++j)
                jump_per_face[q][i] += (gradients[q][i][j] - neighbor_gradients[q][i][j]) *
                  (fe_face_values.normal_vector(q)[j]);
            jump_val += scalar_product(jump_per_face[q],jump_per_face[q])*JxW_values[q];
          }


          unsigned int min_fe_degree = std::min(cell->get_fe().degree,
              cell->neighbor(face_number)->get_fe().degree);

          term3 += (cell->face(face_number)->diameter())/(2.0*min_fe_degree)*jump_val;
        }
      }
    }
    Jump_est_per_cell(cell_index) = term3;
    est_per_cell(cell_index) = sqrt(Jump_est_per_cell(cell_index)+res_est_per_cell(cell_index));
  }
}



template <int dim>
double StokesProblem<dim>::pressure_mean_value()
{
  hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection, update_values|update_JxW_values);
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<double> values;
  double domain_mean_val_p(0.);
  double measure_domain(0.);
  for (auto cell : dof_handler.active_cell_iterators())
  {
    hp_fe_values.reinit(cell);
    const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
    const std::vector<double> &JxW_values = fe_values.get_JxW_values();
    const unsigned int n_q_points = fe_values.n_quadrature_points;
    values.resize(n_q_points);
    fe_values[pressure].get_function_values(solution, values);
    for (unsigned int q=0; q<n_q_points; ++q)
    {
      domain_mean_val_p += values[q]*JxW_values[q];
      measure_domain += JxW_values[q];
    }
  }

  return domain_mean_val_p/measure_domain;
}

// get_patch_around_cell() : For each cell returns a vector of cells
// which are located around that cell
template <int dim>
std::vector<typename StokesProblem<dim>::DoFHandler_active_cell_iterator>
StokesProblem<dim>::get_patch_around_cell(const DoFHandler_active_cell_iterator &cell,
    const unsigned int n_layers)
{
  std::vector<DoFHandler_active_cell_iterator> patch;
  std::set<DoFHandler_active_cell_iterator> cells_done;

  patch.push_back (cell);
  cells_done.insert(cell);
  for (unsigned int i=0; i<n_layers; ++i)
  {
    const unsigned int patch_size = patch.size();
    for (unsigned int j=0; j<patch_size; ++j)
    {
      for (unsigned int face_number=0; face_number< GeometryInfo<dim>::faces_per_cell; 
          ++face_number)
      {
        if (patch[j]->face(face_number)->at_boundary()==false)
        {
          if (patch[j]->face(face_number)->has_children() == false)
          {
            DoFHandler_active_cell_iterator cell = patch[j]->neighbor(face_number);
            if (cells_done.count(cell)==0)
            {
              patch.push_back(cell);
              cells_done.insert(cell);
            }
          }
          else
          {
            for (unsigned int subface=0; subface< patch[j]->face(face_number)->n_children();
                ++subface)
            {
              DoFHandler_active_cell_iterator child_cell =
                patch[j]->neighbor_child_on_subface (face_number, subface);
              if (cells_done.count(child_cell)==0)
              {
                patch.push_back(child_cell);
                cells_done.insert(child_cell);
              }
            }
          }
        }
      }
    }
  }

  return patch;
}

// In order to build triangulation from patch cells, first we need to
// return the vector of the coarsest common level cells as follows
template <int dim>
std::vector<typename StokesProblem<dim>::DoFHandler_cell_iterator>
StokesProblem<dim>::get_cells_at_coarsest_common_level (
  const std::vector<DoFHandler_active_cell_iterator> &patch)
{
  Assert(patch.size()>0, 
      ExcMessage("vector containing patch cells should not be an empty vector!"));
  unsigned int min_level = static_cast<unsigned int> (patch[0]->level());
  unsigned int max_level = static_cast<unsigned int> (patch[0]->level());
  for (unsigned int i=0; i<patch.size(); ++i)
    {
      min_level = std::min (min_level, static_cast<unsigned int> (patch[i]->level()) );
      max_level = std::max (max_level, static_cast<unsigned int> (patch[i]->level()) );
    }

  std::set<DoFHandler_cell_iterator>  uniform_cells;

  typename std::vector<DoFHandler_active_cell_iterator>::const_iterator  patch_cell;

  for (patch_cell=patch.begin(); patch_cell!=patch.end () ; ++patch_cell)
    {
      if (static_cast<unsigned int>((*patch_cell)->level()) == min_level)
        uniform_cells.insert (*patch_cell);
      else
        {
          DoFHandler_cell_iterator parent = *patch_cell;

          while (static_cast<unsigned int> (parent->level()) > min_level)
            parent = parent-> parent();
          uniform_cells.insert (parent);
        }
    }

  return std::vector<DoFHandler_cell_iterator> (uniform_cells.begin(),
                                                uniform_cells.end());
}

// To solve local variational problems on each patch, we have to build the
//  corresponding triangulation on each patch
template <int dim>
void StokesProblem<dim>::build_triangulation_from_patch(
  const std::vector<DoFHandler_active_cell_iterator> &patch,
  Triangulation<dim> &local_triangulation, unsigned int &level_h_refine,
  unsigned int &level_p_refine, std::map<Triangulation_active_cell_iterator,
  DoFHandler_active_cell_iterator> &patch_to_global_tria_map)
{
  std::vector<DoFHandler_cell_iterator> uniform_cells =
    get_cells_at_coarsest_common_level (patch);

  level_h_refine=static_cast<unsigned int> (patch[0]->level());
  level_p_refine=static_cast<unsigned int> (patch[0]->active_fe_index());

  local_triangulation.clear();
  std::vector<Point<dim> > vertices;
  const unsigned int n_uniform_cells=uniform_cells.size();
  std::vector<CellData<dim> > cells(n_uniform_cells);
  unsigned int k=0;// for enumerating cells
  unsigned int i=0;// for enumerating vertices

  typename std::vector<DoFHandler_cell_iterator>::const_iterator uniform_c;

  for (uniform_c=uniform_cells.begin(); uniform_c!=uniform_cells.end(); ++uniform_c)
    {
      bool repeat_vertex;
      for (unsigned int j=0;  j< GeometryInfo<dim>::vertices_per_cell; ++j)
        {
          Point<dim> position=(*uniform_c)->vertex (j);
          repeat_vertex=false;

          for (unsigned int m=0; m<i; ++m)
            {
              if (position == vertices[m])
                {
                  repeat_vertex=true;
                  cells[k].vertices[j]=m ;
                  break;
                }
            }
          if (repeat_vertex==false)
            {
              vertices.push_back(position);
              cells[k].vertices[j]=i;
              i=i+1;
            }

        }//for vertices_per_cell
      k=k+1;
    }

  local_triangulation.create_triangulation(vertices,cells,SubCellData());
  Assert (local_triangulation.n_active_cells() == uniform_cells.size(), ExcInternalError());

  local_triangulation.clear_user_flags ();
  unsigned int index=0;
  std::map<Triangulation_cell_iterator, DoFHandler_cell_iterator> patch_to_global_tria_map_tmp;
  for (Triangulation_cell_iterator coarse_cell_t = local_triangulation.begin();
       coarse_cell_t != local_triangulation.end(); ++coarse_cell_t, ++index)
    {
      patch_to_global_tria_map_tmp.insert (std::make_pair(coarse_cell_t, uniform_cells[index]));
      Assert ((std::fabs(coarse_cell_t->center()(0) - uniform_cells[index]->center()(0))<1e-16 &&
                    std::fabs(coarse_cell_t->center()(1) - uniform_cells[index]->center()(1)) <1e-16),
                   ExcInternalError());
    }

  bool refinement_necessary;
  do
    {
      refinement_necessary = false;
      for (Triangulation_active_cell_iterator cell_tt = local_triangulation.begin_active();
           cell_tt != local_triangulation.end(); ++cell_tt)
        {
          if (patch_to_global_tria_map_tmp[cell_tt]->has_children())
            {
              cell_tt -> set_refine_flag();
              refinement_necessary = true;
            }
          else for (unsigned int i=0; i<patch.size(); ++i)
              {
                if (patch_to_global_tria_map_tmp[cell_tt]==patch[i])
                  {
                    cell_tt->set_user_flag();
                    break;
                  }
              }
        }

      if (refinement_necessary)
        {
          local_triangulation.execute_coarsening_and_refinement ();

          for (Triangulation_cell_iterator cell_ttt = local_triangulation.begin();
              cell_ttt != local_triangulation.end(); ++cell_ttt)
          {

            if (patch_to_global_tria_map_tmp.find(cell_ttt)!=patch_to_global_tria_map_tmp.end())
            {
              if (cell_ttt-> has_children())
              {
                // Note: Since the cell got children, then it should not be in the map anymore...
                // children may be added into the map, instead

                // these children may not yet be in the map
                for (unsigned int c=0; c< cell_ttt ->n_children(); ++c)
                {
                  if (patch_to_global_tria_map_tmp.find(cell_ttt->child(c)) ==
                      patch_to_global_tria_map_tmp.end())
                  {
                    patch_to_global_tria_map_tmp.insert (std::make_pair(cell_ttt ->child(c),
                          patch_to_global_tria_map_tmp[cell_ttt]->child(c)));

                    Assert((std::fabs (cell_ttt ->child(c)->center()(0) -
                            patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(0)) < 1e-16 &&
                          std::fabs (cell_ttt ->child(c)->center()(1) -
                            patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(1)) < 1e-16),
                        ExcInternalError());
                  }
                }
                patch_to_global_tria_map_tmp.erase(cell_ttt);
              }
            }
          }
        }

    }
  while (refinement_necessary);

  typename std::map<Triangulation_cell_iterator,DoFHandler_cell_iterator>::iterator map_tmp_it =
    patch_to_global_tria_map_tmp.begin(),map_tmp_end = patch_to_global_tria_map_tmp.end();

  for (; map_tmp_it!=map_tmp_end; ++map_tmp_it)
    patch_to_global_tria_map[map_tmp_it->first] = map_tmp_it->second;
}


// mark cells in the "local_triangulation" that exist in the patch: go
// through all cells in the "local_triangulation" and see whether the
// corresponding cell in the global triangulation is part of
// the 'patch' list of cells
template <int dim>
void StokesProblem<dim>::set_active_fe_indices(hp::DoFHandler<dim> &local_dof_handler,
                                               std::map<Triangulation_active_cell_iterator, 
                                               DoFHandler_active_cell_iterator>
                                               &patch_to_global_tria_map)
{
  DoFHandler_active_cell_iterator patch_cell = local_dof_handler.begin_active(),
                                  end_patch_cell = local_dof_handler.end();
  for (; patch_cell!=end_patch_cell; ++patch_cell)
    {
      if (patch_cell->user_flag_set()==true)
        {
          DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[patch_cell];

          patch_cell->set_active_fe_index(global_cell->active_fe_index());
        }
      else if (patch_cell->user_flag_set()==false)
        {
          // which assigns FE_Nothing for the cells out of patch
          patch_cell->set_active_fe_index(max_degree);
        }
      else
        Assert (false, ExcNotImplemented());
    }
}


template <int dim>
void StokesProblem<dim>::patch_output (unsigned int patch_number,
                                       const unsigned int cycle, 
                                       hp::DoFHandler<dim> &local_dof_handler,
                                       BlockVector<double> &local_solu)
{
  std::vector<std::string> solution_names;
  solution_names.push_back("patch_x_velocity");
  solution_names.push_back("patch_y_velocity");
  if (dim==3)
    solution_names.push_back("patch_z_velocity");
  solution_names.push_back("patch_pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  DataOut<dim,hp::DoFHandler<dim> > patch_data_out;
  patch_data_out.attach_dof_handler (local_dof_handler);


  patch_data_out.add_data_vector(local_solu, solution_names,
                                 DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, data_component_interpretation);
  patch_data_out.build_patches();

  std::string filename = "patch_solution-" +
                         Utilities::int_to_string (cycle, 2) +
                         +"-"+Utilities::int_to_string (patch_number, 2) +".vtu";
  std::ofstream output (filename.c_str());
  patch_data_out.write_vtu (output);
}


template <int dim>
void StokesProblem<dim>::p_refinement(hp::DoFHandler<dim> &local_dof_handler,
    std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator>
    &patch_to_global_tria_map, unsigned int level_p_refine, BlockVector<double> &local_solu,
    PRIMALDUAL primal_dual)
{
  bool need_to_refine = false;

  DoFHandler_active_cell_iterator active_patch_cell = local_dof_handler.begin_active(),
                                  active_end_patch_cell = local_dof_handler.end();
  for (; active_patch_cell!=active_end_patch_cell; ++active_patch_cell)
  {
    DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[active_patch_cell];

    // Refine the cell if the refinement level is below the refinement of the
    // ``main'' cell + 1 and below the maximum level of refinement.
    if ((global_cell->active_fe_index()+1) < (fe_collection.size()-1) &&
        global_cell->active_fe_index() < (level_p_refine+1))
    {
      need_to_refine = true;
      break;
    }

    // Throw an exception if we cannot increase the polynomial order when doing
    // a goal-oriented simulation. The goal-oriented algorithm requires the
    // error to be computed using higher polynomial order.
    if ((primal_dual == dual) && 
        (global_cell->active_fe_index()+1 == fe_collection.size()-1))
      AssertThrow(false, ExcMessage("Reached the maximum polynomial order."));
  }

  if (need_to_refine==true)
  {
    SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim> > solution_transfer(local_dof_handler);
    solution_transfer.prepare_for_pure_refinement();
    active_patch_cell = local_dof_handler.begin_active();
    for (; active_patch_cell!=active_end_patch_cell; ++active_patch_cell)
    {
      DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[active_patch_cell];
      if ((global_cell->active_fe_index()+1) < (fe_collection.size()-1) &&
          global_cell->active_fe_index() < (level_p_refine+1))
      {
        DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[active_patch_cell];
        active_patch_cell->set_active_fe_index(global_cell->active_fe_index()+1);
      }
    }

    local_dof_handler.distribute_dofs (fe_collection);
    std::vector<unsigned int> block_component_patch(dim+1, 0);
    block_component_patch[dim] = 1;
    DoFRenumbering::component_wise(local_dof_handler, block_component_patch);
    std::vector<types::global_dof_index> dofs_per_block_patch (2);
    DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
    BlockVector<double> temp (dofs_per_block_patch);
    solution_transfer.refine_interpolate(local_solu , temp);
    local_solu = temp;
  }
}


template <int dim>
void StokesProblem<dim>::h_refinement(Triangulation<dim> &local_triangulation,
                                      hp::DoFHandler<dim> &local_dof_handler,
                                      unsigned int level_h_refine,
                                      BlockVector<double> &local_solu)
{
  bool need_to_refine = false;
  DoFHandler_active_cell_iterator active_patch_cell = local_dof_handler.begin_active(),
                                  active_end_patch_cell = local_dof_handler.end();
  for (; active_patch_cell!=active_end_patch_cell; ++active_patch_cell)
    {
      if (static_cast<unsigned int> (active_patch_cell->level()) < (level_h_refine+1))
        {
          need_to_refine = true;
          active_patch_cell->set_refine_flag();
        }
    }

  if (need_to_refine==true)
    {
      // user flags will be overwritten by
      // execute_coarsening_and_refinement. save their values into
      // the material_id, since that one not only survives
      // refinement but is also inherited to the children
      active_patch_cell = local_dof_handler.begin_active();
      for (; active_patch_cell!=active_end_patch_cell; ++active_patch_cell)
        {
          if (active_patch_cell->user_flag_set())
            active_patch_cell->set_material_id (1);
          else
            active_patch_cell->set_material_id (0);
        }

      local_triangulation.prepare_coarsening_and_refinement ();

      SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim>> solution_transfer(local_dof_handler);
      solution_transfer.prepare_for_pure_refinement();

      local_triangulation.execute_coarsening_and_refinement ();

      // get user flags back out of the material_id field
      for (DoFHandler_cell_iterator patch_cell = local_dof_handler.begin ();
           patch_cell!=local_dof_handler.end(); ++patch_cell)
        {
          if (patch_cell->material_id() == 1)
            patch_cell->set_user_flag();
          else
            patch_cell->clear_user_flag();
        }

      local_dof_handler.distribute_dofs(fe_collection);
      std::vector<unsigned int> block_component_patch(dim+1, 0);
      block_component_patch[dim] = 1;
      DoFRenumbering::component_wise(local_dof_handler, block_component_patch);
      std::vector<types::global_dof_index> dofs_per_block_patch(2);
      DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
      // resize the vector temp to the correct size
      BlockVector<double> tmp (dofs_per_block_patch);
      solution_transfer.refine_interpolate(local_solu, tmp);
      local_solu = tmp;
    }
}


template <int dim>
void StokesProblem<dim>::patch_assemble_system(hp::DoFHandler<dim> const &local_dof_handler,
    ConstraintMatrix const &constraints_patch, BlockVector<double> const &local_solu,
    BlockSparseMatrix<double> &patch_system, BlockVector<double> &patch_rhs,
    PRIMALDUAL primal_dual)
{
  hp::FEValues<dim> hp_fe_values(hp::FEValues<dim> (fe_collection, quadrature_collection,
        update_values|update_quadrature_points|update_JxW_values|update_gradients|
        update_hessians));

  FullMatrix<double> local_matrix_patch;
  Vector<double> local_rhs_patch;
  std::vector<types::global_dof_index> local_dof_indices;

  std::vector<Vector<double>> rhs_values;

  FEValuesExtractors::Vector velocities (0);
  FEValuesExtractors::Scalar pressure (dim);

  std::vector<Tensor<2,dim>> grad_phi_u;
  std::vector<Tensor<1,dim>> phi_u;
  std::vector<double> phi_p;

  std::vector<Tensor<1,dim>> gradients_p;
  std::vector<double> divergences;
  std::vector<Tensor<1,dim>> laplacians;

  std::vector<double> values;
  std::vector<Tensor<2,dim>> gradients;

  DoFHandler_active_cell_iterator act_patch_cell = local_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      const unsigned int dofs_per_cell = act_patch_cell->get_fe().dofs_per_cell;
      if (dofs_per_cell!=0)
        {
          local_matrix_patch.reinit (dofs_per_cell, dofs_per_cell);
          local_rhs_patch.reinit (dofs_per_cell);

          hp_fe_values.reinit (act_patch_cell);
          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
          const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
          const unsigned int n_q_points = fe_values.n_quadrature_points;

          rhs_values.resize(n_q_points, Vector<double>(dim+1));
          if (primal_dual==primal)
            rhs_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);
          else
            dual_source->vector_value_list(fe_values.get_quadrature_points(), rhs_values);

          grad_phi_u.resize(dofs_per_cell);
          phi_u.resize (dofs_per_cell);
          phi_p.resize(dofs_per_cell);

          divergences.resize(n_q_points);
          gradients_p.resize(n_q_points);
          laplacians.resize(n_q_points);

          fe_values[pressure].get_function_gradients(local_solu, gradients_p);
          fe_values[velocities].get_function_divergences(local_solu, divergences);
          fe_values[velocities].get_function_laplacians(local_solu, laplacians);


          for (unsigned int q=0; q<n_q_points; ++q)
            {
              // Optimizations
              const double JxW_val(JxW_values[q]);
              const double div(divergences[q]);
              std::vector<double> alpha(dim,0.);
              for (unsigned int d=0; d<dim; ++d)
                alpha[d] = rhs_values[q][d] + laplacians[q][d]-gradients_p[q][d];

              for (unsigned int k=0; k<dofs_per_cell; ++k)
                {
                  grad_phi_u[k] = fe_values[velocities].gradient (k, q);
                  phi_u[k] = fe_values[velocities].value (k, q);
                  phi_p[k] = fe_values[pressure].value (k, q);
                }


              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    local_matrix_patch(i,j) +=
                      (scalar_product(grad_phi_u[i],grad_phi_u[j]) +
                       (phi_p[i]*phi_p[j])) * JxW_val;

                  double local_rhs1 = 0;
                  for (unsigned int d=0; d<dim; ++d)
                    local_rhs1 += alpha[d]*phi_u[i][d];
                  local_rhs1 *= JxW_val;

                  const double local_rhs2 =  phi_p[i]*div*JxW_val;
                  local_rhs_patch[i] += local_rhs1+local_rhs2;
                }
            }

          local_dof_indices.resize(dofs_per_cell);
          act_patch_cell->get_dof_indices (local_dof_indices);
          constraints_patch.distribute_local_to_global(local_matrix_patch, local_rhs_patch,
                                                       local_dof_indices, patch_system, patch_rhs);
        }
    }
}

template <int dim>
void StokesProblem<dim>::patch_solve(hp::DoFHandler<dim> &local_dof_handler,
    unsigned int patch_number, unsigned int cycle, BlockVector<double> &local_solu,
    double &conv_est, double &workload_num)
{
  // setup_patch_system and patch_rhs
  workload_num = local_dof_handler.n_dofs();
  if (verbose==true)
    patch_output(patch_number, cycle, local_dof_handler, local_solu);
  ConstraintMatrix constraints_patch;
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure (dim);
  DoFTools::make_hanging_node_constraints(local_dof_handler, constraints_patch);

  // Zero_Bdry_Condition_on_Patch
  DoFHandler_active_cell_iterator act_patch_cell = local_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
  {
    std::vector<types::global_dof_index> local_face_dof_indices((act_patch_cell->get_fe()).dofs_per_face);
    if (act_patch_cell->user_flag_set() == true)
    {
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        bool face_is_on_patch_Bdry = false;
        if (act_patch_cell->face(f)->at_boundary())
          face_is_on_patch_Bdry = true;
        else
        {
          if (act_patch_cell->neighbor(f)->has_children() == true)
          {
            for (unsigned int sf=0; sf<act_patch_cell->face(f)->n_children(); ++sf)
              if (act_patch_cell->neighbor_child_on_subface(f, sf)->user_flag_set() == false)
              {
                face_is_on_patch_Bdry = true;
                break;
              }
          }
          else
          {
            if (act_patch_cell->neighbor(f)->user_flag_set() == false)
              face_is_on_patch_Bdry = true;
          }
        }

        if (face_is_on_patch_Bdry)
        {
          act_patch_cell->face(f)->get_dof_indices (local_face_dof_indices,
              act_patch_cell->active_fe_index());
          for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
            // system_to_component_index: "Compute vector component and index
            // of this shape function within the shape functions corresponding
            // to this component from the index of a shape function within this
            // finite element"
            if ((act_patch_cell->get_fe()).face_system_to_component_index(i).first < dim)
              constraints_patch.add_line (local_face_dof_indices[i]);
        }
      }
    }
  }
  constraints_patch.close();

  std::vector<unsigned int> block_component_patch (dim+1, 0);
  block_component_patch[dim]=1;
  DoFRenumbering::component_wise(local_dof_handler, block_component_patch);
  std::vector<types::global_dof_index> dofs_per_block_patch (2);
  DoFTools::count_dofs_per_block(local_dof_handler, dofs_per_block_patch, block_component_patch);

  BlockDynamicSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);
  BlockSparsityPattern sparsity_pattern_patch;
  DoFTools::make_sparsity_pattern (local_dof_handler, csp, constraints_patch, false);
  sparsity_pattern_patch.copy_from(csp);

  BlockSparseMatrix<double> patch_system (sparsity_pattern_patch);
  BlockVector<double> patch_solution (dofs_per_block_patch);
  BlockVector<double> patch_rhs (dofs_per_block_patch);

  // assemble patch_system and patch_rhs
  patch_assemble_system(local_dof_handler, constraints_patch, local_solu, patch_system,
                        patch_rhs, primal);

  // direct solver
  SparseDirectUMFPACK A_inverse_stiffness;
  A_inverse_stiffness.initialize (patch_system.block(0,0),
                                  SparseDirectUMFPACK::AdditionalData());
  A_inverse_stiffness.vmult (patch_solution.block(0), patch_rhs.block(0));
  SparseDirectUMFPACK A_inverse_mass;
  A_inverse_mass.initialize (patch_system.block(1,1),
                             SparseDirectUMFPACK::AdditionalData());
  A_inverse_mass.vmult (patch_solution.block(1), patch_rhs.block(1));
  constraints_patch.distribute (patch_solution);

  // get the L2 norm of the gradient of velocity solution and pressure value
  double pressure_val=0;
  double grad_u_val=0;
  double solu_norm_per_patch=0.;
  std::vector<double> values;
  std::vector<Tensor<2,dim>> gradients;
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection,
      update_values|update_quadrature_points|update_JxW_values|update_gradients|
      update_hessians);

  act_patch_cell = local_dof_handler.begin_active();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      const unsigned int dofs_per_cel = act_patch_cell->get_fe().dofs_per_cell;
      if (dofs_per_cel!=0)
        {
          hp_fe_values.reinit (act_patch_cell);
          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
          const std::vector<double> &JxW_values = fe_values.get_JxW_values();
          const unsigned int n_q_points = fe_values.n_quadrature_points;

          gradients.resize(n_q_points);
          values.resize(n_q_points);

          fe_values[velocities].get_function_gradients(patch_solution, gradients);
          fe_values[pressure].get_function_values(patch_solution, values);


          for (unsigned int q=0; q<n_q_points; ++q)
            {
              pressure_val += values[q]*values[q]*JxW_values[q];
              grad_u_val += scalar_product(gradients[q],gradients[q]) * JxW_values[q];
            }
          solu_norm_per_patch += pressure_val + grad_u_val;
        }
    }
  conv_est = std::sqrt(solu_norm_per_patch);
}

// For both p-refinement() and h-refinement() computes the
// error reduction as we named it as patch_convergence_estimator()
template <int dim>
void StokesProblem<dim>::patch_convergence_estimator(const unsigned int cycle,
    SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
    std::vector<unsigned int>::iterator>> const &synch_iterator,
    ScratchData &scratch_data, CopyData<dim> &copy_data)
{
  // Silence a warning
  (void) scratch_data;

  DoFHandler_active_cell_iterator cell = std::get<0>(synch_iterator.iterators);
  const unsigned int patch_number(*std::get<1>(synch_iterator.iterators));
  copy_data.global_cell = cell;
  copy_data.cell_index = patch_number;

  Triangulation<dim> local_triangulation;
  unsigned int level_h_refine(0);
  unsigned int level_p_refine(0);
  std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator> patch_to_global_tria_map;
  std::vector<DoFHandler_active_cell_iterator> patch = get_patch_around_cell(cell, 1);
  build_triangulation_from_patch (patch, local_triangulation, level_h_refine, 
      level_p_refine, patch_to_global_tria_map);
  hp::DoFHandler<dim> local_dof_handler(local_triangulation);
  set_active_fe_indices(local_dof_handler,patch_to_global_tria_map);
  local_dof_handler.distribute_dofs(fe_collection);
  std::vector<unsigned int> block_component_patch(dim+1, 0);
  block_component_patch[dim]=1;
  DoFRenumbering::component_wise(local_dof_handler, block_component_patch);

  std::vector<types::global_dof_index> dofs_per_block_patch (2);
  DoFTools::count_dofs_per_block(local_dof_handler, dofs_per_block_patch, block_component_patch);

  BlockVector<double> local_solu(dofs_per_block_patch);

  // Here we are trying to project the values of the global vector "solution"
  // into vector "local_solu" which is solution over patch cells corresponding
  // to cell "K".
  DoFHandler_active_cell_iterator act_patch_cell = local_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      const unsigned int dofs_per_cell = act_patch_cell->get_fe().dofs_per_cell;
      // we check if the corresponding finite element for this cell is not 'FE_Nothing!'
      // and it takes usual finite element.
      if (dofs_per_cell!=0)
        {
          Vector<double> local_solution_values(dofs_per_cell);
          DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[act_patch_cell];

          global_cell->get_dof_values(solution,local_solution_values);
          act_patch_cell->set_dof_values(local_solution_values, local_solu);
        }
    }

  BlockVector<double> h_local_solu(local_solu);
  // The local triangulation is not modified by the h-refinement so, we do
  // p-refinement first.
  p_refinement(local_dof_handler, patch_to_global_tria_map, level_p_refine, local_solu);
  patch_solve(local_dof_handler, patch_number, cycle, local_solu, copy_data.p_conv_est_per_cell,
              copy_data.p_workload);
  // Reset the local_dof_handler to what it was before p_refinement was called
  set_active_fe_indices(local_dof_handler,patch_to_global_tria_map);
  local_dof_handler.distribute_dofs(fe_collection);
  DoFRenumbering::component_wise(local_dof_handler, block_component_patch);
  h_refinement(local_triangulation, local_dof_handler, level_h_refine, h_local_solu);
  patch_solve(local_dof_handler, patch_number, cycle, h_local_solu, copy_data.h_conv_est_per_cell,
              copy_data.h_workload);
}


// The marking strategy for cells to be selected for h- or p-refinement
template <int dim>
void StokesProblem<dim>::copy_to_refinement_maps(CopyData<dim> const &copy_data)
{
  h_Conv_Est[copy_data.cell_index] = copy_data.h_conv_est_per_cell;
  p_Conv_Est[copy_data.cell_index] = copy_data.p_conv_est_per_cell;

  const double h_improv = copy_data.h_conv_est_per_cell/est_per_cell[copy_data.cell_index];
  const double p_improv = copy_data.p_conv_est_per_cell/est_per_cell[copy_data.cell_index];

  double h_ratio(h_improv/copy_data.h_workload);
  double p_ratio(p_improv/copy_data.p_workload);

  if (refinement==h_refine)
    p_ratio = 0.;
  if (refinement==p_refine)
    h_ratio = 0.;

  double indicator_per_cell(0.);
  if (h_ratio > p_ratio)
    {
      convergence_est_per_cell[copy_data.cell_index] = h_improv;
      indicator_per_cell = est_per_cell[copy_data.cell_index];
      hp_Conv_Est[copy_data.cell_index] = indicator_per_cell;
      p_ref_map[copy_data.global_cell] = false;
    }
  else
    {
      convergence_est_per_cell[copy_data.cell_index] = p_improv;
      indicator_per_cell = est_per_cell[copy_data.cell_index];
      hp_Conv_Est[copy_data.cell_index] = indicator_per_cell;
      p_ref_map[copy_data.global_cell] = true;
    }

  primal_indicator_cell.push_back(std::make_pair(indicator_per_cell, copy_data.global_cell));
}


template <int dim>
void StokesProblem<dim>::compute_local_dual_ritz_residual(
    hp::DoFHandler<dim> &local_dual_dof_handler, 
    unsigned int level_p_refine,
    std::vector<unsigned int> const &block_component_patch,
    std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator>
      &patch_to_global_tria_map,
    BlockVector<double> &local_dual_solu,
    BlockVector<double> &dual_residual_solu)
{
  // Do a p-refinement of the local_dual_dof_handler and project the local
  // solution to the refined mesh.
  std::vector<types::global_dof_index> dofs_per_block_patch(2);
  DoFTools::count_dofs_per_block(local_dual_dof_handler, dofs_per_block_patch, block_component_patch);
  local_dual_solu.reinit(dofs_per_block_patch);

  // Here we are projecting the values of the global vector "dual_solution"
  // into vector "local_dual_solu" which is solution over patch cells corresponding
  // to cell "K".
  DoFHandler_active_cell_iterator act_patch_cell = local_dual_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dual_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      const unsigned int dofs_per_cell = act_patch_cell->get_fe().dofs_per_cell;
      // we check if the corresponding finite element for this cell is not 'FE_Nothing'
      // and it takes usual finite element.
      if (dofs_per_cell!=0)
        {
          Vector<double> local_dual_solution_values(dofs_per_cell);
          DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[act_patch_cell];

          global_cell->get_dof_values(dual_solution, local_dual_solution_values);
          act_patch_cell->set_dof_values(local_dual_solution_values, local_dual_solu);
        }
    }
  p_refinement(local_dual_dof_handler, patch_to_global_tria_map, level_p_refine, 
      local_dual_solu);

  // Build the constraint to force zero boundary condition on patch
  ConstraintMatrix constraints_patch;
  DoFTools::make_hanging_node_constraints(local_dual_dof_handler, constraints_patch);
  DoFHandler_active_cell_iterator dual_cell = local_dual_dof_handler.begin_active(),
                                  dual_endc = local_dual_dof_handler.end();
  for (; dual_cell!=dual_endc; ++dual_cell)
  {
    std::vector<types::global_dof_index> local_face_dof_indices((dual_cell->get_fe()).dofs_per_face);
    if (dual_cell->user_flag_set() == true)
    {
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        bool face_is_on_patch_bdry = false;
        if (dual_cell->face(f)->at_boundary())
          face_is_on_patch_bdry = true;
        else
        {
          if (dual_cell->neighbor(f)->has_children() == true)
          {
            for (unsigned int sf=0; sf<dual_cell->face(f)->n_children(); ++sf)
              if (dual_cell->neighbor_child_on_subface(f,sf)->user_flag_set() == false)
              {
                face_is_on_patch_bdry = true;
                break;
              }
          }
          else
          {
            if (dual_cell->neighbor(f)->user_flag_set() == false)
              face_is_on_patch_bdry = true;
          }
        }

        if (face_is_on_patch_bdry)
        {
          dual_cell->face(f)->get_dof_indices (local_face_dof_indices,
              dual_cell->active_fe_index());
          for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
            // system_to_component_index: "Compute vector component and index
            // of this shape function within the shape functions corresponding
            // to this component from the index of a shape function within this
            // finite element"
            if ((dual_cell->get_fe()).face_system_to_component_index(i).first < dim)
              constraints_patch.add_line(local_face_dof_indices[i]);
        }
      }
    }
  }
  constraints_patch.close();

  // After the p-refinement, dofs_per_block_patch has been changed.
  DoFTools::count_dofs_per_block(local_dual_dof_handler, dofs_per_block_patch, block_component_patch);
  BlockDynamicSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);
  BlockSparsityPattern sparsity_pattern_patch;
  DoFTools::make_sparsity_pattern (local_dual_dof_handler, csp, constraints_patch, false);
  sparsity_pattern_patch.copy_from(csp);

  BlockSparseMatrix<double> patch_dual_system(sparsity_pattern_patch);
  BlockVector<double> patch_dual_rhs(dofs_per_block_patch);
  dual_residual_solu.reinit(dofs_per_block_patch);
  
  // Assemble the system on the patch.
  patch_assemble_system(local_dual_dof_handler, constraints_patch, local_dual_solu,
      patch_dual_system, patch_dual_rhs, dual);

  // direct solver
  SparseDirectUMFPACK A_inverse_stiffness;
  A_inverse_stiffness.initialize (patch_dual_system.block(0,0),
                                  SparseDirectUMFPACK::AdditionalData());
  A_inverse_stiffness.vmult (dual_residual_solu.block(0), patch_dual_rhs.block(0));
  SparseDirectUMFPACK A_inverse_mass;
  A_inverse_mass.initialize (patch_dual_system.block(1,1),
                             SparseDirectUMFPACK::AdditionalData());
  A_inverse_mass.vmult(dual_residual_solu.block(1), patch_dual_rhs.block(1));
  constraints_patch.distribute(dual_residual_solu);
}


template <int dim>
double StokesProblem<dim>::squared_error_estimator_on_two_layer_patch(
    hp::DoFHandler<dim> const &local_dual_dof_handler,
    BlockVector<double> const &dual_solution,
    BlockVector<double> const &dual_residual_solution)
{
  double square_error_estimator = 0.;

  hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection, 
      update_values|update_gradients|update_hessians|update_quadrature_points
      |update_JxW_values);
  hp::FEFaceValues<dim> hp_fe_face_values(fe_collection, face_quadrature_collection,
      update_JxW_values|update_gradients|update_normal_vectors);
  hp::FEFaceValues<dim> hp_neighbor_face_values(fe_collection,
      face_quadrature_collection, update_gradients);
  hp::FESubfaceValues<dim> hp_subface_values(fe_collection, face_quadrature_collection,
      update_JxW_values|update_gradients|update_normal_vectors);
  hp::FESubfaceValues<dim> hp_neighbor_subface_values(fe_collection,
      face_quadrature_collection, update_gradients);

  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure (dim);

  // Use in the residual term
  std::vector<Vector<double>> dual_source_values;
  std::vector<Tensor<1,dim>> dual_solution_gradients_p;
  std::vector<double> dual_solution_divergences_v;
  std::vector<Tensor<1,dim>> dual_solution_laplacians_v;
  std::vector<double> dual_residual_values_p;
  std::vector<Tensor<1,dim>> dual_residual_laplacians_v;

  // Use in the jump term
  std::vector<Tensor<2,dim>> dual_solution_gradients_v;
  std::vector<Tensor<2,dim>> dual_solution_neighbor_gradients_v;
  std::vector<Tensor<2,dim>> dual_residual_gradients_v;
  std::vector<Tensor<2,dim>> dual_residual_neighbor_gradients_v;

  for (auto cell : local_dual_dof_handler.active_cell_iterators())
  {
    hp_fe_values.reinit(cell);
    FEValues<dim> const &fe_values = hp_fe_values.get_present_fe_values();
    std::vector<double> const &JxW_values = fe_values.get_JxW_values();
    unsigned int const n_q_points = fe_values.n_quadrature_points;

    dual_source_values.resize(n_q_points, Vector<double>(dim+1));
    dual_source->vector_value_list(fe_values.get_quadrature_points(), 
        dual_source_values);

    dual_solution_gradients_p.resize(n_q_points);
    dual_solution_divergences_v.resize(n_q_points);
    dual_solution_laplacians_v.resize(n_q_points);
    dual_residual_values_p.resize(n_q_points);
    dual_residual_laplacians_v.resize(n_q_points);

    fe_values[pressure].get_function_gradients(dual_solution, 
        dual_solution_gradients_p);
    fe_values[velocities].get_function_divergences(dual_solution,
        dual_solution_divergences_v);
    fe_values[velocities].get_function_laplacians(dual_solution,
        dual_solution_laplacians_v);
    fe_values[pressure].get_function_values(dual_residual_solution,
        dual_residual_values_p);
    fe_values[velocities].get_function_laplacians(dual_residual_solution,
        dual_residual_laplacians_v);
    
    // Compute the residual-based term.
    double residual_term = 0.;
    double const scaling = pow((cell->diameter())/(cell->get_fe().degree), 2.0);
    for (unsigned int q=0; q<n_q_points; ++q)
    {
      Tensor<1,dim> term_1;
      for (unsigned int i=0; i<dim; ++i)
        term_1[i] = dual_source_values[q][i] + dual_solution_laplacians_v[q][i] -
          dual_solution_gradients_p[q][i] + dual_residual_laplacians_v[q][i];

      residual_term += scaling * (scalar_product(term_1,term_1)*JxW_values[q]);

      residual_term += pow((dual_solution_divergences_v[q]-dual_residual_values_p[q]),
         2) * JxW_values[q]; 
    }

    // Compute the jump-based term.
    double jump_term = 0.;
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    {
      // Check that the face is not at the boundary
      if (cell->face(face)->at_boundary() == false)
      {
        // If the neighbor cell is active and on the same level
        if ((cell->neighbor(face)->active()) &&
            (cell->neighbor(face)->level() == cell->level()))
        {
          unsigned int const q_index = std::max(cell->active_fe_index(),
              cell->neighbor(face)->active_fe_index());

          hp_fe_face_values.reinit(cell, face, q_index);
          hp_neighbor_face_values.reinit(cell->neighbor(face),
              cell->neighbor_of_neighbor(face), q_index);

          FEFaceValues<dim> const &fe_face_values = 
            hp_fe_face_values.get_present_fe_values();
          FEFaceValues<dim> const &neighbor_face_values = 
            hp_neighbor_face_values.get_present_fe_values();
          std::vector<double> const &JxW_values = fe_face_values.get_JxW_values();
          unsigned int const n_face_q_points = fe_face_values.n_quadrature_points;

          dual_solution_gradients_v.resize(n_face_q_points);
          dual_solution_neighbor_gradients_v.resize(n_face_q_points);
          dual_residual_gradients_v.resize(n_face_q_points);
          dual_residual_neighbor_gradients_v.resize(n_face_q_points);

          fe_face_values[velocities].get_function_gradients(dual_solution,
              dual_solution_gradients_v);
          neighbor_face_values[velocities].get_function_gradients(dual_solution,
              dual_residual_neighbor_gradients_v);
          fe_face_values[velocities].get_function_gradients(dual_residual_solution,
              dual_residual_gradients_v);
          neighbor_face_values[velocities].get_function_gradients(dual_residual_solution,
              dual_residual_neighbor_gradients_v);
        
          std::vector<Tensor<1,dim>> jump_per_face;
          jump_per_face.resize(n_face_q_points);
          double jump_val =0.;
          for (unsigned int q=0; q<n_face_q_points; ++q)
          {
            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int j=0; j<dim; ++j)
                jump_per_face[q][i] += (dual_solution_gradients_v[q][i][j] +
                    dual_residual_gradients_v[q][i][j] - 
                    dual_solution_neighbor_gradients_v[q][i][j] -
                    dual_residual_neighbor_gradients_v[q][i][j]) *
                  fe_face_values.normal_vector(q)[j];
            jump_val += scalar_product(jump_per_face[q], jump_per_face[q]) *
              JxW_values[q];
          }
          unsigned int min_fe_degree = std::min(cell->get_fe().degree,
              cell->neighbor(face)->get_fe().degree);

          jump_term += (cell->face(face)->diameter()/(2.*min_fe_degree)) * 
            jump_val;
        }
      }
      // If the neighbor has children
      else if (cell->neighbor(face)->has_children() == true)
      {
        for (unsigned int subface=0; subface<cell->face(face)->n_children();
            ++subface)
        {
          unsigned int const q_index = std::max(cell->neighbor_child_on_subface(
                face, subface)->active_fe_index(), cell->active_fe_index());

          hp_neighbor_face_values.reinit(cell->neighbor_child_on_subface(
                face, subface), cell->neighbor_of_neighbor(face), q_index);
          hp_subface_values.reinit(cell, face, subface, q_index);

          FEFaceValues<dim> const &neighbor_face_values = 
            hp_neighbor_face_values.get_present_fe_values();
          FESubfaceValues<dim> const &fe_subface_values =
            hp_subface_values.get_present_fe_values();
          std::vector<double> const &JxW_values = fe_subface_values.get_JxW_values();
          unsigned int const n_subface_q_points = fe_subface_values.n_quadrature_points;

          dual_solution_gradients_v.resize(n_subface_q_points);
          dual_solution_neighbor_gradients_v.resize(n_subface_q_points);
          dual_residual_gradients_v.resize(n_subface_q_points);
          dual_residual_neighbor_gradients_v.resize(n_subface_q_points);

          fe_subface_values[velocities].get_function_gradients(dual_solution,
              dual_solution_gradients_v);
          neighbor_face_values[velocities].get_function_gradients(dual_solution,
              dual_solution_neighbor_gradients_v);
          fe_subface_values[velocities].get_function_gradients(dual_residual_solution,
              dual_residual_gradients_v);
          neighbor_face_values[velocities].get_function_gradients(dual_residual_solution,
              dual_residual_neighbor_gradients_v);

          std::vector<Tensor<1,dim>> jump_per_subface;
          jump_per_subface.resize(n_subface_q_points);

          double jump_val = 0.;
          for (unsigned int q=0; q<n_subface_q_points; ++q)
          {
            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int j=0; j<dim; ++j)
                jump_per_subface[q][j] += (dual_solution_gradients_v[q][i][j] +
                    dual_residual_gradients_v[q][i][j] - 
                    dual_solution_neighbor_gradients_v[q][i][j] -
                    dual_residual_neighbor_gradients_v[q][i][j]) *
                  fe_subface_values.normal_vector(q)[j];
            jump_val += scalar_product(jump_per_subface[q], jump_per_subface[q]) *
              JxW_values[q];
          }
          unsigned int min_fe_degree = std::min(cell->get_fe().degree, 
              cell->neighbor_child_on_subface(face, subface)->get_fe().degree);

          jump_term += (cell->face(face)->child(subface)->diameter()/(2*
                min_fe_degree)) * jump_val;
        }
      }
      // If the neighbor is coarser
      else
      {
        Assert(cell->neighbor_is_coarser(face),
            ExcMessage("Problem while computing the error estimator"));
        unsigned int const q_index = std::max(cell->active_fe_index(),
            cell->neighbor(face)->active_fe_index());
        hp_fe_face_values.reinit(cell, face, q_index);
        hp_neighbor_subface_values.reinit(cell->neighbor(face),
            cell->neighbor_of_coarser_neighbor(face).first,
            cell->neighbor_of_coarser_neighbor(face).second,
            q_index);

        FEFaceValues<dim> const &fe_face_values = 
          hp_fe_face_values.get_present_fe_values();
        FESubfaceValues<dim> const &neighbor_subface_values = 
          hp_neighbor_subface_values.get_present_fe_values();

        std::vector<double> const &JxW_values = fe_face_values.get_JxW_values();
        unsigned int const n_face_q_points = fe_face_values.n_quadrature_points;

        dual_solution_gradients_v.resize(n_face_q_points);
        dual_solution_neighbor_gradients_v.resize(n_face_q_points);
        dual_residual_gradients_v.resize(n_face_q_points);
        dual_residual_neighbor_gradients_v.resize(n_face_q_points);

        fe_face_values[velocities].get_function_gradients(dual_solution, 
            dual_solution_gradients_v);
        neighbor_subface_values[velocities].get_function_gradients(dual_solution,
            dual_solution_neighbor_gradients_v);
        fe_face_values[velocities].get_function_gradients(dual_residual_solution, 
            dual_residual_gradients_v);
        neighbor_subface_values[velocities].get_function_gradients(dual_residual_solution,
            dual_residual_neighbor_gradients_v);
                                 
        std::vector<Tensor<1,dim>> jump_per_face;
        jump_per_face.resize(n_face_q_points);
        double jump_val = 0;
        for (unsigned int q=0; q<n_face_q_points; ++q)
        {
          for (unsigned int i=0; i<dim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              jump_per_face[q][i] += (dual_solution_gradients_v[q][i][j] +
                    dual_residual_gradients_v[q][i][j] - 
                    dual_solution_neighbor_gradients_v[q][i][j] -
                    dual_residual_neighbor_gradients_v[q][i][j]) *
                fe_face_values.normal_vector(q)[j];
          jump_val += scalar_product(jump_per_face[q], jump_per_face[q]) * JxW_values[q];
        }

        unsigned int min_fe_degree = std::min(cell->get_fe().degree,
            cell->neighbor(face)->get_fe().degree);
        jump_term += (cell->face(face)->diameter()/(2.*min_fe_degree)) * jump_val;
      }
    }

    square_error_estimator += residual_term + jump_term;
  }

  return square_error_estimator;
}


template <int dim>
double StokesProblem<dim>::compute_local_go_error_estimator_square(
    hp::DoFHandler<dim> const &local_dual_dof_handler,
    BlockVector<double> const &local_dual_solution,
    BlockVector<double> const &dual_residual_solution)
{
  double go_error_estimator_square = 0.;
  
  return go_error_estimator_square;
}


template <int dim>
std::vector<std::pair<double, typename hp::DoFHandler<dim>::active_cell_iterator>>
StokesProblem<dim>::compute_goal_oriented_error_estimator()
{
  std::vector<std::pair<double, DoFHandler_active_cell_iterator>>
    go_error_estimator_square(triangulation.n_active_cells());
  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  endc = dof_handler.end();
  for (unsigned int i=0; cell!=endc; ++cell, ++i)
  {
    // Create the local triangulation associated to a patch. Unlike the primal
    // case, for the dual problem the patch is composed of two layers of cells
    // and not one.
    Triangulation<dim> local_triangulation;
    unsigned int level_h_refine(0);
    unsigned int level_p_refine(0);
    std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator> 
      patch_to_global_tria_map;
    std::vector<DoFHandler_active_cell_iterator> patch = get_patch_around_cell(cell, 2);
    build_triangulation_from_patch(patch, local_triangulation, level_h_refine, level_p_refine, 
        patch_to_global_tria_map);

    // Create the local dof_handler associated to the local triangulation and
    // the primal space. This is necessary so that the solution of the dual
    // problem can be refined to p+1 and then, projected to local_dual_dof_handler.
    hp::DoFHandler<dim> local_dual_dof_handler(local_triangulation);
    set_active_fe_indices(local_dual_dof_handler, patch_to_global_tria_map);
    local_dual_dof_handler.distribute_dofs(fe_collection);
    std::vector<unsigned int> block_component_patch (dim+1, 0);
    block_component_patch[dim]=1;
    DoFRenumbering::component_wise(local_dual_dof_handler, block_component_patch);

    // Compute the Ritz representation of dual residual on the patch:
    // local_dual_solution is projected on the p-refined space and
    // dual_residual_solution is computed.
    BlockVector<double> local_dual_solution, local_dual_residual_solution;
    compute_local_dual_ritz_residual(local_dual_dof_handler, level_p_refine, 
        block_component_patch, patch_to_global_tria_map,
        local_dual_solution, local_dual_residual_solution);

    // Compute the square of the a posteriori error estimator on the two-layer patch. 
    double squared_two_layer_error_estimator = squared_error_estimator_on_two_layer_patch(
        local_dual_dof_handler, local_dual_solution, local_dual_residual_solution);
    
    // Compute the square of the goal-oriented a posteriori error estimator
    go_error_estimator_square.push_back(std::pair<double, DoFHandler_active_cell_iterator>
        (compute_local_go_error_estimator_square(local_dual_dof_handler, 
                                                 local_dual_solution,
                                                 local_dual_residual_solution),
         cell));
  }

  return go_error_estimator_square;
}


template <int dim>
void StokesProblem<dim>::mark_cells_goal_oriented(const unsigned int cycle, 
    const double theta, std::vector<std::pair<double, DoFHandler_active_cell_iterator>> 
    const &go_error_estimator_square)
{
  primal_indicator_cell.clear();
  candidate_cell_set.clear();
  marked_cells.reinit(triangulation.n_active_cells());

  // this vector "convergence_est_per_cell" will be finalized after checking out
  // which h- or p- refinement are going to be chosen for each cell
  convergence_est_per_cell.reinit(triangulation.n_active_cells());

  h_Conv_Est.reinit(triangulation.n_active_cells());
  p_Conv_Est.reinit(triangulation.n_active_cells());
  hp_Conv_Est.reinit(triangulation.n_active_cells());

  // Create the synchronous iterators used by WorkStream
  std::vector<unsigned int> cell_indices(
      dof_handler.get_triangulation().n_active_cells(),0);
  for (unsigned int i=0; i<dof_handler.get_triangulation().n_active_cells(); ++i)
    cell_indices[i] = i;
  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  end_cell = dof_handler.end();
  SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
    std::vector<unsigned int>::iterator>> synch_iter(std::tuple<DoFHandler_active_cell_iterator,
        std::vector<unsigned int>::iterator> (cell, cell_indices.begin()));
  SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
    std::vector<unsigned int>::iterator>> end_synch_iter(std::tuple<DoFHandler_active_cell_iterator,
        std::vector<unsigned int>::iterator> (end_cell, cell_indices.end()));
  // Compute the convergence estimator and then mark the cell for h- or p-refinement
  WorkStream::run(synch_iter, end_synch_iter,
      std::bind(&StokesProblem<dim>::patch_convergence_estimator,this, cycle,
        std::placeholders::_1, std::placeholders::_2, std::placeholders::_3),
      std::bind(&StokesProblem<dim>::copy_to_refinement_maps,this,
        std::placeholders::_1),
      ScratchData(),CopyData<dim>());

  // Compute go_error_estimator, i.e. the square root of the product of the
  // square of the goal-oriented error estimator and the square of the primal
  // error estimator.
  std::vector<std::pair<double, DoFHandler_active_cell_iterator>> go_error_estimator(
      go_error_estimator_square);
  for (unsigned int i=0; i<triangulation.n_active_cells(); ++i)
  {
    go_error_estimator[i].first += std::pow(est_per_cell[i],2);
    go_error_estimator[i].first = std::sqrt(go_error_estimator[i].first);
  }

  // Sort the go_error_estimator from the largest to smallest value.
  std::sort(go_error_estimator.begin(), go_error_estimator.end(),
      std::bind(&StokesProblem<dim>::sort_decreasing_order, this,
        std::placeholders::_1, std::placeholders::_2));

  // Mark for refinement the cells with the largest go_error_estimator.
  double L2_norm(0.);
  for (auto const &error_estimator : go_error_estimator)
    L2_norm += std::pow(error_estimator.first, 2.);
  L2_norm = std::sqrt(L2_norm);
  double sum(0.);
  for (auto const &error_estimator : go_error_estimator)
  {
    sum += std::pow(error_estimator.first, 2);
    candidate_cell_set.push_back(error_estimator.second);
    // if theta is one, refine every cell
    if ((theta<1.) && (sum>=std::pow(theta*L2_norm,2)))
      break;
  }
  cell = dof_handler.begin_active();
  for (unsigned int i=0; cell!=end_cell; ++cell,++i)
  {
    typename std::vector<DoFHandler_active_cell_iterator>::iterator  mark_candidate;
    for (mark_candidate=candidate_cell_set.begin();
        mark_candidate!=candidate_cell_set.end(); ++mark_candidate)
      if (cell == *mark_candidate)
        marked_cells[i]=1;
  }
}


// This function controls the portion of cells for refinement for the primal
// problem.
// The parameter theta here plays an important role in this step
template <int dim>
void StokesProblem<dim>::mark_cells(const unsigned int cycle, const double theta)
{
  primal_indicator_cell.clear();
  candidate_cell_set.clear();
  marked_cells.reinit(triangulation.n_active_cells());

  // this vector "convergence_est_per_cell" will be finalized after checking out
  // which h- or p- refinement are going to be chosen for each cell
  convergence_est_per_cell.reinit(triangulation.n_active_cells());

  h_Conv_Est.reinit(triangulation.n_active_cells());
  p_Conv_Est.reinit(triangulation.n_active_cells());
  hp_Conv_Est.reinit(triangulation.n_active_cells());

  // Create the synchronous iterators used by WorkStream
  std::vector<unsigned int> cell_indices(
      dof_handler.get_triangulation().n_active_cells(),0);
  for (unsigned int i=0; i<dof_handler.get_triangulation().n_active_cells(); ++i)
    cell_indices[i] = i;
  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  end_cell = dof_handler.end();
  SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
    std::vector<unsigned int>::iterator>> synch_iter(std::tuple<DoFHandler_active_cell_iterator,
        std::vector<unsigned int>::iterator> (cell, cell_indices.begin()));
  SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
    std::vector<unsigned int>::iterator>> end_synch_iter(std::tuple<DoFHandler_active_cell_iterator,
        std::vector<unsigned int>::iterator> (end_cell, cell_indices.end()));
  WorkStream::run(synch_iter, end_synch_iter,
      std_cxx11::bind(&StokesProblem<dim>::patch_convergence_estimator,this, cycle,
        std_cxx11::_1, std_cxx11::_2, std_cxx11::_3),
      std_cxx11::bind(&StokesProblem<dim>::copy_to_refinement_maps,this,std_cxx11::_1),
      ScratchData(),CopyData<dim>());

  std::sort(primal_indicator_cell.begin(), primal_indicator_cell.end(),
      std::bind(&StokesProblem<dim>::sort_decreasing_order, this,
        std::placeholders::_1, std::placeholders::_2));

  double const L2_norm(est_per_cell.l2_norm());
  double sum(0.);
  for (unsigned int i=0; i<primal_indicator_cell.size(); ++i)
    {
      DoFHandler_active_cell_iterator cell_sort = primal_indicator_cell[i].second;
      sum += (primal_indicator_cell[i].first)*(primal_indicator_cell[i].first);

      candidate_cell_set.push_back(cell_sort);
      // if theta is one, refine every cell
      if ((theta<1.0) && (sum>=std::pow(theta*L2_norm,2)))
        break;
    }

  cell = dof_handler.begin_active();
  for (unsigned int i=0; cell!=end_cell; ++cell,++i)
    {
      typename std::vector<DoFHandler_active_cell_iterator>::iterator  mark_candidate;
      for (mark_candidate=candidate_cell_set.begin();
           mark_candidate!=candidate_cell_set.end(); ++mark_candidate)
        if (cell == *mark_candidate)
          marked_cells[i]=1;
    }
}


template <int dim>
void StokesProblem<dim>::output_results (const unsigned int cycle)
{
  Vector<float> fe_degrees(triangulation.n_active_cells());
  {
    DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                    endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      fe_degrees[index] = fe_collection[cell->active_fe_index()].degree;
  }


  std::vector<std::string> solution_names;
  solution_names.push_back ("x_velocity");
  solution_names.push_back ("y_velocity");
  if (dim==3)
    solution_names.push_back ("z_velocity");
  solution_names.push_back ("pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  if (dim==3)
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);


  DataOut<dim,hp::DoFHandler<dim> > data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector(solution, solution_names,
                           DataOut<dim,hp::DoFHandler<dim>>::type_dof_data, data_component_interpretation);

  data_out.add_data_vector (marked_cells, "marked_cells");
  data_out.add_data_vector (fe_degrees, "fe_degree");
  data_out.add_data_vector (est_per_cell, "Error_Estimator");
  data_out.add_data_vector (error_per_cell, "Error");
  data_out.add_data_vector (Vect_Pressure_Err, "Pressure_Error");
  data_out.add_data_vector (Vect_grad_Velocity_Err, "Grad_Velocity_Error");

  data_out.add_data_vector (h_Conv_Est, "h_refine_Conv_Est");
  data_out.add_data_vector (p_Conv_Est, "p_refine_Conv_Est");
  data_out.add_data_vector (hp_Conv_Est, "hp_refine_Conv_Est");

  data_out.build_patches();
  std::string filename = "solution-" +
                         Utilities::int_to_string (cycle, 2) +".vtu";
  std::ofstream output(filename.c_str());
  data_out.write_vtu (output);
}

// This function goes through all cells and tries to do h- or p- refinement
// on all candidates cells which already have been marked for h- or p-refinment
template <int dim>
void StokesProblem<dim>::refine_in_h_p()
{
  bool need_to_h_refine=false;
  typename std::vector<DoFHandler_active_cell_iterator>::iterator  cell_candidate;
  for (cell_candidate=candidate_cell_set.begin(); cell_candidate!=candidate_cell_set.end();
       ++cell_candidate)
    {
      // Mark the cell for h-refinement
      if (p_ref_map[*cell_candidate]==false)
        {
          need_to_h_refine = true;
          (*cell_candidate)->set_refine_flag();
        }
      // Perform the p-refinement, i.e. increase the fe index of the cell
      else
      {
        // Increase the order of the FE on the cell if we haven't 
        // reached the maximum polynomial order.
        if (((*cell_candidate)->active_fe_index()+1) < (fe_collection.size()-1))
        {
          (*cell_candidate)->set_active_fe_index((*cell_candidate)->active_fe_index()+1);
        }
        // Mark the cell for h-refinement if they were marked for p-refinement
        // but we reached the maximum polynomial order.
        else
        {
          if ((p_ref_map[*cell_candidate]==true) && 
              ((*cell_candidate)->active_fe_index()+1) >= (fe_collection.size()-1))
          {
            need_to_h_refine = true;
            (*cell_candidate)->set_refine_flag();
          }
          else
            Assert (false, ExcNotImplemented()); // should never get here!
        }
      }
    }

  // Clear the p-refinement map
  p_ref_map.clear();

  // Perform the h-refinement
  if (need_to_h_refine==true)
    triangulation.execute_coarsening_and_refinement();

  // Smooth the p-refinement so that we don't get polynomial orders too
  // different in adjacent cells.
  bool cell_changed=false;
  do
    {
      cell_changed=false;
      unsigned int count_cell=0;
      DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();
      for (; cell!=endc; ++cell,++count_cell)
        {
          std::vector<DoFHandler_active_cell_iterator> patch_cells = 
            get_patch_around_cell(cell, 1);

          for (unsigned int i=1; i<patch_cells.size(); ++i)
            {

              // Increase the order of the cells around patch_cells[0] if their
              // order is less than patch_cells[0] order - 1.
              if (patch_cells[i]->active_fe_index()+1 < (patch_cells[0]->active_fe_index()))
                {
                  patch_cells[i]->set_active_fe_index(patch_cells[0]->active_fe_index()-1);
                  cell_changed=true;
                }
              // Decrease the order of the cells around patch_cells[0] if their
              // order is greater than patch_cells[0] order + 1.
              else if (patch_cells[i]->active_fe_index() > (patch_cells[0]->active_fe_index()+1))
                {
                  patch_cells[0]-> set_active_fe_index (patch_cells[i]->active_fe_index()-1);
                  cell_changed=true;
                }
            }
        }

    }
  while (cell_changed==true);
}



//Explicit initialization
template class StokesProblem<2>;
template class StokesProblem<3>;
