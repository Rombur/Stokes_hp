#ifndef _COPYDATA_HH_
#define _COPYDATA_HH_

#include <deal.II/hp/dof_handler.h>

template <int dim>
struct CopyData
{
  CopyData() {};

  CopyData(CopyData<dim> const &data);

  typename dealii::hp::DoFHandler<dim>::active_cell_iterator global_cell;
  unsigned int cell_index;
  double h_conv_est_per_cell;
  double h_workload;
  double p_conv_est_per_cell;
  double p_workload;
};


template <int dim>
CopyData<dim>::CopyData(CopyData<dim> const &data)
  :
  global_cell(data.global_cell),
  cell_index(data.cell_index),
  h_conv_est_per_cell(data.h_conv_est_per_cell),
  h_workload(data.h_workload),
  p_conv_est_per_cell(data.p_conv_est_per_cell),
  p_workload(data.p_workload)
{}

#endif
