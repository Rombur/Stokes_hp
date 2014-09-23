vector<hp::DoFHandler<2>::active_cell_iterator>
get_cells_at_coarsest_common_level (const vector<hp::DoFHandler<2>::active_cell_iterator>  &cells)
{
  Assert (cells.size() > 0, ...);

  unsigned int min_level = cells[0]->level();
  unsigned int max_level = cells[0]->level();
for (cell in cells)
  {
    min_level = std::min (min_level, cell->level());
    max_level = std::max (max_level, cell->level());
  }


if (min_level == max_level)
  return cells;
 else
   {
     set<hp::DoFHandler<2>::active_cell_iterator>  uniform_cells;
     for (cell in cells)
       {
         if (cell->level() == min_level)
           uniform_cells.push_back (cell);
         else
           {
             hp::DoFHandler<2>::active_cell_iterator parent = cell;
             while (parent->level() > min_level)
               parent = parent->parent();

             // add parent cell; if it's already in the set, this does nothing
             uniform_cells.insert (parent);
           }
       }

     return vector<hp::DoFHandler<2>::active_cell_iterator> (uniform_cells.begin(),
                                                             uniform_cells.end());
   }
}




void build_triangulation_from_patch (const vector<hp::DoFHandler<2>::active_cell_iterator>  &patch_cells,
                                     Triangulation<dim> &triangulation)
{
  const vector<hp::DoFHandler<2>::active_cell_iterator>
    uniform_cells = get_cells_at_coarsest_common_level (cells);

  triangulation.clear();
  ...as before, just use uniform_cells...;
  triangulation.create_triangulation(vertices,cells,SubCellData());

  bool need_to_refine = false;
  do
    {
      need_to_refine = false;
      for (cell=triangulation.begin_active()...end())
        if (cell is a parent of any of the cells in 'patch_cells')
          {
            cell->set_refine_flag();
            need_to_refine = true;
          }

      if (need_to_refine == true)
        triangulation.execute_coarsening_and_refinement ();
    }
  while (need_to_refine == true);
}


