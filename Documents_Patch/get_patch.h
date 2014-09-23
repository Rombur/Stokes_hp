
/**
*The type of the Container template argument would be <hp::DoFHandler<2>> . 
*
*This Function, returns a list of the all neighbor cells of the given cell.
*Obviuously the first element of the list inside patch, would be the current cell. After pushing back the cell into vector patch, It loops over all faces of that given cell and checks if that face is not on the *boundary of the domain. Then, one considers two cases, if the neighbor cell does not have any children then the neighbor cell which shares the same face_number with our current cell would be add up to the list *of cells. In the second case, when the neighbor cell has already been refined and therefore ahs children, then it loops over all subfaces of current face and similarly pushing all those children back to the *vector nemed patch.
*/




template<class Container >
std::vector< typename Container::active_cell_iterator> get_patch_around_cell(const typename Container::active_cell_iterator &cell)   
{
   std::vector< typename Container::active_cell_iterator> patch;  
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

//**********************************************************************************************************************************************************//

/**
*The type of the Container template argument would be <hp::DoFHandler<2>> . 
*It Counts how many degrees of freedom out of the total number belong to each patch
*The input argument is the patch corresponding to the given cell and the output would be an insigned integer which is the number of DoFs included in patch.
*
*/

template<class Container >
unsigned int count_dofs_on_patch(const std::vector< typename Container::active_cell_iterator> &patch)
{
   std::set<types::global_dof_index> dofs_on_patch;
   std::vector<types::global_dof_index> local_dof_indices;
   typename Container::active_cell_iterator  cell=patch.begin();
   const typename Container::active_cell_iterator endc=patch.end ();
   for(; cell!=endc ; ++cell)
     {
       
       local_dof_indices.resize (cell->get_fe().dofs_per_cell);
       cell->get_dof_indices (local_dof_indices);
       dofs_on_patch.insert (local_dof_indices.begin(),local_dof_indices.end());
     }
   return dofs_on_patch.size();  
}

//**********************************************************************************************************************************************************//


/**
*The type of the Container template argument would be <hp::DoFHandler<2>> . 
*
*This function gets a patch corresponding to a given cell as input argument and create a mapping from global DoFs index to patch indices. 
*This mapping is used when we want to solve a system of equation just on a patch. In that case, we assemble local system and local rhs on each cell included in corresponding patch, and then using the indeces 
out of this mapping, we can easily distribute local system and rhs to the patch system and patch rhs.
*/

template<class Container >
std::map<types::global_dof_index,unsigned int> map_global_dofs_to_patch_indices (const std::vector<typename Container::active_cell_iterator> &patch)
{
 std::map<types::global_dof_index,unsigned int> dofs_mapping;
  std::vector<types::global_dof_index> local_dof_indices;

   unsigned int next_unused_patch_dof_index = 0;
   typename Container::active_cell_iterator  cell=patch.begin();
   const typename Container::active_cell_iterator endc=patch.end ();
   for(; cell!=endc ; ++cell)
     {
       local_dof_indices.resize (cell->get_fe().dofs_per_cell);
       cell->get_dof_indices(local_dof_indices);

       for (unsigned int i=0; i < cell->get_fe().dofs_per_cell ; ++i)
        if (dofs_mapping.find(local_dof_indices[i]) == dofs_mapping.end())
          {
            dofs_mapping.insert(make_pair(local_dof_indices[i], next_unused_patch_dof_index));
            ++next_unused_patch_dof_index;
          }
     }
   return dof_mapping;
}

//**********************************************************************************************************************************************************//


