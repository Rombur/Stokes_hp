 DoFTools::extract_subset (const std::map<....> &map_global_to_local,ConstraintMatrix &subset)
  {
                    for (unsigned int i=0; i<lines.size(); ++i)
                      if (map_global_to_local.find(lines[i].line) != map_g_t_l.end())
                        {                      
                          const unsigned int local_dof_index = map_global_to_local.find(lines[i].line)->second;
                          bool constraints_all_local = true;
                          for (unsigned int j=0; j<lines[i].entries.size(); ++j)
                            if (map_global_to_local.find(lines[i].entries[j].first) == map_g_t_l.end())
                              {
                                constraits_all_local = false;
                                break;
                              }
                        
                          if (constraints_all_local == true)
                            {
                              subset.add_line (local_dof_index);
                              for (unsigned int j=0; j<lines[i].entries.size(); ++j)
                                subset.add_entry (local_dof_index,
                                              map_global_to_local.find(lines[i].entries[j].first)->second, 
// what if this does not exist?
                                                  lines[i].entries[j].second);
                              subset.add_inhomogeneity (local_dof_index,
                                                        lines[i].inhomogeneity);
                            }
                        }
                  }
                

