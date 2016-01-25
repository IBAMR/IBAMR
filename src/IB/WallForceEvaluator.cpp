#include "WallForceEvaluator.h"

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>


WallForceEvaluator::WallForceEvaluator(Pointer<Database> input_db,
                                       Pointer<CartesianGridGeometry<NDIM> > grid_geometry) :
    d_grid_geometry(grid_geometry)
{
   // register walls, call this during IB Hierarchy Integrator constructor

    // get number of walls:
    int n_walls = 0;
    if (input_db->keyExists("wall_number"))
    {
        n_walls = input_db->getInteger("wall_number");
    }
    else
    {
        pout << "WARNING: wall_number does not exist in IBHierarchyIntegrator database, defaulting to 0 walls." << std::endl;
        n_walls = 0;
    }
    
    // set up padding distance to track particles that might come within
    // the force area soon.  Error if this is not specified in the file.
    double wall_ghost_dist = -1.0; 
    if(input_db->keyExists("regrid_alpha"))
    {
        wall_ghost_dist = 2.0*(input_db->getDouble("regrid_alpha"));
    }
    else
    {
        TBOX_ERROR("ERROR: regrid_alpha not specified in WallForceEvaluator database.");
    }

    //wall parameters with defaults (no walls)
    Pointer<Database> wall_db;

    string wall_str, num_str;
    for(int wall = 0; wall < n_walls; ++ wall)
    {
        stringstream ss;
        ss << wall;
        num_str = ss.str();
        wall_str = "Wall_"+num_str;
        
        wall_db = input_db->getDatabase(wall_str);

        addWall(wall_db,wall_ghost_dist);
    }
    
}// Constructor


void
WallForceEvaluator::registerWallForceFcn(WallForceFcnPtr wall_force_fcn)
{
    // register the given force function pointer with all walls
    for(std::vector<Wall>::iterator wall_it = d_walls_vec.begin();
        wall_it != d_walls_vec.end();
        ++wall_it)
    {// iterate through walls
        wall_it->registerWallForceFcn(wall_force_fcn);
    }
    
    return;
} // registerWallForceFcn

void
WallForceEvaluator::addWall(Pointer<Database>  wall_db, double wall_ghost_dist)
{
    
    // Create the new wall from the wall_db
    Wall new_wall(wall_db,d_grid_geometry, wall_ghost_dist);

    // add wall to vector of present walls
    d_walls_vec.push_back(new_wall);
    
    return;

} // addWall

void
WallForceEvaluator::computeLagrangianForce(Pointer<LData> F_data,
                                           Pointer<LData> X_data,
                                           Pointer<LData> /*U_data*/,
                                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           const int level_number,
                                           const double eval_time,
                                           LDataManager* const l_data_manager)
{


    PetscScalar *force;
    VecGetArray(F_data->getVec(),&force);
    
    // get access to petsc vector of positions.
    PetscScalar *lposition;
    VecGetArray(X_data->getVec(),&lposition);


    const int lag_node_idx_current_idx = l_data_manager->getLNodePatchDescriptorIndex();
    
    for(std::vector<Wall>::iterator wall_it = d_walls_vec.begin();
        wall_it != d_walls_vec.end();
        ++wall_it)
    {// iterate through walls

        // get axis (which direction the wall is normal to)
        int axis = wall_it->getAxis();
        
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {//iterate through patches
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<LNodeSetData> current_idx_data = patch->getPatchData(lag_node_idx_current_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            // get just the area near the wall
            const Box<NDIM> intersect_box = patch_box * (wall_it->getForceArea());
            
            // iterate through cells in relevant area
            for (LNodeSetData::CellIterator scit(intersect_box); scit; scit++)
            {
                // get current nodes in the cell.
                const Index<NDIM>& search_cell_idx = *scit;
                LNodeSet* search_node_set = current_idx_data->getItem(search_cell_idx);

                // if particles exist in this cell, add forces to them
                if(search_node_set)
                {
                    for (LNodeSet::iterator it = search_node_set->begin(); it != search_node_set->end(); ++it)
                    {
                        // get node data
                        LNodeSet::value_type& particle_node_idx = *it;
                        const int particle_petsc_idx = particle_node_idx->getLocalPETScIndex();

                        // get wall distance
                        double wall_distance = lposition[particle_petsc_idx*NDIM + axis] - wall_it->getLocation();
                        force[particle_petsc_idx*NDIM + axis] += wall_it->applyForce(wall_distance,eval_time);
                    } // iterate through nodes
                } // if search_node_set
            } // iterate through cells in wall area
        }// iterate through patches
    } // iterate through walls

    VecRestoreArray(F_data->getVec(),&force);
    
    return;

} //computeLagrangianForce


