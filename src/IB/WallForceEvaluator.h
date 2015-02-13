#ifndef included_WallForceEvaluator
#define included_WallForceEvaluator

// IBAMR INCLUDES
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LNodeSetData.h"
#include "ibamr/IBSpringForceSpec.h"
#include <ibtk/IndexUtilities.h>
#include <ibtk/AppInitializer.h>
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/app_namespaces.h>
#include <ibamr/IBLagrangianForceStrategy.h>
#include "Wall.h"

/////////////////////////////// WallForceEvaluator Class /////////////////////////


class WallForceEvaluator : public IBLagrangianForceStrategy
{

  public:

    // constructor
    WallForceEvaluator(Pointer<Database> input_db,Pointer<CartesianGridGeometry<NDIM> > grid_geometry);

    // register a WallForceFcnPtr to be used by all walls
    void
    registerWallForceFcn(WallForceFcnPtr wall_force_fcn);
    
    // add a wall to the WallForceEvaluator
    void
    addWall(Pointer<Database>  wall_db, double wall_ghost_dist);

    // compute forces from all walls
    void
    computeLagrangianForce(Pointer<LData> F_data,
                           Pointer<LData> X_data,
                           Pointer<LData> /*U_data*/,
                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                           const int level_number,
                           const double /*data_time*/,
                           LDataManager* const l_data_manager);

  private:

    // vector of walls.
    std::vector<Wall> d_walls_vec;
        
    // grid geometry, used when making walls
    Pointer<CartesianGridGeometry<NDIM> > d_grid_geometry;


};// class WallForceEvaluator


#endif
