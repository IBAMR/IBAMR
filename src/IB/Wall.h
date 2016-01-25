#ifndef included_Wall
#define included_Wall

#include "tbox/Array.h"
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

// IBTK THIRD-PARTY INCLUDES
#include <muParser.h>


// typedef for Wall Force Function Pointer
typedef double (*WallForceFcnPtr)(double D,
                                  const SAMRAI::tbox::Array<double> params);


/////////////////////////////// Wall Class //////////////////////////////////

class Wall 
{
    // wall class to store information about a wall specified in the input file.
    // this wall only interacts directly with particles in a specified distance from it,
    // and only applies forces to particles on one side (a 2-sided wall can be made by
    // specifying 2 separate walls, one with side = 0, one with side = 1).
    
  public:
    
    Wall(Pointer<Database> wall_db,
         Pointer<CartesianGridGeometry<NDIM> > grid_geometry,
         double wall_ghost_dist);

    // copy constructor
    Wall(const Wall& other_wall);

    // copy assignment operator
    Wall& operator=(const Wall& other_wall);
    
    // register the inline function for wall forces.
    void
    registerWallForceFcn(WallForceFcnPtr wall_force_fcn);

    // function to get the force area box of the wall
    Box<NDIM>
    getForceArea(void);

    // function to get axis
    int
    getAxis(void);

    // function to get location
    int
    getLocation(void);

    // function that applies force based on distance
    double
    applyForce(double wall_distance, double eval_time); 
    
    
  private:

    // information about where the wall is, which way it points, and
    // how far it affects particles
    int d_axis, d_side;
    int d_location;
    double d_force_distance;

    // double* of parameters used in the wall force function
    SAMRAI::tbox::Array<double> d_parameters;
    
    // just need a box here, to intersect.  This is the area of the domain that
    // the wall influences.
    Box<NDIM> d_force_area;

    // inline function for applying wall force
    WallForceFcnPtr d_wall_force_fcn;

};

#endif
