// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_Wall
#define included_IBAMR_Wall

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/IBLagrangianForceStrategy.h>

#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIArray.h>
#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAIDatabase.h>
#include <SAMRAIPointer.h>

namespace SAMRAI
{
namespace geom
{
template <int DIM>
class CartesianGridGeometry;
} // namespace geom
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
// Simple wall class to store information about a wall specified in the input file.
//
// This wall only interacts directly with particles in a specified distance from it,
// and only applies forces to particles on one side (a 2-sided wall can be made by
// specifying 2 separate walls, one with side = 0, one with side = 1).
class Wall
{
public:
    // typedef for Wall Force Function Pointer
    using WallForceFcnPtr = double (*)(double D, const SAMRAIArray<double> params);

    // Constructor.
    Wall(SAMRAIPointer<SAMRAIDatabase> wall_db,
         SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry,
         double wall_ghost_dist);

    // Copy constructor.
    Wall(const Wall& other_wall);

    // Assignment operator.
    Wall& operator=(const Wall& other_wall);

    // Register the function for evaluating wall forces.
    void registerWallForceFcn(WallForceFcnPtr wall_force_fcn);

    // Get the force area box of the wall.
    SAMRAIBox getForceArea();

    // Function to get the wall axis.
    int getAxis();

    // Function to get the wall location.
    int getLocation();

    // Apply force based on distance
    double applyForce(double wall_distance, double eval_time);

private:
    // Information about where the wall is, which way it points, and
    // how far it affects particles:
    int d_axis, d_side;
    int d_location;
    double d_force_distance;

    // Parameters used in the wall force function:
    SAMRAIArray<double> d_parameters;

    // The area of the domain that the wall influences:
    SAMRAIBox d_force_area;

    // Function for applying wall force:
    WallForceFcnPtr d_wall_force_fcn;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_Wall
