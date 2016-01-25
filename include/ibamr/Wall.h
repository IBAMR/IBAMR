// Filename: Wall.h
// Created by Steven Delong
//
// Copyright (c) 2002-2014, Steven Delong, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_Wall
#define included_Wall

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "tbox/Array.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibamr/IBLagrangianForceStrategy.h"
#include "muParser.h"

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
    typedef double (*WallForceFcnPtr)(double D, const SAMRAI::tbox::Array<double> params);

    // Constructor.
    Wall(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> wall_db,
         SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry,
         double wall_ghost_dist);

    // Copy constructor.
    Wall(const Wall& other_wall);

    // Assignment operator.
    Wall& operator=(const Wall& other_wall);

    // Register the function for evaluating wall forces.
    void registerWallForceFcn(WallForceFcnPtr wall_force_fcn);

    // Get the force area box of the wall.
    SAMRAI::hier::Box<NDIM> getForceArea();

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
    SAMRAI::tbox::Array<double> d_parameters;

    // The area of the domain that the wall influences:
    SAMRAI::hier::Box<NDIM> d_force_area;

    // Function for applying wall force:
    WallForceFcnPtr d_wall_force_fcn;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_Wall
