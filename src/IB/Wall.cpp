// Filename: Wall.cpp
// Created by Steven DeLong
//
// Copyright (c) 2002-2015, Steven DeLong and Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/Wall.h"
#include "ibamr/namespaces.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

Wall::Wall(Pointer<Database> wall_db, Pointer<CartesianGridGeometry<NDIM> > grid_geometry, double wall_ghost_dist)
{
    // get geometry for the wall force box and scaling from cells to numerical postions
    const Box<NDIM> domain_box = grid_geometry->getPhysicalDomain()[0];
    const double* dx = grid_geometry->getDx();

    // get direction the wall is facing, dim
    if (wall_db->keyExists("dim"))
    {
        d_axis = wall_db->getInteger("dim");
    }
    else
    {
        TBOX_ERROR("walls must have a dim parameter");
    }

    // get the location of the wall
    if (wall_db->keyExists("location"))
    {
        // location gives the numerical location of the wall, but is read in
        // as a cell number for box creation.  The convention is that the cell
        // indicated is the first cell the wall applies forces to
        // e.g. if location is 1, and side = 1, the wall will be at location h = cell width.
        //      if location is 1, and side = 0, the wall will be at location 2h,
        //      where h = cell width
        d_location = wall_db->getInteger("location");
    }
    else
    {
        TBOX_ERROR("walls must have a location parameter");
    }

    // which way is the wall 'pushing'.  0 - pushes negative direction
    //                                   1 - pushes positive direction
    if (wall_db->keyExists("side"))
    {
        d_side = wall_db->getInteger("side");
    }
    else
    {
        TBOX_ERROR("walls must have a side parameter");
    }

    // how far to search for particles to apply forces to.  This is padded by lagrangian ghost cell width.
    if (wall_db->keyExists("force_distance"))
    {
        d_force_distance = wall_db->getDouble("force_distance");
    }
    else
    {
        TBOX_ERROR("walls must have a force_distance parameter");
    }

    Box<NDIM> force_box = domain_box;

    if (d_side == 0)
    {
        force_box.upper(d_axis) = d_location + 1;
        force_box.lower(d_axis) = d_location + 1 - d_force_distance - wall_ghost_dist;
    }
    else if (d_side == 1)
    {
        force_box.lower(d_axis) = d_location;
        force_box.upper(d_axis) = d_location + d_force_distance + wall_ghost_dist;
    }
    else
    {
        TBOX_ERROR("Side must be 0 or 1");
    }

    d_location = (d_location + 1 - d_side) * dx[d_axis];

    d_force_area = force_box;
    d_force_distance = d_force_distance * dx[d_axis];

    // get parameters for wall force
    d_parameters = wall_db->getDoubleArray("parameters");
}

Wall::Wall(const Wall& other_wall)
    : d_axis(other_wall.d_axis),
      d_side(other_wall.d_side),
      d_location(other_wall.d_location),
      d_force_distance(other_wall.d_force_distance),
      d_parameters(other_wall.d_parameters),
      d_force_area(other_wall.d_force_area)
{
    // this body intentionally blank
    return;
}

Wall& Wall::operator=(const Wall& other_wall)
{
    if (&other_wall != this)
    {
        d_axis = other_wall.d_axis;
        d_side = other_wall.d_side;
        d_location = other_wall.d_location;
        d_force_distance = other_wall.d_force_distance;
        d_force_area = other_wall.d_force_area;
        d_parameters = other_wall.d_parameters;
    }
    return *this;
}

void
Wall::registerWallForceFcn(WallForceFcnPtr wall_force_fcn)
{
    // set the local d_wall_force_fcn to the provided
    // function pointer
    d_wall_force_fcn = wall_force_fcn;
    return;
} // registerWallForceFcn

Box<NDIM>
Wall::getForceArea()
{
    return d_force_area;
} // getForceArea

int
Wall::getAxis()
{
    return d_axis;
} // getAxis

int
Wall::getLocation()
{
    return d_location;
} // getLocation

double
Wall::applyForce(double wall_distance, double /*eval_time*/)
{
    // apply forces from this wall.  for now time dependence is not implemented.
    int sgn = (wall_distance > 0.0) ? 1 : ((wall_distance < 0.0) ? -1 : 0);
    if (sgn == 0)
    {
        TBOX_ERROR("Particle is on the wall, must be inside the domain");
        return 0.0;
    }
    else
    {
        // make sure we don't apply forces to particles that moved out of the
        // wall area, but haven't been regridded yet.
        if (std::abs(wall_distance) < d_force_distance)
        {
            return sgn * (d_wall_force_fcn)(std::abs(wall_distance), d_parameters);
        }
        else
        {
            return 0.0;
        }
    }
} // applyForce

/////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
