// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "LevelSetInitialCondition.h"

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <SAMRAI_config.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialCondition::LevelSetInitialCondition(const std::string& object_name,
                                                   const Pointer<CartesianGridGeometryNd> grid_geom,
                                                   const double radius,
                                                   const IBTK::VectorNd& origin,
                                                   const bool fluid_is_interior_to_cylinder)
    : d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_radius(radius),
      d_origin(origin),
      d_fluid_is_interior_to_cylinder(fluid_is_interior_to_cylinder)
{
    // intentionally blank
    return;
} // LevelSetInitialCondition

bool
LevelSetInitialCondition::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LevelSetInitialCondition::setDataOnPatch(const int data_idx,
                                         Pointer<VariableNd> /*var*/,
                                         Pointer<PatchNd> patch,
                                         const double /*data_time*/,
                                         const bool initial_time,
                                         Pointer<PatchLevelNd> patch_level)
{
    // Set the level set function throughout the domain
    if (initial_time)
    {
        const BoxNd& patch_box = patch->getBox();
        Pointer<CellDataNd<double> > D_data = patch->getPatchData(data_idx);

        Pointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double* const grid_x_lower = d_grid_geom->getXLower();
        IntVectorNd ratio = patch_level->getRatio();
        const SAMRAI::hier::BoxNd domain_box = SAMRAI::hier::BoxNd::refine(d_grid_geom->getPhysicalDomain()[0], ratio);
        const hier::IndexNd& grid_lower_idx = domain_box.lower();

        for (BoxNd::Iterator it(patch_box); it; it++)
        {
            CellIndexNd ci(it());

            // Get physical coordinates
            IBTK::Vector coord = IBTK::Vector::Zero();
            for (int d = 0; d < NDIM; ++d)
                coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - grid_lower_idx(d)) + 0.5);

            double distance = std::numeric_limits<double>::quiet_NaN();
#if (NDIM == 2)
            distance = std::sqrt((coord[0] - d_origin[0]) * (coord[0] - d_origin[0]) +
                                 (coord[1] - d_origin[1]) * (coord[1] - d_origin[1]));
#endif
#if (NDIM == 3)
            distance = std::sqrt((coord[0] - d_origin[0]) * (coord[0] - d_origin[0]) +  ( coord[1] - d_origin[1]) * (coord[1] - d_origin[1]) + (coord[2] - d_origin[2]) * (coord[2] - d_origin[2]);
#endif
            double sign = 0.0;
            if (d_fluid_is_interior_to_cylinder) sign =  distance < d_radius ? 1.0 : -1.0;
            else sign =  distance < d_radius ? -1.0 : 1.0;
            (*D_data)(ci) = sign * std::abs(distance - d_radius);
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
