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

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialCondition::LevelSetInitialCondition(const std::string& object_name,
                                                   const Pointer<CartesianGridGeometry<NDIM> > grid_geom,
                                                   const IBTK::VectorNd& interface_loc,
                                                   const bool left_side)
    : d_object_name(object_name), d_grid_geom(grid_geom), d_interface_loc(interface_loc), d_left_side(left_side)
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
                                         Pointer<Variable<NDIM> > /*var*/,
                                         Pointer<Patch<NDIM> > patch,
                                         const double /*data_time*/,
                                         const bool initial_time,
                                         Pointer<PatchLevel<NDIM> > patch_level)
{
    // Set the level set function throughout the domain
    if (initial_time)
    {
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);

        Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double* const grid_x_lower = d_grid_geom->getXLower();
        IntVector<NDIM> ratio = patch_level->getRatio();
        const SAMRAI::hier::Box<NDIM> domain_box =
            SAMRAI::hier::Box<NDIM>::refine(d_grid_geom->getPhysicalDomain()[0], ratio);
        const hier::Index<NDIM>& grid_lower_idx = domain_box.lower();

        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());

            // Get physical coordinates
            IBTK::Vector coord = IBTK::Vector::Zero();
            for (int d = 0; d < NDIM; ++d)
                coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - grid_lower_idx(d)) + 0.5);

            double distance = std::numeric_limits<double>::quiet_NaN();
            distance = d_interface_loc[0] - coord[0];
            double sign = 0.0;
            if (d_left_side)
                sign = distance < 0.0 ? 1.0 : -1.0;
            else
                sign = distance < 0.0 ? -1.0 : 1.0;
            (*D_data)(ci) = sign * std::abs(distance);
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
