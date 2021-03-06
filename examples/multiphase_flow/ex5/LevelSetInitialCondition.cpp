// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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

LevelSetInitialCondition::LevelSetInitialCondition(const std::string& object_name, ColumnInterface init_column)
    : d_object_name(object_name), d_init_column(init_column)
{
    // intentionally blank
    return;
} // LevelSetInitialCondition

LevelSetInitialCondition::~LevelSetInitialCondition()
{
    // intentionally blank
    return;
} // ~LevelSetInitialCondition

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
                                         Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    // Set the level set function throughout the domain
    if (initial_time)
    {
        // Get the parameters for the interface
        const IBTK::Vector& X_UR = d_init_column.X_UR;

        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);
        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());

            // Get physical coordinates
            IBTK::Vector coord = IBTK::Vector::Zero();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();
            for (int d = 0; d < NDIM; ++d)
            {
                coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
            }

            // Check if the coordinate is inside the interface
            const bool inside_interface = (coord[0] <= X_UR[0]) && (coord[1] <= X_UR[1])
#if (NDIM == 3)
                                          && (coord[2] <= X_UR[2])
#endif
                ;
            if (inside_interface)
            {
                // If inside the interface, simply set the distance to be the minimum distance from all faces of the
                // column
                double abs_dist[NDIM];
                for (int d = 0; d < NDIM; ++d)
                {
                    abs_dist[d] = std::abs(coord[d] - X_UR[d]);
                }
                (*D_data)(ci) = -(*std::min_element(abs_dist, abs_dist + NDIM));
            }
            else
            {
                // If outside the interface, figure out the closest face and figure out the distance from that.
                // Note that this will make a slight error in distances near the corner, but likely does not matter.
                if (coord[0] >= X_UR[0])
                    (*D_data)(ci) = std::abs(coord[0] - X_UR[0]);
                else if (coord[1] >= X_UR[1])
                    (*D_data)(ci) = std::abs(coord[1] - X_UR[1]);
#if (NDIM == 3)
                else if (coord[2] >= X_UR[2])
                    (*D_data)(ci) = std::abs(coord[2] - X_UR[2]);
#endif
                else
                    TBOX_ERROR("This statement should not be reached");
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
