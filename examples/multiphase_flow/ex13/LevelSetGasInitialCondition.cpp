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

#include "LevelSetGasInitialCondition.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetGasInitialCondition::LevelSetGasInitialCondition(const std::string& object_name,
                                                         const double greater_x_column,
                                                         const double less_z_column)
    : d_object_name(object_name), d_greater_x_column(greater_x_column), d_less_z_column(less_z_column)
{
    // intentionally blank
    return;
} // LevelSetGasInitialCondition

LevelSetGasInitialCondition::~LevelSetGasInitialCondition()
{
    // intentionally blank
    return;
} // ~LevelSetGasInitialCondition

bool
LevelSetGasInitialCondition::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LevelSetGasInitialCondition::setDataOnPatch(const int data_idx,
                                            Pointer<Variable<NDIM> > /*var*/,
                                            Pointer<Patch<NDIM> > patch,
                                            const double /*data_time*/,
                                            const bool initial_time,
                                            Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    // Set the level set function throughout the domain
    if (initial_time)
    {
        if (NDIM == 2)
        {
            TBOX_ERROR("Presently not implemented for NDIM != 3");
        }

        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
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
            const bool inside_interface = (coord[0] >= d_greater_x_column) && (coord[2] < d_less_z_column);

            // Initialize the locator data to be zero on the interface,
            // negative inside gas region, and positive outside gas region

            // If outside the interface, figure out the closest face and figure out the distance from that.
            // Note that this will make a slight error in distances near the corner, but likely does not matter.
            if (!inside_interface)
            {
                if (coord[0] < d_greater_x_column)
                {
                    (*D_data)(ci) = -std::abs(coord[0] - d_greater_x_column);
                }
                else if (coord[2] >= d_less_z_column)
                {
                    (*D_data)(ci) = -std::abs(coord[2] - d_less_z_column);
                }
                else
                {
                    TBOX_ERROR("This statement should not be reached");
                }
            }
            else
            {
                double abs_dist_x = std::abs(coord[0] - d_greater_x_column);
                double abs_dist_z = std::abs(coord[2] - d_less_z_column);
                (*D_data)(ci) = std::min(abs_dist_x, abs_dist_z);
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
