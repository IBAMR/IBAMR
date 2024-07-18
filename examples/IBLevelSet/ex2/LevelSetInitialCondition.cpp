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

#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialCondition::LevelSetInitialCondition(const std::string& object_name, CircularInterface init_circle)
    : d_object_name(object_name), d_init_circle(init_circle)
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
                                         SAMRAIPointer<VariableNd> /*var*/,
                                         SAMRAIPointer<PatchNd> patch,
                                         const double /*data_time*/,
                                         const bool initial_time,
                                         SAMRAIPointer<PatchLevelNd> /*patch_level*/)
{
    // Set the level set function throughout the domain
    if (initial_time)
    {
        // Get the parameters for the interface
        const double& R = d_init_circle.R;
        const Eigen::Vector3d& X0 = d_init_circle.X0;

        const BoxNd& patch_box = patch->getBox();
        SAMRAIPointer<CellDataNd<double> > D_data = patch->getPatchData(data_idx);
        for (BoxNd::Iterator it(patch_box); it; it++)
        {
            CellIndexNd ci(it());

            // Get physical coordinates
            IBTK::Vector coord = IBTK::Vector::Zero();
            SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::IndexNd& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();
            for (int d = 0; d < NDIM; ++d)
            {
                coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
            }

            // Distance from the circle
            const double distance_circle =
                std::sqrt(std::pow((coord[0] - X0(0)), 2.0) + std::pow((coord[1] - X0(1)), 2.0)
#if (NDIM == 3)
                          + std::pow((coord[2] - X0(2)), 2.0)
#endif
                              ) -
                R;

            (*D_data)(ci) = distance_circle;
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
