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

#include "LevelSetInitialConditionTorus.h"

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <SAMRAI_config.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialConditionTorus::LevelSetInitialConditionTorus(const std::string& object_name,
                                                             const SAMRAIPointer<CartesianGridGeometryNd> grid_geom,
                                                             const IBTK::VectorNd& origin,
                                                             const IBTK::Vector2d& t)
    : d_object_name(object_name), d_grid_geom(grid_geom), d_origin(origin), d_t(t)
{
    // intentionally blank
    return;
} // LevelSetInitialConditionTorus

bool
LevelSetInitialConditionTorus::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LevelSetInitialConditionTorus::setDataOnPatch(const int data_idx,
                                              SAMRAIPointer<VariableNd> /*var*/,
                                              SAMRAIPointer<PatchNd> patch,
                                              const double /*data_time*/,
                                              const bool initial_time,
                                              SAMRAIPointer<PatchLevelNd> patch_level)
{
    // Set the level set function throughout the domain
    if (initial_time)
    {
        const BoxNd& patch_box = patch->getBox();
        SAMRAIPointer<CellDataNd<double> > D_data = patch->getPatchData(data_idx);

        // Get physical coordinates
        IBTK::VectorNd coord = IBTK::Vector::Zero();
        IBTK::VectorNd p = IBTK::Vector::Zero();
        IBTK::Vector2d q(0.0, 0.0);

        SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double* const grid_x_lower = d_grid_geom->getXLower();
        IntVectorNd ratio = patch_level->getRatio();
        const SAMRAI::hier::BoxNd domain_box = SAMRAI::hier::BoxNd::refine(d_grid_geom->getPhysicalDomain()[0], ratio);
        const hier::IndexNd& grid_lower_idx = domain_box.lower();

        for (BoxNd::Iterator it(patch_box); it; it++)
        {
            CellIndexNd ci(it());

            for (int d = 0; d < NDIM; ++d)
                coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - grid_lower_idx(d)) + 0.5);

            p = coord - d_origin;
            IBTK::Vector2d p_xz(p[0], p[2]);
            q << p_xz.norm() - d_t[0], p[1];
            (*D_data)(ci) = (q.norm() - d_t[1]);
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
