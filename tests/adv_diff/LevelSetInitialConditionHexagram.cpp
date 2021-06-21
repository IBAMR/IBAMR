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

#include "LevelSetInitialConditionHexagram.h"

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "tbox/Pointer.h"

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialConditionHexagram::LevelSetInitialConditionHexagram(
    const std::string& object_name,
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom,
    const IBTK::VectorNd& origin)
    : d_object_name(object_name), d_grid_geom(grid_geom), d_origin(origin)
{
    // intentionally blank
    return;
} // LevelSetInitialCondition

bool
LevelSetInitialConditionHexagram::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LevelSetInitialConditionHexagram::setDataOnPatch(const int data_idx,
                                                 Pointer<SAMRAI::hier::Variable<NDIM> > /*var*/,
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
        std::vector<double> k{ -0.5, 0.8660254038, 0.5773502692, 1.7320508076 };
        IBTK::VectorNd kxy{ k[0], k[1] };
        IBTK::VectorNd kyx{ k[1], k[0] };

        // Get physical coordinates
        IBTK::VectorNd coord = IBTK::Vector::Zero();
        IBTK::VectorNd p = IBTK::Vector::Zero();
        IBTK::VectorNd q = IBTK::Vector::Zero();

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

            for (int d = 0; d < NDIM; ++d)
                coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - grid_lower_idx(d)) + 0.5);

            p = coord - d_origin;
            q = p.cwiseAbs();
            q -= 2.0 * std::min(kxy.dot(q), 0.0) * kxy;
            q -= 2.0 * std::min(kyx.dot(q), 0.0) * kyx;
            const double width = 1.0;
            const double min_val = width * k[2];
            const double max_value = width * k[3];
            const double clamp = std::min(std::max(q[0], min_val), max_value);
            IBTK::VectorNd mm(clamp, width);
            q = q - mm;
            const double sign = q[1] > 0 ? 1.0 : ((q[1] < 0) ? -1 : 0);
            (*D_data)(ci) = q.norm() * sign;
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
