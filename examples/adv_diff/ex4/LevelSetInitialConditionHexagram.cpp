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
#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIndex.h>
#include <SAMRAIHierarchyDataOpsManager.h>
#include <SAMRAIIndex.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAIVariable.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialConditionHexagram::LevelSetInitialConditionHexagram(
    const std::string& object_name,
    const SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geom,
    const IBTK::VectorNd& origin)
    : d_object_name(object_name), d_grid_geom(grid_geom), d_origin(origin)
{
    // intentionally blank
    return;
} // LevelSetInitialConditionHexagram

bool
LevelSetInitialConditionHexagram::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LevelSetInitialConditionHexagram::setDataOnPatch(const int data_idx,
                                                 Pointer<SAMRAIVariable> /*var*/,
                                                 Pointer<SAMRAIPatch> patch,
                                                 const double /*data_time*/,
                                                 const bool initial_time,
                                                 Pointer<SAMRAIPatchLevel> patch_level)
{
    // Set the level set function throughout the domain
    if (initial_time)
    {
        const SAMRAIBox& patch_box = patch->getBox();
        Pointer<SAMRAICellData<double>> D_data = patch->getPatchData(data_idx);

        Pointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double* const grid_x_lower = d_grid_geom->getXLower();
        SAMRAIIntVector ratio = patch_level->getRatio();
        const SAMRAIBox domain_box = SAMRAIBox::refine(d_grid_geom->getPhysicalDomain()[0], ratio);
        const SAMRAIIndex& grid_lower_idx = domain_box.lower();

        std::vector<double> k{ -0.5, 0.8660254038, 0.5773502692, 1.7320508076 };
        IBTK::VectorNd kxy{ k[0], k[1] };
        IBTK::VectorNd kyx{ k[1], k[0] };
        // Get physical coordinates
        IBTK::VectorNd coord = IBTK::Vector::Zero();
        IBTK::VectorNd p = IBTK::Vector::Zero();
        IBTK::VectorNd q = IBTK::Vector::Zero();

        for (SAMRAIBox::Iterator it(patch_box); it; it++)
        {
            SAMRAICellIndex ci(it());

            for (int d = 0; d < NDIM; ++d)
                coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - grid_lower_idx(d)) + 0.5);

            p = coord - d_origin;
            q = p.cwiseAbs();
            q -= 2.0 * std::min(kxy.dot(q), 0.0) * kxy;
            q -= 2.0 * std::min(kyx.dot(q), 0.0) * kyx;
            const double width = 1.5;
            const double min_val = width * k[2];
            const double max_value = width * k[3];
            const double clamp = std::min(std::max(q[0], min_val), max_value);
            IBTK::VectorNd mm(clamp, width);
            q = q - mm;
            const double sign = q[1] < 0 ? 1.0 : ((q[1] > 0) ? -1 : 0);
            (*D_data)(ci) = q.norm() * sign;
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
