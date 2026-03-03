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

#include <ibtk/samrai_compatibility_names.h>

#include "LevelSetInitialConditionEgg.h"

#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIndex.h>
#include <SAMRAIIndex.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIVariable.h>

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <SAMRAI_config.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialConditionEgg::LevelSetInitialConditionEgg(const std::string& object_name,
                                                         const Pointer<SAMRAICartesianGridGeometry> grid_geom,
                                                         const IBTK::VectorNd& origin)
    : d_object_name(object_name), d_grid_geom(grid_geom), d_origin(origin)
{
    // intentionally blank
    return;
} // LevelSetInitialConditionEgg

bool
LevelSetInitialConditionEgg::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LevelSetInitialConditionEgg::setDataOnPatch(const int data_idx,
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

        // Get physical coordinates
        IBTK::VectorNd coord = IBTK::Vector::Zero();
        IBTK::VectorNd p = IBTK::Vector::Zero();

        for (SAMRAIBox::Iterator it(patch_box); it; it++)
        {
            SAMRAICellIndex ci(it());

            for (int d = 0; d < NDIM; ++d)
                coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - grid_lower_idx(d)) + 0.5);

            p = coord - d_origin;
            const double k = std::sqrt(3.0);
            p[0] = std::abs(p[0]);
            const double ra = 1.2;
            const double rb = 0.2;
            const double r = ra - rb;

            IBTK::Vector2d p_xy(p(0), p(1));
            double q = 0.0;
            if (p(1) < 0.0)
            {
                q = p_xy.norm() - r;
            }
            else
            {
                if (k * (p(0) + r) < p(1))
                {
                    p_xy[1] = p[1] - k * r;
                    q = p_xy.norm();
                }
                else
                {
                    p_xy[0] = p[0] + r;
                    q = p_xy.norm() - 2.0 * r;
                }
            }
            q = q - rb;
            (*D_data)(ci) = q;
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
