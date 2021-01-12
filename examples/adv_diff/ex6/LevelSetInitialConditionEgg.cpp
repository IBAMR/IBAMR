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

#include "LevelSetInitialConditionEgg.h"

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <SAMRAI_config.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialConditionEgg::LevelSetInitialConditionEgg(const std::string& object_name, const IBTK::VectorNd& origin)
    : d_object_name(object_name), d_origin(origin)
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
                                            Pointer<Variable<NDIM> > /*var*/,
                                            Pointer<Patch<NDIM> > patch,
                                            const double /*data_time*/,
                                            const bool initial_time,
                                            Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    // Set the level set function throughout the domain
    if (initial_time)
    {
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);

        // Get physical coordinates
        IBTK::VectorNd coord = IBTK::Vector::Zero();
        IBTK::VectorNd p = IBTK::Vector::Zero();

        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            for (int d = 0; d < NDIM; ++d)
            {
                coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
            }
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
