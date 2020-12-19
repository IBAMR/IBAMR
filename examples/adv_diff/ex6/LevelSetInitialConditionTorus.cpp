// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
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
                                                             const IBTK::VectorNd& origin,
                                                             const IBTK::Vector2d& t)
    : d_object_name(object_name), d_origin(origin), d_t(t)
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
        IBTK::Vector2d q(0.0, 0.0);
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
            IBTK::Vector2d p_xz(p[0], p[2]);
            q << p_xz.norm() - d_t[0], p[1];
            (*D_data)(ci) = (q.norm() - d_t[1]);
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
