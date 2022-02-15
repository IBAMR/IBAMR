// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "LevelSetInitialCondition.h"

#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialCondition::LevelSetInitialCondition(const std::string& object_name, RectangularInterface rectangle)
    : d_object_name(object_name), d_rectangle(rectangle)
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
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);

        Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_X_lower = patch_geom->getXLower();
        const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
        const double* const patch_dx = patch_geom->getDx();

        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());

            // Get physical coordinates
            IBTK::VectorNd coord = IBTK::VectorNd::Zero();
            for (int d = 0; d < NDIM; ++d)
            {
                coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
            }

            // rectangle distance
            const IBTK::VectorNd& origin = d_rectangle.X0;
            IBTK::VectorNd p = coord - origin;
            const IBTK::VectorNd& b = d_rectangle.S;
            IBTK::VectorNd d = p.cwiseAbs() - b;
            const double mm = std::max(d[0], d[1]);
            d[0] = std::max(d[0], 0.0);
            d[1] = std::max(d[1], 0.0);
            const double distance = (d.norm() + std::min(0.0, mm));
            (*D_data)(ci) = distance;
        }
    }

    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
