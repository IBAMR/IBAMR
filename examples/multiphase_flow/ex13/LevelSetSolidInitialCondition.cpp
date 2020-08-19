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

#include "LevelSetSolidInitialCondition.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetSolidInitialCondition::LevelSetSolidInitialCondition(const std::string& object_name,
                                                             RectangleInterface init_rectangle)
    : d_object_name(object_name), d_init_rectangle(init_rectangle)
{
    // intentionally blank
    return;
} // LevelSetSolidInitialCondition

LevelSetSolidInitialCondition::~LevelSetSolidInitialCondition()
{
    // intentionally blank
    return;
} // ~LevelSetSolidInitialCondition

bool
LevelSetSolidInitialCondition::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LevelSetSolidInitialCondition::setDataOnPatch(const int data_idx,
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

        const double x0 = d_init_rectangle.X0(0);
        const double y0 = d_init_rectangle.X0(1);
        const double z0 = d_init_rectangle.X0(2);
        const double l = d_init_rectangle.length;
        const double w = d_init_rectangle.width;
        const double h = d_init_rectangle.height;

        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());

            // Get physical coordinates
            IBTK::Vector X = IBTK::Vector::Zero();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const SAMRAI::hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
            }

            const double relx = X[0] - x0;
            const double rely = X[1] - y0;
            const double relz = X[2] - z0;

            bool inside = (-0.5 * l < relx && relx < 0.5 * l && -0.5 * w < rely && rely < 0.5 * w && -0.5 * h < relz &&
                           relz < 0.5 * h);

            if (!inside)
            {
                const double dx = std::max(abs(relx) - 0.5 * l, 0.0);
                const double dy = std::max(abs(rely) - 0.5 * w, 0.0);
                const double dz = std::max(abs(relz) - 0.5 * h, 0.0);
                (*D_data)(ci) = std::sqrt(dx * dx + dy * dy + dz * dz);
            }
            else
            {
                double dx_min = std::min(abs(relx + 0.5 * l), abs(0.5 * l - relx));
                double dy_min = std::min(abs(rely + 0.5 * w), abs(0.5 * w - rely));
                double dz_min = std::min(abs(relz + 0.5 * h), abs(0.5 * h - relz));
                (*D_data)(ci) = -std::min(dx_min, std::min(dy_min, dz_min));
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
