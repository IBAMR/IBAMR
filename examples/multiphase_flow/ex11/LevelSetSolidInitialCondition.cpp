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

LevelSetSolidInitialCondition::LevelSetSolidInitialCondition(const std::string& object_name, BargeInterface* init_barge)
    : d_object_name(object_name), d_init_barge(init_barge)
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
        // Set the initial condition for locating the interface
        double Xcom = d_init_barge->COM(0);
        double Ycom = d_init_barge->COM(1);
        double L = d_init_barge->length;
        double W = d_init_barge->width;
        double theta = d_init_barge->theta;

        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);
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

#if (NDIM == 2)
            const double relx = X[0] - Xcom;
            const double rely = X[1] - Ycom;
            const double rotx = relx * std::cos(-theta) - rely * std::sin(-theta);
            const double roty = relx * std::sin(-theta) + rely * std::cos(-theta);

            bool inside = (-0.5 * L < rotx && rotx < 0.5 * L && -0.5 * W < roty && roty < 0.5 * W);
            if (!inside)
            {
                const double dx = std::max(abs(rotx) - 0.5 * L, 0.0);
                const double dy = std::max(abs(roty) - 0.5 * W, 0.0);
                (*D_data)(ci) = std::sqrt(dx * dx + dy * dy);
            }
            else
            {
                double dx_min = std::min(abs(rotx + 0.5 * L), abs(0.5 * L - rotx));
                double dy_min = std::min(abs(roty + 0.5 * W), abs(0.5 * W - roty));
                (*D_data)(ci) = -std::min(dx_min, dy_min);
            }

#endif

#if (NDIM == 3)
            TBOX_ERROR("Presently not implemented");
#endif
        }
    }
    return;
} // setDataOnPatch
