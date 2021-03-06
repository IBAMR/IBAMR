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

LevelSetSolidInitialCondition::LevelSetSolidInitialCondition(const std::string& object_name, WedgeInterface* init_wedge)
    : d_object_name(object_name), d_init_wedge(init_wedge)
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
    if (initial_time && d_init_wedge->wedge_locate_method == "GEOMETRY_METHOD")
    {
        // Set the initial condition for locating the interface
        static const double m = tan(d_init_wedge->wedge_angle);
        static const double wedge_height = 0.5 * (d_init_wedge->wedge_length) * m;

#if (NDIM == 3)
        // Normal of right plane
        static const double ar = -sin(d_init_wedge->wedge_angle);
        static const double br = 0.0;
        static const double cr = cos(d_init_wedge->wedge_angle);

        // Normal of left plane
        static const double al = sin(d_init_wedge->wedge_angle);
        static const double bl = 0.0;
        static const double cl = cos(d_init_wedge->wedge_angle);
#endif

        // Set the initial condition for locating the interface
        IBTK::Vector& X0 = d_init_wedge->X0;

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

            double distance[3];

            // a) Top plane
            const double H = wedge_height + X0(1);
            distance[0] = X(1) - H;

            // b) Right line of the form y = slope*x + c
            const double cr = X0(1) - m * X0(0);
            distance[1] = (m * X(0) - X(1) + cr) / sqrt(1.0 + m * m);

            // c) Left line of the form y = slope*x + c
            const double cl = X0(1) + m * X0(0);
            distance[2] = (-m * X(0) - X(1) + cl) / sqrt(1.0 + m * m);

            (*D_data)(ci) = *std::max_element(distance, distance + 3);

#endif

#if (NDIM == 3)

            double distance[5];

            // a) Top plane
            const double H = wedge_height + X0(2);
            distance[0] = X(2) - H;

            // b) Right plane of the form: ax + by + cz + d = 0
            const double dr = 0.0 - ar * X0(0) - br * X0(1) - cr * X0(2);
            distance[1] = -dr - ar * X(0) - br * X(1) - cr * X(2);

            // c) Left plane  of the form ax + by + cz + d = 0
            const double dl = 0.0 - al * X0(0) - bl * X0(1) - cl * X0(2);
            distance[2] = -dl - al * X(0) - bl * X(1) - cl * X(2);

            // d) Front plane
            const double Yf = X0(1) - (d_init_wedge->wedge_length) / 2.0;
            distance[3] = Yf - X(1);

            // d) Back plane
            const double Yb = X0(1) + (d_init_wedge->wedge_length) / 2.0;
            distance[4] = X(1) - Yb;

            (*D_data)(ci) = *std::max_element(distance, distance + 5);
#endif
        }
    }
    return;
} // setDataOnPatch
