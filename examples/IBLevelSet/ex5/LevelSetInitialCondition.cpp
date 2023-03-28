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

#include <HierarchyDataOpsManager.h>

// Application includes
#include "LevelSetInitialCondition.h"

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetInitialCondition::LevelSetInitialCondition(const std::string& object_name, FoilInterface init_foil)
    : d_object_name(object_name), d_init_foil(init_foil)
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
        // Get the parameters for the interface
        const double& R = d_init_foil.R;
        const Eigen::Vector3d& X0 = d_init_foil.X0;
        const Eigen::Vector3d& X1 = d_init_foil.X1;
        const Eigen::Vector3d& X2 = d_init_foil.X2;
        const Eigen::Vector3d& X3 = d_init_foil.X3;
        const Eigen::Vector3d& X_T = (X1 + X2 + X3) / 3.0;
        Eigen::Vector3d check1, check2, check3, check4, check5, check6;

        const double slope1 = (X3[1] - X1[1]) / (X3[0] - X1[0]);
        const double slope2 = (X3[1] - X2[1]) / (X3[0] - X2[0]);
        const double slope3 = (X2[0] - X1[0]) / (X2[1] - X1[1]);

        const double y_intercept1 = X1[1] - slope1 * X1[0];
        const double y_intercept2 = X2[1] - slope2 * X2[0];
        const double x_intercept3 = X1[0] - slope3 * X1[1];

        double distance1[2], distance2[3]; // Foil has three surfaces and 1 surface for circle.

        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(data_idx);
        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());

            // Get physical coordinates
            IBTK::Vector3d coord = IBTK::Vector3d::Zero();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();
            for (int d = 0; d < NDIM; ++d)
            {
                coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
            }

            // Distance from the semi-circle
            distance1[0] = std::sqrt(std::pow((coord[0] - X0[0]), 2.0) + std::pow((coord[1] - X0[1]), 2.0)) - R;

            // Distance from top triangle surface.
            distance2[0] = std::abs(coord[1] - slope1 * coord[0] - y_intercept1) / std::sqrt(1.0 + slope1 * slope1);

            check1 = (X1 - X3).cross(X_T - X3);
            check2 = (X1 - X3).cross(coord - X3);

            distance2[0] *= (-sgn(check1[2]) * sgn(check2[2]));

            // Distance from bottom triangle surface.
            distance2[1] = std::abs(coord[1] - slope2 * coord[0] - y_intercept2) / std::sqrt(1.0 + slope2 * slope2);

            check3 = (X2 - X3).cross(X_T - X3);
            check4 = (X2 - X3).cross(coord - X3);

            distance2[1] *= (-sgn(check3[2]) * sgn(check4[2]));

            // Distance from base of triangle.
            distance2[2] = std::abs(coord[0] - slope3 * coord[1] - x_intercept3) / std::sqrt(1.0 + slope3 * slope3);

            check5 = (X2 - X1).cross(X_T - X1);
            check6 = (X2 - X1).cross(coord - X1);

            distance2[2] *= (-sgn(check5[2]) * sgn(check6[2]));

            distance1[1] =
                std::max({ distance2[0], distance2[1], distance2[2] }); // intersection of three traingle surfaces

            (*D_data)(ci) = std::min({ distance1[0], distance1[1] }); // union of circle and triangle
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
