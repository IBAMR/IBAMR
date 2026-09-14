// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/CartGridPatchwiseFunction.h>
#include <ibtk/CartGridPointwiseFunction.h>

#include <CartesianPatchGeometry.h>
#include <CellIndex.h>
#include <Patch.h>
#include <SideData.h>
#include <SideGeometry.h>

#include <cmath>

#include "InitialConditions.h"

#include <ibamr/app_namespaces.h>

namespace
{
double
sphere_distance(const VectorNd& X, const VectorNd& center, const double radius)
{
    return std::sqrt(std::pow(X[0] - center[0], 2.0) + std::pow(X[1] - center[1], 2.0)
#if (NDIM == 3)
                     + std::pow(X[2] - center[2], 2.0)
#endif
                         ) -
           radius;
}
} // namespace

namespace MultiphaseEx3
{
Pointer<CartGridFunction>
make_sphere_initial_condition(const std::string& object_name,
                              Pointer<CellVariable<NDIM, double>> var,
                              const VectorNd& center,
                              const double radius)
{
    return make_cart_grid_pointwise_function<double>(object_name,
                                                     var,
                                                     [center, radius](const VectorNd& X, double, int, int)
                                                     { return sphere_distance(X, center, radius); });
}

Pointer<CartGridFunction>
make_velocity_initial_condition(const std::string& object_name,
                                const VectorNd& center,
                                const double radius,
                                const double num_interface_cells,
                                const std::vector<double>& inside_velocity,
                                const std::vector<double>& outside_velocity)
{
    return make_cart_grid_patchwise_function(
        object_name,
        [center, radius, num_interface_cells, inside_velocity, outside_velocity](const int data_idx,
                                                                                 Pointer<Variable<NDIM>>,
                                                                                 Pointer<Patch<NDIM>> patch,
                                                                                 double,
                                                                                 const bool initial_time,
                                                                                 Pointer<PatchLevel<NDIM>>)
        {
            if (!initial_time)
            {
                return;
            }

            Pointer<SideData<NDIM, double>> U_data = patch->getPatchData(data_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_X_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                vol_cell *= patch_dx[d];
            }
            const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const double u_inside = inside_velocity[axis], u_outside = outside_velocity[axis];
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    const SideIndex<NDIM> side(it(), axis, SideIndex<NDIM>::Lower);
                    const CellIndex<NDIM> c_l = side.toCell(0), c_u = side.toCell(1);
                    VectorNd coord_lower = VectorNd::Zero(), coord_upper = VectorNd::Zero();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        coord_lower[d] =
                            patch_X_lower[d] + patch_dx[d] * (static_cast<double>(c_l(d) - patch_lower_idx(d)) + 0.5);
                        coord_upper[d] =
                            patch_X_lower[d] + patch_dx[d] * (static_cast<double>(c_u(d) - patch_lower_idx(d)) + 0.5);
                    }
                    // Average adjacent cell distances before applying the smoothed
                    // Heaviside function.
                    const double phi = 0.5 * (sphere_distance(coord_lower, center, radius) +
                                              sphere_distance(coord_upper, center, radius));
                    double h;
                    if (phi < -alpha)
                    {
                        h = 0.0;
                    }
                    else if (std::abs(phi) <= alpha)
                    {
                        h = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                    }
                    else
                    {
                        h = 1.0;
                    }
                    (*U_data)(side) = (u_outside - u_inside) * h + u_inside;
                }
            }
        });
}
} // namespace MultiphaseEx3
