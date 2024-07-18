// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/HierarchyMathOps.h>

#include "LSLocateInterface.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

// Initialize the neighborhood of a circular interface.
void
circular_interface_neighborhood(int D_idx,
                                SAMRAIPointer<HierarchyMathOps> hier_math_ops,
                                double /*time*/,
                                bool /*initial_time*/,
                                void* /*ctx*/)
{
    SAMRAIPointer<PatchHierarchyNd> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            SAMRAIPointer<CellDataNd<double> > D_data = patch->getPatchData(D_idx);
            for (BoxNd::Iterator it(patch_box); it; it++)
            {
                CellIndexNd ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const hier::IndexNd& patch_lower_idx = patch_box.lower();
                const double* const patch_dx = patch_geom->getDx();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }
                const double x = coord[0];
                const double y = coord[1];
#if (NDIM == 3)
                const double z = coord[2];
#endif
                (*D_data)(ci) = 1.0 *
                                (std::pow(x - 1.0, 2.0) + std::pow(y - 1.0, 2.0)
#if (NDIM == 3)
                                 + std::pow(z - 1.0, 2.0)
#endif
                                 + 0.1) *
                                (std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)
#if (NDIM == 3)
                                           + std::pow(z, 2.0)
#endif
                                               ) -
                                 1.0);
            }
        }
    }
    return;
} // circular_interface_neighborhood
