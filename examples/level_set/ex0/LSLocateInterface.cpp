// Filename LSLocateInterface.cpp
// Created on Oct 1, 2017 by Amneet Bhalla

#include "LSLocateInterface.h"

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

static const double init_positive = 1e8;
static const double init_negative = -1e8;
static const double init_interface = 0.0;

// Initialize the neighborhood of a circular interface.
void
circular_interface_neighborhood(int D_idx,
                                Pointer<HierarchyMathOps> hier_math_ops,
                                double /*time*/,
                                bool /*initial_time*/,
                                void* ctx)
{
    CircularInterface* circle = static_cast<CircularInterface*>(ctx);
    const double& R = circle->R;
    const IBTK::Vector& X0 = circle->X0;

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(D_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
                const double* const patch_dx = patch_geom->getDx();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

                const double distance = std::sqrt(std::pow((coord[0] - X0(0)), 2.0) + std::pow((coord[1] - X0(1)), 2.0)
#if (NDIM == 3)
                                                  + std::pow((coord[2] - X0(2)), 2.0)
#endif
                );
                if (distance < R - 2.0 * patch_dx[0])
                {
                    (*D_data)(ci) = init_negative;
                }
                else if (distance > R + 2.0 * patch_dx[0])
                {
                    (*D_data)(ci) = init_positive;
                }
                else
                {
                    (*D_data)(ci) = init_interface;
                }
            }
        }
    }
    return;
} // circular_interface_neighborhood
