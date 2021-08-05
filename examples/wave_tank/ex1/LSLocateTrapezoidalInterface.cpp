// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#include "LSLocateTrapezoidalInterface.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateTrapezoidalInterfaceCallbackFunction(int D_idx,
                                                 Pointer<HierarchyMathOps> hier_math_ops,
                                                 double time,
                                                 bool initial_time,
                                                 void* ctx)
{
    // Set the level set information
    static LSLocateTrapezoidalInterface* ptr_LSLocateTrapezoidalInterface =
        static_cast<LSLocateTrapezoidalInterface*>(ctx);
    ptr_LSLocateTrapezoidalInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateTrapezoidalInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

LSLocateTrapezoidalInterface::LSLocateTrapezoidalInterface(const std::string& object_name,
                                                           Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                           Pointer<CellVariable<NDIM, double> > ls_var,
                                                           TrapezoidalInterface* trapezoid)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_trapezoid(trapezoid)
{
    // intentionally left blank
    return;
} // LSLocateTrapezoidalInterface

void
LSLocateTrapezoidalInterface::setLevelSetPatchData(int D_idx,
                                                   Pointer<HierarchyMathOps> hier_math_ops,
                                                   double /*time*/,
                                                   bool /*initial_time*/)
{
#if (NDIM == 3)
    TBOX_ERROR("LSLocateTrapezoidalInterface is only implemented for NDIM = 2");
#endif
    // Note that this class assumes the object remains stationary
    // Also note that I did not test this algorithm on very many trapezoids, so this may need to be modified accordingly
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Set the initial condition for locating the interface
    const IBTK::Vector& BL = d_trapezoid->BL;
    const IBTK::Vector& BR = d_trapezoid->BR;
    const IBTK::Vector& TL = d_trapezoid->TL;
    const IBTK::Vector& TR = d_trapezoid->TR;

    // Slopes and intercepts for the four lines that make up the trapezoid
    const double m0 = (TL(1) - BL(1)) / (TL(0) - BL(0));
    const double b0 = TL(1) - m0 * TL(0);

    const double m1 = (TR(1) - TL(1)) / (TR(0) - TL(0));
    const double b1 = TR(1) - m1 * TR(0);

    const double m2 = (BR(1) - TR(1)) / (BR(0) - TR(0));
    const double b2 = BR(1) - m2 * BR(0);

    const double m3 = (BL(1) - BR(1)) / (BL(0) - BR(0));
    const double b3 = BL(1) - m3 * BL(0);

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
                const double x = coord[0];
                const double y = coord[1];

                // Whether or not the Eulerian cell is inside the trapezoid
                const double x_lower = std::max(BL(0), TL(0));
                // const double x_upper = std::min(BR(0), TR(0));
                const double y_lower = std::max(BL(1), BR(1));
                const double y_upper = std::min(TL(1), TR(1));
                const bool inside = (y <= m0 * x + b0 && y <= m1 * x + b1 && y <= m2 * x + b2 && y >= m3 * x + b3);
                double signed_distance = std::numeric_limits<double>::max();

                if (inside)
                {
                    // Min distance between point and lines
                    const double d1 = std::min(pointToLineDistance(coord, BL, TL), pointToLineDistance(coord, TL, TR));
                    const double d2 = std::min(pointToLineDistance(coord, TR, BR), pointToLineDistance(coord, BR, BL));
                    signed_distance = -std::min(d1, d2);
                }
                else
                {
                    if (y >= y_upper && TL(0) <= x && x <= TR(0))
                    {
                        signed_distance = pointToLineDistance(coord, TL, TR);
                    }
                    else if (y <= y_lower && BL(0) <= x && x <= BR(0))
                    {
                        signed_distance = pointToLineDistance(coord, BL, BR);
                    }
                    else if (x <= x_lower)
                    {
                        signed_distance = pointToLineDistance(coord, BL, TL);
                    }
                    else
                    {
                        signed_distance = pointToLineDistance(coord, TR, BR);
                    }
                }
                (*D_data)(ci) = signed_distance;
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////

double
LSLocateTrapezoidalInterface::pointToLineDistance(IBTK::Vector X0, IBTK::Vector P1, IBTK::Vector P2)
{
    const double x0 = X0(0), y0 = X0(1);
    const double x1 = P1(0), y1 = P1(1);
    const double x2 = P2(0), y2 = P2(1);

    return std::abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) /
           std::sqrt(std::pow(y2 - y1, 2.0) + std::pow(x2 - x1, 2.0));
}
