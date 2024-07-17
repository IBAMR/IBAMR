// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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

#include "LSLocateCircularInterface.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateCircularInterfaceCallbackFunction(int D_idx,
                                              Pointer<HierarchyMathOps> hier_math_ops,
                                              double time,
                                              bool initial_time,
                                              void* ctx)
{
    // Set the level set information
    static LSLocateCircularInterface* ptr_LSLocateCircularInterface = static_cast<LSLocateCircularInterface*>(ctx);
    ptr_LSLocateCircularInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateCircularInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateCircularInterface::LSLocateCircularInterface(const std::string& object_name,
                                                     Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                     Pointer<CellVariableNd<double> > ls_var,
                                                     CircularInterface* circle)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_circle(circle)
{
    // intentionally left blank
    return;
} // LSLocateCircularInterface

LSLocateCircularInterface::~LSLocateCircularInterface()
{
    // intentionally left blank
    return;
}

void
LSLocateCircularInterface::setLevelSetPatchData(int D_idx,
                                                Pointer<HierarchyMathOps> hier_math_ops,
                                                double /*time*/,
                                                bool /*initial_time*/)
{
    // In this version of this class, the initial level set location is set to be
    // exact since we always know the radius of the ball

    Pointer<PatchHierarchyNd> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Set the initial condition for locating the interface
    const double& R = d_circle->R;
    const IBTK::Vector& X0 = d_circle->X0;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            Pointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            Pointer<CellDataNd<double> > D_data = patch->getPatchData(D_idx);
            for (BoxNd::Iterator it(patch_box); it; it++)
            {
                CellIndexNd ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                Pointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const hier::IndexNd& patch_lower_idx = patch_box.lower();
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

                (*D_data)(ci) = distance - R;
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
