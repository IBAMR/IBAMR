// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/app_namespaces.h>

#include <ibtk/HierarchyMathOps.h>

#include "LSLocateCircularInterface.h"

#include <CartesianGridGeometry.h>

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
                                                     Pointer<CellVariable<NDIM, double> > ls_var,
                                                     CircularInterface init_circle)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_init_circle(init_circle)
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
                                                bool initial_time)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // If not the intial time, set the level set to the current value maintained by the integrator
    if (!initial_time)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int ls_current_idx =
            var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        hier_cc_data_ops.copyData(D_idx, ls_current_idx);

        return;
    }

    // Set the initial condition for locating the interface
    const double& R = d_init_circle.R;
    const IBTK::Vector& X0 = d_init_circle.X0;

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

                // Initialize the locator data
                (*D_data)(ci) = distance - R;
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
