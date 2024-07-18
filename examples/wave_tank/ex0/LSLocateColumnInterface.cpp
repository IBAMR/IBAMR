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

#include "LSLocateColumnInterface.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateColumnInterfaceCallbackFunction(int D_idx,
                                            SAMRAIPointer<HierarchyMathOps> hier_math_ops,
                                            double time,
                                            bool initial_time,
                                            void* ctx)
{
    // Set the level set information
    static LSLocateColumnInterface* ptr_LSLocateColumnInterface = static_cast<LSLocateColumnInterface*>(ctx);
    ptr_LSLocateColumnInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateColumnInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateColumnInterface::LSLocateColumnInterface(const std::string& object_name,
                                                 SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                 SAMRAIPointer<CellVariableNd<double> > ls_var,
                                                 ColumnInterface init_column)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_init_column(init_column)
{
    // intentionally left blank
    return;
} // LSLocateColumnInterface

void
LSLocateColumnInterface::setLevelSetPatchData(int D_idx,
                                              SAMRAIPointer<HierarchyMathOps> hier_math_ops,
                                              double /*time*/,
                                              bool initial_time)
{
    SAMRAIPointer<PatchHierarchyNd> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // If not the intial time, set the level set to the current value maintained by the integrator
    if (!initial_time)
    {
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        const int ls_current_idx =
            var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
        HierarchyCellDataOpsRealNd<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        hier_cc_data_ops.copyData(D_idx, ls_current_idx);

        return;
    }

    // Set the initial condition for locating the interface
    const double& depth = d_init_column.DEPTH;
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
                (*D_data)(ci) = coord[NDIM - 1] - depth;
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
