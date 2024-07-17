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

#include "LSLocateGasInterface.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateGasInterfaceCallbackFunction(int D_idx,
                                         Pointer<HierarchyMathOps> hier_math_ops,
                                         double time,
                                         bool initial_time,
                                         void* ctx)
{
    // Set the level set information
    static LSLocateGasInterface* ptr_LSLocateGasInterface = static_cast<LSLocateGasInterface*>(ctx);
    ptr_LSLocateGasInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // callLSLocateGasInterfaceCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateGasInterface::LSLocateGasInterface(const std::string& object_name,
                                           Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                           Pointer<CellVariableNd<double> > ls_var,
                                           const double greater_x_column,
                                           const double less_z_column)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(ls_var),
      d_greater_x_column(greater_x_column),
      d_less_z_column(less_z_column)
{
    // intentionally left blank
    return;
} // LSLocateGasInterface

LSLocateGasInterface::~LSLocateGasInterface()
{
    // intentionally left blank
    return;
}

void
LSLocateGasInterface::setLevelSetPatchData(int D_idx,
                                           Pointer<HierarchyMathOps> hier_math_ops,
                                           double /*time*/,
                                           bool initial_time)
{
    if (NDIM == 2)
    {
        TBOX_ERROR("Presently not implemented for NDIM != 3");
    }

    Pointer<PatchHierarchyNd> patch_hierarchy = hier_math_ops->getPatchHierarchy();
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

                // Check if the coordinate is inside the interface
                const bool inside_interface = (coord[0] >= d_greater_x_column) && (coord[2] < d_less_z_column);

                // Initialize the locator data to be zero on the interface,
                // negative inside gas region, and positive outside gas region

                // If outside the interface, figure out the closest face and figure out the distance from that.
                // Note that this will make a slight error in distances near the corner, but likely does not matter.
                if (!inside_interface)
                {
                    if (coord[0] < d_greater_x_column)
                    {
                        (*D_data)(ci) = -std::abs(coord[0] - d_greater_x_column);
                    }
                    else if (coord[2] >= d_less_z_column)
                    {
                        (*D_data)(ci) = -std::abs(coord[2] - d_less_z_column);
                    }
                    else
                    {
                        TBOX_ERROR("This statement should not be reached");
                    }
                }
                else
                {
                    double abs_dist_x = std::abs(coord[0] - d_greater_x_column);
                    double abs_dist_z = std::abs(coord[2] - d_less_z_column);
                    (*D_data)(ci) = std::min(abs_dist_x, abs_dist_z);
                }
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
