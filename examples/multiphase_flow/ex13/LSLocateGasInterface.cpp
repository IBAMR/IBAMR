// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
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
                                           Pointer<CellVariable<NDIM, double> > ls_var,
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
