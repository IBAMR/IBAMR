// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// SAMRAI INCLUDES
#include "ibtk/samrai_compatibility_names.h"
#include <ibtk/HierarchyMathOps.h>

#include "LSLocateGasInterface.h"
#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIHierarchyCellDataOpsReal.h"
#include "SAMRAIIndex.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIVariableDatabase.h"

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateGasInterfaceCallbackFunction(int D_idx,
                                         SAMRAIPointer<HierarchyMathOps> hier_math_ops,
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
                                           SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                           SAMRAIPointer<SAMRAICellVariable<double>> ls_var,
                                           const double init_height)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(ls_var), d_init_height(init_height)
{
    // intentionally left blank
    return;
} // LSLocateGasInterface

void
LSLocateGasInterface::setLevelSetPatchData(int D_idx,
                                           SAMRAIPointer<HierarchyMathOps> hier_math_ops,
                                           double /*time*/,
                                           bool initial_time)
{
    SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // If not the intial time, set the level set to the current value maintained by the integrator
    if (!initial_time)
    {
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        const int ls_current_idx =
            var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
        SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        hier_cc_data_ops.copyData(D_idx, ls_current_idx);

        return;
    }

    // Set the initial condition for locating the interface
    const double H = d_init_height;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        for (SAMRAIPatchLevel::Iterator p(level); p; p++)
        {
            SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
            const SAMRAIBox& patch_box = patch->getBox();
            SAMRAIPointer<SAMRAICellData<double>> D_data = patch->getPatchData(D_idx);
            for (SAMRAIBox::Iterator it(patch_box); it; it++)
            {
                SAMRAICellIndex ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const SAMRAIIndex& patch_lower_idx = patch_box.lower();
                const double* const patch_dx = patch_geom->getDx();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

                const double distance = NDIM < 3 ? coord[1] - H : coord[2] - H;

                // Initialize the locator data to be zero on the interface,
                // positive inside, and negative outside
                (*D_data)(ci) = -distance;
            }
        }
    }
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
