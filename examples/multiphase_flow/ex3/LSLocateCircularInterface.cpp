// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2021 by the IBAMR developers
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
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/samrai_compatibility_names.h>

#include "LSLocateCircularInterface.h"

#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIndex.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIHierarchyCellDataOpsReal.h>
#include <SAMRAIIndex.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAIVariableDatabase.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callLSLocateCircularInterfaceCallbackFunction(int D_idx,
                                              SAMRAIPointer<HierarchyMathOps> hier_math_ops,
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
                                                     SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                     SAMRAIPointer<SAMRAICellVariable<double>> ls_var,
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
    const double& R = d_init_circle.R;
    const IBTK::Vector& X0 = d_init_circle.X0;

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
