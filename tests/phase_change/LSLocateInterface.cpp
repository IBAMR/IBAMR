// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
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

#include <HierarchyCellDataOpsReal.h>
#include <VariableDatabase.h>

#include "LSLocateInterface.h"

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
call_ls_locate_interface_callback(int D_idx,
                                  Pointer<HierarchyMathOps> hier_math_ops,
                                  double time,
                                  bool initial_time,
                                  void* ctx)
{
    // Set the level set information
    LSLocateInterface* ptr_LSLocateInterface = static_cast<LSLocateInterface*>(ctx);
    ptr_LSLocateInterface->setLevelSetPatchData(D_idx, hier_math_ops, time, initial_time);

    return;
} // call_ls_locate_interface_callback

/////////////////////////////// PUBLIC //////////////////////////////////////
LSLocateInterface::LSLocateInterface(const std::string& object_name,
                                     Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                     Pointer<CellVariable<NDIM, double>> ls_var,
                                     Pointer<CartGridFunction> initial_condition)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(ls_var),
      d_initial_condition(initial_condition)
{
    // intentionally left blank
    return;
} // LSLocateInterface

void
LSLocateInterface::setLevelSetPatchData(int D_idx,
                                        Pointer<HierarchyMathOps> hier_math_ops,
                                        const double time,
                                        bool initial_time)
{
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // If not the initial time, set the level set to the current value maintained
    // by the integrator
    if (!initial_time)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int ls_current_idx =
            var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        hier_cc_data_ops.copyData(D_idx, ls_current_idx);

        return;
    }

    d_initial_condition->setDataOnPatchHierarchy(D_idx, d_ls_var, patch_hierarchy, time, true, coarsest_ln, finest_ln);
    return;
} // setLevelSetPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
