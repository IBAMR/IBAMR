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

void
call_locate_interface(const int data_idx,
                      Pointer<HierarchyMathOps> hier_math_ops,
                      const double time,
                      const bool initial_time,
                      void* ctx)
{
    auto* locator = static_cast<LSLocateInterface*>(ctx);
    locator->setLevelSetPatchData(data_idx, hier_math_ops, time, initial_time);
}

LSLocateInterface::LSLocateInterface(Pointer<AdvDiffHierarchyIntegrator> integrator,
                                     Pointer<CellVariable<NDIM, double>> ls_var,
                                     Pointer<CartGridFunction> initial_conditions)
    : d_integrator(integrator), d_ls_var(ls_var), d_initial_conditions(initial_conditions)
{
}

void
LSLocateInterface::setLevelSetPatchData(const int data_idx,
                                        Pointer<HierarchyMathOps> hier_math_ops,
                                        const double time,
                                        const bool initial_time)
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = hier_math_ops->getPatchHierarchy();
    if (initial_time)
    {
        d_initial_conditions->setDataOnPatchHierarchy(data_idx, d_ls_var, hierarchy, time, true);
    }
    else
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int current_idx = var_db->mapVariableAndContextToIndex(d_ls_var, d_integrator->getCurrentContext());
        HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, hierarchy->getFinestLevelNumber());
        data_ops.copyData(data_idx, current_idx);
    }
}
