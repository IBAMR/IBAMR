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

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/HeavisideForcingFunction.h"

#include <SAMRAI_config.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

HeavisideForcingFunction::HeavisideForcingFunction(const std::string& object_name,
                                                   const Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                   const Pointer<CellVariable<NDIM, double> > H_var,
                                                   const Pointer<FaceVariable<NDIM, double> > U_adv_var)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver), d_H_var(H_var), d_U_adv_var(U_adv_var)
{
    // intentionally blank
    return;
} // HeavisideForcingFunction

bool
HeavisideForcingFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
HeavisideForcingFunction::setDataOnPatchHierarchy(const int data_idx,
                                                  Pointer<Variable<NDIM> > var,
                                                  Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                  const double data_time,
                                                  const bool initial_time,
                                                  const int coarsest_ln_in,
                                                  const int finest_ln_in)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    d_hier_math_ops = new HierarchyMathOps("HierarchyMathOps", hierarchy, coarsest_ln, finest_ln);
    d_hier_cc_data_ops = new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, coarsest_ln, finest_ln);

    // Registering temporary variable div U.
    Pointer<CellVariable<NDIM, double> > div_U_var = d_H_var;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_div_U_scratch_idx = var_db->registerVariableAndContext(div_U_var, var_db->getContext(d_object_name + "::div_U"));

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_div_U_scratch_idx, data_time);
    }

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    // Deallocate and remove scratch/smooth phi.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_div_U_scratch_idx);
    }
    var_db->removePatchDataIndex(d_div_U_scratch_idx);

    return;
}

void
HeavisideForcingFunction::setDataOnPatch(const int data_idx,
                                         Pointer<Variable<NDIM> > /*var*/,
                                         Pointer<Patch<NDIM> > /*patch*/,
                                         const double data_time,
                                         const bool initial_time,
                                         Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    if (initial_time) return;

    // Compute H*div U which is to be added in Heaviside transport equation.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_new_idx = var_db->mapVariableAndContextToIndex(d_U_adv_var, d_adv_diff_solver->getNewContext());
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());

    d_hier_math_ops->div(d_div_U_scratch_idx,
                         d_H_var,
                         1.0,
                         U_new_idx,
                         d_U_adv_var,
                         nullptr,
                         data_time,
                         /*synch_cf_bdry*/ false);

    // Multiply H with div U.
    d_hier_cc_data_ops->multiply(d_div_U_scratch_idx, H_new_idx, d_div_U_scratch_idx);
    d_hier_cc_data_ops->copyData(data_idx, d_div_U_scratch_idx);

    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
