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
                                                  Pointer<Variable<NDIM> > /*var*/,
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
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);

    // NOTE: At the initial time we take div u = 0 for the lack of knowledge of the
    // velocity field.
    if (initial_time)
    {
        hier_cc_data_ops.setToScalar(data_idx, 0.0);
        return;
    }

    // Compute H*div U which is to be added in Heaviside transport equation.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int U_new_idx = var_db->mapVariableAndContextToIndex(d_U_adv_var, d_adv_diff_solver->getNewContext());
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());

    HierarchyMathOps hier_math_ops("HierarchyMathOps", hierarchy, coarsest_ln, finest_ln);
    hier_math_ops.div(data_idx,
                      d_H_var,
                      1.0,
                      U_new_idx,
                      d_U_adv_var,
                      nullptr,
                      data_time,
                      /*synch_cf_bdry*/ false);

    // Multiply H with div U.
    hier_cc_data_ops.multiply(data_idx, H_new_idx, data_idx);

    return;
} // setDataOnPatchHierarchy

void
HeavisideForcingFunction::setDataOnPatch(const int /*data_idx*/,
                                         Pointer<Variable<NDIM> > /*var*/,
                                         Pointer<Patch<NDIM> > /*patch*/,
                                         const double /*data_time*/,
                                         const bool /*initial_time*/,
                                         Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    // As we directly compute and set data on the patch hierarchy don't do anything over here.
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
