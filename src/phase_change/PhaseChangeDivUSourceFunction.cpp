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

#include "ibamr/PhaseChangeDivUSourceFunction.h"
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

PhaseChangeDivUSourceFunction::PhaseChangeDivUSourceFunction(
    const std::string& object_name,
    const Pointer<PhaseChangeHierarchyIntegrator> pc_hier_integrator)
    : d_object_name(object_name), d_pc_hier_integrator(pc_hier_integrator)
{
    // intentionally blank
    return;
} // PhaseChangeDivUSourceFunction

bool
PhaseChangeDivUSourceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
PhaseChangeDivUSourceFunction::setDataOnPatchHierarchy(const int data_idx,
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
    d_hier_cc_data_ops = new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, coarsest_ln, finest_ln);

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);
    return;
}

void
PhaseChangeDivUSourceFunction::setDataOnPatch(const int data_idx,
                                              Pointer<Variable<NDIM> > var,
                                              Pointer<Patch<NDIM> > patch,
                                              const double data_time,
                                              const bool initial_time,
                                              Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    if (initial_time) return;

    // set -Div U = -S where source term S is computed from PhaseChangeHierarchyIntegrator.
    const int S_idx = d_pc_hier_integrator->getDivergenceVelocitySourceTermIndex();

    // Set Div_U_F_idx = - S_idx.
    d_hier_cc_data_ops->scale(data_idx, -1.0, S_idx);

    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
