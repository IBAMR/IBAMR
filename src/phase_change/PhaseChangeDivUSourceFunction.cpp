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
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/PhaseChangeDivUSourceFunction.h"
#include "ibamr/PhaseChangeHierarchyIntegrator.h"

#include "SAMRAIVectorReal.h"

#include "ibamr/app_namespaces.h"

namespace IBAMR
{
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

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);
    return;
}

void
PhaseChangeDivUSourceFunction::setDataOnPatch(const int data_idx,
                                              Pointer<Variable<NDIM> > var,
                                              Pointer<Patch<NDIM> > patch,
                                              const double /*data_time*/,
                                              const bool initial_time,
                                              Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
#if !defined(NDEBUG)
    Pointer<CellVariable<NDIM, double> > cc_var = var;
    TBOX_ASSERT(cc_var);
#else
    NULL_USE(var);
#endif

    Pointer<CellData<NDIM, double> > div_u_cc_data = patch->getPatchData(data_idx);
    if (initial_time)
    {
        div_u_cc_data->fill(0.0);
        return;
    }

    // Set Div U = S where source term S is computed from PhaseChangeHierarchyIntegrator.
    const int S_idx = d_pc_hier_integrator->getVelocityDivergencePatchDataIndex();
    Pointer<CellData<NDIM, double> > S_cc_data = patch->getPatchData(S_idx);
    PatchCellDataOpsReal<NDIM, double> patch_cc_data_ops;
    patch_cc_data_ops.copyData(div_u_cc_data, S_cc_data, patch->getBox());

    return;
} // setDataOnPatch

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
