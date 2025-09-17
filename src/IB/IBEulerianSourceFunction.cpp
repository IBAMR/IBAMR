// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include "ibamr/IBHierarchyIntegrator.h"

#include "ibtk/CartGridFunction.h"

#include "CellData.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataBasicOps.h"
#include "tbox/Pointer.h"

#include <string>

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
template <int DIM>
class PatchLevel;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBHierarchyIntegrator::IBEulerianSourceFunction::IBEulerianSourceFunction(IBHierarchyIntegrator* ib_solver)
    : CartGridFunction(ib_solver->getName() + "::IBEulerianSourceFunction"), d_ib_solver(ib_solver)
{
    // intentionally blank
    return;
} // IBEulerianSourceFunction

bool
IBHierarchyIntegrator::IBEulerianSourceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
IBHierarchyIntegrator::IBEulerianSourceFunction::setDataOnPatchHierarchy(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > /*var*/,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const double data_time,
    const bool initial_time,
    const int coarsest_ln_in,
    const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == IBTK::invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        TBOX_ASSERT(level->checkAllocated(data_idx));
    }

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_ops(hierarchy, coarsest_ln, finest_ln);
    hier_cc_ops.setToScalar(data_idx, 0.0, false);

    // At the initial time the structures may not yet be set up, so we cannot
    // call computeFluidSources()
    if (initial_time) return;

    d_ib_solver->computeFluidSources(data_idx, data_time);

    return;
}

void
IBHierarchyIntegrator::IBEulerianSourceFunction::setDataOnPatch(const int data_idx,
                                                                Pointer<Variable<NDIM> > /*var*/,
                                                                Pointer<Patch<NDIM> > patch,
                                                                const double /*data_time*/,
                                                                const bool initial_time,
                                                                Pointer<PatchLevel<NDIM> > /*level*/)
{
    // This function is called during initialization, but at that point we are
    // not guaranteed that Lagrangian data is set up in a way that we can
    // actually compute fluid sources: in that case zero everything out.
    if (initial_time)
    {
        Pointer<CellData<NDIM, double> > q_cc_data = patch->getPatchData(data_idx);
        TBOX_ASSERT(q_cc_data);
        q_cc_data->fillAll(0.0);
    }
    else
    {
        TBOX_ERROR("Use setDataOnPatchHierarchy() with this class instead.");
    }
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
