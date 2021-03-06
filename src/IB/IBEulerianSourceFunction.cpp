// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

IBHierarchyIntegrator::IBEulerianSourceFunction::IBEulerianSourceFunction(const IBHierarchyIntegrator* const ib_solver)
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
IBHierarchyIntegrator::IBEulerianSourceFunction::setDataOnPatch(const int data_idx,
                                                                Pointer<Variable<NDIM> > /*var*/,
                                                                Pointer<Patch<NDIM> > patch,
                                                                const double /*data_time*/,
                                                                const bool initial_time,
                                                                Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<CellData<NDIM, double> > q_cc_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(q_cc_data);
#endif
    q_cc_data->fillAll(0.0);
    if (initial_time) return;
    Pointer<CellData<NDIM, double> > q_ib_cc_data = patch->getPatchData(d_ib_solver->d_q_idx);
    PatchCellDataBasicOps<NDIM, double> patch_ops;
    patch_ops.add(q_cc_data, q_cc_data, q_ib_cc_data, patch->getBox());
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
