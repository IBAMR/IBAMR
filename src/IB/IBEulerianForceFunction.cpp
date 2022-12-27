// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataBasicOps.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataBasicOps.h"
#include "SideData.h"
#include "Variable.h"
#include "tbox/Pointer.h"

#include <string>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBHierarchyIntegrator::IBEulerianForceFunction::IBEulerianForceFunction(const IBHierarchyIntegrator* const ib_solver)
    : CartGridFunction(ib_solver->getName() + "::IBEulerianForceFunction"), d_ib_solver(ib_solver)
{
    // intentionally blank
    return;
} // IBEulerianForceFunction

bool
IBHierarchyIntegrator::IBEulerianForceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
IBHierarchyIntegrator::IBEulerianForceFunction::setDataOnPatchHierarchy(const int data_idx,
                                                                        Pointer<Variable<NDIM> > var,
                                                                        Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                        const double data_time,
                                                                        const bool initial_time,
                                                                        const int coarsest_ln_in,
                                                                        const int finest_ln_in)
{
    if (initial_time)
    {
        d_ib_solver->d_hier_velocity_data_ops->setToScalar(data_idx, 0.0);
        return;
    }
    if (d_ib_solver->d_body_force_fcn)
    {
        d_ib_solver->d_body_force_fcn->setDataOnPatchHierarchy(
            data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);
    }
    else
    {
        d_ib_solver->d_hier_velocity_data_ops->setToScalar(data_idx, 0.0);
    }
    const int coarsest_ln = (coarsest_ln_in == IBTK::invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    return;
} // setDataOnPatchHierarchy

void
IBHierarchyIntegrator::IBEulerianForceFunction::setDataOnPatch(const int data_idx,
                                                               Pointer<Variable<NDIM> > /*var*/,
                                                               Pointer<Patch<NDIM> > patch,
                                                               const double /*data_time*/,
                                                               const bool initial_time,
                                                               Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<PatchData<NDIM> > f_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(f_data);
#endif
    Pointer<CellData<NDIM, double> > f_cc_data = f_data;
    Pointer<SideData<NDIM, double> > f_sc_data = f_data;
#if !defined(NDEBUG)
    TBOX_ASSERT(f_cc_data || f_sc_data);
#endif
    if (initial_time)
    {
        if (f_cc_data) f_cc_data->fillAll(0.0);
        if (f_sc_data) f_sc_data->fillAll(0.0);
        return;
    }
    Pointer<PatchData<NDIM> > f_ib_data = patch->getPatchData(d_ib_solver->d_f_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(f_ib_data);
#endif
    Pointer<CellData<NDIM, double> > f_ib_cc_data = f_ib_data;
    Pointer<SideData<NDIM, double> > f_ib_sc_data = f_ib_data;
#if !defined(NDEBUG)
    TBOX_ASSERT(f_ib_cc_data || f_ib_sc_data);
    TBOX_ASSERT((f_ib_cc_data && f_cc_data) || (f_ib_sc_data && f_sc_data));
#endif
    if (f_cc_data)
    {
        PatchCellDataBasicOps<NDIM, double> patch_ops;
        patch_ops.add(f_cc_data, f_cc_data, f_ib_cc_data, patch->getBox());
    }
    if (f_sc_data)
    {
        PatchSideDataBasicOps<NDIM, double> patch_ops;
        patch_ops.add(f_sc_data, f_sc_data, f_ib_sc_data, patch->getBox());
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
