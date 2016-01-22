// Filename: IBEulerianForceFunction.cpp
// Created on 28 Sep 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>
#include <string>

#include "CellData.h"
#include "HierarchyDataOpsReal.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataBasicOps.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataBasicOps.h"
#include "SideData.h"
#include "Variable.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

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

IBHierarchyIntegrator::IBEulerianForceFunction::~IBEulerianForceFunction()
{
    // intentionally blank
    return;
} // ~IBEulerianForceFunction

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
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
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
