// Filename: IBEulerianSourceFunction.cpp
// Created on 18 Jun 2005 by Boyce Griffith
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
#include "ibamr/IBHierarchyIntegrator.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataBasicOps.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

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

IBHierarchyIntegrator::IBEulerianSourceFunction::~IBEulerianSourceFunction()
{
    // intentionally blank
    return;
} // ~IBEulerianSourceFunction

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
