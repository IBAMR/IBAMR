// Filename: CartGridBodyForce.cpp
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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
// POSSIBILITY OF SUCH DAMAGE

#include "CartGridBodyForce.h"

////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

#include "CellData.h"
#include "CellVariable.h"
#include "SideData.h"
#include "SideVariable.h"
#include "VariableDatabase.h"
#include "ibamr/namespaces.h"
#include "tbox/Utilities.h"

namespace IBTK
{
CartGridBodyForce::CartGridBodyForce(const int body_force_idx) : d_body_force_idx(body_force_idx)
{
    // this is intentionally left blank
    return;
} // CartGridBodyForce

bool
CartGridBodyForce::isTimeDependent() const
{
    return true;

} // isTimeDependent

void
CartGridBodyForce::setDataOnPatch(const int data_idx,
                                  Pointer<Variable<NDIM> > var,
                                  Pointer<Patch<NDIM> > patch,
                                  const double /*data_time*/,
                                  const bool /*initial_time*/,
                                  Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(!var.isNull());
#endif

    Pointer<SideVariable<NDIM, double> > copy_to_sc_var = var;
    Pointer<CellVariable<NDIM, double> > copy_to_cc_var = var;

    if (!copy_to_sc_var.isNull())
    {
        Pointer<SideData<NDIM, double> > copy_to_sc_data = patch->getPatchData(data_idx);
        const Pointer<SideData<NDIM, double> > copy_from_sc_data = patch->getPatchData(d_body_force_idx);

#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!copy_to_sc_data.isNull());
        TBOX_ASSERT(!copy_from_sc_data.isNull());
#endif

        copy_to_sc_data->copy(*copy_from_sc_data);
    }
    else if (!copy_to_cc_var.isNull())
    {
        Pointer<CellData<NDIM, double> > copy_to_cc_data = patch->getPatchData(data_idx);
        const Pointer<CellData<NDIM, double> > copy_from_cc_data = patch->getPatchData(d_body_force_idx);

#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!copy_to_cc_data.isNull());
        TBOX_ASSERT(!copy_from_cc_data.isNull());
#endif

        copy_to_cc_data->copy(*copy_from_cc_data);
    }
    else
    {
        TBOX_ERROR("CartGridBodyForce::setDataOnPatch() "
                   << "UNKNOWN DATA TYPE ENCOUNTERED"
                   << std::endl);
    }

    return;

} // setDataonPatch

} // namespace IBTK
