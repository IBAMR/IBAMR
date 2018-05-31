// Filename: RelaxationLSBcCoefs.cpp
// Created on 25 Feb 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2017, Amneet Bhalla
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

#include "ibamr/RelaxationLSBcCoefs.h"
#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "VariableDatabase.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RelaxationLSBcCoefs::RelaxationLSBcCoefs(const std::string& object_name) : d_object_name(object_name), d_phi_idx(-1)
{
    // intentionally blank
    return;
} // RelaxationLSBcCoefs

RelaxationLSBcCoefs::~RelaxationLSBcCoefs()
{
    // intentionally blank
    return;
} // ~RelaxationLSBcCoefs

void
RelaxationLSBcCoefs::setLSPatchDataIndex(int phi_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(phi_idx >= 0);
#endif
    d_phi_idx = phi_idx;
    return;
} // setLSPatchDataIndex

void
RelaxationLSBcCoefs::resetLSPatchDataIndex()
{
    d_phi_idx = -1;
    return;
} // resetLSPatchDataIndex

void
RelaxationLSBcCoefs::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                const Pointer<Variable<NDIM> >& /*variable*/,
                                const Patch<NDIM>& patch,
                                const BoundaryBox<NDIM>& bdry_box,
                                double /*fill_time*/) const
{
    Pointer<CellData<NDIM, double> > phi_data = patch.getPatchData(d_phi_idx);

    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const bool is_upper = (location_index % 2 == 1);
#if !defined(NDEBUG)
    TBOX_ASSERT(!acoef_data.isNull());
#endif
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
#endif
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const Index<NDIM>& i = bc();
        Index<NDIM> i_bdry = bc();
        Index<NDIM> i_intr = bc();

        if (is_upper)
        {
            i_bdry(axis) -= 1;
            i_intr(axis) -= 2;
        }
        else
        {
            i_intr(axis) += 1;
        }
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);

        a = 1.0;
        b = 0.0;
        g = 1.5 * (*phi_data)(i_bdry, 0) - 0.5 * (*phi_data)(i_intr, 0);
    }
    return;
} // setBcCoefs

IntVector<NDIM>
RelaxationLSBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
