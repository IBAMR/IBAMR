// Filename: LagMarkerCoarsen.C
// Created on 30 Sep 2006 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
// POSSIBILITY OF SUCH DAMAGE.

#include "LagMarkerCoarsen.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/LagMarker.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <CellGeometry.h>
#include <IndexData.h>
#include <IndexVariable.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string LagMarkerCoarsen::s_op_name = "LAG_MARKER_COARSEN";

namespace
{
static const int COARSEN_OP_PRIORITY = 0;
static const int COARSEN_OP_STENCIL_WIDTH = 0;

inline int
coarsen(
    const int index,
    const int ratio)
{
    return (index < 0 ? (index+1)/ratio-1 : index/ratio);
}// coarsen

inline Index<NDIM>
coarsen_index(
    const Index<NDIM>& i,
    const IntVector<NDIM>& ratio)
{
    Index<NDIM> coarse_i;
    for (int d = 0; d < NDIM; ++d)
    {
        coarse_i(d) = coarsen(i(d), ratio(d));
    }
    return coarse_i;
}// coarsen_index
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LagMarkerCoarsen::LagMarkerCoarsen()
{
    // intentionally blank
    return;
}// LagMarkerCoarsen

LagMarkerCoarsen::~LagMarkerCoarsen()
{
    // intentionally blank
    return;
}// ~LagMarkerCoarsen

bool
LagMarkerCoarsen::findCoarsenOperator(
    const Pointer<Variable<NDIM> >& var,
    const std::string &op_name) const
{
    Pointer<IndexVariable<NDIM,LagMarker,CellGeometry<NDIM> > > mark_var = var;
    return (!mark_var.isNull() && op_name == s_op_name);
}// findCoarsenOperator

const std::string&
LagMarkerCoarsen::getOperatorName() const
{
    return s_op_name;
}// getOperatorName

int
LagMarkerCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
}// getOperatorPriority

IntVector<NDIM>
LagMarkerCoarsen::getStencilWidth() const
{
    return COARSEN_OP_STENCIL_WIDTH;
}// getStencilWidth

void
LagMarkerCoarsen::coarsen(
    Patch<NDIM>& coarse,
    const Patch<NDIM>& fine,
    const int dst_component,
    const int src_component,
    const Box<NDIM>& coarse_box,
    const IntVector<NDIM>& ratio) const
{
    Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > dst_mark_data = coarse.getPatchData(dst_component);
    Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > src_mark_data = fine  .getPatchData(src_component);

    const Box<NDIM> fine_box = Box<NDIM>::refine(coarse_box,ratio);
    for (IndexData<NDIM,LagMarker,CellGeometry<NDIM> >::Iterator it(*src_mark_data); it; it++)
    {
        const Index<NDIM>& fine_i = it.getIndex();
        const Index<NDIM> coarse_i = coarsen_index(fine_i,ratio);
        if (fine_box.contains(fine_i) && coarse_box.contains(coarse_i))
        {
            if (!dst_mark_data->isElement(coarse_i))
            {
                dst_mark_data->appendItem(coarse_i, LagMarker());
            }
            LagMarker& coarse_mark = *(dst_mark_data->getItem(coarse_i));
            const LagMarker& fine_mark = it();
            coarse_mark.addMarker(fine_mark);
        }
    }
    return;
}// coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::LagMarkerCoarsen>;

//////////////////////////////////////////////////////////////////////////////
