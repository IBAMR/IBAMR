// Filename: LagMarkerRefine.C
// Created on 04 Oct 2007 by Boyce Griffith
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

#include "LagMarkerRefine.h"

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
#include <ibtk/IndexUtilities.h>
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

const std::string LagMarkerRefine::s_op_name = "LAG_MARKER_REFINE";

namespace
{
static const int REFINE_OP_PRIORITY = 0;
static const int REFINE_OP_STENCIL_WIDTH = 0;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LagMarkerRefine::LagMarkerRefine()
{
    // intentionally blank
    return;
}// LagMarkerRefine

LagMarkerRefine::~LagMarkerRefine()
{
    // intentionally blank
    return;
}// ~LagMarkerRefine

bool
LagMarkerRefine::findRefineOperator(
    const Pointer<Variable<NDIM> >& var,
    const std::string& op_name) const
{
    Pointer<IndexVariable<NDIM,LagMarker,CellGeometry<NDIM> > > mark_var = var;
    return (!mark_var.isNull() && op_name == s_op_name);
}// findRefineOperator

const std::string&
LagMarkerRefine::getOperatorName() const
{
    return s_op_name;
}// getOperatorName

int
LagMarkerRefine::getOperatorPriority() const
{
    return REFINE_OP_PRIORITY;
}// getOperatorPriority

IntVector<NDIM>
LagMarkerRefine::getStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
}// getStencilWidth

void
LagMarkerRefine::refine(
    Patch<NDIM>& fine,
    const Patch<NDIM>& coarse,
    const int dst_component,
    const int src_component,
    const Box<NDIM>& fine_box,
    const IntVector<NDIM>& ratio) const
{
    Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > dst_mark_data = fine  .getPatchData(dst_component);
    Pointer<IndexData<NDIM,LagMarker,CellGeometry<NDIM> > > src_mark_data = coarse.getPatchData(src_component);

    const Box<NDIM>& fine_patch_box = fine.getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > fine_patch_geom = fine.getPatchGeometry();
    const Index<NDIM>& fine_patch_lower = fine_patch_box.lower();
    const Index<NDIM>& fine_patch_upper = fine_patch_box.upper();
    const double* const fine_patchXLower = fine_patch_geom->getXLower();
    const double* const fine_patchXUpper = fine_patch_geom->getXUpper();
    const double* const fine_patchDx = fine_patch_geom->getDx();

    //const Box<NDIM>& coarse_patch_box = coarse.getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > coarse_patch_geom = coarse.getPatchGeometry();
    //const Index<NDIM>& coarse_patch_lower = coarse_patch_box.lower();
    //const Index<NDIM>& coarse_patch_upper = coarse_patch_box.upper();
    //const double* const coarse_patchXLower = coarse_patch_geom->getXLower();
    //const double* const coarse_patchXUpper = coarse_patch_geom->getXUpper();
    const double* const coarse_patchDx = coarse_patch_geom->getDx();

    const Box<NDIM> coarse_box = Box<NDIM>::coarsen(fine_box,ratio);
    for (IndexData<NDIM,LagMarker,CellGeometry<NDIM> >::Iterator it(*src_mark_data); it; it++)
    {
        const Index<NDIM>& coarse_i = it.getIndex();
        if (coarse_box.contains(coarse_i))
        {
            const LagMarker& coarse_mark = it();
            const std::vector<double>& coarse_X = coarse_mark.getPositions();
            const std::vector<double>& coarse_U = coarse_mark.getVelocities();
            const std::vector<int>& coarse_idx = coarse_mark.getIndices();
            const IntVector<NDIM>& coarse_offset = coarse_mark.getPeriodicOffset();
            double X_shifted[NDIM];
            for (int k = 0; k < coarse_mark.getNumberOfMarkers(); ++k)
            {
                const double* const X = &coarse_X[NDIM*k];
                const double* const U = &coarse_U[NDIM*k];
                const int& idx = coarse_idx[k];

                for (int d = 0; d < NDIM; ++d)
                {
                    X_shifted[d] = X[d] + double(coarse_offset(d))*coarse_patchDx[d];
                }

                const Index<NDIM> fine_i = IndexUtilities::getCellIndex(X_shifted, fine_patchXLower, fine_patchXUpper, fine_patchDx, fine_patch_lower, fine_patch_upper);
                if (fine_box.contains(fine_i))
                {
                    if (!dst_mark_data->isElement(fine_i))
                    {
                        dst_mark_data->appendItem(fine_i, LagMarker());
                    }
                    LagMarker& fine_mark = *(dst_mark_data->getItem(fine_i));
                    std::vector<double>& fine_X = fine_mark.getPositions();
                    std::vector<double>& fine_U = fine_mark.getVelocities();
                    std::vector<int>& fine_idx = fine_mark.getIndices();

                    fine_X.insert(fine_X.end(),X_shifted,X_shifted+NDIM);
                    fine_U.insert(fine_U.end(),U,U+NDIM);
                    fine_idx.push_back(idx);
                }
            }
        }
    }
    return;
}// refine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::LagMarkerRefine>;

//////////////////////////////////////////////////////////////////////////////
