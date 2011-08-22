// Filename: VecCellRefineAdapter.C
// Created on 09 Apr 2010 by Boyce Griffith
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

#include "VecCellRefineAdapter.h"

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
#include <ibtk/VecCellData.h>
#include <ibtk/VecCellVariable.h>
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VecCellRefineAdapter::VecCellRefineAdapter(
    Pointer<RefineOperator<NDIM> > cell_refine_op)
    : RefineOperator<NDIM>(),
      d_cell_refine_op(cell_refine_op)
{
    // intentionally blank
    return;
}// VecCellRefineAdapter

VecCellRefineAdapter::~VecCellRefineAdapter()
{
    // intentionally blank
    return;
}// ~VecCellRefineAdapter

bool
VecCellRefineAdapter::findRefineOperator(
    const Pointer<Variable<NDIM> >& var,
    const std::string &op_name) const
{
    const Pointer<VecCellVariable<double> > cast_var(var);
    if (!cast_var.isNull() && (op_name == getOperatorName()))
    {
        return true;
    }
    else
    {
        return false;
    }
}// findRefineOperator

const std::string&
VecCellRefineAdapter::getOperatorName() const
{
    return d_cell_refine_op->getOperatorName();
}// getOperatorName

int
VecCellRefineAdapter::getOperatorPriority() const
{
    return d_cell_refine_op->getOperatorPriority();
}// getOperatorPriority

IntVector<NDIM>
VecCellRefineAdapter::getStencilWidth() const
{
    return d_cell_refine_op->getStencilWidth();
}// getStencilWidth

void
VecCellRefineAdapter::refine(
    Patch<NDIM>& fine,
    const Patch<NDIM>& coarse,
    const int dst_component,
    const int src_component,
    const Box<NDIM>& fine_box,
    const IntVector<NDIM>& ratio) const
{
    Pointer<VecCellData<double> > dst_data = fine  .getPatchData(dst_component);
    Pointer<VecCellData<double> > src_data = coarse.getPatchData(src_component);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!dst_data.isNull());
    TBOX_ASSERT(!src_data.isNull());
#endif
    // Create "dummy" patches.
    Patch<NDIM> fine_cell(fine.getBox(), fine.getPatchDescriptor());
    fine_cell.setPatchGeometry(fine.getPatchGeometry());
    fine_cell.setPatchInHierarchy(fine.inHierarchy());
    fine_cell.setPatchLevelNumber(fine.getPatchLevelNumber());
    fine_cell.setPatchNumber(fine.getPatchNumber());

    Patch<NDIM> coarse_cell(coarse.getBox(), coarse.getPatchDescriptor());
    coarse_cell.setPatchGeometry(coarse.getPatchGeometry());
    coarse_cell.setPatchInHierarchy(coarse.inHierarchy());
    coarse_cell.setPatchLevelNumber(coarse.getPatchLevelNumber());
    coarse_cell.setPatchNumber(coarse.getPatchNumber());

    // Make copies of the dst and src data.
    CellData<NDIM,double> dst_cell_data(dst_data->getBox(), dst_data->getDepth(), dst_data->getGhostCellWidth());
    dst_data->copy2(dst_cell_data);
    fine_cell.allocatePatchData(dst_component);
    fine_cell.setPatchData(dst_component,Pointer<PatchData<NDIM> >(&dst_cell_data,false));

    CellData<NDIM,double> src_cell_data(src_data->getBox(), src_data->getDepth(), src_data->getGhostCellWidth());
    src_data->copy2(src_cell_data);
    coarse_cell.allocatePatchData(src_component);
    coarse_cell.setPatchData(src_component,Pointer<PatchData<NDIM> >(&src_cell_data,false));

    // Refine data from the coarse dummy patch to the fine dummy patch.
    d_cell_refine_op->refine(fine_cell, coarse_cell, dst_component, src_component, fine_box, ratio);

    // Copy the result into the VecCellData.
    dst_data->copy(dst_cell_data);
    return;
}// refine

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
