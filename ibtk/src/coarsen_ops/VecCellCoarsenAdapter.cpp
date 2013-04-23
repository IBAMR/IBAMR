// Filename: VecCellCoarsenAdapter.cpp
// Created on 09 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>

#include "CellData.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchDescriptor.h"
#include "PatchGeometry.h"
#include "SAMRAI_config.h"
#include "VecCellCoarsenAdapter.h"
#include "ibtk/VecCellData.h"
#include "ibtk/VecCellData-inl.h"
#include "ibtk/VecCellVariable.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Utilities.h"

namespace SAMRAI {
namespace hier {
template <int DIM> class Variable;
}  // namespace hier
}  // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VecCellCoarsenAdapter::VecCellCoarsenAdapter(
    Pointer<CoarsenOperator<NDIM> > cell_coarsen_op)
    : CoarsenOperator<NDIM>(),
      d_cell_coarsen_op(cell_coarsen_op)
{
    // intentionally blank
    return;
}// VecCellCoarsenAdapter

VecCellCoarsenAdapter::~VecCellCoarsenAdapter()
{
    // intentionally blank
    return;
}// ~VecCellCoarsenAdapter

bool
VecCellCoarsenAdapter::findCoarsenOperator(
    const Pointer<Variable<NDIM> >& var,
    const std::string &op_name) const
{
    const Pointer<VecCellVariable<double> > cast_var(var);
    if (cast_var && (op_name == getOperatorName()))
    {
        return true;
    }
    else
    {
        return false;
    }
}// findCoarsenOperator

const std::string&
VecCellCoarsenAdapter::getOperatorName() const
{
    return d_cell_coarsen_op->getOperatorName();
}// getOperatorName

int
VecCellCoarsenAdapter::getOperatorPriority() const
{
    return d_cell_coarsen_op->getOperatorPriority();
}// getOperatorPriority

IntVector<NDIM>
VecCellCoarsenAdapter::getStencilWidth() const
{
    return d_cell_coarsen_op->getStencilWidth();
}// getStencilWidth

void
VecCellCoarsenAdapter::coarsen(
    Patch<NDIM>& coarse,
    const Patch<NDIM>& fine,
    const int dst_component,
    const int src_component,
    const Box<NDIM>& coarse_box,
    const IntVector<NDIM>& ratio) const
{
    Pointer<VecCellData<double> > dst_data = coarse.getPatchData(dst_component);
    Pointer<VecCellData<double> > src_data = fine  .getPatchData(src_component);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst_data);
    TBOX_ASSERT(src_data);
#endif
    // Create "dummy" patches.
    Patch<NDIM> coarse_cell(coarse.getBox(), coarse.getPatchDescriptor());
    coarse_cell.setPatchGeometry(coarse.getPatchGeometry());
    coarse_cell.setPatchInHierarchy(coarse.inHierarchy());
    coarse_cell.setPatchLevelNumber(coarse.getPatchLevelNumber());
    coarse_cell.setPatchNumber(coarse.getPatchNumber());

    Patch<NDIM> fine_cell(fine.getBox(), fine.getPatchDescriptor());
    fine_cell.setPatchGeometry(fine.getPatchGeometry());
    fine_cell.setPatchInHierarchy(fine.inHierarchy());
    fine_cell.setPatchLevelNumber(fine.getPatchLevelNumber());
    fine_cell.setPatchNumber(fine.getPatchNumber());

    // Make copies of the dst and src data.
    CellData<NDIM,double> dst_cell_data(dst_data->getBox(), dst_data->getDepth(), dst_data->getGhostCellWidth());
    dst_data->copy2(dst_cell_data);
    coarse_cell.allocatePatchData(dst_component);
    coarse_cell.setPatchData(dst_component,Pointer<PatchData<NDIM> >(&dst_cell_data,false));

    CellData<NDIM,double> src_cell_data(src_data->getBox(), src_data->getDepth(), src_data->getGhostCellWidth());
    src_data->copy2(src_cell_data);
    fine_cell.allocatePatchData(src_component);
    fine_cell.setPatchData(src_component,Pointer<PatchData<NDIM> >(&src_cell_data,false));

    // Coarsen data from the coarse dummy patch to the fine dummy patch.
    d_cell_coarsen_op->coarsen(coarse_cell, fine_cell, dst_component, src_component, coarse_box, ratio);

    // Copy the result into the VecCellData.
    dst_data->copy(dst_cell_data);
    return;
}// coarsen

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
