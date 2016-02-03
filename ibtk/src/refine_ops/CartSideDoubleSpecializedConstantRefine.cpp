// Filename: CartSideDoubleSpecializedConstantRefine.cpp
// Created on 17 Sep 2011 by Boyce Griffith
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

#include "Box.h"
#include "IBTK_config.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "SideData.h"
#include "SideVariable.h"
#include "ibtk/CartSideDoubleSpecializedConstantRefine.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CART_SIDE_SPECIALIZED_CONSTANT_REFINE_FC                                                                       \
    IBTK_FC_FUNC(cart_side_specialized_constant_refine2d, CART_SIDE_SPECIALIZED_CONSTANT_REFINE2D)
#endif
#if (NDIM == 3)
#define CART_SIDE_SPECIALIZED_CONSTANT_REFINE_FC                                                                       \
    IBTK_FC_FUNC(cart_side_specialized_constant_refine3d, CART_SIDE_SPECIALIZED_CONSTANT_REFINE3D)
#endif

// Function interfaces
extern "C" {
void CART_SIDE_SPECIALIZED_CONSTANT_REFINE_FC(
#if (NDIM == 2)
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
#endif
#if (NDIM == 3)
    double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
#endif
    const int*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string CartSideDoubleSpecializedConstantRefine::s_op_name = "SPECIALIZED_CONSTANT_REFINE";

namespace
{
static const int REFINE_OP_PRIORITY = 0;
static const int REFINE_OP_STENCIL_WIDTH = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleSpecializedConstantRefine::CartSideDoubleSpecializedConstantRefine()
{
    // intentionally blank
    return;
} // CartSideDoubleSpecializedConstantRefine

CartSideDoubleSpecializedConstantRefine::~CartSideDoubleSpecializedConstantRefine()
{
    // intentionally blank
    return;
} // ~CartSideDoubleSpecializedConstantRefine

bool
CartSideDoubleSpecializedConstantRefine::findRefineOperator(const Pointer<Variable<NDIM> >& var,
                                                            const std::string& op_name) const
{
    const Pointer<SideVariable<NDIM, double> > sc_var = var;
    return (sc_var && op_name == s_op_name);
} // findRefineOperator

const std::string&
CartSideDoubleSpecializedConstantRefine::getOperatorName() const
{
    return s_op_name;
} // getOperatorName

int
CartSideDoubleSpecializedConstantRefine::getOperatorPriority() const
{
    return REFINE_OP_PRIORITY;
} // getOperatorPriority

IntVector<NDIM>
CartSideDoubleSpecializedConstantRefine::getStencilWidth() const
{
    return REFINE_OP_STENCIL_WIDTH;
} // getStencilWidth

void
CartSideDoubleSpecializedConstantRefine::refine(Patch<NDIM>& fine,
                                                const Patch<NDIM>& coarse,
                                                const int dst_component,
                                                const int src_component,
                                                const Box<NDIM>& fine_box,
                                                const IntVector<NDIM>& ratio) const
{
    // Get the patch data.
    Pointer<SideData<NDIM, double> > fdata = fine.getPatchData(dst_component);
    Pointer<SideData<NDIM, double> > cdata = coarse.getPatchData(src_component);
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata);
    TBOX_ASSERT(cdata);
    TBOX_ASSERT(fdata->getDepth() == cdata->getDepth());
#endif
    const int data_depth = fdata->getDepth();

    const Box<NDIM>& fdata_box = fdata->getBox();
    const int fdata_gcw = fdata->getGhostCellWidth().max();
#if !defined(NDEBUG)
    TBOX_ASSERT(fdata_gcw == fdata->getGhostCellWidth().min());
#endif

    const Box<NDIM>& cdata_box = cdata->getBox();
    const int cdata_gcw = cdata->getGhostCellWidth().max();
#if !defined(NDEBUG)
    TBOX_ASSERT(cdata_gcw == cdata->getGhostCellWidth().min());
#endif

    // Refine the data.
    const Box<NDIM> fill_box = Box<NDIM>::refine(Box<NDIM>::coarsen(fine_box, ratio), ratio);
    for (int depth = 0; depth < data_depth; ++depth)
    {
        CART_SIDE_SPECIALIZED_CONSTANT_REFINE_FC(
#if (NDIM == 2)
            fdata->getPointer(0, depth),
            fdata->getPointer(1, depth),
            fdata_gcw,
            fdata_box.lower()(0),
            fdata_box.upper()(0),
            fdata_box.lower()(1),
            fdata_box.upper()(1),
            cdata->getPointer(0, depth),
            cdata->getPointer(1, depth),
            cdata_gcw,
            cdata_box.lower()(0),
            cdata_box.upper()(0),
            cdata_box.lower()(1),
            cdata_box.upper()(1),
            fine_box.lower()(0),
            fine_box.upper()(0),
            fine_box.lower()(1),
            fine_box.upper()(1),
            fill_box.lower()(0),
            fill_box.upper()(0),
            fill_box.lower()(1),
            fill_box.upper()(1),
#endif
#if (NDIM == 3)
            fdata->getPointer(0, depth),
            fdata->getPointer(1, depth),
            fdata->getPointer(2, depth),
            fdata_gcw,
            fdata_box.lower()(0),
            fdata_box.upper()(0),
            fdata_box.lower()(1),
            fdata_box.upper()(1),
            fdata_box.lower()(2),
            fdata_box.upper()(2),
            cdata->getPointer(0, depth),
            cdata->getPointer(1, depth),
            cdata->getPointer(2, depth),
            cdata_gcw,
            cdata_box.lower()(0),
            cdata_box.upper()(0),
            cdata_box.lower()(1),
            cdata_box.upper()(1),
            cdata_box.lower()(2),
            cdata_box.upper()(2),
            fine_box.lower()(0),
            fine_box.upper()(0),
            fine_box.lower()(1),
            fine_box.upper()(1),
            fine_box.lower()(2),
            fine_box.upper()(2),
            fill_box.lower()(0),
            fill_box.upper()(0),
            fill_box.lower()(1),
            fill_box.upper()(1),
            fill_box.lower()(2),
            fill_box.upper()(2),
#endif
            ratio);
    }
    return;
} // refine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
