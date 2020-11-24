// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#include "ibtk/LMarkerCoarsen.h"
#include "ibtk/LMarkerSet.h"
#include "ibtk/LMarkerSetData.h"
#include "ibtk/LMarkerSetVariable.h"

#include "Box.h"
#include "Patch.h"
#include "tbox/Pointer.h"

#include <algorithm>
#include <string>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string LMarkerCoarsen::s_op_name = "LMARKER_COARSEN";

namespace
{
static const int COARSEN_OP_PRIORITY = 0;
static const int COARSEN_OP_STENCIL_WIDTH = 0;

inline int
coarsen(const int index, const int ratio)
{
    return (index < 0 ? (index + 1) / ratio - 1 : index / ratio);
} // coarsen

inline hier::Index<NDIM>
coarsen_index(const hier::Index<NDIM>& i, const IntVector<NDIM>& ratio)
{
    hier::Index<NDIM> coarse_i;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        coarse_i(d) = coarsen(i(d), ratio(d));
    }
    return coarse_i;
} // coarsen_index
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
LMarkerCoarsen::findCoarsenOperator(const Pointer<Variable<NDIM> >& var, const std::string& op_name) const
{
    Pointer<LMarkerSetVariable> mark_var = var;
    return (mark_var && op_name == s_op_name);
} // findCoarsenOperator

const std::string&
LMarkerCoarsen::getOperatorName() const
{
    return s_op_name;
} // getOperatorName

int
LMarkerCoarsen::getOperatorPriority() const
{
    return COARSEN_OP_PRIORITY;
} // getOperatorPriority

IntVector<NDIM>
LMarkerCoarsen::getStencilWidth() const
{
    return COARSEN_OP_STENCIL_WIDTH;
} // getStencilWidth

void
LMarkerCoarsen::coarsen(Patch<NDIM>& coarse,
                        const Patch<NDIM>& fine,
                        const int dst_component,
                        const int src_component,
                        const Box<NDIM>& coarse_box,
                        const IntVector<NDIM>& ratio) const
{
    Pointer<LMarkerSetData> dst_mark_data = coarse.getPatchData(dst_component);
    Pointer<LMarkerSetData> src_mark_data = fine.getPatchData(src_component);

    const Box<NDIM> fine_box = Box<NDIM>::refine(coarse_box, ratio);
    for (LMarkerSetData::SetIterator it(*src_mark_data); it; it++)
    {
        const hier::Index<NDIM>& fine_i = it.getIndex();
        const hier::Index<NDIM> coarse_i = coarsen_index(fine_i, ratio);
        if (fine_box.contains(fine_i) && coarse_box.contains(coarse_i))
        {
            const LMarkerSet& fine_mark_set = it();
            if (!dst_mark_data->isElement(coarse_i))
            {
                dst_mark_data->appendItemPointer(coarse_i, new LMarkerSet());
            }
            LMarkerSet& coarse_mark_set = *(dst_mark_data->getItem(coarse_i));
            coarse_mark_set.insert(coarse_mark_set.end(), fine_mark_set.begin(), fine_mark_set.end());
        }
    }
    return;
} // coarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
