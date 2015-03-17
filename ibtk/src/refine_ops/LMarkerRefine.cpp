// Filename: LMarkerRefine.cpp
// Created on 04 Oct 2007 by Boyce Griffith
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

#include <string>
#include <vector>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellOverlap.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "boost/array.hpp"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LMarker.h"
#include "ibtk/LMarkerRefine.h"
#include "ibtk/LMarkerSet.h"
#include "ibtk/LMarkerSetData.h"
#include "ibtk/LMarkerSetVariable.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{

class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string LMarkerRefine::OP_NAME = "LMARKER_REFINE";

namespace
{
static const int REFINE_OP_PRIORITY = 0;
static const int REFINE_OP_STENCIL_WIDTH = 0;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LMarkerRefine::LMarkerRefine() : RefineOperator(OP_NAME)
{
    // intentionally blank
    return;
}

LMarkerRefine::~LMarkerRefine()
{
    // intentionally blank
    return;
}

int LMarkerRefine::getOperatorPriority() const
{
    return REFINE_OP_PRIORITY;
}

IntVector LMarkerRefine::getStencilWidth(const Dimension& dim) const
{
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
}

void LMarkerRefine::refine(Patch& fine,
                           const Patch& coarse,
                           const int dst_component,
                           const int src_component,
                           const BoxOverlap& fine_overlap,
                           const IntVector& ratio) const
{
    auto dst_mark_data = BOOST_CAST<LMarkerSetData>(fine.getPatchData(dst_component));
    auto src_mark_data = BOOST_CAST<LMarkerSetData>(coarse.getPatchData(src_component));

    const Box& fine_patch_box = fine.getBox();
    auto fine_patch_geom = BOOST_CAST<CartesianPatchGeometry>(fine.getPatchGeometry());
    const Index& fine_patch_lower = fine_patch_box.lower();
    const Index& fine_patch_upper = fine_patch_box.upper();
    const double* const fine_patch_x_lower = fine_patch_geom->getXLower();
    const double* const fine_patch_x_upper = fine_patch_geom->getXUpper();
    const double* const fine_patch_dx = fine_patch_geom->getDx();

    auto coarse_patch_geom = BOOST_CAST<CartesianPatchGeometry>(coarse.getPatchGeometry());
    const double* const coarse_patch_dx = coarse_patch_geom->getDx();

    auto fine_cell_overlap = CPP_CAST<const CellOverlap*>(&fine_overlap);
    const BoxContainer& fine_boxes = fine_cell_overlap->getDestinationBoxContainer();
    for (auto bl = fine_boxes.begin(), e = fine_boxes.end(); bl != e; ++bl)
    {
        const Box& fine_box = *bl;
        const Box coarse_box = Box::coarsen(fine_box, ratio);
        const Box fill_box = Box::refine(Box::coarsen(fine_box, ratio), ratio);
        for (LMarkerSetData::SetIterator it(*src_mark_data, /*begin*/ true), e(*src_mark_data, /*begin*/ false); it != e; ++it)
        {
            const Index& coarse_i = it.getIndex();
            if (coarse_box.contains(coarse_i))
            {
                const LMarkerSet& coarse_mark_set = *it;
                for (auto cit = coarse_mark_set.begin(); cit != coarse_mark_set.end(); ++cit)
                {
                    const LMarkerSet::value_type& coarse_mark = *cit;
                    const Point& X = coarse_mark->getPosition();
                    const IntVector& offset = coarse_mark->getPeriodicOffset();
                    boost::array<double, NDIM> X_shifted;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X_shifted[d] = X[d] + static_cast<double>(offset(d)) * coarse_patch_dx[d];
                    }

                    const Index fine_i =
                        IndexUtilities::getCellIndex(X_shifted, fine_patch_x_lower, fine_patch_x_upper, fine_patch_dx,
                                                     fine_patch_lower, fine_patch_upper);
                    if (fine_box.contains(fine_i))
                    {
                        if (!dst_mark_data->isElement(fine_i))
                        {
                            dst_mark_data->appendItemPointer(fine_i, new LMarkerSet());
                        }
                        LMarkerSet& fine_mark_set = *(dst_mark_data->getItem(fine_i));
                        fine_mark_set.push_back(coarse_mark);
                    }
                }
            }
        }
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
