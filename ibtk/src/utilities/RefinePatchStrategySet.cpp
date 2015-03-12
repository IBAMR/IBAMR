// Filename: RefinePatchStrategySet.cpp
// Created on 11 Sep 2006 by Boyce Griffith
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

#include "ibtk/RefinePatchStrategySet.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{

class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RefinePatchStrategySet::~RefinePatchStrategySet()
{
    // intentionally blank
    return;
}

void RefinePatchStrategySet::setPhysicalBoundaryConditions(Patch& patch,
                                                           const double fill_time,
                                                           const IntVector& ghost_width_to_fill)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    }
    return;
}

IntVector RefinePatchStrategySet::getRefineOpStencilWidth(const Dimension& dim) const
{
    IntVector width = IntVector::getZero(dim);
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        width = IntVector::max(width, (*it)->getRefineOpStencilWidth(dim));
    }
    return width;
}

void
RefinePatchStrategySet::preprocessRefine(Patch& fine, const Patch& coarse, const Box& fine_box, const IntVector& ratio)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->preprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
}

void
RefinePatchStrategySet::postprocessRefine(Patch& fine, const Patch& coarse, const Box& fine_box, const IntVector& ratio)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->postprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
}

void RefinePatchStrategySet::preprocessRefineBoxes(Patch& fine,
                                                   const Patch& coarse,
                                                   const BoxContainer& fine_boxes,
                                                   const IntVector& ratio)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->preprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
}

void RefinePatchStrategySet::postprocessRefineBoxes(Patch& fine,
                                                    const Patch& coarse,
                                                    const BoxContainer& fine_boxes,
                                                    const IntVector& ratio)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->postprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
