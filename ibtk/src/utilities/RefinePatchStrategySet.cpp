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
#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
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
    if (d_managed)
    {
        typedef std::vector<RefinePatchStrategy<NDIM>*> refine_strategy_set;
        for (refine_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
        {
            delete (*it);
        }
    }
    return;
} // ~RefinePatchStrategySet

void
RefinePatchStrategySet::setPhysicalBoundaryConditions(Patch<NDIM>& patch,
                                                      const double fill_time,
                                                      const IntVector<NDIM>& ghost_width_to_fill)
{
    typedef std::vector<RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    }
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
RefinePatchStrategySet::getRefineOpStencilWidth() const
{
    IntVector<NDIM> width = 0;
    typedef std::vector<RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::const_iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        width = IntVector<NDIM>::max(width, (*it)->getRefineOpStencilWidth());
    }
    return width;
} // getRefineOpStencilWidth()

void
RefinePatchStrategySet::preprocessRefine(Patch<NDIM>& fine,
                                         const Patch<NDIM>& coarse,
                                         const Box<NDIM>& fine_box,
                                         const IntVector<NDIM>& ratio)
{
    typedef std::vector<RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->preprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
} // preprocessRefine

void
RefinePatchStrategySet::postprocessRefine(Patch<NDIM>& fine,
                                          const Patch<NDIM>& coarse,
                                          const Box<NDIM>& fine_box,
                                          const IntVector<NDIM>& ratio)
{
    typedef std::vector<RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->postprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
} // postprocessRefine

void
RefinePatchStrategySet::preprocessRefineBoxes(Patch<NDIM>& fine,
                                              const Patch<NDIM>& coarse,
                                              const BoxList<NDIM>& fine_boxes,
                                              const IntVector<NDIM>& ratio)
{
    typedef std::vector<RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->preprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
} // preprocessRefineBoxes

void
RefinePatchStrategySet::postprocessRefineBoxes(Patch<NDIM>& fine,
                                               const Patch<NDIM>& coarse,
                                               const BoxList<NDIM>& fine_boxes,
                                               const IntVector<NDIM>& ratio)
{
    typedef std::vector<RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->postprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
} // postprocessRefineBoxes

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
