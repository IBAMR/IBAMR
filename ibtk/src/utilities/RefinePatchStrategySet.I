// Filename: RefinePatchStrategySet.I
// Created on 11 Sep 2006 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
RefinePatchStrategySet::~RefinePatchStrategySet()
{
    if (d_managed)
    {
        typedef std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> refine_strategy_set;
        for (refine_strategy_set::iterator it = d_strategy_set.begin();
             it != d_strategy_set.end(); ++it)
        {
            delete (*it);
        }
    }
    return;
}// ~RefinePatchStrategySet

inline void
RefinePatchStrategySet::setPhysicalBoundaryConditions(
    SAMRAI::hier::Patch<NDIM>& patch,
    const double fill_time,
    const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill)
{
    typedef std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->setPhysicalBoundaryConditions(
            patch, fill_time, ghost_width_to_fill);
    }
    return;
}// setPhysicalBoundaryConditions

inline SAMRAI::hier::IntVector<NDIM>
RefinePatchStrategySet::getRefineOpStencilWidth() const
{
    SAMRAI::hier::IntVector<NDIM> width = 0;
    typedef std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::const_iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        width = SAMRAI::hier::IntVector<NDIM>::max(
            width, (*it)->getRefineOpStencilWidth());
    }
    return width;
}// getRefineOpStencilWidth()

inline void
RefinePatchStrategySet::preprocessRefine(
    SAMRAI::hier::Patch<NDIM>& fine,
    const SAMRAI::hier::Patch<NDIM>& coarse,
    const SAMRAI::hier::Box<NDIM>& fine_box,
    const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    typedef std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->preprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
}// preprocessRefine

inline void
RefinePatchStrategySet::postprocessRefine(
    SAMRAI::hier::Patch<NDIM>& fine,
    const SAMRAI::hier::Patch<NDIM>& coarse,
    const SAMRAI::hier::Box<NDIM>& fine_box,
    const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    typedef std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->postprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
}// postprocessRefine

inline void
RefinePatchStrategySet::preprocessRefineBoxes(
    SAMRAI::hier::Patch<NDIM>& fine,
    const SAMRAI::hier::Patch<NDIM>& coarse,
    const SAMRAI::hier::BoxList<NDIM>& fine_boxes,
    const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    typedef std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->preprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
}// preprocessRefineBoxes

inline void
RefinePatchStrategySet::postprocessRefineBoxes(
    SAMRAI::hier::Patch<NDIM>& fine,
    const SAMRAI::hier::Patch<NDIM>& coarse,
    const SAMRAI::hier::BoxList<NDIM>& fine_boxes,
    const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    typedef std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> refine_strategy_set;
    for (refine_strategy_set::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->postprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
}// postprocessRefineBoxes

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
