// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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
        for (const auto& strategy : d_strategy_set)
        {
            delete strategy;
        }
    }
    return;
} // ~RefinePatchStrategySet

void
RefinePatchStrategySet::setPhysicalBoundaryConditions(PatchNd& patch,
                                                      const double fill_time,
                                                      const IntVectorNd& ghost_width_to_fill)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    }
    return;
} // setPhysicalBoundaryConditions

IntVectorNd
RefinePatchStrategySet::getRefineOpStencilWidth() const
{
    IntVectorNd width = 0;
    for (const auto& strategy : d_strategy_set)
    {
        width = IntVectorNd::max(width, strategy->getRefineOpStencilWidth());
    }
    return width;
} // getRefineOpStencilWidth()

void
RefinePatchStrategySet::preprocessRefine(PatchNd& fine,
                                         const PatchNd& coarse,
                                         const BoxNd& fine_box,
                                         const IntVectorNd& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->preprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
} // preprocessRefine

void
RefinePatchStrategySet::postprocessRefine(PatchNd& fine,
                                          const PatchNd& coarse,
                                          const BoxNd& fine_box,
                                          const IntVectorNd& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->postprocessRefine(fine, coarse, fine_box, ratio);
    }
    return;
} // postprocessRefine

void
RefinePatchStrategySet::preprocessRefineBoxes(PatchNd& fine,
                                              const PatchNd& coarse,
                                              const BoxListNd& fine_boxes,
                                              const IntVectorNd& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->preprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
} // preprocessRefineBoxes

void
RefinePatchStrategySet::postprocessRefineBoxes(PatchNd& fine,
                                               const PatchNd& coarse,
                                               const BoxListNd& fine_boxes,
                                               const IntVectorNd& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->postprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }
    return;
} // postprocessRefineBoxes

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
