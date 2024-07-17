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

#include "ibtk/CoarsenPatchStrategySet.h"

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

CoarsenPatchStrategySet::~CoarsenPatchStrategySet()
{
    if (d_managed)
    {
        for (const auto& strategy : d_strategy_set)
        {
            delete strategy;
        }
    }
    return;
} // ~CoarsenPatchStrategySet

IntVectorNd
CoarsenPatchStrategySet::getCoarsenOpStencilWidth() const
{
    IntVectorNd width = 0;
    for (const auto& strategy : d_strategy_set)
    {
        width = IntVectorNd::max(width, strategy->getCoarsenOpStencilWidth());
    }
    return width;
} // getCoarsenOpStencilWidth

void
CoarsenPatchStrategySet::preprocessCoarsen(PatchNd& coarse,
                                           const PatchNd& fine,
                                           const BoxNd& coarse_box,
                                           const IntVectorNd& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->preprocessCoarsen(coarse, fine, coarse_box, ratio);
    }
    return;
} // preprocessCoarsen

void
CoarsenPatchStrategySet::postprocessCoarsen(PatchNd& coarse,
                                            const PatchNd& fine,
                                            const BoxNd& coarse_box,
                                            const IntVectorNd& ratio)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->postprocessCoarsen(coarse, fine, coarse_box, ratio);
    }
    return;
} // postprocessCoarsen

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
