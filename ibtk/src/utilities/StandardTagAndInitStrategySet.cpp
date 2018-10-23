// Filename: StandardTagAndInitStrategySet.cpp
// Created on 26 Jun 2007 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include <algorithm>
#include <limits>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "StandardTagAndInitStrategy.h"
#include "ibtk/StandardTagAndInitStrategySet.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StandardTagAndInitStrategySet::~StandardTagAndInitStrategySet()
{
    if (d_managed)
    {
        for (const auto& strategy : d_strategy_set)
        {
            delete strategy;
        }
    }
    return;
} // ~StandardTagAndInitStrategySet

double
StandardTagAndInitStrategySet::getLevelDt(const Pointer<BasePatchLevel<NDIM> > level,
                                          const double dt_time,
                                          const bool initial_time)
{
    double dt = std::numeric_limits<double>::max();
    for (const auto& strategy : d_strategy_set)
    {
        dt = std::min(dt, strategy->getLevelDt(level, dt_time, initial_time));
    }
    return dt;
} // getLevelDt

double
StandardTagAndInitStrategySet::advanceLevel(const Pointer<BasePatchLevel<NDIM> > level,
                                            const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                            const double current_time,
                                            const double new_time,
                                            const bool first_step,
                                            const bool last_step,
                                            const bool regrid_advance)
{
    double dt = std::numeric_limits<double>::max();
    for (const auto& strategy : d_strategy_set)
    {
        dt = std::min(
            dt,
            strategy->advanceLevel(level, hierarchy, current_time, new_time, first_step, last_step, regrid_advance));
    }
    return dt;
} // advanceLevel

void
StandardTagAndInitStrategySet::resetTimeDependentData(const Pointer<BasePatchLevel<NDIM> > level,
                                                      const double new_time,
                                                      const bool can_be_refined)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->resetTimeDependentData(level, new_time, can_be_refined);
    }
    return;
} // resetTimeDependentData

void
StandardTagAndInitStrategySet::resetDataToPreadvanceState(const Pointer<BasePatchLevel<NDIM> > level)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->resetDataToPreadvanceState(level);
    }
    return;
} // resetDataToPreadvanceState

void
StandardTagAndInitStrategySet::initializeLevelData(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                   const int level_number,
                                                   const double init_data_time,
                                                   const bool can_be_refined,
                                                   const bool initial_time,
                                                   const Pointer<BasePatchLevel<NDIM> > old_level,
                                                   const bool allocate_data)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->initializeLevelData(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    }
    return;
} // initializeLevelData

void
StandardTagAndInitStrategySet::resetHierarchyConfiguration(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                           const int coarsest_level,
                                                           const int finest_level)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    }
    return;
} // resetHierarchyConfiguration

void
StandardTagAndInitStrategySet::applyGradientDetector(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                     const int level_number,
                                                     const double error_data_time,
                                                     const int tag_index,
                                                     const bool initial_time,
                                                     const bool uses_richardson_extrapolation_too)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }
    return;
} // applyGradientDetector

void
StandardTagAndInitStrategySet::applyRichardsonExtrapolation(const Pointer<PatchLevel<NDIM> > level,
                                                            const double error_data_time,
                                                            const int tag_index,
                                                            const double deltat,
                                                            const int error_coarsen_ratio,
                                                            const bool initial_time,
                                                            const bool uses_gradient_detector_too)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->applyRichardsonExtrapolation(
            level, error_data_time, tag_index, deltat, error_coarsen_ratio, initial_time, uses_gradient_detector_too);
    }
    return;
} // applyRichardsonExtrapolation

void
StandardTagAndInitStrategySet::coarsenDataForRichardsonExtrapolation(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                     const int level_number,
                                                                     const Pointer<PatchLevel<NDIM> > coarser_level,
                                                                     const double coarsen_data_time,
                                                                     const bool before_advance)
{
    for (const auto& strategy : d_strategy_set)
    {
        strategy->coarsenDataForRichardsonExtrapolation(
            hierarchy, level_number, coarser_level, coarsen_data_time, before_advance);
    }
    return;
} // coarsenDataForRichardsonExtrapolation

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
