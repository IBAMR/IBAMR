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

#include "ibtk/StandardTagAndInitStrategySet.h"

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "StandardTagAndInitStrategy.h"
#include "tbox/Pointer.h"

#include <algorithm>
#include <limits>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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
