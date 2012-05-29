// Filename: StandardTagAndInitStrategySet.I
// Created on 26 Jun 2007 by Boyce Griffith
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

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
StandardTagAndInitStrategySet::~StandardTagAndInitStrategySet()
{
    if (d_managed)
    {
        typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> coarsen_strategy_set;
        for (coarsen_strategy_set::const_iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
        {
            delete (*it);
        }
    }
    return;
}// ~StandardTagAndInitStrategySet

inline double
StandardTagAndInitStrategySet::getLevelDt(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > level,
    const double dt_time,
    const bool initial_time)
{
    double dt = std::numeric_limits<double>::max();
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        dt = std::min(dt,(*it)->getLevelDt(level,dt_time,initial_time));
    }
    return dt;
}// getLevelDt

inline double
StandardTagAndInitStrategySet::advanceLevel(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > level,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const double current_time,
    const double new_time,
    const bool first_step,
    const bool last_step,
    const bool regrid_advance)
{
    double dt = std::numeric_limits<double>::max();
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        dt = std::min(dt,(*it)->advanceLevel(level,hierarchy,current_time,new_time,first_step,last_step,regrid_advance));
    }
    return dt;
}// advanceLevel

inline void
StandardTagAndInitStrategySet::resetTimeDependentData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > level,
    const double new_time,
    const bool can_be_refined)
{
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->resetTimeDependentData(level,new_time,can_be_refined);
    }
    return;
}// resetTimeDependentData

inline void
StandardTagAndInitStrategySet::resetDataToPreadvanceState(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > level)
{
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->resetDataToPreadvanceState(level);
    }
    return;
}// resetDataToPreadvanceState

inline void
StandardTagAndInitStrategySet::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    }
    return;
}// initializeLevelData

inline void
StandardTagAndInitStrategySet::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    }
    return;
}// resetHierarchyConfiguration

inline void
StandardTagAndInitStrategySet::applyGradientDetector(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }
    return;
}// applyGradientDetector

inline void
StandardTagAndInitStrategySet::applyRichardsonExtrapolation(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
    const double error_data_time,
    const int tag_index,
    const double deltat,
    const int error_coarsen_ratio,
    const bool initial_time,
    const bool uses_gradient_detector_too)
{
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->applyRichardsonExtrapolation(level, error_data_time, tag_index, deltat, error_coarsen_ratio, initial_time, uses_gradient_detector_too);
    }
    return;
}// applyRichardsonExtrapolation

inline void
StandardTagAndInitStrategySet::coarsenDataForRichardsonExtrapolation(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level,
    const double coarsen_data_time,
    const bool before_advance)
{
    typedef std::vector<SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>*> tag_and_init_strategy_set;
    for (tag_and_init_strategy_set::iterator it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->coarsenDataForRichardsonExtrapolation(hierarchy, level_number, coarser_level, coarsen_data_time, before_advance);
    }
    return;
}// coarsenDataForRichardsonExtrapolation

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
