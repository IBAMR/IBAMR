// Filename: IBStrategy.C
// Created on 21 Sep 2011 by Boyce Griffith
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

#include "IBStrategy.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStrategy::IBStrategy()
{
    // intentionally blank
    return;
}// IBStrategy

IBStrategy::~IBStrategy()
{
    // intentionally blank
    return;
}// ~IBStrategy

void
IBStrategy::registerIBHierarchyIntegrator(
    IBHierarchyIntegrator* ib_solver)
{
    d_ib_solver = ib_solver;
    return;
}// registerIBHierarchyIntegrator

void
IBStrategy::preprocessIntegrateData(
    double /*current_time*/,
    double /*new_time*/,
    int /*num_cycles*/)
{
    // intentionally blank
    return;
}// preprocessIntegrateData

void
IBStrategy::postprocessIntegrateData(
    double /*current_time*/,
    double /*new_time*/,
    int /*num_cycles*/)
{
    // intentionally blank
    return;
}// postprocessIntegrateData

bool
IBStrategy::hasFluidSources() const
{
    return false;
}// hasFluidSources

void
IBStrategy::computeLagrangianFluidSource(
    double /*data_time*/)
{
    // intentionally blank
    return;
}// computeLagrangianFluidSource

void
IBStrategy::spreadFluidSource(
    int /*q_data_idx*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*q_prolongation_scheds*/,
    double /*data_time*/)
{
    // intentionally blank
    return;
}// spreadFluidSource

void
IBStrategy::interpolatePressure(
    int /*p_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*p_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*p_ghost_fill_scheds*/,
    double /*data_time*/)
{
    // intentionally blank
    return;
}// interpolatePressure

void
IBStrategy::postprocessData(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/)
{
    // intentionally blank
    return;
}// postprocessData

void
IBStrategy::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/,
    int /*u_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    int /*integrator_step*/,
    double /*init_data_time*/,
    bool /*initial_time*/)
{
    // intentionally blank
    return;
}// initializePatchHierarchy

void
IBStrategy::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > /*load_balancer*/,
    int /*workload_data_idx*/)
{
    // intentionally blank
    return;
}// registerLoadBalancer

void
IBStrategy::updateWorkloadEstimates(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    int /*workload_data_idx*/)
{
    // intentionally blank
    return;
}// updateWorkloadEstimates

void
IBStrategy::beginDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
}// beginDataRedistribution

void
IBStrategy::endDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
}// endDataRedistribution

void
IBStrategy::initializeLevelData(
    Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
    int /*level_number*/,
    double /*init_data_time*/,
    bool /*can_be_refined*/,
    bool /*initial_time*/,
    Pointer<BasePatchLevel<NDIM> > /*old_level*/,
    bool /*allocate_data*/)
{
    // intentionally blank
    return;
}// initializeLevelData

void
IBStrategy::resetHierarchyConfiguration(
    Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
    int /*coarsest_level*/,
    int /*finest_level*/)
{
    // intentionally blank
    return;
}// resetHierarchyConfiguration

void
IBStrategy::applyGradientDetector(
    Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
    int /*level_number*/,
    double /*error_data_time*/,
    int /*tag_index*/,
    bool /*initial_time*/,
    bool /*uses_richardson_extrapolation_too*/)
{
    // intentionally blank
    return;
}// applyGradientDetector

void
IBStrategy::putToDatabase(
    Pointer<Database> /*db*/)
{
    // intentionally blank
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
