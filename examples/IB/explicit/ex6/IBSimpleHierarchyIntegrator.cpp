// Filename: IBSimpleHierarchyIntegrator.cpp
// Created on 04 Apr 2013 by Boyce Griffith
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

#include "IBSimpleHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>

// IBTK INCLUDES
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSimpleHierarchyIntegrator::IBSimpleHierarchyIntegrator(const std::string& object_name,
                                                         Pointer<Database> input_db,
                                                         Pointer<IBMethod> ib_method_ops,
                                                         Pointer<INSHierarchyIntegrator> ins_hier_integrator)
    : IBHierarchyIntegrator(object_name, input_db, ib_method_ops, ins_hier_integrator, /*register_for_restart*/ false)
{
    // intentionally blank
    return;
} // IBSimpleHierarchyIntegrator

IBSimpleHierarchyIntegrator::~IBSimpleHierarchyIntegrator()
{
    // intentionally blank
    return;
} // ~IBSimpleHierarchyIntegrator

void
IBSimpleHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                          const double new_time,
                                                          const int num_cycles)
{
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_level_num = 0;
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();

    // Allocate Eulerian scratch and new data.
    for (int level_num = coarsest_level_num; level_num <= finest_level_num; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Initialize the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, ins_num_cycles);

    // Initialize IB data.
    //
    // NOTE: We assume here that all IB data are assigned to the finest level of
    // the AMR patch hierarchy.
    d_X_current_data = l_data_manager->getLData(LDataManager::POSN_DATA_NAME, finest_level_num);
    d_X_new_data = l_data_manager->createLData("X_new", finest_level_num, NDIM);
    d_U_data = l_data_manager->getLData(LDataManager::VEL_DATA_NAME, finest_level_num);
    d_F_data = l_data_manager->createLData("F", finest_level_num, NDIM);
    return;
} // preprocessIntegrateHierarchy

void
IBSimpleHierarchyIntegrator::integrateHierarchy(const double current_time, const double new_time, const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    PetscErrorCode ierr;
    const double dt = new_time - current_time;
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();

    // Here we implement a simple time integration scheme:
    //
    // (1) Compute Lagrangian forces and spread those forces to the grid: f =
    //     S[X^{n}] F^{n}.
    //
    // (2) Solve the fluid equations using the fluid solver registered with this
    //     class using the forcing computed in (1).
    //
    // (3) Interpolate the Eulerian velocity and update the positions of the
    //     Lagrangian points: X^{n+1} = X^{n} + dt*S^{*}[X^{n}] u^{n+1}.
    //
    // NOTE: We assume here that all IB data are assigned to the finest level of
    // the AMR patch hierarchy.

    ///////////////////////////////////////////////////////////////////////////
    // (1) Compute Lagrangian forces and spread those forces to the grid: f =
    //     S[X^{n}] F^{n}..

    // For simplicity, set the Lagrangian force density to equal zero.
    ierr = VecZeroEntries(d_F_data->getVec());
    IBTK_CHKERRQ(ierr);

    // Spread the forces to the grid.  We use the "current" Lagrangian position
    // data to define the locations from where the forces are spread.
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
    l_data_manager->spread(d_f_idx,
                           d_F_data,
                           d_X_current_data,
                           d_u_phys_bdry_op,
                           finest_level_num,
                           getProlongRefineSchedules(d_object_name + "::f"),
                           /*F_needs_ghost_fill*/ true,
                           /*X_needs_ghost_fill*/ true);

    // NOTE: Any additional Eulerian forcing should be computed here and added
    // to the data associated with d_f_idx.

    ///////////////////////////////////////////////////////////////////////////
    // (2) Solve the fluid equations using the fluid solver registered with this
    // class using the forcing computed in (1).
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    for (int ins_cycle_num = 0; ins_cycle_num < ins_num_cycles; ++ins_cycle_num)
    {
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time, ins_cycle_num);
    }

    ///////////////////////////////////////////////////////////////////////////
    // (3) Interpolate the Eulerian velocity and update the positions of the
    // Lagrangian points: X^{n+1} = X^{n} + dt*S^{*}[X^{n}] u^{n+1}.
    //
    // NOTE: We use the "new" velocity data (defined at time level n+1) to
    // determine the velocity of the Lagrangian nodes.  We could also use the
    // "current" data (defined at time level n) or some other velocity field
    // here.  We use the "current" Lagrangian position data to define the
    // locations to where the velocities are interpolated.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());
    d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    l_data_manager->interp(d_u_idx,
                           d_U_data,
                           d_X_current_data,
                           finest_level_num,
                           getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                           getGhostfillRefineSchedules(d_object_name + "::u"),
                           current_time);
    ierr = VecWAXPY(d_X_new_data->getVec(), dt, d_U_data->getVec(), d_X_current_data->getVec());
    IBTK_CHKERRQ(ierr);
    return;
} // integrateHierarchy

void
IBSimpleHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                           const double new_time,
                                                           const bool skip_synchronize_new_state_data,
                                                           const int num_cycles)
{
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    const int coarsest_level_num = 0;
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    PetscErrorCode ierr;

    // Deallocate Eulerian scratch data.
    for (int level_num = coarsest_level_num; level_num <= finest_level_num; ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_f_idx);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

    // Deallocate the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, ins_num_cycles);

    // Reset and deallocate IB data.
    //
    // NOTE: We assume here that all IB data are assigned to the finest level of
    // the AMR patch hierarchy.
    ierr = VecSwap(d_X_current_data->getVec(), d_X_new_data->getVec());
    IBTK_CHKERRQ(ierr);
    d_X_current_data = NULL;
    d_X_new_data = NULL;
    d_U_data = NULL;
    d_F_data = NULL;
    return;
} // postprocessIntegrateHierarchy

void
IBSimpleHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                           Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Setup the fluid solver for explicit coupling.
    //
    // NOTE: This will use the data associated with d_f_idx to provide forcing
    // for the fluid equations.
    d_ins_hier_integrator->registerBodyForceFunction(new IBEulerianForceFunction(this));

    // NOTE: Any additional implementation-specific initialization should be
    // performed here.

    // Finish initializing the hierarchy integrator.  This function call should
    // come at the end of this function.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
    return;
} // initializeHierarchyIntegrator

//////////////////////////////////////////////////////////////////////////////
