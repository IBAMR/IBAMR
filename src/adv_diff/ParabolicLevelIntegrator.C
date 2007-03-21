// Filename: ParabolicLevelIntegrator.C
// Last modified: <21.Mar.2007 01:05:29 griffith@box221.cims.nyu.edu>
// Created on 09 Jan 2007 by Boyce Griffith (boyce@box221.cims.nyu.edu)

// NOTE: This implementation is directly derived from the SAMRAI
// HyperbolicLevelIntegrator source code.

#include "ParabolicLevelIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

/*
*************************************************************************
*                                                                       *
* This constructor sets the ParabolicPatchStrategy pointer and          *
* initializes integration parameters to default values.  Communication  *
* algorithms are created here too.  Other data members are read in      *
* from the input database or from the restart database corresponding    *
* to the specified object_name.                                         *
*                                                                       *
*************************************************************************
*/
ParabolicLevelIntegrator::ParabolicLevelIntegrator(
    const string& object_name,
    tbox::Pointer<tbox::Database> input_db,
    ParabolicPatchStrategy<NDIM>* patch_strategy,
    bool register_for_restart,
    bool using_time_refinement)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
    assert(patch_strategy != ((ParabolicPatchStrategy<NDIM>*)NULL));
#endif

    d_object_name = object_name;
    d_using_time_refinement = using_time_refinement;
    d_registered_for_restart = register_for_restart;

    if (d_registered_for_restart)
    {
        tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    d_patch_strategy = patch_strategy;

    /*
     * Default parameter values.
     */
    d_number_time_data_levels = 2;

    d_flux_is_face = true;
    d_flux_face_registered = false;
    d_flux_side_registered = false;

    d_have_flux_on_level_zero = false;

    d_lag_dt_computation = true;
    d_using_ghosts_for_dt = false;
    d_distinguish_mpi_reduction_costs = false;

    d_cfl = tbox::IEEE::getSignalingNaN();
    d_cfl_init = tbox::IEEE::getSignalingNaN();

    /*
     * Communication algorithms.
     */
    d_bdry_fill_advance     = new xfer::RefineAlgorithm<NDIM>();
    d_bdry_fill_advance_new = new xfer::RefineAlgorithm<NDIM>();
    d_bdry_fill_advance_old = new xfer::RefineAlgorithm<NDIM>();
    d_fill_new_level        = new xfer::RefineAlgorithm<NDIM>();
    d_coarsen_fluxsum       = new xfer::CoarsenAlgorithm<NDIM>();
    d_coarsen_sync_data     = new xfer::CoarsenAlgorithm<NDIM>();
    d_sync_initial_data     = new xfer::CoarsenAlgorithm<NDIM>();

    d_coarsen_rich_extrap_init  = new xfer::CoarsenAlgorithm<NDIM>();
    d_coarsen_rich_extrap_final = new xfer::CoarsenAlgorithm<NDIM>();

    /*
     * hier::Variable contexts used in algorithm.  Note that "OLD" context
     * is only created and used in the case of Richardson extrapolation
     * and a refinement ratio of 3 (see registerModelVariables()).
     */
    d_scratch = hier::VariableDatabase<NDIM>::getDatabase()->getContext("SCRATCH");
    d_current = hier::VariableDatabase<NDIM>::getDatabase()->getContext("CURRENT");
    d_new     = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");

    d_plot_context = d_current;

    /*
     * Timers:  one for each of the communication algorithms ("create"
     * indicates schedule creation, "comm" indicates communication)
     */
    t_advance_bdry_fill_comm = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::advance_bdry_fill_comm");
    t_error_bdry_fill_create = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::error_bdry_fill_create");
    t_error_bdry_fill_comm = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::error_bdry_fill_comm");
    t_mpi_reductions = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::mpi_reductions");
    t_initialize_level_data = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::initializeLevelData()");
    t_fill_new_level_create = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::fill_new_level_create");
    t_fill_new_level_comm = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::fill_new_level_comm");
    t_advance_bdry_fill_create =tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::advance_bdry_fill_create");
    t_new_advance_bdry_fill_create = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::new_advance_bdry_fill_create");
    t_apply_gradient_detector = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::applyGradientDetector()");
    t_coarsen_rich_extrap = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::coarsen_rich_extrap");
    t_get_level_dt = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::getLevelDt()");
    t_get_level_dt_sync = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::getLevelDt()_sync");
    t_advance_level = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::advanceLevel()");
    t_new_advance_bdry_fill_comm = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::new_advance_bdry_fill_comm");
    t_patch_num_kernel = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::patch_numerical_kernels");
    t_advance_level_sync = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::advanceLevel()_sync");
    t_std_level_sync = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::standardLevelSynchronization()");
    t_sync_new_levels = tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::synchronizeNewLevels()");

    /*
     * Initialize object with data read from the input and restart databases.
     */

    bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);
    return;
}// ParabolicLevelIntegrator

/*
*************************************************************************
*                                                                       *
* Destructor tells the tbox::RestartManager to remove this object from   *
* the list of restart items.                                            *
*                                                                       *
*************************************************************************
*/
ParabolicLevelIntegrator::~ParabolicLevelIntegrator()
{
    if (d_registered_for_restart)
    {
        tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    return;
}// ~ParabolicLevelIntegrator

/*
*************************************************************************
*                                                                       *
* Initialize level integrator by:                                       *
*                                                                       *
*   (1) Setting the number of time data levels based on needs of        *
*       the gridding algorithm                                          *
*   (2) Invoking variable registration in patch strategy.               *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::initializeLevelIntegrator(
    tbox::Pointer<mesh::BaseGriddingAlgorithm<NDIM> > base_gridding_alg)
{
    tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_alg = base_gridding_alg;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!gridding_alg.isNull());
#endif

    d_number_time_data_levels = 2;

    if ((gridding_alg->getErrorCoarsenRatio() < 1) ||
        (gridding_alg->getErrorCoarsenRatio() > 3))
    {
        TBOX_ERROR("ParabolicLevelIntegrator::initializeLevelIntegrator "
                   << "error...\n" << "   object name = " << d_object_name
                   << "   gridding algorithm has bad error coarsen ratio" << endl);
    }

    if ((gridding_alg->errorEstimationUsesTimeIntegration()) &&
        (gridding_alg->getErrorCoarsenRatio() == 3))
    {
        d_number_time_data_levels = 3;
        d_old = hier::VariableDatabase<NDIM>::getDatabase()->getContext("OLD");
    }

    d_patch_strategy->registerModelVariables(this);

    d_patch_strategy->setupLoadBalancer(this, gridding_alg.getPointer());

    return;
}// initializeLevelIntegrator

/*
*************************************************************************
*                                                                       *
* Invoke dt calculation routines in patch strategy and take a min       *
* over all patches on the level.  The result will be the max of the     *
* next timestep on the level. If the boolean recompute_dt is true,      *
* the max timestep on the level will be computed.  If it is false,      *
* the method will simply access the latest dt stored in the time        *
* refinement integrator.                                                *
*                                                                       *
*************************************************************************
*/
double
ParabolicLevelIntegrator::getLevelDt(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double dt_time,
    const bool initial_time)
{
    tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
    assert(!patch_level.isNull());
#endif
    t_get_level_dt->start();

    double dt = tbox::IEEE::getDBL_MAX();

    if (!d_using_ghosts_for_dt)
    {

        d_patch_strategy->setDataContext(d_current);
        for (hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            tbox::Pointer<hier::Patch<NDIM> > patch = patch_level->getPatch(p());

            patch->allocatePatchData(d_temp_var_scratch_data, dt_time);

            double patch_dt;
            patch_dt = d_patch_strategy->
                computeStableDtOnPatch(*patch,
                                       initial_time,
                                       dt_time);

            dt = tbox::Utilities::dmin(dt, patch_dt);

            patch->deallocatePatchData(d_temp_var_scratch_data);
        }

        d_patch_strategy->clearDataContext();

    }
    else
    {
        patch_level->allocatePatchData(d_saved_var_scratch_data, dt_time);

        d_patch_strategy->setDataContext(d_scratch);

        t_advance_bdry_fill_comm->start();
        d_bdry_sched_advance[patch_level->getLevelNumber()]->fillData(dt_time);
        t_advance_bdry_fill_comm->stop();

        for (hier::PatchLevel<NDIM>::Iterator ip(patch_level); ip; ip++)
        {
            tbox::Pointer<hier::Patch<NDIM> > patch = patch_level->getPatch(ip());

            patch->allocatePatchData(d_temp_var_scratch_data, dt_time);

            double patch_dt;
            patch_dt = d_patch_strategy->
                computeStableDtOnPatch(*patch,
                                       initial_time,
                                       dt_time);

            dt = tbox::Utilities::dmin(dt, patch_dt);

            patch->deallocatePatchData(d_temp_var_scratch_data);
        }

        d_patch_strategy->clearDataContext();

        /*
         * Copy data from scratch to current and de-allocate scratch storage.
         * This may be excessive here, but seems necessary if the
         * computation of dt affects the state of the problem solution.
         * Also, this getLevelDt() routine is called at initialization only
         * in most cases.
         */
        copyTimeDependentData(patch_level, d_scratch, d_current);

        patch_level->deallocatePatchData(d_saved_var_scratch_data);
    }

    t_get_level_dt_sync->start();

    if (d_distinguish_mpi_reduction_costs)
    {
        tbox::MPI::barrier();
        t_get_level_dt_sync->stop();
        t_mpi_reductions->start();
    }

    /*
     * The level time increment is a global min over all patches.
     */
    double global_dt = tbox::MPI::minReduction(dt)
        * tbox::Utilities::dmin(d_cfl_init, d_cfl);

    if (d_distinguish_mpi_reduction_costs)
    {
        t_mpi_reductions->stop();
    }
    else
    {
        t_get_level_dt_sync->stop();
    }

    t_get_level_dt->stop();

    return global_dt;
}// getLevelDt

/*
*************************************************************************
*                                                                       *
* For the standard explicit integration algorithm for hyperbolic        *
* conservation laws, the fine time increment is the coarse increment    *
* divided by the maximum mesh ratio (independent of level number).      *
*                                                                       *
*************************************************************************
*/
double
ParabolicLevelIntegrator::getMaxFinerLevelDt(
    const int finer_level_number,
    const double coarse_dt,
    const hier::IntVector<NDIM>& ratio)
{
    NULL_USE(finer_level_number);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int id = 0; id < DIM; id++)
    {
        assert(ratio(id) > 0);
    }
#endif
    return coarse_dt/static_cast<double>(ratio.max());
}// getMaxFinerLevelDt

/*
*************************************************************************
*                                                                       *
* Integrate data on all patches in patch level from current time        *
* to new time (new_time) using a single time step.  Before the advance  *
* can occur, proper ghost cell information is obtained for all patches  *
* on the level.  Then, local patches are advanced sequentially in the   *
* loop over patches.  The details of the routine are as follows:        *
*                                                                       *
*  0) Allocate storage for new time level data. Also, allocate          *
*     necessary FLUX and flux integral storage if needed                *
*     (i.e., if regrid_advance is false, first_step is true, and        *
*     coarser or finer level than current level exists in hierarchy.)   *
*                                                                       *
*  1) Scratch space is filled so that, for each patch, interior data    *
*     and ghost cell bdry data correspond to specified time.            *
*                                                                       *
*  1a) Call user routines to pre-process advance data, if needed.       *
*                                                                       *
*  2) Compute explicit fluxes in scratch space using data on            *
*     patch + ghosts at given time.                                     *
*                                                                       *
*  3) Apply conservative difference in scratch space to advance patch   *
*     interior data to time = new_time.                                 *
*                                                                       *
*  3a) Call user routines to post-process advance data, if needed.      *
*                                                                       *
*  4) Compute next stable time increment for subsequent level advances: *
*                                                                       *
*     4a) If (d_lag_dt_computation == true)                             *
*         {                                                             *
*            DO NOT RECOMPUTE characteristic data after advancing       *
*            data on patch. Use characteristic data corresponding       *
*            to current time level, computed prior to flux computation, *
*            in dt calculation.                                         *
*            If (d_using_ghosts_for_dt == true)                         *
*               - Compute dt using data on patch+ghosts at time.        *
*            Else                                                       *
*               - Compute dt using data on patch interior ONLY.         *
*         }                                                             *
*                                                                       *
*     4b) Copy data from scratch space patch interior to new data       *
*         storage for patch (i.e., at time = new_time).                 *
*                                                                       *
*     4a) If (d_lag_dt_computation == false)                            *
*         {                                                             *
*            RECOMPUTE characteristic data after advancing data on      *
*            patch. Use characteristic data corresponding to new time   *
*            level in dt calculation.                                   *
*            If (d_using_ghosts_for_dt == true)                         *
*               - Refill scratch space with new interior patch data     *
*                 and ghost cell bdry data correspond to new time.      *
*                 (NOTE: This requires a new boundary schedule.)        *
*               - Compute dt using data on patch+ghosts at new_time.    *
*            Else                                                       *
*               - Compute dt using data on patch interior ONLY.         *
*                 (using patch interior data at new_time)               *
*         }                                                             *
*                                                                       *
*  5) If (ln > 0), update flux integrals by adding patch bdry FLUXes    *
*     to flux sums.                                                     *
*                                                                       *
* Important Notes:                                                      *
*    1) In order to advance finer levels (if they exist), both old      *
*       and new data for each patch on the level must be maintained.    *
*    2) If the timestep is the first in the timestep loop on the level  *
*       (indicated by first_step), then time interpolation is           *
*       is unnecessary to fill ghost cells from the next coarser level. *
*    3) The new dt is not calculated if regrid_advance is true.         *
*       If this is the case, it is assumed that the results of the      *
*       advance and the timestep calculation will be discarded          *
*       (e.g., during regridding, or initialization).  Also, allocation *
*       and post-processing of FLUX/flux integral data is not performed *
*       in this case.                                                   *
*                                                                       *
*************************************************************************
*/
double
ParabolicLevelIntegrator::advanceLevel(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const double current_time,
    const double new_time,
    const bool first_step,
    const bool last_step,
    const bool regrid_advance)
{
    tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!patch_level.isNull());
    assert(!hierarchy.isNull());
    assert(current_time <= new_time);
#endif

    t_advance_level->start();

    const int level_number = patch_level->getLevelNumber();
    const double dt = new_time - current_time;

    /*
     * (1) Allocate data needed for advancing level.
     * (2) Generate temporary communication schedule to fill ghost
     *     cells, if needed.
     * (3) Fill ghost cell data.
     * (4) Process flux storage before the advance.
     */

    patch_level->allocatePatchData(d_new_time_dep_data, new_time);
    patch_level->allocatePatchData(d_saved_var_scratch_data, current_time);

    tbox::Pointer<xfer::RefineSchedule<NDIM> > fill_schedule;
    if (!patch_level->inHierarchy())
    {
        t_error_bdry_fill_create->start();
        if (d_number_time_data_levels == 3)
        {
            fill_schedule = d_bdry_fill_advance_old->
                createSchedule(patch_level,
                               patch_level->
                               getNextCoarserHierarchyLevelNumber(),
                               hierarchy,
                               d_patch_strategy);
        }
        else
        {
            fill_schedule = d_bdry_fill_advance->
                createSchedule(patch_level,
                               patch_level->
                               getNextCoarserHierarchyLevelNumber(),
                               hierarchy,
                               d_patch_strategy);
        }
        t_error_bdry_fill_create->stop();
    }
    else
    {
        fill_schedule = d_bdry_sched_advance[level_number];
    }

    d_patch_strategy->setDataContext(d_scratch);
    if (regrid_advance)
    {
        t_error_bdry_fill_comm->start();
    }
    else
    {
        t_advance_bdry_fill_comm->start();
    }
    fill_schedule->fillData(current_time);
    if (regrid_advance)
    {
        t_error_bdry_fill_comm->stop();
    }
    else
    {
        t_advance_bdry_fill_comm->stop();
    }

    d_patch_strategy->clearDataContext();
    fill_schedule.setNull();

    preprocessFluxData(patch_level,
                       current_time,
                       new_time,
                       regrid_advance,
                       first_step,
                       last_step);

    /*
     * (5) Call user-routine to pre-process state data, if needed.
     * (6) Advance solution on all level patches (scratch storage).
     * (7) Copy new solution to from scratch to new storage.
     * (8) Call user-routine to post-process state data, if needed.
     */
    t_patch_num_kernel->start();
    d_patch_strategy->preprocessAdvanceLevelState(patch_level,
                                                  current_time,
                                                  dt,
                                                  first_step,
                                                  last_step,
                                                  regrid_advance);
    t_patch_num_kernel->stop();

    d_patch_strategy->setDataContext(d_scratch);
    for (hier::PatchLevel<NDIM>::Iterator ip(patch_level); ip; ip++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = patch_level->getPatch(ip());

        patch->allocatePatchData(d_temp_var_scratch_data, current_time);

        t_patch_num_kernel->start();
        d_patch_strategy->computeFluxesOnPatch(*patch,
                                               current_time,
                                               dt);
        t_patch_num_kernel->stop();

        bool at_syncronization = false;

        t_patch_num_kernel->start();
        d_patch_strategy->conservativeDifferenceOnPatch(*patch,
                                                        current_time,
                                                        dt,
                                                        at_syncronization);
        t_patch_num_kernel->stop();


        patch->deallocatePatchData(d_temp_var_scratch_data);
    }
    d_patch_strategy->clearDataContext();

    patch_level->setTime(new_time, d_saved_var_scratch_data);
    patch_level->setTime(new_time, d_flux_var_data);

    copyTimeDependentData(patch_level, d_scratch, d_new);

    t_patch_num_kernel->start();
    d_patch_strategy->postprocessAdvanceLevelState(patch_level,
                                                   current_time,
                                                   dt,
                                                   first_step,
                                                   last_step,
                                                   regrid_advance);
    t_patch_num_kernel->stop();

    /*
     * (9) If the level advance is for regridding, we compute the next timestep:
     *
     * (a) If the dt computation is lagged (i.e., we use pre-advance data
     *     to compute timestep), we reset scratch space on patch interiors
     *     if needed.  Then, we set the strategy context to current or scratch
     *     depending on whether ghost values are used to compute dt.
     * (b) If the dt computation is not lagged (i.e., we use advanced data
     *     to compute timestep), we refill scratch space, including ghost
     *     data with new solution values if needed.  Then, we set the strategy
     *     context to new or scratch depending on whether ghost values are
     *     used to compute dt.
     * (c) Then, we loop over patches and compute the dt on each patch.
     */

    double dt_next = tbox::IEEE::getDBL_MAX();

    if (!regrid_advance)
    {

        if (d_lag_dt_computation)
        {

            if (d_using_ghosts_for_dt)
            {
                d_patch_strategy->setDataContext(d_scratch);
                copyTimeDependentData(patch_level, d_current, d_scratch);
            }
            else
            {
                d_patch_strategy->setDataContext(d_current);
            }
        }
        else
        {

            if (d_using_ghosts_for_dt)
            {

                if (d_bdry_sched_advance_new[level_number].isNull())
                {
                    TBOX_ERROR(d_object_name << ":  "
                               << "Attempt to fill new ghost data for timestep"
                               << "computation, but schedule not defined." << endl);
                }

                d_patch_strategy->setDataContext(d_scratch);
                t_new_advance_bdry_fill_comm->start();
                d_bdry_sched_advance_new[level_number]->fillData(new_time);
                t_new_advance_bdry_fill_comm->stop();

            }
            else
            {
                d_patch_strategy->setDataContext(d_new);
            }

        }

        for (hier::PatchLevel<NDIM>::Iterator ip(patch_level); ip; ip++)
        {
            tbox::Pointer<hier::Patch<NDIM> > patch = patch_level->getPatch(ip());

            patch->allocatePatchData(d_temp_var_scratch_data, new_time);
            // "false" argument indicates "initial_time" is false.
            t_patch_num_kernel->start();
            double patch_dt =
                d_patch_strategy->computeStableDtOnPatch(*patch,
                                                         false,
                                                         new_time);
            t_patch_num_kernel->stop();

            dt_next = tbox::Utilities::dmin(dt_next, patch_dt);

            patch->deallocatePatchData(d_temp_var_scratch_data);

        }
        d_patch_strategy->clearDataContext();

    } // !regrid_advance

    patch_level->deallocatePatchData(d_saved_var_scratch_data);

    postprocessFluxData(patch_level,
                        regrid_advance,
                        first_step,
                        last_step);

    t_advance_level->stop();

    t_advance_level_sync->start();

    if (d_distinguish_mpi_reduction_costs)
    {
        tbox::MPI::barrier();
        t_advance_level_sync->stop();
        t_mpi_reductions->start();
    }

    double next_dt = tbox::MPI::minReduction(dt_next) * d_cfl;

    if (d_distinguish_mpi_reduction_costs)
    {
        t_mpi_reductions->stop();
    }
    else
    {
        t_advance_level_sync->stop();
    }

    returnnext_dt;
}// advanceLevel

/*
*************************************************************************
*                                                                       *
* Synchronize data between patch levels according to the standard       *
* hyperbolic flux correction algorithm.                                 *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::standardLevelSynchronization(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const double old_time)
{
    tbox::Array<double> old_times(finest_level - coarsest_level + 1);
    for (int i=coarsest_level; i <= finest_level; i++)
    {
        old_times[i] = old_time;
    }
    standardLevelSynchronization(hierarchy, coarsest_level, finest_level,
                                 sync_time, old_times);
    return;
}// standardLevelSynchronization

void
ParabolicLevelIntegrator::standardLevelSynchronization(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const tbox::Array<double>& old_times)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
            && (coarsest_level < finest_level)
            && (finest_level <= hierarchy->getFinestLevelNumber()));
    assert(old_times.getSize() >= finest_level);
    for (int ln = coarsest_level; ln < finest_level; ln++)
    {
        assert(!(hierarchy->getPatchLevel(ln)).isNull());
        assert(sync_time >= old_times[ln]);
    }
    assert(!(hierarchy->getPatchLevel(finest_level)).isNull());
#endif
    t_std_level_sync->start();

    for (int fine_ln = finest_level; fine_ln > coarsest_level; fine_ln--)
    {
        const int coarse_ln = fine_ln - 1;

#ifdef DEBUG_CHECK_ASSERTIONS
        assert(sync_time >= old_times[coarse_ln]);
#endif

        tbox::Pointer<hier::PatchLevel<NDIM> >
            fine_level = hierarchy->getPatchLevel(fine_ln);
        tbox::Pointer<hier::PatchLevel<NDIM> >
            coarse_level = hierarchy->getPatchLevel(coarse_ln);

        synchronizeLevelWithCoarser(fine_level,
                                    coarse_level,
                                    sync_time,
                                    old_times[coarse_ln]);

        fine_level->deallocatePatchData(d_fluxsum_data);
        fine_level->deallocatePatchData(d_flux_var_data);

        if (coarse_ln > coarsest_level)
        {
            coarse_level->deallocatePatchData(d_flux_var_data);
        }
        else
        {
            if (coarsest_level == 0)
            {
                coarse_level->deallocatePatchData(d_flux_var_data);
                d_have_flux_on_level_zero = false;
            }
        }

    }

    t_std_level_sync->stop();
    return;
}// standardLevelSynchronization

/*
*************************************************************************
*                                                                       *
* Coarsen current solution data from finest hierarchy level specified   *
* down through the coarsest hierarchy level specified, if initial_time  *
* is true (i.e., hierarchy is being constructed at initial simulation   *
* time).  After data is coarsened, the user's initialization routine    *
* is called to reset data (as needed by the application) before         *
* that solution is further coarsened to the next coarser level in the   *
* hierarchy.  If initial_time is false, then this routine does nothing  *
* In that case, interpolation of data from coarser levels is sufficient *
* to set data on new levels in the hierarchy during regridding.         *
*                                                                       *
* NOTE: The fact that this routine does nothing when called at any      *
*       time later than when the AMR hierarchy is constructed initially *
*        may need to change at some point based on application needs.   *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::synchronizeNewLevels(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
            && (coarsest_level < finest_level)
            && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = coarsest_level; ln <= finest_level; ln++)
    {
        assert(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif

    tbox::Pointer<tbox::Timer> t_sync_initial_create =
        tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::sync_initial_create");
    tbox::Pointer<tbox::Timer> t_sync_initial_comm =
        tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::sync_initial_comm");

    t_sync_new_levels->start();

    if (initial_time)
    {

        d_patch_strategy->setDataContext(d_current);

        for (int fine_ln = finest_level; fine_ln > coarsest_level; fine_ln--)
        {
            const int coarse_ln = fine_ln - 1;

            tbox::Pointer<hier::PatchLevel<NDIM> > fine_level =
                hierarchy->getPatchLevel(fine_ln);

            tbox::Pointer<hier::PatchLevel<NDIM> > coarse_level =
                hierarchy->getPatchLevel(coarse_ln);

            t_sync_initial_create->start();
            tbox::Pointer<xfer::CoarsenSchedule<NDIM> > sched =
                d_sync_initial_data->createSchedule(coarse_level,
                                                    fine_level,
                                                    d_patch_strategy);
            t_sync_initial_create->stop();

            t_sync_initial_comm->start();
            sched->coarsenData();
            t_sync_initial_comm->stop();


            for (hier::PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = coarse_level->getPatch(p());

                patch->allocatePatchData(d_temp_var_scratch_data, sync_time);

                d_patch_strategy->initializeDataOnPatch(*patch,
                                                        sync_time,
                                                        initial_time);
                patch->deallocatePatchData(d_temp_var_scratch_data);
            }
        }

        d_patch_strategy->clearDataContext();

    } // if (initial_time)

    t_sync_new_levels->stop();
    return;
}// synchronizeNewLevels

/*
*************************************************************************
*                                                                       *
* Reset time-dependent data on patch level by replacing current data    *
* with new.  The boolean argument is used for odd refinement ratios     *
* (in particular 3 used in certain applications).                       *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::resetTimeDependentData(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double new_time,
    const bool can_be_refined)
{
    tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!patch_level.isNull());
#endif

    hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();

    double cur_time = 0.;
    for (hier::PatchLevel<NDIM>::Iterator ip(patch_level); ip; ip++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = patch_level->getPatch(ip());

        tbox::List<tbox::Pointer<hier::Variable<NDIM> > >::Iterator
            time_dep_var = d_time_dep_variables.listStart();
        while(time_dep_var)
        {

            int cur_indx =
                variable_db->mapVariableAndContextToIndex(time_dep_var(),
                                                          d_current);
            int new_indx =
                variable_db->mapVariableAndContextToIndex(time_dep_var(),
                                                          d_new);

            cur_time = patch->getPatchData(cur_indx)->getTime();

            if (can_be_refined && d_number_time_data_levels == 3)
            {

                int old_indx =
                    variable_db->mapVariableAndContextToIndex(time_dep_var(),
                                                              d_old);

                patch->setPatchData(old_indx, patch->getPatchData(cur_indx));

                patch->setPatchData(cur_indx, patch->getPatchData(new_indx));


            }
            else
            {

                if (d_number_time_data_levels == 3)
                {

                    int old_indx =
                        variable_db->mapVariableAndContextToIndex(time_dep_var(),
                                                                  d_old);

                    patch->setPatchData(old_indx, patch->getPatchData(cur_indx));

                }

                patch->setPatchData(cur_indx, patch->getPatchData(new_indx));

            }

            patch->deallocatePatchData(new_indx);

            time_dep_var++;

        }

    }

    patch_level->setTime(new_time, d_new_patch_init_data);

    if (d_number_time_data_levels == 3)
    {
        patch_level->setTime(cur_time, d_old_time_dep_data);
    }

    return;
}// resetTimeDependentData

/*
*************************************************************************
*                                                                       *
* Discard new data on level.  This is used primarily to reset patch     *
* data after error estimation (e.g., Richardson extrapolation.)         *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::resetDataToPreadvanceState(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
#endif

    /*
     * De-allocate new context
     */
    level->deallocatePatchData(d_new_time_dep_data);
    return;
}// resetDataToPreadvanceState

/*
*************************************************************************
*                                                                       *
* Initialize integration data on all patches on level.  This process    *
* is used at the start of the simulation to set the initial hierarchy   *
* data and after adaptive regridding.  In the second case, the old      *
* level pointer points to the level that existed in the hierarchy       *
* before regridding.  This pointer may be null, in which case it is     *
* ignored.  If it is non-null, then data is copied from the old level   *
* to the new level before the old level is discarded.                   *
*                                                                       *
* Note that we also allocate flux storage for the coarsest AMR          *
* hierarchy level here (i.e., level 0).  The time step sequence on      *
* level 0 is dictated by the user code; so to avoid any memory          *
* management errors, flux storage on level 0 persists as long as the    *
* level does.
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::initializeLevelData(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
            && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        assert(level_number == old_level->getLevelNumber());
    }
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    t_initialize_level_data->start();

    tbox::Pointer<hier::PatchLevel<NDIM> > level =
        hierarchy->getPatchLevel(level_number);

    /*
     * Allocate storage needed to initialize level and fill data
     * from coarser levels in AMR hierarchy, potentially. Since
     * time gets set when we allocate data, re-stamp it to current
     * time if we don't need to allocate.
     */
    if (allocate_data)
    {
        level->allocatePatchData(d_new_patch_init_data, init_data_time);
        level->allocatePatchData(d_old_time_dep_data, init_data_time);
    }
    else
    {
        level->setTime(init_data_time, d_new_patch_init_data);
    }

    /*
     * Create schedules for filling new level and fill data.
     */

    if ((level_number > 0) || !old_level.isNull())
    {
        t_fill_new_level_create->start();

        tbox::Pointer<xfer::RefineSchedule<NDIM> > sched =
            d_fill_new_level->createSchedule(level,
                                             old_level,
                                             level_number-1,
                                             hierarchy,
                                             d_patch_strategy);
        t_fill_new_level_create->stop();

        d_patch_strategy->setDataContext(d_scratch);

        t_fill_new_level_comm->start();
        sched->fillData(init_data_time);
        t_fill_new_level_comm->stop();

        d_patch_strategy->clearDataContext();
    }


    if ((d_number_time_data_levels == 3) && can_be_refined)
    {

        hier::VariableDatabase<NDIM>* variable_db =
            hier::VariableDatabase<NDIM>::getDatabase();

        for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(ip());

            tbox::List<tbox::Pointer<hier::Variable<NDIM> > >::Iterator
                time_dep_var = d_time_dep_variables.listStart();
            while(time_dep_var)
            {
                int old_indx =
                    variable_db->mapVariableAndContextToIndex(time_dep_var(),
                                                              d_old);
                int cur_indx =
                    variable_db->mapVariableAndContextToIndex(time_dep_var(),
                                                              d_current);

                patch->setPatchData(old_indx, patch->getPatchData(cur_indx));

                time_dep_var++;
            }

        }

    }

    /*
     * Initialize data on patch interiors.
     */
    d_patch_strategy->setDataContext(d_current);
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());

        patch->allocatePatchData(d_temp_var_scratch_data, init_data_time);

        d_patch_strategy->initializeDataOnPatch(*patch,
                                                init_data_time,
                                                initial_time);

        patch->deallocatePatchData(d_temp_var_scratch_data);
    }
    d_patch_strategy->clearDataContext();

    t_initialize_level_data->stop();

    return;
}// initializeLevelData

/*
*************************************************************************
*                                                                       *
* Reset hierarchy configuration information where the range of new      *
* hierarchy levels is specified.   The information updated involves     *
* the cached communication schedules maintained by the algorithm.       *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::resetHierarchyConfiguration(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
            && (coarsest_level <= finest_level)
            && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        assert(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    int finest_hiera_level = hierarchy->getFinestLevelNumber();

    d_bdry_sched_advance.resizeArray(finest_hiera_level+1);
    d_bdry_sched_advance_new.resizeArray(finest_hiera_level+1);

    for (int ln = coarsest_level; ln <= finest_hiera_level; ln++)
    {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        t_advance_bdry_fill_create->start();
        d_bdry_sched_advance[ln] =
            d_bdry_fill_advance->createSchedule(level,
                                                ln-1,
                                                hierarchy,
                                                d_patch_strategy);
        t_advance_bdry_fill_create->stop();

        if (!d_lag_dt_computation && d_using_ghosts_for_dt)
        {
            t_new_advance_bdry_fill_create->start();
            d_bdry_sched_advance_new[ln] =
                d_bdry_fill_advance_new->createSchedule(level,
                                                        ln-1,
                                                        hierarchy,
                                                        d_patch_strategy);
            t_new_advance_bdry_fill_create->stop();
        }

    }
    return;
}// resetHierarchyConfiguration

/*
*************************************************************************
*                                                                       *
* Call patch routines to tag cells near large gradients.                *
* These cells will be refined.                                          *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::applyGradientDetector(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
            && (level_number <= hierarchy->getFinestLevelNumber()));
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    t_apply_gradient_detector->start();

    tbox::Pointer<hier::PatchLevel<NDIM> > level =
        hierarchy->getPatchLevel(level_number);

    level->allocatePatchData(d_saved_var_scratch_data, error_data_time);

    d_patch_strategy->setDataContext(d_scratch);

    t_error_bdry_fill_comm->start();
    d_bdry_sched_advance[level_number]->fillData(error_data_time);
    t_error_bdry_fill_comm->stop();

    for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(ip());
        d_patch_strategy->
            tagGradientDetectorCells(*patch,
                                     error_data_time,
                                     initial_time,
                                     tag_index,
                                     uses_richardson_extrapolation_too);
    }

    d_patch_strategy->clearDataContext();

    level->deallocatePatchData(d_saved_var_scratch_data);

    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

/*
*************************************************************************
*                                                                       *
* Call patch routines to tag cells for refinement using Richardson      *
* extrapolation.    Richardson extrapolation requires two copies of     *
* the solution to compare.  The NEW context holds the solution          *
* computed on the fine level and coarsened, whereas the CURRENT         *
* context holds the solution integrated on the coarse level after       *
* coarsening the initial data from the fine level.                      *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::applyRichardsonExtrapolation(
    const tbox::Pointer<hier::PatchLevel<NDIM> > level,
    const double error_data_time,
    const int tag_index,
    const double deltat,
    const int error_coarsen_ratio,
    const bool initial_time,
    const bool uses_gradient_detector_too)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
#endif
    /*
     * Compare solutions computed on level (stored in NEW context) and on
     * the coarser level (stored in CURR context) on the patches of the
     * coarser level.  The patch strategy implements the compare operations
     * performed on each patch.
     */

    int error_level_number =
        level->getNextCoarserHierarchyLevelNumber() + 1;

    for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(ip());

        d_patch_strategy->
            tagRichardsonExtrapolationCells(*patch,
                                            error_level_number,
                                            d_new,     //  finer context
                                            d_current, //  coarser context
                                            error_data_time,
                                            deltat,
                                            error_coarsen_ratio,
                                            initial_time,
                                            tag_index,
                                            uses_gradient_detector_too);
    }
    return;
}// applyRichardsonExtrapolation

/*
*************************************************************************
*                                                                       *
* The Richardson extrapolation algorithm requires a coarsened version   *
* of the level on which error estiamtion is performed.  This routine    *
* is used to coarsen data from a level in the AMR hierarchy to some     *
* coarsened version of it.  Note that this routine will be called twice *
* The init_coarse_level boolean argument indicates whether data is    *
* set on the coarse level by coarsening the "old" time level solution   *
* or by coarsening the "new" solution on the fine level (i.e., after    *
* it has been advanced).                                                *
*                                                                       *
* The contexts used for coarsening old data depends on the number of    *
* time levels.  We always want to use data at the oldest time on the    *
* fine level, coarsened to the CURRENT context on the coarse level.     *
* Thus, if the problem uses two time levels, we coarsen data from       *
* CURRENT on fine level (since CURRENT is the oldest time maintained)   *
* to CURRENT on the coarse level.  If the problem uses three time       *
* levels, we coarsen from OLD on the fine level (since OLD is the       *
* time maintained) to CURRENT on the coarse level.                      *
*                                                                       *
* When the boolean is false, indicating we are operating at the new     *
* time, we coarsen the time advanced solution at the NEW context on     *
* the fine level to the NEW context on the coarse level so that they    *
* may be compared later.                                                *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::coarsenDataForRichardsonExtrapolation(
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const tbox::Pointer<hier::PatchLevel<NDIM> > coarse_level,
    const double coarsen_data_time,
    const bool before_advance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
            && (level_number <= hierarchy->getFinestLevelNumber()));
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
    assert(!coarse_level.isNull());
#endif
    t_coarsen_rich_extrap->start();

    tbox::Pointer<hier::PatchLevel<NDIM> > hier_level =
        hierarchy->getPatchLevel(level_number);

    if (before_advance)
    {

        coarse_level->allocatePatchData(d_new_patch_init_data,
                                        coarsen_data_time);

        if (d_number_time_data_levels == 3)
        {
            d_patch_strategy->setDataContext(d_old);
        }
        else
        {
            d_patch_strategy->setDataContext(d_current);
        }

        d_coarsen_rich_extrap_init->
            createSchedule(coarse_level, hier_level, d_patch_strategy)->
            coarsenData();

        d_patch_strategy->clearDataContext();

    }
    else
    {

        coarse_level->allocatePatchData(d_new_time_dep_data,
                                        coarsen_data_time);

        d_patch_strategy->setDataContext(d_new);

        d_coarsen_rich_extrap_final->
            createSchedule(coarse_level, hier_level, d_patch_strategy)->
            coarsenData();

        d_patch_strategy->clearDataContext();

    }

    t_coarsen_rich_extrap->stop();
    return;
}// coarsenDataForRichardsonExtrapolation

/*
*************************************************************************
*                                                                       *
* Register given variable with algorithm according to specified         *
* algorithm role (i.e., PAR_VAR_TYPE).  Assignment of descriptor        *
* indices to variable lists, component selectors, and communication     *
* algorithms takes place here.  The different cases are:                *
*                                                                       *
* TIME_DEP:                                                             *
*            The number of factories depends on the number of time      *
*            levels of data that must be stored on patches to satisfy   *
*            regridding reqs.  Currently, there are two possibilities:  *
*                                                                       *
*            (1) If the coarsen ratios between levels are even, the     *
*                error coarsening ratio will be two and so only two     *
*                time levels of data must be maintained on every level  *
*                but the finest as usual.                               *
*                                                                       *
*            (2) If the coarsen ratios between levels are three, and    *
*                time integration is used during regridding (e.g., Rich-*
*                ardson extrapolation), then three time levels of data  *
*                must be maintained on every level but the finest so    *
*                that error estimation can be executed properly.        *
*                                                                       *
*            In case (1), three factories are needed:                   *
*                         SCRATCH, CURRENT, NEW.                        *
*            In case (2), four factories are needed:                    *
*                         SCRATCH, OLD, CURRENT, NEW.                   *
*                                                                       *
*            SCRATCH index is added to d_saved_var_scratch_data.        *
*            CURRENT index is added to d_new_patch_init_data.           *
*            NEW index is added to d_new_time_dep_data.                 *
*                                                                       *
* INPUT:                                                                *
*            Only one time level of data is maintained and once values  *
*            are set on patches, they do not change in time.            *
*                                                                       *
*            Two factories are needed: SCRATCH, CURRENT.                *
*                                                                       *
*            SCRATCH index is added to d_saved_var_scratch_data.        *
*            CURRENT index is added to d_new_patch_init_data.           *
*                                                                       *
* NO_FILL:                                                              *
*            Only one time level of data is stored and no scratch space *
*            is used.  Data may be set and manipulated at will in user  *
*            routines.  Data (including ghost values) is never touched  *
*            outside of user routines.                                  *
*                                                                       *
*            Two factories are needed: CURRENT, SCRATCH.                *
*                                                                       *
*            CURRENT index is added to d_new_patch_init_data.           *
*            SCRATCH index is needed only for temporary work space to   *
*            fill new patch levels.                                     *
*                                                                       *
* FLUX:                                                                 *
*            One factory is needed: SCRATCH.                            *
*                                                                       *
*            SCRATCH index is added to d_flux_var_data.                 *
*                                                                       *
*            Additionally, a variable for flux integral data is created *
*            for each FLUX variable. It has a single factory, SCRATCH,  *
*            which is added to d_fluxsum_data.                          *
*                                                                       *
* TEMPORARY:                                                            *
*            One factory needed: SCRATCH.                               *
*            SCRATCH index is added to d_temp_var_scratch_data.         *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::registerVariable(
    const tbox::Pointer<hier::Variable<NDIM> > var,
    const hier::IntVector<NDIM> ghosts,
    const PAR_VAR_TYPE h_v_type,
    const tbox::Pointer<xfer::Geometry<NDIM> > transfer_geom,
    const string& coarsen_name,
    const string& refine_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!var.isNull());
    assert(!transfer_geom.isNull());
#endif

    hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();

    const hier::IntVector<NDIM> zero_ghosts(0);

    d_all_variables.appendItem(var);

    switch(h_v_type)
    {
        case TIME_DEP:
        {
            d_time_dep_variables.appendItem(var);

            int cur_id = variable_db->registerVariableAndContext(var,
                                                                 d_current,
                                                                 zero_ghosts);
            int new_id = variable_db->registerVariableAndContext(var,
                                                                 d_new,
                                                                 zero_ghosts);
            int scr_id = variable_db->registerVariableAndContext(var,
                                                                 d_scratch,
                                                                 ghosts);

            d_saved_var_scratch_data.setFlag(scr_id);

            d_new_patch_init_data.setFlag(cur_id);

            d_new_time_dep_data.setFlag(new_id);

            /*
             * Register variable and context needed for restart.
             */
            hier::VariableDatabase<NDIM>::getDatabase()->
                registerPatchDataForRestart(cur_id);

            /*
             * Set boundary fill schedules for time-dependent variable.
             * If time interpolation operator is non-NULL, regular advance
             * bdry fill algorithm will time interpolate between current and
             * new data on coarser levels, and fill from current data on
             * same level.  New advance bdry fill algorithm will time interpolate
             * between current and new data on coarser levels, and fill from new
             * data on same level.  If time interpolation operator is NULL,
             * regular and new bdry fill algorithms will use current and new
             * data, respectively.
             */

            tbox::Pointer<xfer::RefineOperator<NDIM> > refine_op =
                transfer_geom->lookupRefineOperator(var, refine_name);

            tbox::Pointer<xfer::TimeInterpolateOperator<NDIM> > time_int =
                transfer_geom->lookupTimeInterpolateOperator(var);

            d_bdry_fill_advance->registerRefine(
                scr_id, cur_id, cur_id, new_id, scr_id, refine_op, time_int);
            d_bdry_fill_advance_new->registerRefine(
                scr_id, new_id, cur_id, new_id, scr_id, refine_op, time_int);
            d_fill_new_level->registerRefine(
                cur_id, cur_id, cur_id, new_id, scr_id, refine_op, time_int);

            /*
             * For data synchronization between levels, the coarsen algorithm
             * will coarsen new data on finer level to new data on coarser.
             * Recall that coarser level data pointers will not be reset until
             * after synchronization so we always coarsen to new
             * (see synchronizeLevelWithCoarser routine).
             */
            tbox::Pointer<xfer::CoarsenOperator<NDIM> > coarsen_op =
                transfer_geom->lookupCoarsenOperator(var, coarsen_name);

            d_coarsen_sync_data->registerCoarsen(new_id, new_id, coarsen_op);

            d_sync_initial_data->registerCoarsen(cur_id, cur_id, coarsen_op);

            /*
             * Coarsen operations used in Richardson extrapolation.  The init
             * initializes data on coarser level, before the coarse level
             * advance.  If two time levels are used, coarsening occurs between
             * the CURRENT context on both levels.  If three levels are used,
             * coarsening occurs between the OLD context on the fine level and
             * the CURRENT context on the coarse level.  The final coarsen
             * algorithm coarsens data after it has been advanced on the fine
             * level to the NEW context on the coarser level.
             */
            if (d_number_time_data_levels == 3)
            {
                int old_id = variable_db->registerVariableAndContext(var,
                                                                     d_old,
                                                                     zero_ghosts);
                d_old_time_dep_data.setFlag(old_id);

                d_bdry_fill_advance_old->registerRefine(
                    scr_id, cur_id, old_id, new_id, scr_id, refine_op, time_int);

                d_coarsen_rich_extrap_init->
                    registerCoarsen(cur_id, old_id, coarsen_op);
            }
            else
            {
                d_coarsen_rich_extrap_init->
                    registerCoarsen(cur_id, cur_id, coarsen_op);
            }

            d_coarsen_rich_extrap_final->
                registerCoarsen(new_id, new_id, coarsen_op);

            break;
        }
        case INPUT:
        {
            int cur_id = variable_db->registerVariableAndContext(var,
                                                                 d_current,
                                                                 zero_ghosts);
            int scr_id = variable_db->registerVariableAndContext(var,
                                                                 d_scratch,
                                                                 ghosts);

            d_saved_var_scratch_data.setFlag(scr_id);

            d_new_patch_init_data.setFlag(cur_id);

            /*
             * Register variable and context needed for restart.
             */
            hier::VariableDatabase<NDIM>::getDatabase()->
                registerPatchDataForRestart(cur_id);

            /*
             * Bdry algorithms for input variables will fill from current only.
             */
            tbox::Pointer<xfer::RefineOperator<NDIM> > refine_op =
                transfer_geom->lookupRefineOperator(var, refine_name);

            d_bdry_fill_advance->registerRefine(
                scr_id, cur_id, scr_id, refine_op);
            d_bdry_fill_advance_new->registerRefine(
                scr_id, cur_id, scr_id, refine_op);
            d_fill_new_level->registerRefine(
                cur_id, cur_id, scr_id, refine_op);

            /*
             * At initialization, it may be necessary to coarsen INPUT data
             * up through the hierarchy so that all levels are consistent.
             */
            tbox::Pointer<xfer::CoarsenOperator<NDIM> > coarsen_op =
                transfer_geom->lookupCoarsenOperator(var, coarsen_name);

            d_sync_initial_data->registerCoarsen(cur_id, cur_id, coarsen_op);

            /*
             * Coarsen operation for setting initial data on coarser level
             * in the Richardson extrapolation algorithm.
             */
            d_coarsen_rich_extrap_init->
                registerCoarsen(cur_id, cur_id, coarsen_op);

            break;
        }
        case NO_FILL:
        {
            int cur_id = variable_db->registerVariableAndContext(var,
                                                                 d_current,
                                                                 ghosts);

            int scr_id = variable_db->registerVariableAndContext(var,
                                                                 d_scratch,
                                                                 ghosts);

            d_new_patch_init_data.setFlag(cur_id);

            /*
             * Register variable and context needed for restart.
             */
            hier::VariableDatabase<NDIM>::getDatabase()->
                registerPatchDataForRestart(cur_id);

            tbox::Pointer<xfer::RefineOperator<NDIM> > refine_op =
                transfer_geom->lookupRefineOperator(var, refine_name);

            d_fill_new_level->registerRefine(
                cur_id, cur_id, scr_id, refine_op);

            /*
             * Coarsen operation for setting initial data on coarser level
             * in the Richardson extrapolation algorithm.
             */

            tbox::Pointer<xfer::CoarsenOperator<NDIM> > coarsen_op =
                transfer_geom->lookupCoarsenOperator(var, coarsen_name);

            d_coarsen_rich_extrap_init->
                registerCoarsen(cur_id, cur_id, coarsen_op);

            break;
        }
        case FLUX:
        {
            /*
             * Note that we force all flux variables to hold double
             * precision data and be face- or side-centered.  Also,
             * for each flux variable, a corresponding "fluxsum"
             * variable is created to manage synchronization of data
             * betweeen patch levels in the hierarchy.
             */
            const tbox::Pointer<pdat::FaceVariable<DIM,double> > face_var(var);
            const tbox::Pointer<pdat::SideVariable<DIM,double> > side_var(var);

            if (!face_var.isNull())
            {
                if (d_flux_side_registered)
                {
                    TBOX_ERROR(d_object_name << ":  "
                               << "Attempt to register FaceVariable when "
                               << "SideVariable already registered." << endl);
                }

                d_flux_is_face = true;

            }
            else if (!side_var.isNull())
            {
                if (d_flux_face_registered)
                {
                    TBOX_ERROR(d_object_name << ":  "
                               << "Attempt to register SideVariable when "
                               << "FaceVariable already registered." << endl);
                }

                d_flux_is_face = false;

            }
            else
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "Flux is neither face- or side-centered." << endl);
            }

            d_flux_variables.appendItem(var);

            int scr_id = variable_db->registerVariableAndContext(var,
                                                                 d_scratch,
                                                                 ghosts);

            d_flux_var_data.setFlag(scr_id);

            string var_name = var->getName();
            string fs_suffix = "_fluxsum";
            string fsum_name = var_name;
            fsum_name += fs_suffix;

            tbox::Pointer<hier::Variable<NDIM> > fluxsum;

            if (d_flux_is_face)
            {
                fluxsum = new pdat::OuterfaceVariable<DIM,double>(
                    fsum_name,
                    ((tbox::Pointer<pdat::FaceDataFactory<DIM,double> >)
                     var->getPatchDataFactory())->getDefaultDepth());
                d_flux_face_registered = true;
            }
            else
            {
                fluxsum = new pdat::OutersideVariable<DIM,double>(
                    fsum_name,
                    ((tbox::Pointer<pdat::SideDataFactory<DIM,double> >)
                     var->getPatchDataFactory())->getDefaultDepth());
                d_flux_side_registered = true;
            }

            d_fluxsum_variables.appendItem(fluxsum);

            int fs_id = variable_db->registerVariableAndContext(fluxsum,
                                                                d_scratch,
                                                                zero_ghosts);

            d_fluxsum_data.setFlag(fs_id);

            tbox::Pointer<xfer::CoarsenOperator<NDIM> > coarsen_op =
                transfer_geom->lookupCoarsenOperator(fluxsum, coarsen_name);

            d_coarsen_fluxsum->registerCoarsen(scr_id, fs_id, coarsen_op);

            break;
        }
        case TEMPORARY:
        {
            int scr_id = variable_db->registerVariableAndContext(var,
                                                                 d_scratch,
                                                                 ghosts);

            d_temp_var_scratch_data.setFlag(scr_id);

            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "unknown PAR_VAR_TYPE = " << h_v_type << endl);
        }
    }
    return;
}// registerVariable

/*
*************************************************************************
*                                                                       *
* Print all class data for ParabolicLevelIntegrator object.             *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::printClassData(ostream& os) const
{
    os << "\nParabolicLevelIntegrator::printClassData..." << endl;
    os << "ParabolicLevelIntegrator: this = "
       << (ParabolicLevelIntegrator*)this << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_cfl = " << d_cfl << "\n"
       << "d_cfl_init = " << d_cfl_init << endl;
    os << "d_lag_dt_computation = " << d_lag_dt_computation << "\n"
       << "d_using_ghosts_for_dt = "
       << d_using_ghosts_for_dt << endl;
    os << "d_patch_strategy = "
       << (ParabolicPatchStrategy<NDIM>*)d_patch_strategy << endl;
    os << "NOTE: Not printing variable arrays, ComponentSelectors, communication schedules, etc." << endl;
    return;
}// printClassData

/*
*************************************************************************
*                                                                       *
* Writes out the class version number, d_cfl, d_cfl_init, 		*
* d_lag_dt_computation, and d_using_ghosts_for_dt to the database.	*
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::putToDatabase(
    tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

    db->putInteger("PARABOLIC_LEVEL_INTEGRATOR_VERSION",
                   PARABOLIC_LEVEL_INTEGRATOR_VERSION);

    db->putDouble("d_cfl", d_cfl);
    db->putDouble("d_cfl_init", d_cfl_init);
    db->putBool("d_lag_dt_computation", d_lag_dt_computation);
    db->putBool("d_using_ghosts_for_dt", d_using_ghosts_for_dt);
    return;
}// putToDatabase

/*
*************************************************************************
*                                                                       *
* Utility routines to retrieve variable contexts used by integrator.    *
*                                                                       *
*************************************************************************
*/
tbox::Pointer<hier::VariableContext>
ParabolicLevelIntegrator::getCurrentContext() const
{
    return d_current;
}// getCurrentContext

tbox::Pointer<hier::VariableContext>
ParabolicLevelIntegrator::getNewContext() const
{
    return d_new;
}// getNewContext

tbox::Pointer<hier::VariableContext>
ParabolicLevelIntegrator::getOldContext() const
{
    return d_old;
}// getOldContext

tbox::Pointer<hier::VariableContext>
ParabolicLevelIntegrator::getScratchContext() const
{
    return d_scratch;
}// getScratchContext

tbox::Pointer<hier::VariableContext>
ParabolicLevelIntegrator::getPlotContext() const
{
    return d_plot_context;
}// getPlotContext

bool
ParabolicLevelIntegrator::usingRefinedTimestepping() const
{
    return d_using_time_refinement;
}// usingRefinedTimestepping

/////////////////////////////// PROTECTED ////////////////////////////////////

/*
*************************************************************************
*                                                                       *
* Reads in cfl, cfl_init, lag_dt_computation, and 			*
* using_ghosts_to_compute_dt from the input database.  			*
* Note all restart values are overriden with values from the input	*
* database.								*
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::getFromInput(
    tbox::Pointer<tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

    if (db->keyExists("cfl"))
    {
        d_cfl = db->getDouble("cfl");
    }
    else
    {
        if (!is_from_restart)
        {
            d_cfl = db->getDoubleWithDefault("cfl", d_cfl);
        }
    }

    if (db->keyExists("cfl_init"))
    {
        d_cfl_init = db->getDouble("cfl_init");
    }
    else
    {
        if (!is_from_restart)
        {
            d_cfl_init = db->getDoubleWithDefault("cfl_init", d_cfl_init);
        }
    }


    if (db->keyExists("lag_dt_computation"))
    {
        d_lag_dt_computation = db->getBool("lag_dt_computation");
    }
    else
    {
        if (!is_from_restart)
        {
            d_lag_dt_computation =
                db->getDoubleWithDefault("lag_dt_computation",
                                         d_lag_dt_computation);
        }
    }


    if (db->keyExists("using_ghosts_to_compute_dt"))
    {
        d_using_ghosts_for_dt = db->getBool("using_ghosts_to_compute_dt");
    }
    else
    {
        if (!is_from_restart)
        {
            d_using_ghosts_for_dt =
                db->getDoubleWithDefault("using_ghosts_for_dt",
                                         d_using_ghosts_for_dt);
            TBOX_WARNING(d_object_name << ":  "
                         << "Key data `using_ghosts_to_compute_dt' not found in input."
                         << "  Using default value " << d_using_ghosts_for_dt << endl);
        }
    }

    if (db->keyExists("distinguish_mpi_reduction_costs"))
    {
        d_distinguish_mpi_reduction_costs =
            db->getBool("distinguish_mpi_reduction_costs");
    }

    return;
}// getFromInput

/*
*************************************************************************
*                                                                       *
* First, gets the database corresponding to the object_name from the    *
* restart file.   If this database exists, this method checks to make   *
* sure that the version number of the class matches the version number  *
* of the restart file.  If they match, then d_cfl, d_cfl_init,          *
* d_lag_dt_computation, and d_using_ghosts_to_compute_dt are read from  *
* restart database.  		                                        *
* Note all restart values can be overriden with values from the input	*
* database.								*
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::getFromRestart()
{

    tbox::Pointer<tbox::Database> root_db =
        tbox::RestartManager::getManager()->getRootDatabase();

    tbox::Pointer<tbox::Database> db;
    if (root_db->isDatabase(d_object_name))
    {
        db = root_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name << " not found in restart file" << endl);
    }

    int ver = db->getInteger("PARABOLIC_LEVEL_INTEGRATOR_VERSION");
    if (ver != PARABOLIC_LEVEL_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different "
                   << "than class version." << endl);
    }

    d_cfl = db->getDouble("d_cfl");
    d_cfl_init = db->getDouble("d_cfl_init");
    d_lag_dt_computation = db->getBool("d_lag_dt_computation");
    d_using_ghosts_for_dt = db->getBool("d_using_ghosts_for_dt");

    return;
}// getFromRestart

/*
*************************************************************************
*                                                                       *
* Process FLUX and FLUX INTEGRAL data before integration on the level.  *
*                                                                       *
* We allocate FLUX storage if appropriate.                              *
*                                                                       *
* If the advance is not temporary, we also zero out the FLUX INTEGRALS  *
* on the first step of any level finer than level zero.                 *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::preprocessFluxData(
    const tbox::Pointer<hier::PatchLevel<NDIM> > level,
    const double cur_time,
    const double new_time,
    const bool regrid_advance,
    const bool first_step,
    const bool last_step)
{
    NULL_USE(last_step);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
    assert(cur_time <= new_time);
#endif

    hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();

    const int level_number = level->getLevelNumber();

    if (!regrid_advance)
    {
        if (((level_number > 0) && first_step) ||
             ((level_number == 0) && !d_have_flux_on_level_zero))
        {
            level->allocatePatchData(d_flux_var_data, new_time);
            if (level_number == 0)
            {
                d_have_flux_on_level_zero = true;
            }
        }
    }
    else
    {
        if (first_step)
        {
            level->allocatePatchData(d_flux_var_data, new_time);
        }
    }

    if (!regrid_advance && (level_number > 0))
    {

        if (first_step)
        {

            level->allocatePatchData(d_fluxsum_data, new_time);

            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());

                tbox::List<tbox::Pointer<hier::Variable<NDIM> > >::Iterator
                    fs_var = d_fluxsum_variables.listStart();

                while (fs_var)
                {
                    int fsum_id =
                        variable_db->mapVariableAndContextToIndex(fs_var(),
                                                                  d_scratch);

                    if (d_flux_is_face)
                    {
                        tbox::Pointer<pdat::OuterfaceData<DIM,double> > fsum_data =
                            patch->getPatchData(fsum_id);


#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(!fsum_data.isNull());
#endif
                        fsum_data->fillAll(0.0);
                    }
                    else
                    {
                        tbox::Pointer<pdat::OutersideData<DIM,double> > fsum_data =
                            patch->getPatchData(fsum_id);

#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(!fsum_data.isNull());
#endif
                        fsum_data->fillAll(0.0);
                    }

                    fs_var++;
                }
            }
        }
        else
        {
            level->setTime(new_time, d_fluxsum_data);
        }
    } // if (!regrid_advance && (level_number > 0))
    return;
}// preprocessFluxData

/*
*************************************************************************
*                                                                       *
* Process FLUX and FLUX INTEGRAL data after advancing the solution on   *
* the level.  During normal integration steps, the flux integrals are   *
* updated for subsequent synchronization by adding FLUX values to       *
* flux integrals.                                                       *
*                                                                       *
* If the advance is not temporary (regular integration step):           *
* 1) If the level is the finest in the hierarchy, FLUX data is          *
*    deallocated.  It is not used during synchronization, and is only   *
*    maintained if needed for the advance.                              *
*                                                                       *
* 2) If the level is not the coarsest in the hierarchy, update the      *
*    flux integrals for later synchronization by adding FLUX values to  *
*    flux integrals.                                                    *
*                                                                       *
* If the advance is temporary, deallocate the flux data if first step.  *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::postprocessFluxData(
    const tbox::Pointer<hier::PatchLevel<NDIM> > level,
    const bool regrid_advance,
    const bool first_step,
    const bool last_step)
{
    NULL_USE(last_step);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
#endif

    if (regrid_advance && first_step)
    {
        level->deallocatePatchData(d_flux_var_data);
    }

    if (!regrid_advance && (level->getLevelNumber() > 0))
    {

        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());

            tbox::List<tbox::Pointer<hier::Variable<NDIM> > >::Iterator
                flux_var = d_flux_variables.listStart();
            tbox::List<tbox::Pointer<hier::Variable<NDIM> > >::Iterator
                fluxsum_var = d_fluxsum_variables.listStart();

            const hier::Index<NDIM>& ilo = patch->getBox().lower();
            const hier::Index<NDIM>& ihi = patch->getBox().upper();

            while (flux_var)
            {

                tbox::Pointer<hier::PatchData<NDIM> > flux_data =
                    patch->getPatchData(flux_var(), d_scratch);
                tbox::Pointer<hier::PatchData<NDIM> > fsum_data =
                    patch->getPatchData(fluxsum_var(), d_scratch);

                tbox::Pointer<pdat::FaceData<DIM,double> > fflux_data;
                tbox::Pointer<pdat::OuterfaceData<DIM,double> > ffsum_data;

                tbox::Pointer<pdat::SideData<DIM,double> > sflux_data;
                tbox::Pointer<pdat::OutersideData<DIM,double> > sfsum_data;

                int ddepth;
                hier::IntVector<NDIM> flux_ghosts;

                if (d_flux_is_face)
                {
                    fflux_data = flux_data;
                    ffsum_data = fsum_data;

#ifdef DEBUG_CHECK_ASSERTIONS
                    assert(!fflux_data.isNull() && !ffsum_data.isNull());
                    assert(fflux_data->getDepth() == ffsum_data->getDepth());
#endif
                    ddepth = fflux_data->getDepth();
                    flux_ghosts = fflux_data->getGhostCellWidth();
                }
                else
                {
                    sflux_data = flux_data;
                    sfsum_data = fsum_data;

#ifdef DEBUG_CHECK_ASSERTIONS
                    assert(!sflux_data.isNull() && !sfsum_data.isNull());
                    assert(sflux_data->getDepth() == sfsum_data->getDepth());
#endif
                    ddepth = sflux_data->getDepth();
                    flux_ghosts = sflux_data->getGhostCellWidth();
                }

                for (int d = 0; d < ddepth; d++)
                {
                    // loop over lower and upper parts of outer face/side arrays
                    for (int ifs = 0; ifs < 2; ifs++)
                    {
                        if (DIM == 1)
                        {
                            upfluxsum1d_(ilo(0),ihi(0),
                                         flux_ghosts(0),
                                         ifs,
                                         fflux_data->getPointer(0,d),
                                         ffsum_data->getPointer(0,ifs,d));
                        }
                        else
                        {

                            if (d_flux_is_face)
                            {
                                if (DIM == 2)
                                {
                                    upfluxsumface2d0_(ilo(0),ilo(1),ihi(0),ihi(1),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      ifs,
                                                      fflux_data->getPointer(0,d),
                                                      ffsum_data->getPointer(0,ifs,d));
                                    upfluxsumface2d1_(ilo(0),ilo(1),ihi(0),ihi(1),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      ifs,
                                                      fflux_data->getPointer(1,d),
                                                      ffsum_data->getPointer(1,ifs,d));
                                }

                                if (DIM == 3)
                                {
                                    upfluxsumface3d0_(ilo(0),ilo(1),ilo(2),
                                                      ihi(0),ihi(1),ihi(2),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      flux_ghosts(2),
                                                      ifs,
                                                      fflux_data->getPointer(0,d),
                                                      ffsum_data->getPointer(0,ifs,d));
                                    upfluxsumface3d1_(ilo(0),ilo(1),ilo(2),
                                                      ihi(0),ihi(1),ihi(2),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      flux_ghosts(2),
                                                      ifs,
                                                      fflux_data->getPointer(1,d),
                                                      ffsum_data->getPointer(1,ifs,d));
                                    upfluxsumface3d2_(ilo(0),ilo(1),ilo(2),
                                                      ihi(0),ihi(1),ihi(2),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      flux_ghosts(2),
                                                      ifs,
                                                      fflux_data->getPointer(2,d),
                                                      ffsum_data->getPointer(2,ifs,d));
                                }
                            }
                            else
                            {
                                if (DIM == 2)
                                {
                                    upfluxsumside2d0_(ilo(0),ilo(1),ihi(0),ihi(1),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      ifs,
                                                      sflux_data->getPointer(0,d),
                                                      sfsum_data->getPointer(0,ifs,d));
                                    upfluxsumside2d1_(ilo(0),ilo(1),ihi(0),ihi(1),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      ifs,
                                                      sflux_data->getPointer(1,d),
                                                      sfsum_data->getPointer(1,ifs,d));
                                }
                                if (DIM == 3)
                                {
                                    upfluxsumside3d0_(ilo(0),ilo(1),ilo(2),
                                                      ihi(0),ihi(1),ihi(2),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      flux_ghosts(2),
                                                      ifs,
                                                      sflux_data->getPointer(0,d),
                                                      sfsum_data->getPointer(0,ifs,d));
                                    upfluxsumside3d1_(ilo(0),ilo(1),ilo(2),
                                                      ihi(0),ihi(1),ihi(2),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      flux_ghosts(2),
                                                      ifs,
                                                      sflux_data->getPointer(1,d),
                                                      sfsum_data->getPointer(1,ifs,d));
                                    upfluxsumside3d2_(ilo(0),ilo(1),ilo(2),
                                                      ihi(0),ihi(1),ihi(2),
                                                      flux_ghosts(0),
                                                      flux_ghosts(1),
                                                      flux_ghosts(2),
                                                      ifs,
                                                      sflux_data->getPointer(2,d),
                                                      sfsum_data->getPointer(2,ifs,d));
                                }
                            }  // if face operations vs. side operations
                        }
                    }  // loop over lower and upper sides/faces
                }  // loop over depth

                flux_var++;
                fluxsum_var++;

            }  // loop over flux variables
        }  // loop over patches
    }  // if !regrid_advance and level number > 0 ....
    return;
}// postprocessFluxData

/*
*************************************************************************
*                                                                       *
* Copy time-dependent data from source to destination on level.         *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::copyTimeDependentData(
    const tbox::Pointer<hier::PatchLevel<NDIM> > level,
    const tbox::Pointer<hier::VariableContext> src_context,
    const tbox::Pointer<hier::VariableContext> dst_context)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!level.isNull());
    assert(!src_context.isNull());
    assert(!src_context.isNull());
#endif

    for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(ip());

        tbox::List<tbox::Pointer<hier::Variable<NDIM> > >::Iterator
            time_dep_var = d_time_dep_variables.listStart();
        while (time_dep_var)
        {
            tbox::Pointer<hier::PatchData<NDIM> > src_data =
                patch->getPatchData(time_dep_var(), src_context);
            tbox::Pointer<hier::PatchData<NDIM> > dst_data =
                patch->getPatchData(time_dep_var(), dst_context);

            dst_data->copy(*src_data);
            time_dep_var++;
        }

    }
    return;
}// copyTimeDependentData

/*
*************************************************************************
*                                                                       *
* Synchronize data between coarse and fine patch levels according to    *
* the standard hyperbolic flux correction algorithm.  The steps of      *
* the algorithm are:                                                    *
*                                                                       *
*    (1) Replace coarse time-space flux integrals at coarse-fine        *
*        boundaries with time-space flux integrals computed on fine     *
*        level.                                                         *
*    (2) Repeat conservative difference on coarse level with corrected  *
*        fluxes.                                                        *
*    (3) Conservatively coarsen solution on interior of fine level to   *
*        coarse level.                                                  *
*                                                                       *
*************************************************************************
*/
void
ParabolicLevelIntegrator::synchronizeLevelWithCoarser(
    const tbox::Pointer<hier::PatchLevel<NDIM> > fine_level,
    const tbox::Pointer<hier::PatchLevel<NDIM> > coarse_level,
    const double sync_time,
    const double coarse_sim_time)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!fine_level.isNull());
    assert(!coarse_level.isNull());
    assert(coarse_level->getLevelNumber() == (fine_level->getLevelNumber()-1));
#endif

    tbox::Pointer<tbox::Timer> t_coarsen_fluxsum_create =
        tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::coarsen_fluxsum_create");
    tbox::Pointer<tbox::Timer> t_coarsen_fluxsum_comm =
        tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::coarsen_fluxsum_comm");
    tbox::Pointer<tbox::Timer> t_coarsen_sync_create =
        tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::coarsen_sync_create");
    tbox::Pointer<tbox::Timer> t_coarsen_sync_comm =
        tbox::TimerManager::getManager()->
        getTimer("ParabolicLevelIntegrator::coarsen_sync_comm");

    /*
     * Coarsen flux integrals around fine patch boundaries to coarser level
     * and replace coarse flux information where appropriate.  NULL patch
     * model is passed in to avoid over complicating coarsen process;
     * i.e. patch model is not needed in coarsening of flux integrals.
     */

    t_coarsen_fluxsum_create->start();
    tbox::Pointer<xfer::CoarsenSchedule<NDIM> > sched =
        d_coarsen_fluxsum->createSchedule(coarse_level, fine_level, NULL);
    t_coarsen_fluxsum_create->stop();

    t_coarsen_fluxsum_comm->start();
    sched->coarsenData();
    t_coarsen_fluxsum_comm->stop();

    /*
     * Repeat conservative difference on coarser level.
     */
    coarse_level->allocatePatchData(d_saved_var_scratch_data, coarse_sim_time);
    coarse_level->setTime(coarse_sim_time, d_flux_var_data);

    d_patch_strategy->setDataContext(d_scratch);
    t_advance_bdry_fill_comm->start();
    d_bdry_sched_advance[coarse_level->getLevelNumber()]->
        fillData(coarse_sim_time);
    t_advance_bdry_fill_comm->stop();

    const double reflux_dt = sync_time - coarse_sim_time;

    for (hier::PatchLevel<NDIM>::Iterator ip(coarse_level); ip; ip++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = coarse_level->getPatch(ip());

        patch->allocatePatchData(d_temp_var_scratch_data, coarse_sim_time);

        bool at_syncronization = true;
        d_patch_strategy->conservativeDifferenceOnPatch(*patch,
                                                        coarse_sim_time,
                                                        reflux_dt,
                                                        at_syncronization);
        patch->deallocatePatchData(d_temp_var_scratch_data);
    }

    d_patch_strategy->clearDataContext();

    copyTimeDependentData(coarse_level, d_scratch, d_new);

    coarse_level->deallocatePatchData(d_saved_var_scratch_data);

    /*
     * Coarsen time-dependent data from fine patch interiors to coarse patches.
     */

    t_coarsen_sync_create->start();
    sched = d_coarsen_sync_data->createSchedule(coarse_level,
                                                fine_level,
                                                d_patch_strategy);
    t_coarsen_sync_create->stop();

    d_patch_strategy->setDataContext(d_new);

    t_coarsen_sync_comm->start();
    sched->coarsenData();
    t_coarsen_sync_comm->stop();

    d_patch_strategy->clearDataContext();
    return;
}// synchronizeLevelWithCoarser

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::ParabolicLevelIntegrator>;

//////////////////////////////////////////////////////////////////////////////
