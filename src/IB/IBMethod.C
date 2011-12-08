// Filename: IBMethod.C
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

#include "IBMethod.h"

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
#include <ibamr/IBAnchorPointSpec.h>
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IndexUtilities.h>
#include <ibtk/LSiloDataWriter.h>
#if (NDIM == 3)
#include <ibtk/LM3DDataWriter.h>
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
cos_delta(
    const double x,
    const double eps)
{
    if (std::abs(x) > eps)
    {
        return 0.0;
    }
    else
    {
        return 0.5*(1.0+cos(M_PI*x/eps))/eps;
    }
}// cos_delta

// Version of IBMethod restart file data.
static const int IB_METHOD_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMethod::IBMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    // Ensure all pointers to helper objects are NULL.
    d_l_initializer = NULL;
    d_ib_force_fcn = NULL;
    d_ib_force_fcn_needs_init = true;
    d_ib_source_fcn = NULL;
    d_ib_source_fcn_needs_init = true;
    d_normalize_source_strength = false;
    d_post_processor = NULL;
    d_post_processor_needs_init = true;
    d_silo_writer = NULL;
#if (NDIM == 3)
    d_m3D_writer = NULL;
#endif

    // Set some default values.
    d_interp_delta_fcn = "IB_4";
    d_spread_delta_fcn = "IB_4";
    d_do_log = false;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Check the choices for the delta function.
    if (d_interp_delta_fcn != d_spread_delta_fcn)
    {
        pout << "WARNING: different delta functions are being used for velocity interpolation and force spreading.\n"
             << "         recommended usage is to employ the same delta functions for both interpolation and spreading.\n";
    }

    // Get the Lagrangian Data Manager.
    d_l_data_manager = LDataManager::getManager(d_object_name+"::LDataManager", d_interp_delta_fcn, d_spread_delta_fcn, d_ghosts, d_registered_for_restart);
    d_ghosts = d_l_data_manager->getGhostCellWidth();

    // Create the instrument panel object.
    d_instrument_panel = new IBInstrumentPanel(d_object_name+"::IBInstrumentPanel", (input_db->isDatabase("IBInstrumentPanel") ? input_db->getDatabase("IBInstrumentPanel") : Pointer<Database>(NULL)));

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time     = std::numeric_limits<double>::quiet_NaN();
    d_half_time    = std::numeric_limits<double>::quiet_NaN();

    // Indicate all Lagrangian data needs ghost values to be refilled, and that
    // all intermediate data needs to be initialized.
    d_X_current_needs_ghost_fill = true;
    d_X_new_needs_ghost_fill     = true;
    d_X_half_needs_ghost_fill    = true;
    d_F_current_needs_ghost_fill = true;
    d_F_new_needs_ghost_fill     = true;
    d_F_half_needs_ghost_fill    = true;
    d_X_half_needs_reinit        = true;
    d_U_half_needs_reinit        = true;
    return;
}// IBMethod

IBMethod::~IBMethod()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
}// ~IBMethod

void
IBMethod::registerIBLagrangianForceFunction(
    Pointer<IBLagrangianForceStrategy> ib_force_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!ib_force_fcn.isNull());
#endif
    d_ib_force_fcn = ib_force_fcn;
    return;
}// registerIBLagrangianForceFunction

void
IBMethod::registerIBLagrangianSourceFunction(
    Pointer<IBLagrangianSourceStrategy> ib_source_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!ib_source_fcn.isNull());
#endif
    d_ib_source_fcn = ib_source_fcn;
    return;
}// registerIBLagrangianSourceFunction

void
IBMethod::registerLInitStrategy(
    Pointer<LInitStrategy> l_initializer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!l_initializer.isNull());
#endif
    d_l_initializer = l_initializer;
    d_l_data_manager->registerLInitStrategy(d_l_initializer);
    return;
}// registerLInitStrategy

void
IBMethod::freeLInitStrategy()
{
    d_l_initializer.setNull();
    d_l_data_manager->freeLInitStrategy();
    return;
}// freeLInitStrategy

LDataManager*
IBMethod::getLDataManager() const
{
    return d_l_data_manager;
}// getLDataManager

Pointer<IBInstrumentPanel>
IBMethod::getIBInstrumentPanel() const
{
    return d_instrument_panel;
}// getIBInstrumentPanel

void
IBMethod::registerLSiloDataWriter(
    Pointer<LSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
#endif
    d_silo_writer = silo_writer;
    d_l_data_manager->registerLSiloDataWriter(d_silo_writer);
    return;
}// registerLSiloDataWriter

#if (NDIM == 3)
void
IBMethod::registerLM3DDataWriter(
    Pointer<LM3DDataWriter> m3D_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!m3D_writer.isNull());
#endif
    d_m3D_writer = m3D_writer;
    d_l_data_manager->registerLM3DDataWriter(d_m3D_writer);
    return;
}// registerLM3DDataWriter
#endif

const IntVector<NDIM>&
IBMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
}// getMinimumGhostCellWidth

void
IBMethod::preprocessIntegrateData(
    double current_time,
    double new_time,
    int /*num_cycles*/)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time+0.5*(new_time-current_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double start_time = d_ib_solver->d_start_time;

    if (!d_ib_force_fcn.isNull())
    {
        if (d_ib_force_fcn_needs_init)
        {
            const bool initial_time = MathUtilities<double>::equalEps(current_time, start_time);
            resetLagrangianForceFunction(current_time, initial_time);
            d_ib_force_fcn_needs_init = false;
        }
        d_ib_force_fcn->setTimeInterval(current_time, new_time);
    }
    if (!d_ib_source_fcn.isNull())
    {
        if (d_ib_source_fcn_needs_init)
        {
            const bool initial_time = MathUtilities<double>::equalEps(current_time, start_time);
            resetLagrangianSourceFunction(current_time, initial_time);
            d_ib_source_fcn_needs_init = false;
        }
        d_ib_source_fcn->setTimeInterval(current_time, new_time);
    }

    // Look-up or allocate Lagangian data.
    d_X_current_data.resize(finest_ln+1);
    d_X_new_data    .resize(finest_ln+1);
    d_X_half_data   .resize(finest_ln+1);
    d_U_current_data.resize(finest_ln+1);
    d_U_new_data    .resize(finest_ln+1);
    d_U_half_data   .resize(finest_ln+1);
    d_F_current_data.resize(finest_ln+1);
    d_F_new_data    .resize(finest_ln+1);
    d_F_half_data   .resize(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_X_current_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
        d_X_new_data    [ln] = d_l_data_manager->createLData("X_new" ,ln,NDIM);
        d_X_half_data   [ln] = d_l_data_manager->createLData("X_half",ln,NDIM);
        d_U_current_data[ln] = d_l_data_manager->getLData(LDataManager:: VEL_DATA_NAME,ln);
        d_U_new_data    [ln] = d_l_data_manager->createLData("U_new" ,ln,NDIM);
        d_U_half_data   [ln] = d_l_data_manager->createLData("U_half",ln,NDIM);

        int ierr;

        // Initialize X^{n+1/2} and X^{n+1} to equal X^{n}, and initialize
        // U^{n+1/2} and U^{n+1} to equal U^{n}.
        ierr = VecCopy(d_X_current_data[ln]->getVec(), d_X_half_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_X_current_data[ln]->getVec(), d_X_new_data [ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_U_current_data[ln]->getVec(), d_U_half_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_U_current_data[ln]->getVec(), d_U_new_data [ln]->getVec());  IBTK_CHKERRQ(ierr);
    }

    // Indicate all updated Lagrangian data need ghost values to be refilled,
    // and that all intermediate data need to be reinitialized.
    d_X_new_needs_ghost_fill  = true;
    d_X_half_needs_ghost_fill = true;
    d_X_half_needs_reinit     = false;
    d_U_half_needs_reinit     = false;
    return;
}// preprocessIntegrateData

void
IBMethod::postprocessIntegrateData(
    double current_time,
    double new_time,
    int /*num_cycles*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    const int integrator_step = d_ib_solver->getIntegratorStep();

    // Update the instrumentation data.
    updateIBInstrumentationData(integrator_step+1,new_time);
    if (d_instrument_panel->isInstrumented())
    {
        const std::vector<std::string>& instrument_name = d_instrument_panel->getInstrumentNames();
        const std::vector<double>& flow_data = d_instrument_panel->getFlowValues();
        for (unsigned int m = 0; m < flow_data.size(); ++m)
        {
            // NOTE: Flow volume is calculated in default units.
            d_total_flow_volume[m] += flow_data[m]*dt;
            if (d_do_log) plog << "flow volume through " << instrument_name[m] << ":\t " << d_total_flow_volume[m] << "\n";
        }
    }

    // Reset time-dependent Lagrangian data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        int ierr;
        ierr = VecSwap(d_X_current_data[ln]->getVec(), d_X_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecSwap(d_U_current_data[ln]->getVec(), d_U_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }
    d_X_current_needs_ghost_fill = true;

    // Deallocate Lagrangian scratch data.
    d_X_current_data.clear();
    d_X_new_data    .clear();
    d_X_half_data   .clear();
    d_U_current_data.clear();
    d_U_new_data    .clear();
    d_U_half_data   .clear();
    d_F_current_data.clear();
    d_F_new_data    .clear();
    d_F_half_data   .clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time     = std::numeric_limits<double>::quiet_NaN();
    d_half_time    = std::numeric_limits<double>::quiet_NaN();
    return;
}// postprocessIntegrateData

void
IBMethod::interpolateVelocity(
    const int u_data_idx,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
    const double data_time)
{
    std::vector<Pointer<LData> >* X_data = NULL;
    std::vector<Pointer<LData> >* U_data = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        X_data = &d_X_current_data;
        U_data = &d_U_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        if (d_X_half_needs_reinit)
        {
            reinitMidpointData(d_X_current_data, d_X_new_data, d_X_half_data);
            d_X_half_needs_reinit = false;
            d_X_half_needs_ghost_fill = true;
        }
        X_data = &d_X_half_data;
        U_data = &d_U_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        X_data = &d_X_new_data;
        U_data = &d_U_new_data;
    }
    d_l_data_manager->interp(u_data_idx, *U_data, *X_data, u_synch_scheds, u_ghost_fill_scheds, data_time);
    resetAnchorPointValues(*U_data, /*coarsest_ln*/ 0, /*finest_ln*/ d_hierarchy->getFinestLevelNumber());
    d_U_half_needs_reinit = !MathUtilities<double>::equalEps(data_time, d_half_time);
    return;
}// interpolateVelocity

void
IBMethod::eulerStep(
    const double current_time,
    const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;

    // Update the value of X^{n+1} using forward Euler.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        int ierr = VecWAXPY(d_X_new_data[ln]->getVec(), dt, d_U_current_data[ln]->getVec(), d_X_current_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }
    d_X_new_needs_ghost_fill  = true;
    d_X_half_needs_ghost_fill = true;
    d_X_half_needs_reinit     = true;
    return;
}// eulerStep

void
IBMethod::midpointStep(
    const double current_time,
    const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;

    // Update the value of X^{n+1} using the midpoint rule.
    if (d_U_half_needs_reinit)
    {
        reinitMidpointData(d_U_current_data, d_U_new_data, d_U_half_data);
        d_U_half_needs_reinit = false;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        int ierr = VecWAXPY(d_X_new_data[ln]->getVec(), dt, d_U_half_data[ln]->getVec(), d_X_current_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }
    d_X_new_needs_ghost_fill  = true;
    d_X_half_needs_ghost_fill = true;
    d_X_half_needs_reinit     = true;
    return;
}// midpointStep

bool
IBMethod::hasFluidSources() const
{
    return !d_ib_source_fcn.isNull();
}// hasFluidSources

void
IBMethod::computeLagrangianForce(
    const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    int ierr;
    std::vector<Pointer<LData> >* F_data = NULL;
    std::vector<Pointer<LData> >* X_data = NULL;
    std::vector<Pointer<LData> >* U_data = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            if (d_F_current_data[ln].isNull()) d_F_current_data[ln] = d_l_data_manager->createLData("F",ln,NDIM);
        }
        d_F_current_needs_ghost_fill = true;
        F_data = &d_F_current_data;
        X_data = &d_X_current_data;
        U_data = &d_U_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        if (d_X_half_needs_reinit)
        {
            reinitMidpointData(d_X_current_data, d_X_new_data, d_X_half_data);
            d_X_half_needs_reinit = false;
            d_X_half_needs_ghost_fill = true;
        }
        if (d_U_half_needs_reinit)
        {
            reinitMidpointData(d_U_current_data, d_U_new_data, d_U_half_data);
            d_U_half_needs_reinit = false;
        }
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            if (d_F_half_data[ln].isNull()) d_F_half_data[ln] = d_l_data_manager->createLData("F_half",ln,NDIM);
        }
        d_F_half_needs_ghost_fill = true;
        F_data = &d_F_half_data;
        X_data = &d_X_half_data;
        U_data = &d_U_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            if (d_F_new_data[ln].isNull()) d_F_new_data[ln] = d_l_data_manager->createLData("F_new",ln,NDIM);
        }
        d_F_new_needs_ghost_fill = true;
        F_data = &d_F_new_data;
        X_data = &d_X_new_data;
        U_data = &d_U_new_data;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecSet((*F_data)[ln]->getVec(), 0.0);  IBTK_CHKERRQ(ierr);
        if (d_ib_force_fcn.isNull()) continue;
        d_ib_force_fcn->computeLagrangianForce((*F_data)[ln], (*X_data)[ln], (*U_data)[ln], d_hierarchy, ln, data_time, d_l_data_manager);
    }
    resetAnchorPointValues(*F_data, coarsest_ln, finest_ln);
    return;
}// computeLagrangianForce

void
IBMethod::spreadForce(
    const int f_data_idx,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& f_prolongation_scheds,
    const double data_time)
{
    std::vector<Pointer<LData> >* F_data = NULL;
    std::vector<Pointer<LData> >* X_data = NULL;
    bool* F_needs_ghost_fill = NULL;
    bool* X_needs_ghost_fill = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        F_data = &d_F_current_data;
        X_data = &d_X_current_data;
        F_needs_ghost_fill = &d_F_current_needs_ghost_fill;
        X_needs_ghost_fill = &d_X_current_needs_ghost_fill;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        if (d_X_half_needs_reinit)
        {
            reinitMidpointData(d_X_current_data, d_X_new_data, d_X_half_data);
            d_X_half_needs_reinit = false;
            d_X_half_needs_ghost_fill = true;
        }
        F_data = &d_F_half_data;
        X_data = &d_X_half_data;
        F_needs_ghost_fill = &d_F_half_needs_ghost_fill;
        X_needs_ghost_fill = &d_X_half_needs_ghost_fill;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        F_data = &d_F_new_data;
        X_data = &d_X_new_data;
        F_needs_ghost_fill = &d_F_new_needs_ghost_fill;
        X_needs_ghost_fill = &d_X_new_needs_ghost_fill;
    }
    d_l_data_manager->spread(f_data_idx, *F_data, *X_data, f_prolongation_scheds, *F_needs_ghost_fill, *X_needs_ghost_fill);
    *F_needs_ghost_fill = false;
    *X_needs_ghost_fill = false;
    return;
}// spreadForce

void
IBMethod::computeLagrangianFluidSource(
    const double data_time)
{
    if (d_ib_source_fcn.isNull()) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Reset the values of Q_src.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::fill(d_Q_src[ln].begin(), d_Q_src[ln].end(), 0.0);
    }

    // Compute the source strengths for each of the sources/sinks.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_n_src[ln] == 0) continue;
        d_ib_source_fcn->computeSourceStrengths(d_Q_src[ln], d_hierarchy, ln, data_time, d_l_data_manager);
    }
    return;
}// computeLagrangianFluidSource

void
IBMethod::spreadFluidSource(
    const int q_data_idx,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*q_prolongation_scheds*/,
    const double data_time)
{
    if (d_ib_source_fcn.isNull()) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Get the present source locations.
    std::vector<Pointer<LData> >* X_data = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        X_data = &d_X_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        if (d_X_half_needs_reinit)
        {
            reinitMidpointData(d_X_current_data, d_X_new_data, d_X_half_data);
            d_X_half_needs_reinit = false;
            d_X_half_needs_ghost_fill = true;
        }
        X_data = &d_X_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        X_data = &d_X_new_data;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_n_src[ln] == 0) continue;
        d_ib_source_fcn->getSourceLocations(d_X_src[ln], d_r_src[ln], (*X_data)[ln], d_hierarchy, ln, data_time, d_l_data_manager);
    }

    // Spread the sources/sinks onto the Cartesian grid.
    d_ib_solver->d_hier_pressure_cc_data_ops->setToScalar(q_data_idx, 0.0);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_n_src[ln] == 0) continue;
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(ln == d_hierarchy->getFinestLevelNumber());
#endif
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Index<NDIM>& patch_lower = patch_box.lower();
            const Index<NDIM>& patch_upper = patch_box.upper();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const xLower = pgeom->getXLower();
            const double* const xUpper = pgeom->getXUpper();
            const double* const dx = pgeom->getDx();
            const Pointer<CellData<NDIM,double> > q_data = patch->getPatchData(q_data_idx);
            for (int n = 0; n < d_n_src[ln]; ++n)
            {
                // The source radius must be an integer multiple of the grid
                // spacing.
                blitz::TinyVector<double,NDIM> r;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    r[d] = std::max(std::floor(d_r_src[ln][n]/dx[d]+0.5),2.0)*dx[d];
                }

                // Determine the approximate source stencil box.
                const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[ln][n], xLower, xUpper, dx, patch_lower, patch_upper);
                Box<NDIM> stencil_box(i_center,i_center);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    stencil_box.grow(d, static_cast<int>(r[d]/dx[d])+1);
                }

                // Spread the source strength onto the Cartesian grid.
                for (Box<NDIM>::Iterator b(patch_box*stencil_box); b; b++)
                {
                    const Index<NDIM>& i = b();
                    double wgt = 1.0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const double X_center = xLower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
                        wgt *= cos_delta(X_center - d_X_src[ln][n][d], r[d]);
                    }
                    (*q_data)(i) += d_Q_src[ln][n]*wgt;
                }
            }
        }
    }

    // Compute the net inflow into the computational domain.
    const int wgt_idx = d_ib_solver->d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    PatchCellDataOpsReal<NDIM,double> patch_cc_data_ops;
    double Q_sum = 0.0;
    double Q_max = 0.0;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Q_sum = std::accumulate(d_Q_src[ln].begin(), d_Q_src[ln].end(), Q_sum);
        for (unsigned int k = 0; k < d_Q_src[ln].size(); ++k)
        {
            Q_max = std::max(Q_max,std::abs(d_Q_src[ln][k]));
        }
    }
    const double q_total = d_ib_solver->d_hier_pressure_cc_data_ops->integral(q_data_idx, wgt_idx);
    if (std::abs(q_total-Q_sum)                     > 1.0e-12 &&
        std::abs(q_total-Q_sum)/std::max(Q_max,1.0) > 1.0e-12)
    {
#if (NDIM == 2)
        TBOX_ERROR(d_object_name << "::spreadFluidSource():\n"
                   << "  Lagrangian and Eulerian source/sink strengths are inconsistent:\n"
                   << "    Sum_{i,j} q_{i,j} h^2     = " << q_total << "\n"
                   << "    Sum_{l=1,...,n_src} Q_{l} = " << Q_sum   << "\n");
#endif
#if (NDIM == 3)
        TBOX_ERROR(d_object_name << "::spreadFluidSource():\n"
                   << "  Lagrangian and Eulerian source/sink strengths are inconsistent:\n"
                   << "    Sum_{i,j,k} q_{i,j,k} h^3 = " << q_total << "\n"
                   << "    Sum_{l=1,...,n_src} Q_{l} = " << Q_sum   << "\n");
#endif
    }

    // Balance the net inflow/outflow with outflow/inflow along the upper/lower
    // boundaries of the computational domain (if needed).
    if (d_normalize_source_strength)
    {
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
        const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];
        const double* const dx_coarsest = grid_geom->getDx();
        Box<NDIM> interior_box = domain_box;
        for (unsigned int d = 0; d < NDIM-1; ++d)
        {
            interior_box.grow(d,-1);
        }
        BoxList<NDIM> bdry_boxes;
        bdry_boxes.removeIntersections(domain_box,interior_box);
        double vol = static_cast<double>(bdry_boxes.getTotalSizeOfBoxes());
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            vol *= dx_coarsest[d];
        }
        const double q_norm = -q_total/vol;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            BoxList<NDIM> level_bdry_boxes(bdry_boxes);
            level_bdry_boxes.refine(level->getRatio());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM,double> > q_data = patch->getPatchData(q_data_idx);
                for (BoxList<NDIM>::Iterator blist(level_bdry_boxes); blist; blist++)
                {
                    for (Box<NDIM>::Iterator b(blist()*patch_box); b; b++)
                    {
                        (*q_data)(b()) += q_norm;
                    }
                }
            }
        }
        const double integral_q = d_ib_solver->d_hier_pressure_cc_data_ops->integral(q_data_idx, wgt_idx);
        if (std::abs(integral_q) > 1.0e-10*std::max(1.0,d_ib_solver->d_hier_pressure_cc_data_ops->maxNorm(q_data_idx, wgt_idx)))
        {
            TBOX_ERROR(d_object_name << "::spreadFluidSource():\n"
                       << "  ``external'' source/sink does not correctly offset net inflow/outflow into domain.\n"
                       << "  integral{q} = " << integral_q << " != 0.\n");
        }
    }
    return;
}// spreadFluidSource

void
IBMethod::interpolatePressure(
    int p_data_idx,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*p_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*p_ghost_fill_scheds*/,
    const double data_time)
{
    if (d_ib_source_fcn.isNull()) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Get the present source locations.
    std::vector<Pointer<LData> >* X_data = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        X_data = &d_X_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        if (d_X_half_needs_reinit)
        {
            reinitMidpointData(d_X_current_data, d_X_new_data, d_X_half_data);
            d_X_half_needs_reinit = false;
            d_X_half_needs_ghost_fill = true;
        }
        X_data = &d_X_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        X_data = &d_X_new_data;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_n_src[ln] == 0) continue;
        d_ib_source_fcn->getSourceLocations(d_X_src[ln], d_r_src[ln], (*X_data)[ln], d_hierarchy, ln, data_time, d_l_data_manager);
    }

    // Compute the normalization pressure (if needed).
    double p_norm = 0.0;
    if (d_normalize_source_strength)
    {
        const int wgt_idx = d_ib_solver->d_hier_math_ops->getCellWeightPatchDescriptorIndex();
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
        const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];
        Box<NDIM> interior_box = domain_box;
        for (unsigned int d = 0; d < NDIM-1; ++d)
        {
            interior_box.grow(d,-1);
        }
        BoxList<NDIM> bdry_boxes;
        bdry_boxes.removeIntersections(domain_box,interior_box);
        double vol = 0.0;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            BoxList<NDIM> level_bdry_boxes(bdry_boxes);
            level_bdry_boxes.refine(level->getRatio());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM,double> > p_data = patch->getPatchData(p_data_idx);
                const Pointer<CellData<NDIM,double> > wgt_data = patch->getPatchData(wgt_idx);
                for (BoxList<NDIM>::Iterator blist(level_bdry_boxes); blist; blist++)
                {
                    for (Box<NDIM>::Iterator b(blist()*patch_box); b; b++)
                    {
                        const Index<NDIM>& i = b();
                        p_norm += (*p_data)(i)*(*wgt_data)(i);
                        vol += (*wgt_data)(i);
                    }
                }
            }
        }
        SAMRAI_MPI::sumReduction(&p_norm,1);
        SAMRAI_MPI::sumReduction(&vol,1);
        p_norm /= vol;
    }

    // Reset the values of P_src.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        std::fill(d_P_src[ln].begin(), d_P_src[ln].end(), 0.0);
    }

    // Compute the mean pressure at the sources/sinks associated with each level
    // of the Cartesian grid.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_n_src[ln] == 0) continue;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Index<NDIM>& patch_lower = patch_box.lower();
            const Index<NDIM>& patch_upper = patch_box.upper();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const xLower = pgeom->getXLower();
            const double* const xUpper = pgeom->getXUpper();
            const double* const dx = pgeom->getDx();
            const Pointer<CellData<NDIM,double> > p_data = patch->getPatchData(p_data_idx);
            for (int n = 0; n < d_n_src[ln]; ++n)
            {
                // The source radius must be an integer multiple of the grid
                // spacing.
                blitz::TinyVector<double,NDIM> r;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    r[d] = std::max(std::floor(d_r_src[ln][n]/dx[d]+0.5),2.0)*dx[d];
                }

                // Determine the approximate source stencil box.
                const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[ln][n], xLower, xUpper, dx, patch_lower, patch_upper);
                Box<NDIM> stencil_box(i_center,i_center);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    stencil_box.grow(d, static_cast<int>(r[d]/dx[d])+1);
                }

                // Interpolate the pressure from the Cartesian grid.
                for (Box<NDIM>::Iterator b(patch_box*stencil_box); b; b++)
                {
                    const Index<NDIM>& i = b();
                    double wgt = 1.0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const double X_center = xLower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
                        wgt *= cos_delta(X_center - d_X_src[ln][n][d], r[d])*dx[d];
                    }
                    d_P_src[ln][n] += (*p_data)(i)*wgt;
                }
            }
        }
        SAMRAI_MPI::sumReduction(&d_P_src[ln][0], d_P_src[ln].size());
        std::transform(d_P_src[ln].begin(), d_P_src[ln].end(), d_P_src[ln].begin(),
                       std::bind2nd(std::plus<double>(),-p_norm));

        // Update the pressures stored by the Lagrangian source strategy.
        d_ib_source_fcn->setSourcePressures(d_P_src[ln], d_hierarchy, ln, data_time, d_l_data_manager);
    }
    return;
}// interpolatePressure

void
IBMethod::postprocessData(
    Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    if (d_post_processor.isNull()) return;

    INSHierarchyIntegrator* ins_hier_integrator = d_ib_solver->d_ins_hier_integrator;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(ins_hier_integrator-> getVelocityVariable(), ins_hier_integrator->getCurrentContext());
    const int p_current_idx = var_db->mapVariableAndContextToIndex(ins_hier_integrator-> getPressureVariable(), ins_hier_integrator->getCurrentContext());
    const int f_current_idx = var_db->mapVariableAndContextToIndex(ins_hier_integrator->getBodyForceVariable(), ins_hier_integrator->getCurrentContext());

    const double current_time = d_ib_solver->d_integrator_time;
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();

    // Initialize data on each level of the patch hierarchy.
    std::vector<Pointer<LData> > X_data(finest_ln+1);
    std::vector<Pointer<LData> > F_data(finest_ln+1);
    std::vector<Pointer<LData> > U_data(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        X_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
        U_data[ln] = d_l_data_manager->getLData(LDataManager:: VEL_DATA_NAME,ln);
        F_data[ln] = d_l_data_manager->createLData("F",ln,NDIM);
    }

    // Compute F(n) = F(X(n),U(n),n), the Lagrangian force corresponding to
    // configuration X(n) and velocity U(n) at time t_{n}.
    if (d_ib_force_fcn_needs_init)
    {
        const bool initial_time = MathUtilities<double>::equalEps(current_time,d_ib_solver->d_start_time);
        resetLagrangianForceFunction(current_time, initial_time);
        d_ib_force_fcn_needs_init = false;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        int ierr = VecSet(F_data[ln]->getVec(), 0.0);  IBTK_CHKERRQ(ierr);
        if (!d_ib_force_fcn.isNull())
        {
            d_ib_force_fcn->computeLagrangianForce(F_data[ln], X_data[ln], U_data[ln], hierarchy, ln, current_time, d_l_data_manager);
        }
    }
    resetAnchorPointValues(F_data, coarsest_ln, finest_ln);

    // Perform the user-defined post-processing.
    d_post_processor->postprocessData(u_current_idx, p_current_idx, f_current_idx, F_data, X_data, U_data, hierarchy, coarsest_ln, finest_ln, current_time, d_l_data_manager);
    return;
}// postprocessData

void
IBMethod::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
    int u_data_idx,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
    int integrator_step,
    double init_data_time,
    bool initial_time)
{
    // Cache pointers to the patch hierarchy and gridding algorithm.
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Lookup the range of hierarchy levels.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Initialize various Lagrangian data objects.
    if (initial_time)
    {
        // Initialize the interpolated velocity field.
        std::vector<Pointer<LData> > X_data(finest_ln+1);
        std::vector<Pointer<LData> > U_data(finest_ln+1);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            X_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
            U_data[ln] = d_l_data_manager->getLData(LDataManager:: VEL_DATA_NAME,ln);
        }
        d_l_data_manager->interp(u_data_idx, U_data, X_data, u_synch_scheds, u_ghost_fill_scheds, init_data_time);
        resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

        // Initialize source/sink data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_ib_source_fcn.isNull())
            {
                d_ib_source_fcn->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
                d_n_src[ln] = d_ib_source_fcn->getNumSources(d_hierarchy, ln, init_data_time, d_l_data_manager);
                d_X_src[ln].resize(d_n_src[ln], blitz::TinyVector<double,NDIM>(std::numeric_limits<double>::quiet_NaN()));
                d_r_src[ln].resize(d_n_src[ln], std::numeric_limits<double>::quiet_NaN());
                d_P_src[ln].resize(d_n_src[ln], std::numeric_limits<double>::quiet_NaN());
                d_Q_src[ln].resize(d_n_src[ln], std::numeric_limits<double>::quiet_NaN());
                d_ib_source_fcn->getSourceLocations(d_X_src[ln], d_r_src[ln], X_data[ln], d_hierarchy, ln, init_data_time, d_l_data_manager);
            }
        }
    }

    // Initialize the instrumentation data.
    d_instrument_panel->initializeHierarchyIndependentData(d_hierarchy, d_l_data_manager);
    if (d_instrument_panel->isInstrumented())
    {
        d_instrument_panel->initializeHierarchyDependentData(d_hierarchy, d_l_data_manager, integrator_step, init_data_time);
        if (d_total_flow_volume.empty())
        {
            d_total_flow_volume.resize(d_instrument_panel->getFlowValues().size(),0.0);
        }
    }

    // Indicate that the force and source strategies and the post processor need
    // to be re-initialized.
    d_ib_force_fcn_needs_init   = true;
    d_ib_source_fcn_needs_init  = true;
    d_post_processor_needs_init = true;
    return;
}// initializePatchHierarchy

void
IBMethod::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer,
    int workload_data_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!load_balancer.isNull());
#endif
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;
    d_l_data_manager->registerLoadBalancer(load_balancer, workload_data_idx);
    return;
}// registerLoadBalancer

void
IBMethod::updateWorkloadEstimates(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    int /*workload_data_idx*/)
{
    d_l_data_manager->updateWorkloadEstimates();
    return;
}// updateWorkloadEstimates

void
IBMethod::beginDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    d_l_data_manager->beginDataRedistribution();
    return;
}// beginDataRedistribution

void
IBMethod::endDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    d_l_data_manager->endDataRedistribution();

    // Look up the re-distributed Lagrangian position data.
    std::vector<Pointer<LData> > X_data(hierarchy->getFinestLevelNumber()+1);
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        X_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME,ln);
    }

    // Compute the set of local anchor points.
    static const double eps = 2.0*sqrt(std::numeric_limits<double>::epsilon());
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const grid_xLower = grid_geom->getXLower();
    const double* const grid_xUpper = grid_geom->getXUpper();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        d_anchor_point_local_idxs[ln].clear();
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
             cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const IBAnchorPointSpec* const anchor_point_spec = node_idx->getNodeDataItem<IBAnchorPointSpec>();
            if (anchor_point_spec != NULL)
            {
                d_anchor_point_local_idxs[ln].insert(node_idx->getLocalPETScIndex());
            }
        }

        const blitz::Array<double,2>& X_array = *(X_data[ln]->getLocalFormVecArray());
        for (int i = 0; i < static_cast<int>(X_data[ln]->getLocalNodeCount()); ++i)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                if ((periodic_shift[d] == 0) &&
                    (X_array(i,d) <= grid_xLower[d]+eps || X_array(i,d) >= grid_xUpper[d]-eps))
                {
                    d_anchor_point_local_idxs[ln].insert(i);
                    break;
                }
            }
        }
        X_data[ln]->restoreArrays();
    }

    // Indicate that the force and source strategies and the post processor need
    // to be re-initialized.
    d_ib_force_fcn_needs_init   = true;
    d_ib_source_fcn_needs_init  = true;
    d_post_processor_needs_init = true;
    return;
}// endDataRedistribution

void
IBMethod::initializeLevelData(
    Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    int level_number,
    double init_data_time,
    bool can_be_refined,
    bool initial_time,
    Pointer<BasePatchLevel<NDIM> > old_level,
    bool allocate_data)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->resetLevels(0, finest_hier_level);
    d_l_data_manager->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    if (!d_load_balancer.isNull() && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
        d_l_data_manager->updateWorkloadEstimates(level_number, level_number);
    }
    return;
}// initializeLevelData

void
IBMethod::resetHierarchyConfiguration(
    Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    int coarsest_level,
    int finest_level)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->resetLevels(0, finest_hier_level);
    d_l_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // If we have added or removed a level, resize the anchor point vectors.
    d_anchor_point_local_idxs.clear();
    d_anchor_point_local_idxs.resize(finest_hier_level+1);

    // If we have added or removed a level, resize the source/sink data vectors.
    d_X_src.resize(finest_hier_level+1);
    d_r_src.resize(finest_hier_level+1);
    d_P_src.resize(finest_hier_level+1);
    d_Q_src.resize(finest_hier_level+1);
    d_n_src.resize(finest_hier_level+1,0);
    return;
}// resetHierarchyConfiguration

void
IBMethod::applyGradientDetector(
    Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    int level_number,
    double error_data_time,
    int tag_index,
    bool initial_time,
    bool uses_richardson_extrapolation_too)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Tag cells that contain Lagrangian nodes.
    d_l_data_manager->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

    // Tag cells where the Cartesian source/sink strength is nonzero.
    if (!d_ib_source_fcn.isNull() && !initial_time && hierarchy->finerLevelExists(level_number))
    {
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
        if (!grid_geom->getDomainIsSingleBox()) TBOX_ERROR("physical domain must be a single box...\n");

        const Index<NDIM>& lower = grid_geom->getPhysicalDomain()[0].lower();
        const Index<NDIM>& upper = grid_geom->getPhysicalDomain()[0].upper();
        const double* const xLower = grid_geom->getXLower();
        const double* const xUpper = grid_geom->getXUpper();
        const double* const dx = grid_geom->getDx();

        const int finer_level_number = level_number+1;
        Pointer<PatchLevel<NDIM> > finer_level = hierarchy->getPatchLevel(finer_level_number);
        for (int n = 0; n < d_n_src[finer_level_number]; ++n)
        {
            blitz::TinyVector<double,NDIM> dx_finer;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dx_finer[d] = dx[d]/static_cast<double>(finer_level->getRatio()(d));
            }

            // The source radius must be an integer multiple of the grid
            // spacing.
            blitz::TinyVector<double,NDIM> r;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                r[d] = std::max(std::floor(d_r_src[finer_level_number][n]/dx_finer[d]+0.5),2.0)*dx_finer[d];
            }

            // Determine the approximate source stencil box.
            const Index<NDIM> i_center = IndexUtilities::getCellIndex(d_X_src[finer_level_number][n], xLower, xUpper, dx_finer.data(), lower, upper);
            Box<NDIM> stencil_box(i_center,i_center);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                stencil_box.grow(d, static_cast<int>(r[d]/dx_finer[d])+1);
            }
            const Box<NDIM> coarsened_stencil_box = Box<NDIM>::coarsen(stencil_box, finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
                tags_data->fillAll(1, coarsened_stencil_box);
            }
        }
    }
    return;
}// applyGradientDetector

void
IBMethod::putToDatabase(
    Pointer<Database> db)
{
    db->putInteger("IB_METHOD_VERSION", IB_METHOD_VERSION);
    db->putString("d_interp_delta_fcn", d_interp_delta_fcn);
    db->putString("d_spread_delta_fcn", d_spread_delta_fcn);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    const std::vector<std::string>& instrument_names = IBInstrumentationSpec::getInstrumentNames();
    if (!instrument_names.empty())
    {
        db->putInteger("instrument_names_sz", instrument_names.size());
        db->putStringArray("instrument_names", &instrument_names[0], instrument_names.size());
    }
    db->putInteger("d_total_flow_volume_sz", d_total_flow_volume.size());
    if (!d_total_flow_volume.empty())
    {
        db->putDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    db->putInteger("finest_hier_level", finest_hier_level);
    db->putIntegerArray("d_n_src", &d_n_src[0], finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            db->putDoubleArray("d_X_src_"+id_string, &d_X_src[ln][n][0], NDIM);
            db->putDouble("d_r_src_"+id_string, d_r_src[ln][n]);
            db->putDouble("d_P_src_"+id_string, d_P_src[ln][n]);
            db->putDouble("d_Q_src_"+id_string, d_Q_src[ln][n]);
        }
    }
    db->putBool("d_normalize_source_strength", d_normalize_source_strength);
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBMethod::reinitMidpointData(
    const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& current_data,
    const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& new_data,
    const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& half_data)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    int ierr;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecAXPBYPCZ(half_data[ln]->getVec(), 0.5, 0.5, 0.0, current_data[ln]->getVec(), new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }
    return;
}// reinitMidpointData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBMethod::resetLagrangianForceFunction(
    const double init_data_time,
    const bool initial_time)
{
    if (d_ib_force_fcn.isNull()) return;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_ib_force_fcn->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
    }
    return;
}// resetLagrangianForceFunction

void
IBMethod::resetLagrangianSourceFunction(
    const double init_data_time,
    const bool initial_time)
{
    if (d_ib_source_fcn.isNull()) return;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_ib_source_fcn->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
    }
    return;
}// resetLagrangianSourceFunction

void
IBMethod::resetPostProcessor(
    const double init_data_time,
    const bool initial_time)
{
    if (d_post_processor.isNull()) return;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_post_processor->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
    }
    return;
}// resetPostProcessor

void
IBMethod::updateIBInstrumentationData(
    const int timestep_num,
    const double data_time)
{
    if (!d_instrument_panel->isInstrumented()) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Compute the positions of the flow meter nets.
    d_instrument_panel->initializeHierarchyDependentData(d_hierarchy, d_l_data_manager, timestep_num, data_time);

    // Compute the flow rates and pressures.
    INSHierarchyIntegrator* ins_hier_integrator = d_ib_solver->d_ins_hier_integrator;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(ins_hier_integrator->getVelocityVariable(), ins_hier_integrator->getScratchContext());
    const int p_scratch_idx = var_db->mapVariableAndContextToIndex(ins_hier_integrator->getPressureVariable(), ins_hier_integrator->getScratchContext());

    std::vector<bool> deallocate_u_scratch_data(finest_ln+1,false);
    std::vector<bool> deallocate_p_scratch_data(finest_ln+1,false);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(u_scratch_idx))
        {
            deallocate_u_scratch_data[ln] = true;
            level->allocatePatchData(u_scratch_idx, data_time);
        }
        if (!level->checkAllocated(p_scratch_idx))
        {
            deallocate_p_scratch_data[ln] = true;
            level->allocatePatchData(p_scratch_idx, data_time);
        }
        d_ib_solver->getGhostfillRefineSchedules(d_ib_solver->getName()+"::INSTRUMENTATION_DATA_FILL")[ln]->fillData(data_time);
    }

    d_instrument_panel->readInstrumentData(u_scratch_idx, p_scratch_idx, d_hierarchy, d_l_data_manager, timestep_num, data_time);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (deallocate_u_scratch_data[ln]) level->deallocatePatchData(u_scratch_idx);
        if (deallocate_p_scratch_data[ln]) level->deallocatePatchData(p_scratch_idx);
    }
    return;
}// updateIBInstrumentationData

void
IBMethod::resetAnchorPointValues(
    std::vector<Pointer<LData> > U_data,
    const int coarsest_ln,
    const int finest_ln)
{
    PetscErrorCode ierr;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        const int depth = U_data[ln]->getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(depth == NDIM);
#endif
        Vec U_vec = U_data[ln]->getVec();
        double* U_arr;
        ierr = VecGetArray(U_vec, &U_arr);  IBTK_CHKERRQ(ierr);
        for (std::set<int>::const_iterator cit = d_anchor_point_local_idxs[ln].begin();
             cit != d_anchor_point_local_idxs[ln].end(); ++cit)
        {
            const int& i = *cit;
            for (int d = 0; d < depth; ++d)
            {
                U_arr[depth*i+d] = 0.0;
            }
        }
        ierr = VecRestoreArray(U_vec, &U_arr);  IBTK_CHKERRQ(ierr);
    }
    return;
}// resetAnchorPointValues

void
IBMethod::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->isString("interp_delta_fcn") && db->isString("spread_delta_fcn"))
        {
            d_interp_delta_fcn = db->getStringWithDefault("interp_delta_fcn", d_interp_delta_fcn);
            d_spread_delta_fcn = db->getStringWithDefault("spread_delta_fcn", d_spread_delta_fcn);
        }
        else
        {
            d_interp_delta_fcn = db->getStringWithDefault("delta_fcn", d_interp_delta_fcn);
            d_spread_delta_fcn = db->getStringWithDefault("delta_fcn", d_spread_delta_fcn);
        }
        if (db->isInteger("min_ghost_cell_width"))
        {
            d_ghosts = db->getInteger("min_ghost_cell_width");
        }
        else if (db->isDouble("min_ghost_cell_width"))
        {
            d_ghosts = static_cast<int>(std::ceil(db->getDouble("min_ghost_cell_width")));
        }
        if (db->isBool("normalize_source_strength")) d_normalize_source_strength = db->getBool("normalize_source_strength");
    }
    if      (db->keyExists("do_log"        )) d_do_log = db->getBool("do_log"        );
    else if (db->keyExists("enable_logging")) d_do_log = db->getBool("enable_logging");
    return;
}// getFromInput

void
IBMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IB_METHOD_VERSION");
    if (ver != IB_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_interp_delta_fcn = db->getString("d_interp_delta_fcn");
    d_spread_delta_fcn = db->getString("d_spread_delta_fcn");
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);
    if (db->keyExists("instrument_names"))
    {
        const int sz = db->getInteger("instrument_names_sz");
        std::vector<std::string> instrument_names(sz);
        db->getStringArray("instrument_names", &instrument_names[0], sz);
        IBInstrumentationSpec::setInstrumentNames(instrument_names);
    }
    const int total_flow_volume_sz = db->getInteger("d_total_flow_volume_sz");
    d_total_flow_volume.resize(total_flow_volume_sz, std::numeric_limits<double>::quiet_NaN());
    if (!d_total_flow_volume.empty())
    {
        db->getDoubleArray("d_total_flow_volume", &d_total_flow_volume[0], d_total_flow_volume.size());
    }
    const int finest_hier_level = db->getInteger("finest_hier_level");
    d_X_src.resize(finest_hier_level+1);
    d_r_src.resize(finest_hier_level+1);
    d_P_src.resize(finest_hier_level+1);
    d_Q_src.resize(finest_hier_level+1);
    d_n_src.resize(finest_hier_level+1,0);
    db->getIntegerArray("d_n_src", &d_n_src[0], finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        d_X_src[ln].resize(d_n_src[ln],blitz::TinyVector<double,NDIM>(std::numeric_limits<double>::quiet_NaN()));
        d_r_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        d_P_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        d_Q_src[ln].resize(d_n_src[ln],std::numeric_limits<double>::quiet_NaN());
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            db->getDoubleArray("d_X_src_"+id_string, &d_X_src[ln][n][0], NDIM);
            d_r_src[ln][n] = db->getDouble("d_r_src_"+id_string);
            d_P_src[ln][n] = db->getDouble("d_P_src_"+id_string);
            d_Q_src[ln][n] = db->getDouble("d_Q_src_"+id_string);
        }
    }
    d_normalize_source_strength = db->getBool("d_normalize_source_strength");
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
