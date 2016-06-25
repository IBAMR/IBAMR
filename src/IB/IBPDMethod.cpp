// Filename: IBPDMethod.cpp
// Created on 08 Apr 2016 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

#include <math.h>
#include <stddef.h>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBPDForceGen.h"
#include "ibamr/IBPDMethod.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LInitStrategy.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBPDMethod restart file data.
static const int IB_PD_METHOD_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBPDMethod::IBPDMethod(const std::string& object_name, Pointer<Database> input_db, bool register_for_restart)
    : IBMethod(object_name, input_db, register_for_restart)
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object
    // name.

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Indicate if initilization is needed for the objects.
    d_ib_pd_force_fcn_needs_init = true;

    return;
} // IBPDMethod

IBPDMethod::~IBPDMethod()
{
    // intentionally blank
    //
    // NOTE: Parent class constructor unregisters class with the restart
    // manager.
    return;
} // ~IBPDMethod

void
IBPDMethod::registerIBPDForceGen(Pointer<IBPDForceGen> ib_pd_force_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ib_pd_force_fcn);
#endif
    d_ib_pd_force_fcn = ib_pd_force_fcn;
    return;
} // registerIBPDForceGen

void
IBPDMethod::setIBPDForceGenNeedsInit()
{
    d_ib_pd_force_fcn_needs_init = true;
    return;
} // setIBPDForceGenNeedsInit()

void
IBPDMethod::registerEulerianVariables()
{
    IBMethod::registerEulerianVariables();
    return;
} // registerEulerianVariables

void
IBPDMethod::registerEulerianCommunicationAlgorithms()
{
    IBMethod::registerEulerianCommunicationAlgorithms();
    return;
} // registerEulerianCommunicationAlgorithms

void
IBPDMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);
    const double start_time = d_ib_solver->getStartTime();
    if (d_ib_pd_force_fcn)
    {
        if (d_ib_pd_force_fcn_needs_init)
        {
            const bool initial_time = MathUtilities<double>::equalEps(current_time, start_time);
            resetLagrangianPDForceFunction(current_time, initial_time);
            d_ib_pd_force_fcn_needs_init = false;
        }
    }

    return;
} // preprocessIntegrateData

void
IBPDMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{

    IBMethod::postprocessIntegrateData(current_time, new_time, num_cycles);
    return;
} // postprocessIntegrateData

void
IBPDMethod::interpolateVelocity(const int u_data_idx,
                                const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                const double data_time)
{
    // Interpolate fluid velocity.
    IBMethod::interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);

    return;
} // interpolateVelocity

void
IBPDMethod::eulerStep(const double current_time, const double new_time)
{
    IBMethod::eulerStep(current_time, new_time);

    return;
} // eulerStep

void
IBPDMethod::midpointStep(const double current_time, const double new_time)
{
    IBMethod::midpointStep(current_time, new_time);
    return;

} // midpointStep

void
IBPDMethod::trapezoidalStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR(d_object_name << "::trapezoidalStep():\n"
                             << "  time-stepping type TRAPEZOIDAL_RULE not supported by class IBPDMethod;\n"
                             << "  use MIDPOINT_RULE instead.\n");
    return;
} // trapezoidalStep

void
IBPDMethod::computeLagrangianForce(const double data_time)
{

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> >* F_data = NULL;
    std::vector<Pointer<LData> >* X_data = NULL;
    std::vector<Pointer<LData> >* U_data = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        d_F_current_needs_ghost_fill = true;
        F_data = &d_F_current_data;
        X_data = &d_X_current_data;
        U_data = &d_U_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        d_F_half_needs_ghost_fill = true;
        F_data = &d_F_half_data;
        X_data = &d_X_half_data;
        U_data = &d_U_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        d_F_new_needs_ghost_fill = true;
        F_data = &d_F_new_data;
        X_data = &d_X_new_data;
        U_data = &d_U_new_data;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        if (d_ib_pd_force_fcn)
        {
            d_ib_pd_force_fcn->computeLagrangianForceAndDamage((*F_data)[ln],
                                                               d_l_data_manager->getLData("damage", ln),
                                                               (*X_data)[ln],
                                                               (*U_data)[ln],
                                                               d_hierarchy,
                                                               ln,
                                                               data_time,
                                                               d_l_data_manager);
        }
    }

    return;
} // computeLagrangianForce

void
IBPDMethod::spreadForce(const int f_data_idx,
                        RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                        const std::vector<Pointer<RefineSchedule<NDIM> > >& f_prolongation_scheds,
                        const double data_time)
{
    IBMethod::spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);

    return;
} // spreadForce

void
IBPDMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                     Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                     int u_data_idx,
                                     const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                     const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                     int integrator_step,
                                     double init_data_time,
                                     bool initial_time)
{
    // Initialize various Lagrangian data objects required by the conventional
    // IB method.
    IBMethod::initializePatchHierarchy(hierarchy,
                                       gridding_alg,
                                       u_data_idx,
                                       u_synch_scheds,
                                       u_ghost_fill_scheds,
                                       integrator_step,
                                       init_data_time,
                                       initial_time);

    // Indicate that the force strategy needs to be re-initialized.
    d_ib_pd_force_fcn_needs_init = true;

    const int struct_ln = d_hierarchy->getFinestLevelNumber();
    if (initial_time)
    {
        VecSet(d_l_data_manager->getLData("damage", struct_ln)->getVec(), 0.0);
    }

    // Register plot quantities.
    if (d_silo_writer)
    {
        Pointer<LData> dmg_data = d_l_data_manager->getLData("damage", struct_ln);
        d_silo_writer->registerVariableData("damage", dmg_data, struct_ln);
    }

    return;
} // initializePatchHierarchy

void
IBPDMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                int level_number,
                                double init_data_time,
                                bool can_be_refined,
                                bool initial_time,
                                Pointer<BasePatchLevel<NDIM> > old_level,
                                bool allocate_data)
{
    IBMethod::initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    // Allocate LData corresponding to the Lagrange multiplier.
    if (initial_time && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        // Create Lagrange multiplier and regularization data.
        d_l_data_manager->createLData("damage", level_number, 1, /*manage_data*/ true);
    }

    return;
} // initializeLevelData

void
IBPDMethod::putToDatabase(Pointer<Database> db)
{
    IBMethod::putToDatabase(db);
    db->putInteger("IB_PD_METHOD_VERSION", IB_PD_METHOD_VERSION);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBPDMethod::resetLagrangianPDForceFunction(const double init_data_time, const bool initial_time)
{
    if (!d_ib_pd_force_fcn) return;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_ib_pd_force_fcn->initializeLevelData(d_hierarchy, ln, init_data_time, initial_time, d_l_data_manager);
    }
    return;
} // resetLagrangianPDForceFunction

void
IBPDMethod::getFromInput(Pointer<Database> /*db*/, bool /*is_from_restart*/)
{
    // intentionally blank
    return;
} // getFromInput

void
IBPDMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file."
                                 << std::endl);
    }
    int ver = db->getInteger("IB_PD_METHOD_VERSION");
    if (ver != IB_PD_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
