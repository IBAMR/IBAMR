// Filename: GeneralizedIBMethod.C
// Created on 12 Dec 2011 by Boyce Griffith
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

#include "GeneralizedIBMethod.h"

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
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/LSiloDataWriter.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of GeneralizedIBMethod restart file data.
static const int GENERALIZED_IB_METHOD_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

GeneralizedIBMethod::GeneralizedIBMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : IBMethod(object_name, input_db, register_for_restart)
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object name.

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);
    return;
}// GeneralizedIBMethod

GeneralizedIBMethod::~GeneralizedIBMethod()
{
    // intentionally blank
    //
    // NOTE: Parent class constructor unregisters class with the restart
    // manager.
    return;
}// ~GeneralizedIBMethod

void
GeneralizedIBMethod::registerEulerianVariables()
{
    IBMethod::registerEulerianVariables();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const IntVector<NDIM> ib_ghosts = getMinimumGhostCellWidth();
    const IntVector<NDIM>    ghosts = 1;

    Pointer<Variable<NDIM> > u_var = d_ib_solver->getINSHierarchyIntegrator()->getVelocityVariable();
    Pointer<CellVariable<NDIM,double> > u_cc_var = u_var;
    Pointer<SideVariable<NDIM,double> > u_sc_var = u_var;
    if (!u_cc_var.isNull())
    {
        d_w_var = new CellVariable<NDIM,double>(d_object_name+"::w", NDIM);
        d_n_var = new CellVariable<NDIM,double>(d_object_name+"::n", NDIM);
    }
    else if (!u_sc_var.isNull())
    {
        d_w_var = new SideVariable<NDIM,double>(d_object_name+"::w");
        d_n_var = new SideVariable<NDIM,double>(d_object_name+"::n");
    }
    else
    {
        TBOX_ERROR(d_object_name << "::registerEulerianVariables():\n"
                   << "  unsupported velocity data centering" << std::endl);
    }
    d_w_idx = var_db->registerVariableAndContext(d_w_var, d_ib_solver->getScratchContext(), ib_ghosts);
    d_n_idx = var_db->registerVariableAndContext(d_n_var, d_ib_solver->getScratchContext(),    ghosts);
    return;
}// registerEulerianVariables

void
GeneralizedIBMethod::registerEulerianCommunicationAlgorithms()
{
    IBMethod::registerEulerianCommunicationAlgorithms();

    Pointer<Geometry<NDIM> > grid_geom = d_ib_solver->getPatchHierarchy()->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg;
    Pointer<RefineOperator<NDIM> > refine_op;

    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = NULL;
    refine_alg->registerRefine(d_w_idx, d_w_idx, d_w_idx, refine_op);
    registerGhostfillRefineAlgorithm(d_object_name+"::w", refine_alg);

    refine_alg = new RefineAlgorithm<NDIM>();
    refine_op = NULL;
    refine_alg->registerRefine(d_n_idx, d_n_idx, d_n_idx, refine_op);
    registerGhostfillRefineAlgorithm(d_object_name+"::n", refine_alg);
    return;
}// registerEulerianCommunicationAlgorithms

void
GeneralizedIBMethod::preprocessIntegrateData(
    double current_time,
    double new_time,
    int num_cycles)
{
    IBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Look-up or allocate Lagangian data.
    d_D_current_data.resize(finest_ln+1);
    d_D_new_data    .resize(finest_ln+1);
    d_N_current_data.resize(finest_ln+1);
    d_N_new_data    .resize(finest_ln+1);
    d_W_current_data.resize(finest_ln+1);
    d_W_new_data    .resize(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_D_current_data[ln] = d_l_data_manager->getLData("D",ln);
        d_D_new_data    [ln] = d_l_data_manager->createLData("D_new",ln,NDIM*NDIM);
        d_N_current_data[ln] = d_l_data_manager->createLData("N",ln,NDIM);
        d_N_new_data    [ln] = d_l_data_manager->createLData("N_new",ln,NDIM);
        d_W_current_data[ln] = d_l_data_manager->getLData("W",ln);
        d_W_new_data    [ln] = d_l_data_manager->createLData("W_new",ln,NDIM);

        // Initialize D^{n+1} to equal D^{n}, and initialize W^{n+1} to equal
        // W^{n}.
        int ierr;
        ierr = VecCopy(d_D_current_data[ln]->getVec(), d_D_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_W_current_data[ln]->getVec(), d_W_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }
    return;
}// preprocessIntegrateData

void
GeneralizedIBMethod::postprocessIntegrateData(
    double current_time,
    double new_time,
    int num_cycles)
{
    IBMethod::postprocessIntegrateData(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Reset time-dependent Lagrangian data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        int ierr;
        ierr = VecSwap(d_D_current_data[ln]->getVec(), d_D_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecSwap(d_W_current_data[ln]->getVec(), d_W_new_data[ln]->getVec());  IBTK_CHKERRQ(ierr);
    }

    // Deallocate Lagrangian scratch data.
    d_D_current_data.clear();
    d_D_new_data    .clear();
    d_N_current_data.clear();
    d_N_new_data    .clear();
    d_W_current_data.clear();
    d_W_new_data    .clear();
    return;
}// postprocessIntegrateData

void
GeneralizedIBMethod::eulerStep(
    const double current_time,
    const double new_time)
{
    IBMethod::eulerStep(current_time, new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
//  const double dt = new_time-current_time;

    // XXXX.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
    }
    return;
}// eulerStep

void
GeneralizedIBMethod::midpointStep(
    const double current_time,
    const double new_time)
{
    IBMethod::midpointStep(current_time, new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
//  const double dt = new_time-current_time;

    // XXXX.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
    }
    return;
}// midpointStep

void
GeneralizedIBMethod::computeLagrangianForce(
    const double data_time)
{
    IBMethod::computeLagrangianForce(data_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        }
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        }
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        }
    }
    return;
}// computeLagrangianForce

void
GeneralizedIBMethod::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
    int u_data_idx,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
    int integrator_step,
    double init_data_time,
    bool initial_time)
{
    IBMethod::initializePatchHierarchy(hierarchy, gridding_alg, u_data_idx, u_synch_scheds, u_ghost_fill_scheds, integrator_step, init_data_time, initial_time);

    // Initialize various Lagrangian data objects.
    if (initial_time)
    {
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

            const bool can_be_refined = ln < finest_ln || d_gridding_alg->levelCanBeRefined(ln);

            Pointer<LData> D_data = d_l_data_manager->createLData("D",ln,NDIM*NDIM,/*manage_data*/ true);
            Pointer<LData> W_data = d_l_data_manager->createLData("W",ln,NDIM     ,/*manage_data*/ true);
            static const int global_index_offset = 0;
            static const int local_index_offset = 0;
            d_l_initializer->initializeDirectorDataOnPatchLevel(
                global_index_offset, local_index_offset,
                D_data,
                hierarchy, ln,
                init_data_time, can_be_refined, initial_time, d_l_data_manager);
            if (!d_silo_writer.isNull())
            {
                d_silo_writer->registerVariableData("D1", D_data, 0, 3, ln);
                d_silo_writer->registerVariableData("D2", D_data, 3, 3, ln);
                d_silo_writer->registerVariableData("D3", D_data, 6, 3, ln);
                d_silo_writer->registerVariableData("W", W_data, ln);
            }
        }

        TBOX_ERROR("need to initialize W\n");
    }
    return;
}// initializePatchHierarchy

void
GeneralizedIBMethod::putToDatabase(
    Pointer<Database> db)
{
    IBMethod::putToDatabase(db);
    db->putInteger("GENERALIZED_IB_METHOD_VERSION", GENERALIZED_IB_METHOD_VERSION);
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
GeneralizedIBMethod::getFromInput(
    Pointer<Database> /*db*/,
    bool /*is_from_restart*/)
{
    // intentionally blank
    return;
}// getFromInput

void
GeneralizedIBMethod::getFromRestart()
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
    int ver = db->getInteger("GENERALIZED_IB_METHOD_VERSION");
    if (ver != GENERALIZED_IB_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
