// Filename: PenaltyIBMethod.cpp
// Created on 28 Sep 2011 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>
#include <string>
#include <vector>

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "boost/multi_array.hpp"
#include "ibamr/IBMethod.h"
#include "ibamr/PenaltyIBMethod.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LInitStrategy.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of PenaltyIBMethod restart file data.
static const int PENALTY_IB_METHOD_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

PenaltyIBMethod::PenaltyIBMethod(const std::string& object_name, Pointer<Database> input_db, bool register_for_restart)
    : IBMethod(object_name, input_db, register_for_restart)
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object
    // name.

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);
    return;
} // PenaltyIBMethod

PenaltyIBMethod::~PenaltyIBMethod()
{
    // intentionally blank
    //
    // NOTE: Parent class constructor unregisters class with the restart
    // manager.
    return;
} // ~PenaltyIBMethod

void
PenaltyIBMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Look-up or allocate Lagangian data.
    d_K_data.resize(finest_ln + 1);
    d_M_data.resize(finest_ln + 1);
    d_Y_current_data.resize(finest_ln + 1);
    d_Y_new_data.resize(finest_ln + 1);
    d_V_current_data.resize(finest_ln + 1);
    d_V_new_data.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_K_data[ln] = d_l_data_manager->getLData("K", ln);
        d_M_data[ln] = d_l_data_manager->getLData("M", ln);
        d_Y_current_data[ln] = d_l_data_manager->getLData("Y", ln);
        d_Y_new_data[ln] = d_l_data_manager->createLData("Y_new", ln, NDIM);
        d_V_current_data[ln] = d_l_data_manager->getLData("V", ln);
        d_V_new_data[ln] = d_l_data_manager->createLData("V_new", ln, NDIM);

        // Initialize Y^{n+1} and V^{n+1} to equal Y^{n} and V^{n}.
        int ierr;
        ierr = VecCopy(d_Y_current_data[ln]->getVec(), d_Y_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_V_current_data[ln]->getVec(), d_V_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }
    return;
} // preprocessIntegrateData

void
PenaltyIBMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBMethod::postprocessIntegrateData(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Reset time-dependent Lagrangian data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        int ierr;
        ierr = VecSwap(d_Y_current_data[ln]->getVec(), d_Y_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        ierr = VecSwap(d_V_current_data[ln]->getVec(), d_V_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }

    // Deallocate Lagrangian scratch data.
    d_K_data.clear();
    d_M_data.clear();
    d_Y_current_data.clear();
    d_Y_new_data.clear();
    d_V_current_data.clear();
    d_V_new_data.clear();
    return;
} // postprocessIntegrateData

void
PenaltyIBMethod::forwardEulerStep(const double current_time, const double new_time)
{
    IBMethod::forwardEulerStep(current_time, new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Update the values of Y^{n+1} and V^{n+1} using forward Euler.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        const double* const K = d_K_data[ln]->getLocalFormArray()->data();
        const double* const M = d_M_data[ln]->getLocalFormArray()->data();
        const double* const X = d_X_current_data[ln]->getLocalFormVecArray()->data();
        const double* const Y = d_Y_current_data[ln]->getLocalFormVecArray()->data();
        const double* const V = d_V_current_data[ln]->getLocalFormVecArray()->data();
        double* const Y_new = d_Y_new_data[ln]->getLocalFormVecArray()->data();
        double* const V_new = d_V_new_data[ln]->getLocalFormVecArray()->data();
        const unsigned int n_local = d_X_current_data[ln]->getLocalNodeCount();
        unsigned int i, d;
        for (i = 0; i < n_local; ++i)
        {
            for (d = 0; d < NDIM; ++d)
            {
                Y_new[NDIM * i + d] = Y[NDIM * i + d] + dt * V[NDIM * i + d];
                V_new[NDIM * i + d] =
                    V[NDIM * i + d] +
                    dt * (-K[i] * (Y[NDIM * i + d] - X[NDIM * i + d]) / M[i] + d_gravitational_acceleration[d]);
            }
        }
    }
    return;
} // eulerStep

void
PenaltyIBMethod::midpointStep(const double current_time, const double new_time)
{
    IBMethod::midpointStep(current_time, new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Update the values of Y^{n+1} and V^{n+1} using the midpoint rule.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        const double* const K = d_K_data[ln]->getLocalFormArray()->data();
        const double* const M = d_M_data[ln]->getLocalFormArray()->data();
        const double* const X = d_X_current_data[ln]->getLocalFormVecArray()->data();
        const double* const Y = d_Y_current_data[ln]->getLocalFormVecArray()->data();
        const double* const V = d_V_current_data[ln]->getLocalFormVecArray()->data();
        const double* const X_new = d_X_new_data[ln]->getLocalFormVecArray()->data();
        double* const Y_new = d_Y_new_data[ln]->getLocalFormVecArray()->data();
        double* const V_new = d_V_new_data[ln]->getLocalFormVecArray()->data();
        const unsigned int n_local = d_X_current_data[ln]->getLocalNodeCount();
        unsigned int i, d;
        double X_half, Y_half, V_half;
        for (i = 0; i < n_local; ++i)
        {
            for (d = 0; d < NDIM; ++d)
            {
                X_half = 0.5 * (X[NDIM * i + d] + X_new[NDIM * i + d]);
                Y_half = 0.5 * (Y[NDIM * i + d] + Y_new[NDIM * i + d]);
                V_half = 0.5 * (V[NDIM * i + d] + V_new[NDIM * i + d]);
                Y_new[NDIM * i + d] = Y[NDIM * i + d] + dt * V_half;
                V_new[NDIM * i + d] =
                    V[NDIM * i + d] + dt * (-K[i] * (Y_half - X_half) / M[i] + d_gravitational_acceleration[d]);
            }
        }
    }
    return;
} // midpointStep

void
PenaltyIBMethod::trapezoidalStep(const double current_time, const double new_time)
{
    IBMethod::trapezoidalStep(current_time, new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Update the values of Y^{n+1} and V^{n+1} using the trapezoidal rule.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        const double* const K = d_K_data[ln]->getLocalFormArray()->data();
        const double* const M = d_M_data[ln]->getLocalFormArray()->data();
        const double* const X = d_X_current_data[ln]->getLocalFormVecArray()->data();
        const double* const Y = d_Y_current_data[ln]->getLocalFormVecArray()->data();
        const double* const V = d_V_current_data[ln]->getLocalFormVecArray()->data();
        const double* const X_new = d_X_new_data[ln]->getLocalFormVecArray()->data();
        double* const Y_new = d_Y_new_data[ln]->getLocalFormVecArray()->data();
        double* const V_new = d_V_new_data[ln]->getLocalFormVecArray()->data();
        const unsigned int n_local = d_X_current_data[ln]->getLocalNodeCount();
        unsigned int i, d;
        double X_half, Y_half, V_half;
        for (i = 0; i < n_local; ++i)
        {
            for (d = 0; d < NDIM; ++d)
            {
                X_half = 0.5 * (X[NDIM * i + d] + X_new[NDIM * i + d]);
                Y_half = 0.5 * (Y[NDIM * i + d] + Y_new[NDIM * i + d]);
                V_half = 0.5 * (V[NDIM * i + d] + V_new[NDIM * i + d]);
                Y_new[NDIM * i + d] = Y[NDIM * i + d] + dt * V_half;
                V_new[NDIM * i + d] =
                    V[NDIM * i + d] + dt * (-K[i] * (Y_half - X_half) / M[i] + d_gravitational_acceleration[d]);
            }
        }
    }
    return;
} // trapezoidalStep

void
PenaltyIBMethod::computeLagrangianForce(const double data_time)
{
    IBMethod::computeLagrangianForce(data_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    double max_displacement = 0.0;

    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            double* const F = d_F_current_data[ln]->getLocalFormVecArray()->data();
            const double* const K = d_K_data[ln]->getLocalFormArray()->data();
            const double* const X = d_X_current_data[ln]->getLocalFormVecArray()->data();
            const double* const Y = d_Y_current_data[ln]->getLocalFormVecArray()->data();
            const unsigned int n_local = d_X_current_data[ln]->getLocalNodeCount();
            unsigned int i, d;
            double dX;
            for (i = 0; i < n_local; ++i)
            {
                dX = 0.0;
                for (d = 0; d < NDIM; ++d)
                {
                    F[NDIM * i + d] += K[i] * (Y[NDIM * i + d] - X[NDIM * i + d]);
                    dX += (Y[NDIM * i + d] - X[NDIM * i + d]) * (Y[NDIM * i + d] - X[NDIM * i + d]);
                }
                dX = sqrt(dX);
                max_displacement = std::max(max_displacement, dX);
            }
        }
        d_F_current_needs_ghost_fill = true;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            double* const F = d_F_half_data[ln]->getLocalFormVecArray()->data();
            const double* const K = d_K_data[ln]->getLocalFormArray()->data();
            const double* const X = d_X_current_data[ln]->getLocalFormVecArray()->data();
            const double* const Y = d_Y_current_data[ln]->getLocalFormVecArray()->data();
            const double* const X_new = d_X_new_data[ln]->getLocalFormVecArray()->data();
            const double* const Y_new = d_Y_new_data[ln]->getLocalFormVecArray()->data();
            const unsigned int n_local = d_X_current_data[ln]->getLocalNodeCount();
            unsigned int i, d;
            double X_half, Y_half, dX;
            for (i = 0; i < n_local; ++i)
            {
                dX = 0.0;
                for (d = 0; d < NDIM; ++d)
                {
                    X_half = 0.5 * (X[NDIM * i + d] + X_new[NDIM * i + d]);
                    Y_half = 0.5 * (Y[NDIM * i + d] + Y_new[NDIM * i + d]);
                    F[NDIM * i + d] += K[i] * (Y_half - X_half);
                    dX += (Y_half - X_half) * (Y_half - X_half);
                }
                dX = sqrt(dX);
                max_displacement = std::max(max_displacement, dX);
            }
        }
        d_F_half_needs_ghost_fill = true;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
            double* const F = d_F_new_data[ln]->getLocalFormVecArray()->data();
            const double* const K = d_K_data[ln]->getLocalFormArray()->data();
            const double* const X = d_X_new_data[ln]->getLocalFormVecArray()->data();
            const double* const Y = d_Y_new_data[ln]->getLocalFormVecArray()->data();
            const unsigned int n_local = d_X_current_data[ln]->getLocalNodeCount();
            unsigned int i, d;
            double dX;
            for (i = 0; i < n_local; ++i)
            {
                dX = 0.0;
                for (d = 0; d < NDIM; ++d)
                {
                    F[NDIM * i + d] += K[i] * (Y[NDIM * i + d] - X[NDIM * i + d]);
                    dX += (Y[NDIM * i + d] - X[NDIM * i + d]) * (Y[NDIM * i + d] - X[NDIM * i + d]);
                }
                dX = sqrt(dX);
                max_displacement = std::max(max_displacement, dX);
            }
        }
        d_F_new_needs_ghost_fill = true;
    }

    if (d_do_log)
    {
        max_displacement = SAMRAI_MPI::maxReduction(max_displacement);
        plog << d_object_name << "::computeLagrangianForce(): maximum point displacement: " << max_displacement << "\n";
    }
    return;
} // computeLagrangianForce

void
PenaltyIBMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                          Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                          int u_data_idx,
                                          const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                          const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                          int integrator_step,
                                          double init_data_time,
                                          bool initial_time)
{
    IBMethod::initializePatchHierarchy(hierarchy,
                                       gridding_alg,
                                       u_data_idx,
                                       u_synch_scheds,
                                       u_ghost_fill_scheds,
                                       integrator_step,
                                       init_data_time,
                                       initial_time);

    // Initialize various Lagrangian data objects.
    if (initial_time)
    {
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

            const bool can_be_refined = ln < finest_ln || d_gridding_alg->levelCanBeRefined(ln);

            Pointer<LData> X_data = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
            Pointer<LData> U_data = d_l_data_manager->getLData(LDataManager::VEL_DATA_NAME, ln);
            Pointer<LData> M_data = d_l_data_manager->createLData("M", ln, 1, /*manage_data*/ true);
            Pointer<LData> K_data = d_l_data_manager->createLData("K", ln, 1, /*manage_data*/ true);
            Pointer<LData> Y_data = d_l_data_manager->createLData("Y", ln, NDIM, /*manage_data*/ true);
            Pointer<LData> V_data = d_l_data_manager->createLData("V", ln, NDIM, /*manage_data*/ true);
            static const int global_index_offset = 0;
            static const int local_index_offset = 0;
            d_l_initializer->initializeMassDataOnPatchLevel(global_index_offset,
                                                            local_index_offset,
                                                            M_data,
                                                            K_data,
                                                            d_hierarchy,
                                                            ln,
                                                            init_data_time,
                                                            can_be_refined,
                                                            initial_time,
                                                            d_l_data_manager);
            if (d_silo_writer)
            {
                d_silo_writer->registerVariableData("M", M_data, ln);
                d_silo_writer->registerVariableData("Y", Y_data, ln);
            }

            // Set initial conditions.
            int ierr;
            ierr = VecCopy(X_data->getVec(), Y_data->getVec());
            IBTK_CHKERRQ(ierr);
            ierr = VecCopy(U_data->getVec(), V_data->getVec());
            IBTK_CHKERRQ(ierr);
        }
    }
    return;
} // initializePatchHierarchy

void
PenaltyIBMethod::putToDatabase(Pointer<Database> db)
{
    IBMethod::putToDatabase(db);

    db->putInteger("PENALTY_IB_METHOD_VERSION", PENALTY_IB_METHOD_VERSION);
    db->putDoubleArray("d_gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PenaltyIBMethod::getFromInput(Pointer<Database> db, bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->keyExists("gravitational_acceleration"))
        {
            db->getDoubleArray("gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
        }
        else
        {
            TBOX_WARNING(d_object_name << ":  "
                                       << "Using penalty-IB method but key data "
                                          "`gravitational_acceleration' not found in input.");
        }
    }
    return;
} // getFromInput

void
PenaltyIBMethod::getFromRestart()
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
    int ver = db->getInteger("PENALTY_IB_METHOD_VERSION");
    if (ver != PENALTY_IB_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    db->getDoubleArray("d_gravitational_acceleration", &d_gravitational_acceleration[0], NDIM);
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
