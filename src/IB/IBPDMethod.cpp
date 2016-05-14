// Filename: IBPDMethod.cpp
// Created on 08 Apr 2016 by Amneet Bhalla
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
static const int interior_begin = 0;
static const int interior_end = 39999;
static const int bottom_begin = 40000;
static const int bottom_end = 40617;
static const int top_begin = 40618;
static const int top_end = 41235;
static const int left_begin = 41236;
static const int left_end = 41835;
static const int right_begin = 41836;
static const int right_end = 42435;

static const double dens = 1.0;
static const double DX  = 1.0/199.0;
static const double t_ramp = 5.0;
void
get_bodyforce(double* f_vec, const double* X, const double t)
{
    const double x = X[0];
    const double y = X[1];

    f_vec[0] = -2 * x * (1/cosh(t) * 1/cosh(t)) *(-5*y + 6*y*cosh(2*t) + 4*sinh(2*t))* tanh(t);
    f_vec[1] = -tanh(t)* (11 + 2*y*y*(1/cosh(t) * 1/cosh(t)) + 61*y*tanh(t) +
                         4*(x*x + 14*y*y)* tanh(t)*tanh(t));
    return;
} // get_bodyforce

void
get_trac(Eigen::Vector2d& trac, const Eigen::Vector2d& normal, const double* X_0, const double t)
{
    const double x = X_0[0];
    const double y = X_0[1];

    static const double L = 2.0;
    static const double M = 1.0;
    static const Eigen::Matrix2d II = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d FF;
    FF << y*tanh(t) + 1, x*tanh(t), 0, 2*y*tanh(t) + 1;
    const Eigen::Matrix2d E = 0.5 * (FF.transpose() * FF - II);
    const double trE = E.trace();

    const Eigen::Matrix2d PK1 = L * trE * FF + 2 * M * FF * E;
    trac = PK1 * normal;

} // get_trac
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
    const int struct_ln = d_hierarchy->getFinestLevelNumber();

    const int step_no = d_ib_solver->getIntegratorStep() + 1;
    if (step_no % 500 == 0)
    {
        Pointer<LData> D_LData = d_l_data_manager->getLData("damage", struct_ln);
        Vec D_petsc_vec_parallel = D_LData->getVec();
        Vec D_lag_vec_parallel = NULL;
        Vec D_lag_vec_seq = NULL;
        VecDuplicate(D_petsc_vec_parallel, &D_lag_vec_parallel);
        d_l_data_manager->scatterPETScToLagrangian(D_petsc_vec_parallel, D_lag_vec_parallel, struct_ln);
        d_l_data_manager->scatterToZero(D_lag_vec_parallel, D_lag_vec_seq);

        Pointer<LData> X0_LData = d_l_data_manager->getLData("X0", struct_ln);
        Vec X0_petsc_vec_parallel = X0_LData->getVec();
        Vec X0_lag_vec_parallel = NULL;
        Vec X0_lag_vec_seq = NULL;
        VecDuplicate(X0_petsc_vec_parallel, &X0_lag_vec_parallel);
        d_l_data_manager->scatterPETScToLagrangian(X0_petsc_vec_parallel, X0_lag_vec_parallel, struct_ln);
        d_l_data_manager->scatterToZero(X0_lag_vec_parallel, X0_lag_vec_seq);

        Pointer<LData> X_LData = d_X_new_data[struct_ln];
        Vec X_petsc_vec_parallel = X_LData->getVec();
        Vec X_lag_vec_parallel = NULL;
        Vec X_lag_vec_seq = NULL;
        VecDuplicate(X_petsc_vec_parallel, &X_lag_vec_parallel);
        d_l_data_manager->scatterPETScToLagrangian(X_petsc_vec_parallel, X_lag_vec_parallel, struct_ln);
        d_l_data_manager->scatterToZero(X_lag_vec_parallel, X_lag_vec_seq);

        if (SAMRAI_MPI::getRank() == 0)
        {
            const PetscScalar* D;
            VecGetArrayRead(D_lag_vec_seq, &D);
            int counter_D = -1;

            const PetscScalar* X0;
            VecGetArrayRead(X0_lag_vec_seq, &X0);
            int counter_X0 = -1;

            const PetscScalar* X;
            VecGetArrayRead(X_lag_vec_seq, &X);
            int counter_X = -1;

            std::fstream D_stream;
            std::ostringstream D_sstream;
            D_sstream << "./data/D_" << step_no << "_" << new_time;
            D_stream.open(D_sstream.str().c_str(), std::fstream::out);

            int ib_pts;
            VecGetSize(D_lag_vec_seq, &ib_pts);
            for (int i = 0; i < ib_pts; ++i)
            {
                const double X0_0 = X0[++counter_X0];
                const double X0_1 = X0[++counter_X0];
                const double X_0 = X[++counter_X];
                const double X_1 = X[++counter_X];
                const double dmg = D[++counter_D];
                D_stream << X0_0 << "\t" << X0_1 << "\t" << X_0 - X0_0 << "\t" << X_1 - X0_1 << "\t" << dmg
                         << std::endl;
            }

            VecRestoreArrayRead(D_lag_vec_seq, &D);
            VecRestoreArrayRead(X0_lag_vec_seq, &X0);
            VecRestoreArrayRead(X_lag_vec_seq, &X);
        }
        VecDestroy(&D_lag_vec_parallel);
        VecDestroy(&D_lag_vec_seq);
        VecDestroy(&X0_lag_vec_parallel);
        VecDestroy(&X0_lag_vec_seq);
        VecDestroy(&X_lag_vec_parallel);
        VecDestroy(&X_lag_vec_seq);
    }
    IBMethod::postprocessIntegrateData(current_time, new_time, num_cycles);
    return;
} // postprocessIntegrateData

void
IBPDMethod::interpolateVelocity(const int u_data_idx,
                                const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                const double data_time)
{
    // Interpolate the linear velocities.
    // IBMethod::interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);

    return;
} // interpolateVelocity

void
IBPDMethod::eulerStep(const double current_time, const double new_time)
{
    // IBMethod::eulerStep(current_time, new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        Pointer<LData> X_0_data = d_l_data_manager->getLData("X0", ln);
        Pointer<LData> X_half_data = d_X_half_data[ln];
        Pointer<LData> U_half_data = d_U_half_data[ln];

        boost::multi_array_ref<double, 2>& X_0_data_array = *X_0_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& X_half_data_array = *X_half_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& U_half_data_array = *U_half_data->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();

            const double* X_0 = &X_0_data_array[local_idx][0];
            double* U_half = &U_half_data_array[local_idx][0];
            double* X_half = &X_half_data_array[local_idx][0];

            const double u = X_0[0] * X_0[1] * tanh(new_time);
            const double v = X_0[1] * X_0[1] * tanh(new_time);
            if (lag_idx >= bottom_begin && lag_idx <= bottom_end)
            {
                for (int d = 0; d < NDIM; ++d) U_half[d] = 0.0;
                X_half[0] = X_0[0];
                X_half[1] = X_0[1];
            }
        }
    }
    return;
} // eulerStep

void
IBPDMethod::midpointStep(const double current_time, const double new_time)
{

    const double dt = new_time - current_time;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        Pointer<LData> X_0_data = d_l_data_manager->getLData("X0", ln);
        Pointer<LData> X_current_data = d_X_current_data[ln];
        Pointer<LData> X_new_data = d_X_new_data[ln];

        Pointer<LData> U_current_data = d_U_current_data[ln];
        Pointer<LData> U_new_data = d_U_new_data[ln];

        Pointer<LData> F_half_data = d_F_half_data[ln];

        boost::multi_array_ref<double, 2>& X_0_data_array = *X_0_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& X_current_data_array = *X_current_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& X_new_data_array = *X_new_data->getLocalFormVecArray();

        boost::multi_array_ref<double, 2>& F_half_data_array = *F_half_data->getLocalFormVecArray();

        boost::multi_array_ref<double, 2>& U_current_data_array = *U_current_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& U_new_data_array = *U_new_data->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();

            const double* U_current = &U_current_data_array[local_idx][0];
            double* U_new = &U_new_data_array[local_idx][0];

            const double* X_0 = &X_0_data_array[local_idx][0];
            const double* X_current = &X_current_data_array[local_idx][0];
            double* X_new = &X_new_data_array[local_idx][0];
            const double* F_half = &F_half_data_array[local_idx][0];

            const double u = X_0[0] * X_0[1] * tanh(new_time);
            const double v = X_0[1] * X_0[1] * tanh(new_time);
            if (lag_idx >= interior_begin && lag_idx <= interior_end)
            {
                double bforce[NDIM];
                get_bodyforce(bforce, X_0, current_time);
                for (int d = 0; d < NDIM; ++d)
                {
                    U_new[d] = U_current[d] + (dt / dens) * (1.0*F_half[d] + 1.0*bforce[d]);
                    X_new[d] = X_current[d] + dt * U_new[d];
                }
            }
            else if (lag_idx >= left_begin && lag_idx <= left_end)
            {
                Eigen::Vector2d normal(-1.0, 0.0);
                Eigen::Vector2d trac;

                double bforce[NDIM];
                get_bodyforce(bforce, X_0, current_time);
                get_trac(trac, normal, X_0, current_time);

                for (int d = 0; d < NDIM; ++d)
                {
                    U_new[d] = U_current[d] + (dt / dens) * (1.0 * F_half[d] + 1.0 * bforce[d] + 1.0 * trac[d] / (DX));
                    X_new[d] = X_current[d] + dt * U_new[d];
                }
            }
            else if (lag_idx >= right_begin && lag_idx <= right_end)
            {
                Eigen::Vector2d normal(1.0, 0.0);
                Eigen::Vector2d trac;

                double bforce[NDIM];
                get_bodyforce(bforce, X_0, current_time);
                get_trac(trac, normal, X_0, current_time);

                for (int d = 0; d < NDIM; ++d)
                {
                    U_new[d] = U_current[d] + (dt / dens) * (1.0 * F_half[d] + 1.0 * bforce[d] + 1.0 * trac[d] / (DX));
                    X_new[d] = X_current[d] + dt * U_new[d];
                }
            }
            else if (lag_idx >= top_begin && lag_idx <= top_end)
            {
                Eigen::Vector2d normal(0.0, 1.0);
                Eigen::Vector2d trac;

                double bforce[NDIM];
                get_bodyforce(bforce, X_0, current_time);
                get_trac(trac, normal, X_0, current_time);

                for (int d = 0; d < NDIM; ++d)
                {
                    U_new[d] = U_current[d] + (dt / dens) * (1.0 * F_half[d] + 1.0 * bforce[d] + 1.0 * trac[d] / (DX));
                    X_new[d] = X_current[d] + dt * U_new[d];
                }
            }
            else
            {
                TBOX_ASSERT(lag_idx >= bottom_begin && lag_idx <= bottom_end);
                for (int d = 0; d < NDIM; ++d) U_new[d] = 0.0;
                if (new_time <= t_ramp)
                {
                    X_new[0] = X_0[0] + u * new_time / t_ramp;
                    X_new[1] = X_0[1] + v * new_time / t_ramp;
                }
                else
                {
                    X_new[0] = X_0[0] + u;
                    X_new[1] = X_0[1] + v;
                }
            }
        }

        X_0_data->restoreArrays();
        X_current_data->restoreArrays();
        X_new_data->restoreArrays();
        F_current_data->restoreArrays();
        F_half_data->restoreArrays();
        U_current_data->restoreArrays();
        U_new_data->restoreArrays();
    }

    return;

} // midpointStep

void
IBPDMethod::trapezoidalStep(const double current_time, const double new_time)
{
    TBOX_ERROR(d_object_name << "::trapezoidalStep():\n"
                             << "  time-stepping type TRAPEZOIDAL_RULE not supported by class IBPDMethod;\n"
                             << "  use MIDPOINT_RULE instead.\n");
    return;
} // trapezoidalStep

void
IBPDMethod::computeLagrangianForce(const double data_time)
{
    // IBMethod::computeLagrangianForce(data_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> >* F_data = NULL;
    std::vector<Pointer<LData> >* X_data = NULL;
    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        d_F_current_needs_ghost_fill = true;
        F_data = &d_F_current_data;
        X_data = &d_X_current_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_half_time))
    {
        d_F_half_needs_ghost_fill = true;
        F_data = &d_F_half_data;
        X_data = &d_X_half_data;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        d_F_new_needs_ghost_fill = true;
        F_data = &d_F_new_data;
        X_data = &d_X_new_data;
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        if (d_ib_pd_force_fcn)
        {
            d_ib_pd_force_fcn->computeLagrangianForceAndDamage((*F_data)[ln],
                                                               d_l_data_manager->getLData("damage", ln),
                                                               (*X_data)[ln],
                                                               Pointer<LData>(NULL) /*U*/,
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
    // IBMethod::spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);

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

    if (initial_time)
    {
        const int struct_ln = d_hierarchy->getFinestLevelNumber();
        VecSet(d_l_data_manager->getLData("damage", struct_ln)->getVec(), 0.0);
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
