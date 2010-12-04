// Filename: IBImplicitSJROperator.C
// Created on 30 Aug 2010 by Boyce Griffith
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

#include "IBImplicitSJROperator.h"

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
#include <ibamr/IBImplicitHierarchyIntegrator.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/PETScVecUtilities.h>

// SAMRAI INCLUDES
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_apply;
static Pointer<Timer> t_initialize_operator_state;
static Pointer<Timer> t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitSJROperator::IBImplicitSJROperator(
    IBImplicitHierarchyIntegrator* ib_implicit_integrator)
    : d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_ib_implicit_integrator(ib_implicit_integrator),
      d_SJR_mats(),
      d_d_nnz(),
      d_o_nnz()
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSJROperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSJROperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSJROperator::deallocateOperatorState()");
        timers_need_init = false;
    }
    return;
}// IBImplicitSJROperator

IBImplicitSJROperator::~IBImplicitSJROperator()
{
    deallocateOperatorState();
    return;
}// ~IBImplicitSJROperator

void
IBImplicitSJROperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    return;
}// setTimeInterval

void
IBImplicitSJROperator::formJacobian(
    SAMRAIVectorReal<NDIM,double>& x)
{
    int ierr;

    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();

    // Get the vector components.
    const int u_new_idx  = x.getComponentDescriptorIndex(0);
    const int u_half_ib_idx = d_ib_implicit_integrator->d_u_half_ib_idx;
    const int u_current_idx = d_ib_implicit_integrator->d_u_current_idx;

    std::vector<Pointer<RefineSchedule<NDIM> > >& u_half_ib_rscheds = d_ib_implicit_integrator->d_rscheds["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"];
    Pointer<HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops = d_ib_implicit_integrator->d_hier_sc_data_ops;

    LDataManager* lag_data_manager = d_ib_implicit_integrator->d_lag_data_manager;

    std::vector<Pointer<LNodeLevelData> >& X_data = d_ib_implicit_integrator->d_X_data;
    std::vector<Pointer<LNodeLevelData> >& X_mid_data = d_ib_implicit_integrator->d_X_mid_data;
    std::vector<Pointer<LNodeLevelData> >& X_half_data = d_ib_implicit_integrator->d_X_half_data;
    std::vector<Pointer<LNodeLevelData> >& U_half_data = d_ib_implicit_integrator->d_U_half_data;

    std::vector<Mat>& R_mats = d_ib_implicit_integrator->d_R_mats;

    Pointer<IBLagrangianForceStrategy> lag_force_strategy = d_ib_implicit_integrator->d_lag_force_strategy;

    // Set u(n+1/2) = 0.5*(u(n) + u(n+1)).
    hier_sc_data_ops->linearSum(u_half_ib_idx, 0.5, u_current_idx, 0.5, u_new_idx);

    // Interpolate u(n+1/2) to U(n+1/2).
    lag_data_manager->interp(u_half_ib_idx, U_half_data, X_mid_data, u_half_ib_rscheds, d_current_time);
    d_ib_implicit_integrator->resetAnchorPointValues(U_half_data, coarsest_ln, finest_ln);

    // Set X(n+1/2) = X(n) + 0.5*dt*U(n+1/2).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (lag_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_vec = X_data[ln]->getGlobalVec();
            Vec X_half_vec = X_half_data[ln]->getGlobalVec();
            Vec U_half_vec = U_half_data[ln]->getGlobalVec();
            ierr = VecWAXPY(X_half_vec, 0.5*d_dt, U_half_vec, X_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Compute S dF/dX[X(n+1/2),U(n+1/2)] R.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_SJR_mats[ln] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_SJR_mats[ln]);  IBTK_CHKERRQ(ierr);
        }
        if (lag_data_manager->levelContainsLagrangianData(ln))
        {
            // Compute dF/dX.
            const int num_local_nodes = lag_data_manager->getNumberOfLocalNodes(ln);
            Mat J_mat;
            ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,
                                   NDIM*num_local_nodes, NDIM*num_local_nodes,
                                   PETSC_DETERMINE, PETSC_DETERMINE,
                                   PETSC_DEFAULT, &d_d_nnz[ln][0],
                                   PETSC_DEFAULT, &d_o_nnz[ln][0],
                                   &J_mat);  IBTK_CHKERRQ(ierr);
#ifdef DEBUG_CHECK_ASSERTIONS
            ierr = MatSetOption(J_mat, MAT_NEW_NONZERO_LOCATION_ERR  , PETSC_TRUE); IBTK_CHKERRQ(ierr);
            ierr = MatSetOption(J_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); IBTK_CHKERRQ(ierr);
#endif
            ierr = MatSetBlockSize(J_mat, NDIM);  IBTK_CHKERRQ(ierr);
            lag_force_strategy->computeLagrangianForceJacobian(
                J_mat, MAT_FINAL_ASSEMBLY, -0.25*d_dt, X_half_data[ln], -0.5, U_half_data[ln],
                hierarchy, ln, d_current_time+0.5*d_dt, lag_data_manager);

            // Compute S dF/dX R.
            ierr = MatPtAP(J_mat, R_mats[ln], MAT_INITIAL_MATRIX, 32.0, &d_SJR_mats[ln]); IBTK_CHKERRQ(ierr);
            Pointer<PatchLevel<NDIM> > patch_level = hierarchy->getPatchLevel(ln);
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_level->getGridGeometry();
            const double* const dx0 = grid_geom->getDx();
            const IntVector<NDIM>& ratio = patch_level->getRatio();
            double cell_vol = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                cell_vol *= dx0[d] / double(ratio(d));
            }
            ierr = MatScale(d_SJR_mats[ln], 1.0/cell_vol);  IBTK_CHKERRQ(ierr);

            // Destroy dF/dX.
            ierr = MatDestroy(J_mat);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// formJacobian

Pointer<SAMRAIVectorReal<NDIM,double> >
IBImplicitSJROperator::getBaseVector() const
{
    return NULL; // XXXX
}// getBaseVector

void
IBImplicitSJROperator::apply(
    const bool zero_y_before_spread,
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    t_apply->start();

    int ierr;

    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();

    // Get the vector components.
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int f_idx = y.getComponentDescriptorIndex(0);
    Pointer<SideVariable<NDIM,double> > u_var = x.getComponentVariable(0);
    Pointer<SideVariable<NDIM,double> > f_var = y.getComponentVariable(0);

    LDataManager* lag_data_manager = d_ib_implicit_integrator->d_lag_data_manager;

    // Compute f = S (dF/dX) R u.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (lag_data_manager->levelContainsLagrangianData(ln))
        {
            Pointer<PatchLevel<NDIM> > patch_level = hierarchy->getPatchLevel(ln);
            Vec u_vec = static_cast<Vec>(NULL);
            PETScVecUtilities::constructPatchLevelVec(u_vec, u_idx, u_var, patch_level);
            PETScVecUtilities::copyToPatchLevelVec(u_vec, u_idx, u_var, patch_level);
            Vec f_vec = static_cast<Vec>(NULL);
            PETScVecUtilities::constructPatchLevelVec(f_vec, f_idx, f_var, patch_level);
            if (zero_y_before_spread)
            {
                ierr = MatMult(d_SJR_mats[ln], u_vec, f_vec); IBTK_CHKERRQ(ierr);
            }
            else
            {
                PETScVecUtilities::copyToPatchLevelVec(f_vec, f_idx, f_var, patch_level);
                ierr = MatMultAdd(d_SJR_mats[ln], u_vec, f_vec, f_vec); IBTK_CHKERRQ(ierr);
            }
            PETScVecUtilities::copyFromPatchLevelVec(f_vec, f_idx, f_var, patch_level);
            ierr = VecDestroy(u_vec); IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(f_vec); IBTK_CHKERRQ(ierr);
        }
    }

    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent; // XXXX
    SynchronizationTransactionComponent f_synch_transaction = SynchronizationTransactionComponent(f_idx, "CONSERVATIVE_COARSEN");
    Pointer<SideDataSynchronization> side_synch_op = new SideDataSynchronization();
    side_synch_op->initializeOperatorState(f_synch_transaction, y.getPatchHierarchy());
    side_synch_op->synchronizeData(0.0);

    t_apply->stop();
    return;
}// apply

void
IBImplicitSJROperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    apply(true, x, y);
    return;
}// apply

void
IBImplicitSJROperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    t_initialize_operator_state->start();

    if (d_is_initialized) deallocateOperatorState();

    Pointer<PatchHierarchy<NDIM> > hierarchy = in.getPatchHierarchy();
    const int coarsest_ln = in.getCoarsestLevelNumber();
    const int finest_ln = in.getFinestLevelNumber();

    d_SJR_mats.resize(finest_ln+1, static_cast<Mat>(NULL));
    d_d_nnz.resize(finest_ln+1);
    d_o_nnz.resize(finest_ln+1);

    // Determine the non-zero structure of the force Jacobian matrix.
    LDataManager* lag_data_manager = d_ib_implicit_integrator->d_lag_data_manager;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (lag_data_manager->levelContainsLagrangianData(ln))
        {
            // Compute the nonzero structure for the block matrix representation
            // of dF/dX.
            const int num_local_nodes = lag_data_manager->getNumberOfLocalNodes(ln);
            std::vector<int> d_nnz(num_local_nodes), o_nnz(num_local_nodes);
            Pointer<IBLagrangianForceStrategy> lag_force_strategy = d_ib_implicit_integrator->d_lag_force_strategy;
            lag_force_strategy->computeLagrangianForceJacobianNonzeroStructure(
                d_nnz, o_nnz, hierarchy, ln, d_current_time+0.5*d_dt, lag_data_manager);

            // Convert the nonzero structure data into a non-blocked format.
            d_d_nnz[ln].resize(NDIM*num_local_nodes);
            d_o_nnz[ln].resize(NDIM*num_local_nodes);
            for (int k = 0; k < num_local_nodes; ++k)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    d_d_nnz[ln][NDIM*k+d] = NDIM*d_nnz[k];
                    d_o_nnz[ln][NDIM*k+d] = NDIM*o_nnz[k];
                }
            }
        }
    }

    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
IBImplicitSJROperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    t_deallocate_operator_state->start();

    int ierr;

    for (unsigned k = 0; k < d_SJR_mats.size(); ++k)
    {
        if (d_SJR_mats[k] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_SJR_mats[k]);  IBTK_CHKERRQ(ierr);
        }
    }

    d_SJR_mats.clear();
    d_d_nnz.clear();
    d_o_nnz.clear();

    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
IBImplicitSJROperator::enableLogging(
    bool enabled)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBImplicitSJROperator>;

//////////////////////////////////////////////////////////////////////////////
