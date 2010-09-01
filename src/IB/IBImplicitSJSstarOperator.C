// Filename: IBImplicitSJSstarOperator.C
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

#include "IBImplicitSJSstarOperator.h"

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

IBImplicitSJSstarOperator::IBImplicitSJSstarOperator(
    IBImplicitHierarchyIntegrator* ib_implicit_integrator)
    : d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_ib_implicit_integrator(ib_implicit_integrator)
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSJSstarOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSJSstarOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSJSstarOperator::deallocateOperatorState()");
        timers_need_init = false;
    }
    return;
}// IBImplicitSJSstarOperator

IBImplicitSJSstarOperator::~IBImplicitSJSstarOperator()
{
    deallocateOperatorState();
    return;
}// ~IBImplicitSJSstarOperator

void
IBImplicitSJSstarOperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    return;
}// setTimeInterval

void
IBImplicitSJSstarOperator::formJacobian(
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

    Pointer<IBLagrangianForceStrategy> lag_force_strategy = d_ib_implicit_integrator->d_lag_force_strategy;

    // Set u(n+1/2) = 0.5*(u(n) + u(n+1)).
    hier_sc_data_ops->linearSum(u_half_ib_idx, 0.5, u_current_idx, 0.5, u_new_idx);

    // Interpolate u(n+1/2) to U(n+1/2).
    lag_data_manager->interpolate(u_half_ib_idx, U_half_data, X_mid_data,
                                  u_half_ib_rscheds, d_current_time);
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

    // Compute dF/dX[X(n+1/2),U(n+1/2)].
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (d_dF_dX_mat[ln] != static_cast<Mat>(NULL))
        {
            ierr = MatDestroy(d_dF_dX_mat[ln]);  IBTK_CHKERRQ(ierr);
        }
        if (lag_data_manager->levelContainsLagrangianData(ln))
        {
            // Determine the non-zero structure of the force Jacobian matrix and
            // allocate a block AIJ matrix.
            //
            // XXXX: This should be done elsewhere, e.g. in
            // initializeOperatorState(), so that the sparsity pattern may be
            // re-used.
            const int num_local_nodes = lag_data_manager->getNumberOfLocalNodes(ln);
            std::vector<int> d_nnz(num_local_nodes), o_nnz(num_local_nodes);
            lag_force_strategy->computeLagrangianForceJacobianNonzeroStructure(
                d_nnz, o_nnz, hierarchy, ln, d_current_time+0.5*d_dt, lag_data_manager);
            ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                                    NDIM, NDIM*num_local_nodes, NDIM*num_local_nodes,
                                    PETSC_DETERMINE, PETSC_DETERMINE,
                                    PETSC_DEFAULT, &d_nnz[0],
                                    PETSC_DEFAULT, &o_nnz[0],
                                    &d_dF_dX_mat[ln]);  IBTK_CHKERRQ(ierr);

            // Compute the Jacobian of the force.
            lag_force_strategy->computeLagrangianForceJacobian(
                d_dF_dX_mat[ln], MAT_FINAL_ASSEMBLY, -0.25*d_dt, X_half_data[ln], -0.5, U_half_data[ln],
                hierarchy, ln, d_current_time+0.5*d_dt, lag_data_manager);
        }
    }
    return;
}// formJacobian

Pointer<SAMRAIVectorReal<NDIM,double> >
IBImplicitSJSstarOperator::getBaseVector() const
{
    return NULL; // XXXX
}// getBaseVector

void
IBImplicitSJSstarOperator::apply(
    const bool zero_y_before_spread,
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    t_apply->start();

    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();

    // Get the vector components.
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int f_idx = y.getComponentDescriptorIndex(0);

    const int u_ib_idx = d_ib_implicit_integrator->d_u_half_ib_idx;
    std::vector<Pointer<RefineSchedule<NDIM> > >& u_ib_rscheds = d_ib_implicit_integrator->d_rscheds["V->V::S->S::CONSERVATIVE_LINEAR_REFINE"];
    Pointer<HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops = d_ib_implicit_integrator->d_hier_sc_data_ops;

    LDataManager* lag_data_manager = d_ib_implicit_integrator->d_lag_data_manager;

    std::vector<Pointer<LNodeLevelData> >& X_mid_data = d_ib_implicit_integrator->d_X_mid_data;
    std::vector<Pointer<LNodeLevelData> >& U_data = d_ib_implicit_integrator->d_U_half_data;
    std::vector<Pointer<LNodeLevelData> >& F_data = d_ib_implicit_integrator->d_F_half_data;

    // Interpolate u to U.
    hier_sc_data_ops->copyData(u_ib_idx, u_idx);
    lag_data_manager->interpolate(u_ib_idx, U_data, X_mid_data, u_ib_rscheds, d_current_time);
    d_ib_implicit_integrator->resetAnchorPointValues(U_data, coarsest_ln, finest_ln);

    // Compute F = (dF/dX) U.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (lag_data_manager->levelContainsLagrangianData(ln))
        {
            int ierr;
            Vec U_vec = U_data[ln]->getGlobalVec();
            Vec F_vec = F_data[ln]->getGlobalVec();
            ierr = MatMult(d_dF_dX_mat[ln], U_vec, F_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Spread F(n+1/2) to f(n+1/2).
    d_ib_implicit_integrator->resetAnchorPointValues(F_data, coarsest_ln, finest_ln);
    if (zero_y_before_spread) y.setToScalar(0.0, false);
    lag_data_manager->spread(f_idx, F_data, X_mid_data, true, true);

    t_apply->stop();
    return;
}// apply

void
IBImplicitSJSstarOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    apply(true, x, y);
    return;
}// apply

void
IBImplicitSJSstarOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    t_initialize_operator_state->start();

    if (d_is_initialized) deallocateOperatorState();

    d_dF_dX_mat.resize(in.getPatchHierarchy()->getFinestLevelNumber()+1, static_cast<Mat>(NULL));

    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
IBImplicitSJSstarOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    t_deallocate_operator_state->start();

    for (unsigned k = 0; k < d_dF_dX_mat.size(); ++k)
    {
        if (d_dF_dX_mat[k] != static_cast<Mat>(NULL))
        {
            int ierr = MatDestroy(d_dF_dX_mat[k]);  IBTK_CHKERRQ(ierr);
        }
    }

    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
IBImplicitSJSstarOperator::enableLogging(
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
template class Pointer<IBAMR::IBImplicitSJSstarOperator>;

//////////////////////////////////////////////////////////////////////////////
