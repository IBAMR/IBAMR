// Filename: IBImplicitSFROperator.C
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

#include "IBImplicitSFROperator.h"

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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
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
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitSFROperator::IBImplicitSFROperator(
    IBImplicitHierarchyIntegrator* ib_implicit_integrator)
    : d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_ib_implicit_integrator(ib_implicit_integrator)
{
    // Setup Timers.
    IBAMR_DO_ONCE(
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSFROperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSFROperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitSFROperator::deallocateOperatorState()");
                  );
    return;
}// IBImplicitSFROperator

IBImplicitSFROperator::~IBImplicitSFROperator()
{
    deallocateOperatorState();
    return;
}// ~IBImplicitSFROperator

void
IBImplicitSFROperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    return;
}// setTimeInterval

void
IBImplicitSFROperator::apply(
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
    const int u_new_idx  = x.getComponentDescriptorIndex(0);
    const int f_half_idx = y.getComponentDescriptorIndex(0);
    Pointer<SideVariable<NDIM,double> > f_half_var = y.getComponentVariable(0);
    const int u_current_idx = d_ib_implicit_integrator->d_u_current_idx;

    const int u_half_ib_idx = d_ib_implicit_integrator->d_u_half_ib_idx;
    Pointer<SideVariable<NDIM,double> > u_half_ib_var = d_ib_implicit_integrator->d_u_half_ib_var;

    Pointer<HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops = d_ib_implicit_integrator->d_hier_sc_data_ops;

    LDataManager* l_data_manager = d_ib_implicit_integrator->d_l_data_manager;

    std::vector<Pointer<LData> >& X_data = d_ib_implicit_integrator->d_X_data;
    std::vector<Pointer<LData> >& X_half_data = d_ib_implicit_integrator->d_X_half_data;
    std::vector<Pointer<LData> >& U_half_data = d_ib_implicit_integrator->d_U_half_data;
    std::vector<Pointer<LData> >& F_half_data = d_ib_implicit_integrator->d_F_half_data;

    std::vector<Mat>& R_mats = d_ib_implicit_integrator->d_R_mats;

    Pointer<IBLagrangianForceStrategy> lag_force_strategy = d_ib_implicit_integrator->d_lag_force_strategy;

    // Set u(n+1/2) = 0.5*(u(n) + u(n+1)).
    hier_sc_data_ops->linearSum(u_half_ib_idx, 0.5, u_current_idx, 0.5, u_new_idx);

    // Interpolate u(n+1/2) to U(n+1/2).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec U_half_vec = U_half_data[ln]->getVec();
            Vec u_half_ib_vec = static_cast<Vec>(NULL);
            Pointer<PatchLevel<NDIM> > patch_level = hierarchy->getPatchLevel(ln);
            PETScVecUtilities::constructPatchLevelVec(u_half_ib_vec, u_half_ib_idx, u_half_ib_var, patch_level);
            PETScVecUtilities::copyToPatchLevelVec(u_half_ib_vec, u_half_ib_idx, u_half_ib_var, patch_level);
            ierr = MatMult(R_mats[ln], u_half_ib_vec, U_half_vec); IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(u_half_ib_vec); IBTK_CHKERRQ(ierr);
        }
    }

    // Set X(n+1/2) = X(n) + 0.5*dt*U(n+1/2).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec X_vec = X_data[ln]->getVec();
            Vec X_half_vec = X_half_data[ln]->getVec();
            Vec U_half_vec = U_half_data[ln]->getVec();
            ierr = VecWAXPY(X_half_vec, 0.5*d_dt, U_half_vec, X_vec);  IBTK_CHKERRQ(ierr);
        }
    }

    // Compute F(n+1/2) = F(X(n+1/2),U(n+1/2),t_{n+1/2}).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (l_data_manager->levelContainsLagrangianData(ln))
        {
            Vec F_half_vec = F_half_data[ln]->getVec();
            ierr = VecSet(F_half_vec, 0.0);  IBTK_CHKERRQ(ierr);
            lag_force_strategy->computeLagrangianForce(F_half_data[ln], X_half_data[ln], U_half_data[ln], hierarchy, ln, d_current_time+0.5*d_dt, l_data_manager);
            Pointer<PatchLevel<NDIM> > patch_level = hierarchy->getPatchLevel(ln);
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_level->getGridGeometry();
            const double* const dx0 = grid_geom->getDx();
            const IntVector<NDIM>& ratio = patch_level->getRatio();
            double cell_vol = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                cell_vol *= dx0[d] / double(ratio(d));
            }
            ierr = VecScale(F_half_vec, -1.0/cell_vol);  IBTK_CHKERRQ(ierr);
        }
    }

    // Spread F(n+1/2) to f(n+1/2).
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (l_data_manager->levelContainsLagrangianData(ln))
        {
            Pointer<PatchLevel<NDIM> > patch_level = hierarchy->getPatchLevel(ln);
            Vec F_half_vec = F_half_data[ln]->getVec();
            Vec f_half_vec = static_cast<Vec>(NULL);
            PETScVecUtilities::constructPatchLevelVec(f_half_vec, f_half_idx, f_half_var, patch_level);
            if (zero_y_before_spread)
            {
                ierr = MatMultTranspose(R_mats[ln], F_half_vec, f_half_vec); IBTK_CHKERRQ(ierr);
            }
            else
            {
                PETScVecUtilities::copyToPatchLevelVec(f_half_vec, f_half_idx, f_half_var, patch_level);
                ierr = MatMultTransposeAdd(R_mats[ln], F_half_vec, f_half_vec, f_half_vec); IBTK_CHKERRQ(ierr);
            }
            PETScVecUtilities::copyFromPatchLevelVec(f_half_vec, f_half_idx, f_half_var, patch_level);
            ierr = VecDestroy(f_half_vec); IBTK_CHKERRQ(ierr);
        }
    }

    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent; // XXXX
    SynchronizationTransactionComponent f_half_synch_transaction = SynchronizationTransactionComponent(f_half_idx, "CONSERVATIVE_COARSEN");
    Pointer<SideDataSynchronization> side_synch_op = new SideDataSynchronization();
    side_synch_op->initializeOperatorState(f_half_synch_transaction, y.getPatchHierarchy());
    side_synch_op->synchronizeData(0.0);

    t_apply->stop();
    return;
}// apply

void
IBImplicitSFROperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    apply(true, x, y);
    return;
}// apply

void
IBImplicitSFROperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    t_initialize_operator_state->start();

    if (d_is_initialized) deallocateOperatorState();

    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
IBImplicitSFROperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    t_deallocate_operator_state->start();

    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
IBImplicitSFROperator::enableLogging(
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
template class Pointer<IBAMR::IBImplicitSFROperator>;

//////////////////////////////////////////////////////////////////////////////
