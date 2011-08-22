// Filename: INSStaggeredCenteredConvectiveOperator.C
// Created on 30 Oct 2008 by Boyce Griffith
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

#include "INSStaggeredCenteredConvectiveOperator.h"

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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartExtrapPhysBdryOp.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <FaceData.h>
#include <FaceGeometry.h>
#include <SideData.h>
#include <SideGeometry.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_adv_derivative2d,NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE2D)
#define NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_div_derivative2d,NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE2D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_skew_sym_derivative2d,NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_adv_derivative3d,NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE3D)
#define NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_div_derivative3d,NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE3D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE_FC FC_FUNC_(navier_stokes_staggered_skew_sym_derivative3d,NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE3D)
#endif

extern "C"
{
    void
    NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE_FC(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double* , double* , double*
#endif
                                              );

    void
    NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE_FC(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double* , double* , double*
#endif
                                              );

    void
    NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE_FC(
        const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double* , double* , double*
#endif
                                                   );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int GADVECTG = 1;

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredCenteredConvectiveOperator::INSStaggeredCenteredConvectiveOperator(
    const ConvectiveDifferencingType difference_form)
    : d_is_initialized(false),
      d_difference_form(difference_form),
      d_refine_alg(NULL),
      d_refine_op(NULL),
      d_refine_scheds(),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_U_var(NULL),
      d_U_scratch_idx(-1)
{
    if (d_difference_form != ADVECTIVE &&
        d_difference_form != CONSERVATIVE &&
        d_difference_form != SKEW_SYMMETRIC)
    {
        TBOX_ERROR("INSStaggeredCenteredConvectiveOperator::INSStaggeredCenteredConvectiveOperator():\n"
                   << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form) << " \n"
                   << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSStaggeredCenteredConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSStaggeredCenteredConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var.isNull())
    {
        d_U_var = new SideVariable<NDIM,double>(U_var_name);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector<NDIM>(GADVECTG));
    }
    else
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_apply_convective_operator = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator()");
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredCenteredConvectiveOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredCenteredConvectiveOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredCenteredConvectiveOperator::deallocateOperatorState()");
                  );
    return;
}// INSStaggeredCenteredConvectiveOperator

INSStaggeredCenteredConvectiveOperator::~INSStaggeredCenteredConvectiveOperator()
{
    deallocateOperatorState();
    return;
}// ~INSStaggeredCenteredConvectiveOperator

void
INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator(
    const int U_idx,
    const int N_idx)
{
    IBAMR_TIMER_START(t_apply_convective_operator);

    if (!d_is_initialized)
    {
        TBOX_ERROR("INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }

    // Setup communications schedules.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
    refine_alg->registerRefine(d_U_scratch_idx, // destination
                               U_idx,           // source
                               d_U_scratch_idx, // temporary work space
                               d_refine_op);

    // Compute the convective derivative.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        refine_alg->resetSchedule(d_refine_scheds[ln]);
        d_refine_scheds[ln]->fillData(0.0);
        d_refine_alg->resetSchedule(d_refine_scheds[ln]);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            Pointer<SideData<NDIM,double> > N_data = patch->getPatchData(N_idx);
            Pointer<SideData<NDIM,double> > U_data = patch->getPatchData(d_U_scratch_idx);

            const IntVector<NDIM>& N_ghosts = N_data->getGhostCellWidth();
            const IntVector<NDIM>& U_ghosts = U_data->getGhostCellWidth();

            switch (d_difference_form)
            {
                case CONSERVATIVE:
                    NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE_FC(
                        dx,
#if (NDIM == 2)
                        patch_lower(0), patch_upper(0),
                        patch_lower(1), patch_upper(1),
                        U_ghosts(0), U_ghosts(1),
                        U_data->getPointer(0), U_data->getPointer(1),
                        N_ghosts(0), N_ghosts(1),
                        N_data->getPointer(0), N_data->getPointer(1)
#endif
#if (NDIM == 3)
                        patch_lower(0), patch_upper(0),
                        patch_lower(1), patch_upper(1),
                        patch_lower(2), patch_upper(2),
                        U_ghosts(0), U_ghosts(1), U_ghosts(2),
                        U_data->getPointer(0), U_data->getPointer(1), U_data->getPointer(2),
                        N_ghosts(0), N_ghosts(1), N_ghosts(2),
                        N_data->getPointer(0), N_data->getPointer(1), N_data->getPointer(2)
#endif
                                                              );
                    break;
                case ADVECTIVE:
                    NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE_FC(
                        dx,
#if (NDIM == 2)
                        patch_lower(0), patch_upper(0),
                        patch_lower(1), patch_upper(1),
                        U_ghosts(0), U_ghosts(1),
                        U_data->getPointer(0), U_data->getPointer(1),
                        N_ghosts(0), N_ghosts(1),
                        N_data->getPointer(0), N_data->getPointer(1)
#endif
#if (NDIM == 3)
                        patch_lower(0), patch_upper(0),
                        patch_lower(1), patch_upper(1),
                        patch_lower(2), patch_upper(2),
                        U_ghosts(0), U_ghosts(1), U_ghosts(2),
                        U_data->getPointer(0), U_data->getPointer(1), U_data->getPointer(2),
                        N_ghosts(0), N_ghosts(1), N_ghosts(2),
                        N_data->getPointer(0), N_data->getPointer(1), N_data->getPointer(2)
#endif
                                                              );
                    break;
                case SKEW_SYMMETRIC:
                    NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE_FC(
                        dx,
#if (NDIM == 2)
                        patch_lower(0), patch_upper(0),
                        patch_lower(1), patch_upper(1),
                        U_ghosts(0), U_ghosts(1),
                        U_data->getPointer(0), U_data->getPointer(1),
                        N_ghosts(0), N_ghosts(1),
                        N_data->getPointer(0), N_data->getPointer(1)
#endif
#if (NDIM == 3)
                        patch_lower(0), patch_upper(0),
                        patch_lower(1), patch_upper(1),
                        patch_lower(2), patch_upper(2),
                        U_ghosts(0), U_ghosts(1), U_ghosts(2),
                        U_data->getPointer(0), U_data->getPointer(1), U_data->getPointer(2),
                        N_ghosts(0), N_ghosts(1), N_ghosts(2),
                        N_data->getPointer(0), N_data->getPointer(1), N_data->getPointer(2)
#endif
                                                                   );
                    break;
                default:
                    TBOX_ERROR("INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator():\n"
                               << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form) << " \n"
                               << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
            }
        }
    }

    IBAMR_TIMER_STOP(t_apply_convective_operator);
    return;
}// applyConvectiveOperator

void
INSStaggeredCenteredConvectiveOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    IBAMR_TIMER_START(t_apply);

    // Get the vector components.
    const int U_idx = x.getComponentDescriptorIndex(0);
    const int N_idx = y.getComponentDescriptorIndex(0);

    // Compute the action of the operator.
    applyConvectiveOperator(U_idx, N_idx);

    IBAMR_TIMER_STOP(t_apply);
    return;
}// apply

void
INSStaggeredCenteredConvectiveOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    if (d_is_initialized) deallocateOperatorState();

    // Get the hierarchy configuration.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#else
    NULL_USE(out);
#endif
    // Setup the refine algorithm, operator, patch strategy, and schedules.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    d_refine_op = grid_geom->lookupRefineOperator(d_U_var, "CONSERVATIVE_LINEAR_REFINE");
    d_refine_alg = new RefineAlgorithm<NDIM>();
    d_refine_alg->registerRefine(d_U_scratch_idx,                   // destination
                                 in.getComponentDescriptorIndex(0), // source
                                 d_U_scratch_idx,                   // temporary work space
                                 d_refine_op);
    d_refine_strategy = new CartExtrapPhysBdryOp(d_U_scratch_idx, BDRY_EXTRAP_TYPE);
    d_refine_scheds.resize(d_finest_ln+1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_refine_scheds[ln] = d_refine_alg->createSchedule(level, ln-1, d_hierarchy, d_refine_strategy);
    }

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_U_scratch_idx))
        {
            level->allocatePatchData(d_U_scratch_idx);
        }
    }
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
}// initializeOperatorState

void
INSStaggeredCenteredConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_U_scratch_idx))
        {
            level->deallocatePatchData(d_U_scratch_idx);
        }
    }

    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_refine_op.setNull();
    d_refine_alg.setNull();
    d_refine_strategy.setNull();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_refine_scheds[ln].setNull();
    }
    d_refine_scheds.clear();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
}// deallocateOperatorState

void
INSStaggeredCenteredConvectiveOperator::enableLogging(
    bool /*enabled*/)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
