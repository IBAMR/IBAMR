// Filename: INSCollocatedPPMConvectiveOperator.C
// Created on 24 Aug 2011 by Boyce Griffith
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

#include "INSCollocatedPPMConvectiveOperator.h"

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

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_DERIVATIVE_FC FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define CONVECT_DERIVATIVE_FC FC_FUNC_(convect_derivative2d, CONVECT_DERIVATIVE2D)
#define GODUNOV_EXTRAPOLATE_FC FC_FUNC_(godunov_extrapolate2d, GODUNOV_EXTRAPOLATE2D)
#define SKEW_SYM_DERIVATIVE_FC FC_FUNC_(skew_sym_derivative2d, SKEW_SYM_DERIVATIVE2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_FC FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define CONVECT_DERIVATIVE_FC FC_FUNC_(convect_derivative3d, CONVECT_DERIVATIVE3D)
#define GODUNOV_EXTRAPOLATE_FC FC_FUNC_(godunov_extrapolate3d, GODUNOV_EXTRAPOLATE3D)
#define SKEW_SYM_DERIVATIVE_FC FC_FUNC_(skew_sym_derivative3d, SKEW_SYM_DERIVATIVE3D)
#endif

extern "C"
{
    void
    ADVECT_DERIVATIVE_FC(
        const double*,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const int& , const int& ,
        double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double*
#endif
                         );

    void
    CONVECT_DERIVATIVE_FC(
        const double*,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const int& , const int& ,
        double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double*
#endif
                          );

    void
    GODUNOV_EXTRAPOLATE_FC(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , double* , double* , double* , double* ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , double* , double* , double* , double* , double* ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        double* , double* , double*
#endif
                           );

    void
    SKEW_SYM_DERIVATIVE_FC(
        const double*,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const int& , const int& ,
        double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double*
#endif
                         );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// NOTE: The number of ghost cells required by the Godunov advection scheme
// depends on the order of the reconstruction.  These values were chosen to work
// with xsPPM7 (the modified piecewise parabolic method of Rider, Greenough, and
// Kamm).
static const int GADVECTG = 4;

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSCollocatedPPMConvectiveOperator::INSCollocatedPPMConvectiveOperator(
    const ConvectiveDifferencingType difference_form)
    : ConvectiveOperator(difference_form),
      d_is_initialized(false),
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
        TBOX_ERROR("INSCollocatedPPMConvectiveOperator::INSCollocatedPPMConvectiveOperator():\n"
                   << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form) << " \n"
                   << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSCollocatedPPMConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSCollocatedPPMConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var.isNull())
    {
        d_U_var = new CellVariable<NDIM,double>(U_var_name);
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
        t_apply_convective_operator = TimerManager::getManager()->getTimer("IBAMR::INSCollocatedPPMConvectiveOperator::applyConvectiveOperator()");
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::INSCollocatedPPMConvectiveOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSCollocatedPPMConvectiveOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSCollocatedPPMConvectiveOperator::deallocateOperatorState()");
                  );
    return;
}// INSCollocatedPPMConvectiveOperator

INSCollocatedPPMConvectiveOperator::~INSCollocatedPPMConvectiveOperator()
{
    deallocateOperatorState();
    return;
}// ~INSCollocatedPPMConvectiveOperator

void
INSCollocatedPPMConvectiveOperator::applyConvectiveOperator(
    const int U_idx,
    const int N_idx)
{
    TBOX_ERROR("incomplete implementation.\n");
    return;
}// applyConvectiveOperator

void
INSCollocatedPPMConvectiveOperator::initializeOperatorState(
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
    d_refine_alg->registerRefine(d_U_scratch_idx, in.getComponentDescriptorIndex(0), d_U_scratch_idx, d_refine_op);
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
INSCollocatedPPMConvectiveOperator::deallocateOperatorState()
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
INSCollocatedPPMConvectiveOperator::enableLogging(
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
