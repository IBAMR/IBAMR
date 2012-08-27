// Filename: INSStaggeredPPMConvectiveOperator.C
// Created on 08 May 2008 by Boyce Griffith
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

#include "INSStaggeredPPMConvectiveOperator.h"

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
#include <ibtk/CartSideDoubleSpecializedLinearRefine.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <FaceData.h>
#include <FaceGeometry.h>
#include <SideData.h>
#include <SideGeometry.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec.h>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_DERIVATIVE_FC FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define CONVECT_DERIVATIVE_FC FC_FUNC_(convect_derivative2d, CONVECT_DERIVATIVE2D)
#define GODUNOV_EXTRAPOLATE_FC FC_FUNC_(godunov_extrapolate2d, GODUNOV_EXTRAPOLATE2D)
#define NAVIER_STOKES_INTERP_COMPS_FC FC_FUNC_(navier_stokes_interp_comps2d, NAVIER_STOKES_INTERP_COMPS2D)
#define NAVIER_STOKES_RESET_ADV_VELOCITY_FC FC_FUNC_(navier_stokes_reset_adv_velocity2d, NAVIER_STOKES_RESET_ADV_VELOCITY2D)
#define SKEW_SYM_DERIVATIVE_FC FC_FUNC_(skew_sym_derivative2d, SKEW_SYM_DERIVATIVE2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_FC FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define CONVECT_DERIVATIVE_FC FC_FUNC_(convect_derivative3d, CONVECT_DERIVATIVE3D)
#define GODUNOV_EXTRAPOLATE_FC FC_FUNC_(godunov_extrapolate3d, GODUNOV_EXTRAPOLATE3D)
#define NAVIER_STOKES_INTERP_COMPS_FC FC_FUNC_(navier_stokes_interp_comps3d, NAVIER_STOKES_INTERP_COMPS3D)
#define NAVIER_STOKES_RESET_ADV_VELOCITY_FC FC_FUNC_(navier_stokes_reset_adv_velocity3d, NAVIER_STOKES_RESET_ADV_VELOCITY3D)
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
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
#endif
        double*
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
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
#endif
        double*
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
    NAVIER_STOKES_INTERP_COMPS_FC(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double*
#endif
                                  );

    void
    NAVIER_STOKES_RESET_ADV_VELOCITY_FC(
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& ,
        const double* , const double* ,
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        double* , double* ,
        const int& , const int& ,
        const double* , const double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        double* , double* , double* ,
        const int& , const int& , const int& ,
        const double* , const double* , const double*
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
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
#endif
        double*
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

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredPPMConvectiveOperator::INSStaggeredPPMConvectiveOperator(
    const std::string& object_name,
    const ConvectiveDifferencingType difference_form,
    const std::string& bdry_extrap_type)
    : ConvectiveOperator(object_name, difference_form),
      d_ghostfill_alg(NULL),
      d_ghostfill_scheds(),
      d_bdry_extrap_type(bdry_extrap_type),
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
        TBOX_ERROR("INSStaggeredPPMConvectiveOperator::INSStaggeredPPMConvectiveOperator():\n"
                   << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form) << " \n"
                   << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSStaggeredPPMConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSStaggeredPPMConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (!d_U_var)
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
        t_apply_convective_operator = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredPPMConvectiveOperator::applyConvectiveOperator()");
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredPPMConvectiveOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredPPMConvectiveOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::INSStaggeredPPMConvectiveOperator::deallocateOperatorState()");
                  );
    return;
}// INSStaggeredPPMConvectiveOperator

INSStaggeredPPMConvectiveOperator::~INSStaggeredPPMConvectiveOperator()
{
    deallocateOperatorState();
    return;
}// ~INSStaggeredPPMConvectiveOperator

void
INSStaggeredPPMConvectiveOperator::applyConvectiveOperator(
    const int U_idx,
    const int N_idx)
{
    IBAMR_TIMER_START(t_apply_convective_operator);
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSStaggeredPPMConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }
    TBOX_ASSERT(U_idx == d_u_idx);
#endif

    // Setup communications algorithm.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(d_U_var, "SPECIALIZED_LINEAR_REFINE");
    refine_alg->registerRefine(d_U_scratch_idx, U_idx, d_U_scratch_idx, refine_op);

    // Compute the convective derivative.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        refine_alg->resetSchedule(d_ghostfill_scheds[ln]);
        d_ghostfill_scheds[ln]->fillData(d_solution_time);
        d_ghostfill_alg->resetSchedule(d_ghostfill_scheds[ln]);
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

            const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
            blitz::TinyVector<Box<NDIM>,NDIM> side_boxes;
            blitz::TinyVector<Pointer<FaceData<NDIM,double> >,NDIM>  U_adv_data;
            blitz::TinyVector<Pointer<FaceData<NDIM,double> >,NDIM> U_half_data;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                side_boxes [axis] = SideGeometry<NDIM>::toSideBox(patch_box,axis);
                U_adv_data [axis] = new FaceData<NDIM,double>(side_boxes[axis],1,ghosts);
                U_half_data[axis] = new FaceData<NDIM,double>(side_boxes[axis],1,ghosts);
            }
#if (NDIM == 2)
            NAVIER_STOKES_INTERP_COMPS_FC(
                patch_lower(0), patch_upper(0),
                patch_lower(1), patch_upper(1),
                U_data->getGhostCellWidth()(0),         U_data->getGhostCellWidth()(1),
                U_data->getPointer(0),                  U_data->getPointer(1),
                side_boxes[0].lower(0),                 side_boxes[0].upper(0),
                side_boxes[0].lower(1),                 side_boxes[0].upper(1),
                U_adv_data[0]->getGhostCellWidth()(0),  U_adv_data[0]->getGhostCellWidth()(1),
                U_adv_data[0]->getPointer(0),           U_adv_data[0]->getPointer(1),
                side_boxes[1].lower(0),                 side_boxes[1].upper(0),
                side_boxes[1].lower(1),                 side_boxes[1].upper(1),
                U_adv_data[1]->getGhostCellWidth()(0),  U_adv_data[1]->getGhostCellWidth()(1),
                U_adv_data[1]->getPointer(0),           U_adv_data[1]->getPointer(1));
#endif
#if (NDIM == 3)
            NAVIER_STOKES_INTERP_COMPS_FC(
                patch_lower(0), patch_upper(0),
                patch_lower(1), patch_upper(1),
                patch_lower(2), patch_upper(2),
                U_data->getGhostCellWidth()(0),         U_data->getGhostCellWidth()(1),         U_data->getGhostCellWidth()(2),
                U_data->getPointer(0),                  U_data->getPointer(1),                  U_data->getPointer(2),
                side_boxes[0].lower(0),                 side_boxes[0].upper(0),
                side_boxes[0].lower(1),                 side_boxes[0].upper(1),
                side_boxes[0].lower(2),                 side_boxes[0].upper(2),
                U_adv_data[0]->getGhostCellWidth()(0),  U_adv_data[0]->getGhostCellWidth()(1),  U_adv_data[0]->getGhostCellWidth()(2),
                U_adv_data[0]->getPointer(0),           U_adv_data[0]->getPointer(1),           U_adv_data[0]->getPointer(2),
                side_boxes[1].lower(0),                 side_boxes[1].upper(0),
                side_boxes[1].lower(1),                 side_boxes[1].upper(1),
                side_boxes[1].lower(2),                 side_boxes[1].upper(2),
                U_adv_data[1]->getGhostCellWidth()(0),  U_adv_data[1]->getGhostCellWidth()(1),  U_adv_data[1]->getGhostCellWidth()(2),
                U_adv_data[1]->getPointer(0),           U_adv_data[1]->getPointer(1),           U_adv_data[1]->getPointer(2),
                side_boxes[2].lower(0),                 side_boxes[2].upper(0),
                side_boxes[2].lower(1),                 side_boxes[2].upper(1),
                side_boxes[2].lower(2),                 side_boxes[2].upper(2),
                U_adv_data[2]->getGhostCellWidth()(0),  U_adv_data[2]->getGhostCellWidth()(1),  U_adv_data[2]->getGhostCellWidth()(2),
                U_adv_data[2]->getPointer(0),           U_adv_data[2]->getPointer(1),           U_adv_data[2]->getPointer(2));
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                Pointer<SideData<NDIM,double> > dU_data =
                    new SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                Pointer<SideData<NDIM,double> > U_L_data =
                    new SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                Pointer<SideData<NDIM,double> > U_R_data =
                    new SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
                Pointer<SideData<NDIM,double> > U_scratch1_data =
                    new SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
#if (NDIM == 3)
                Pointer<SideData<NDIM,double> > U_scratch2_data =
                    new SideData<NDIM,double>(U_data->getBox(), U_data->getDepth(), U_data->getGhostCellWidth());
#endif
#if (NDIM == 2)
                GODUNOV_EXTRAPOLATE_FC(
                    side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                    side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                    U_data->getGhostCellWidth()(0), U_data->getGhostCellWidth()(1),
                    U_data       ->getPointer(axis),       U_scratch1_data->getPointer(axis),
                    dU_data      ->getPointer(axis),       U_L_data       ->getPointer(axis),       U_R_data->getPointer(axis),
                    U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1),
                    U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1),
                    U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),
                    U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1));
#endif
#if (NDIM == 3)
                GODUNOV_EXTRAPOLATE_FC(
                    side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                    side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                    side_boxes[axis].lower(2), side_boxes[axis].upper(2),
                    U_data->getGhostCellWidth()(0), U_data->getGhostCellWidth()(1), U_data->getGhostCellWidth()(2),
                    U_data       ->getPointer(axis),       U_scratch1_data->getPointer(axis),       U_scratch2_data->getPointer(axis),
                    dU_data      ->getPointer(axis),       U_L_data       ->getPointer(axis),       U_R_data       ->getPointer(axis),
                    U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1), U_adv_data [axis]->getGhostCellWidth()(2),
                    U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1), U_half_data[axis]->getGhostCellWidth()(2),
                    U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),          U_adv_data [axis]->getPointer(2),
                    U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),          U_half_data[axis]->getPointer(2));
#endif
            }
#if (NDIM == 2)
            NAVIER_STOKES_RESET_ADV_VELOCITY_FC(
                side_boxes[0].lower(0), side_boxes[0].upper(0),
                side_boxes[0].lower(1), side_boxes[0].upper(1),
                U_adv_data [0]->getGhostCellWidth()(0), U_adv_data [0]->getGhostCellWidth()(1),
                U_adv_data [0]->getPointer(0),          U_adv_data [0]->getPointer(1),
                U_half_data[0]->getGhostCellWidth()(0), U_half_data[0]->getGhostCellWidth()(1),
                U_half_data[0]->getPointer(0),          U_half_data[0]->getPointer(1),
                side_boxes[1].lower(0), side_boxes[1].upper(0),
                side_boxes[1].lower(1), side_boxes[1].upper(1),
                U_adv_data [1]->getGhostCellWidth()(0), U_adv_data [1]->getGhostCellWidth()(1),
                U_adv_data [1]->getPointer(0),          U_adv_data [1]->getPointer(1),
                U_half_data[1]->getGhostCellWidth()(0), U_half_data[1]->getGhostCellWidth()(1),
                U_half_data[1]->getPointer(0),          U_half_data[1]->getPointer(1));
#endif
#if (NDIM == 3)
            NAVIER_STOKES_RESET_ADV_VELOCITY_FC(
                side_boxes[0].lower(0), side_boxes[0].upper(0),
                side_boxes[0].lower(1), side_boxes[0].upper(1),
                side_boxes[0].lower(2), side_boxes[0].upper(2),
                U_adv_data [0]->getGhostCellWidth()(0), U_adv_data [0]->getGhostCellWidth()(1), U_adv_data [0]->getGhostCellWidth()(2),
                U_adv_data [0]->getPointer(0),          U_adv_data [0]->getPointer(1),          U_adv_data [0]->getPointer(2),
                U_half_data[0]->getGhostCellWidth()(0), U_half_data[0]->getGhostCellWidth()(1), U_half_data[0]->getGhostCellWidth()(2),
                U_half_data[0]->getPointer(0),          U_half_data[0]->getPointer(1),          U_half_data[0]->getPointer(2),
                side_boxes[1].lower(0), side_boxes[1].upper(0),
                side_boxes[1].lower(1), side_boxes[1].upper(1),
                side_boxes[1].lower(2), side_boxes[1].upper(2),
                U_adv_data [1]->getGhostCellWidth()(0), U_adv_data [1]->getGhostCellWidth()(1), U_adv_data [1]->getGhostCellWidth()(2),
                U_adv_data [1]->getPointer(0),          U_adv_data [1]->getPointer(1),          U_adv_data [1]->getPointer(2),
                U_half_data[1]->getGhostCellWidth()(0), U_half_data[1]->getGhostCellWidth()(1), U_half_data[1]->getGhostCellWidth()(2),
                U_half_data[1]->getPointer(0),          U_half_data[1]->getPointer(1),          U_half_data[1]->getPointer(2),
                side_boxes[2].lower(0), side_boxes[2].upper(0),
                side_boxes[2].lower(1), side_boxes[2].upper(1),
                side_boxes[2].lower(2), side_boxes[2].upper(2),
                U_adv_data [2]->getGhostCellWidth()(0), U_adv_data [2]->getGhostCellWidth()(1), U_adv_data [2]->getGhostCellWidth()(2),
                U_adv_data [2]->getPointer(0),          U_adv_data [2]->getPointer(1),          U_adv_data [2]->getPointer(2),
                U_half_data[2]->getGhostCellWidth()(0), U_half_data[2]->getGhostCellWidth()(1), U_half_data[2]->getGhostCellWidth()(2),
                U_half_data[2]->getPointer(0),          U_half_data[2]->getPointer(1),          U_half_data[2]->getPointer(2));
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                switch (d_difference_form)
                {
                    case CONSERVATIVE:
#if (NDIM == 2)
                        CONVECT_DERIVATIVE_FC(
                            dx,
                            side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                            side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                            U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1),
                            U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1),
                            U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),
                            U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),
                            N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1),
                            N_data->getPointer(axis));
#endif
#if (NDIM == 3)
                        CONVECT_DERIVATIVE_FC(
                            dx,
                            side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                            side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                            side_boxes[axis].lower(2), side_boxes[axis].upper(2),
                            U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1), U_adv_data [axis]->getGhostCellWidth()(2),
                            U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1), U_half_data[axis]->getGhostCellWidth()(2),
                            U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),          U_adv_data [axis]->getPointer(2),
                            U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),          U_half_data[axis]->getPointer(2),
                            N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1), N_data->getGhostCellWidth()(2),
                            N_data->getPointer(axis));
#endif
                        break;
                    case ADVECTIVE:
#if (NDIM == 2)
                        ADVECT_DERIVATIVE_FC(
                            dx,
                            side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                            side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                            U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1),
                            U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1),
                            U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),
                            U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),
                            N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1),
                            N_data->getPointer(axis));
#endif
#if (NDIM == 3)
                        ADVECT_DERIVATIVE_FC(
                            dx,
                            side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                            side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                            side_boxes[axis].lower(2), side_boxes[axis].upper(2),
                            U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1), U_adv_data [axis]->getGhostCellWidth()(2),
                            U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1), U_half_data[axis]->getGhostCellWidth()(2),
                            U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),          U_adv_data [axis]->getPointer(2),
                            U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),          U_half_data[axis]->getPointer(2),
                            N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1), N_data->getGhostCellWidth()(2),
                            N_data->getPointer(axis));
#endif
                        break;
                    case SKEW_SYMMETRIC:
#if (NDIM == 2)
                        SKEW_SYM_DERIVATIVE_FC(
                            dx,
                            side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                            side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                            U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1),
                            U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1),
                            U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),
                            U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),
                            N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1),
                            N_data->getPointer(axis));
#endif
#if (NDIM == 3)
                        SKEW_SYM_DERIVATIVE_FC(
                            dx,
                            side_boxes[axis].lower(0), side_boxes[axis].upper(0),
                            side_boxes[axis].lower(1), side_boxes[axis].upper(1),
                            side_boxes[axis].lower(2), side_boxes[axis].upper(2),
                            U_adv_data [axis]->getGhostCellWidth()(0), U_adv_data [axis]->getGhostCellWidth()(1), U_adv_data [axis]->getGhostCellWidth()(2),
                            U_half_data[axis]->getGhostCellWidth()(0), U_half_data[axis]->getGhostCellWidth()(1), U_half_data[axis]->getGhostCellWidth()(2),
                            U_adv_data [axis]->getPointer(0),          U_adv_data [axis]->getPointer(1),          U_adv_data [axis]->getPointer(2),
                            U_half_data[axis]->getPointer(0),          U_half_data[axis]->getPointer(1),          U_half_data[axis]->getPointer(2),
                            N_data->getGhostCellWidth()(0), N_data->getGhostCellWidth()(1), N_data->getGhostCellWidth()(2),
                            N_data->getPointer(axis));
#endif
                        break;
                    default:
                        TBOX_ERROR("INSStaggeredPPMConvectiveOperator::applyConvectiveOperator():\n"
                                   << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form) << " \n"
                                   << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
                }
            }
        }
    }

    IBAMR_TIMER_STOP(t_apply_convective_operator);
    return;
}// applyConvectiveOperator

void
INSStaggeredPPMConvectiveOperator::initializeOperatorState(
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
    grid_geom->addSpatialRefineOperator(new CartSideDoubleSpecializedLinearRefine());
    Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(d_U_var, "SPECIALIZED_LINEAR_REFINE");
    d_ghostfill_alg = new RefineAlgorithm<NDIM>();
    d_ghostfill_alg->registerRefine(d_U_scratch_idx, in.getComponentDescriptorIndex(0), d_U_scratch_idx, refine_op);
    d_ghostfill_strategy = new CartExtrapPhysBdryOp(d_U_scratch_idx, d_bdry_extrap_type);
    d_ghostfill_scheds.resize(d_finest_ln+1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_ghostfill_scheds[ln] = d_ghostfill_alg->createSchedule(level, ln-1, d_hierarchy, d_ghostfill_strategy);
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
INSStaggeredPPMConvectiveOperator::deallocateOperatorState()
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
    d_ghostfill_alg.setNull();
    d_ghostfill_strategy.setNull();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_ghostfill_scheds[ln].setNull();
    }
    d_ghostfill_scheds.clear();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
}// deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
