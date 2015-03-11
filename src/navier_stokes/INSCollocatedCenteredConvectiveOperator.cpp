// Filename: INSCollocatedCenteredConvectiveOperator.cpp
// Created on 24 Aug 2011 by Boyce Griffith
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

#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "IBAMR_config.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSCollocatedCenteredConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartExtrapPhysBdryOp.h"
#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{

class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux2d, ADVECT_FLUX2D)
#define C_TO_F_CWISE_INTERP_2ND_FC IBAMR_FC_FUNC_(ctofcwiseinterp2nd2d, CTOFINTERP2ND2D)
#define F_TO_C_DIV_FC IBAMR_FC_FUNC_(ftocdiv2d, FTOCDIV2D)
#define F_TO_C_DIV_ADD_FC IBAMR_FC_FUNC_(ftocdivadd2d, FTOCDIVADD2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux3d, ADVECT_FLUX3D)
#define C_TO_F_CWISE_INTERP_2ND_FC IBAMR_FC_FUNC_(ctofcwiseinterp2nd3d, CTOFINTERP2ND3D)
#define F_TO_C_DIV_FC IBAMR_FC_FUNC_(ftocdiv3d, FTOCDIV3D)
#define F_TO_C_DIV_ADD_FC IBAMR_FC_FUNC_(ftocdivadd3d, FTOCDIVADD3D)
#endif

extern "C" {
void ADVECT_DERIVATIVE_FC(const double*,
#if (NDIM == 2)
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const double*,
                          const double*,
                          const double*,
                          const double*,
                          const int&,
                          const int&,
#endif
#if (NDIM == 3)
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
                          const double*,
                          const double*,
                          const double*,
                          const double*,
                          const double*,
                          const double*,
                          const int&,
                          const int&,
                          const int&,
#endif
                          double*);

void ADVECT_FLUX_FC(const double&,
#if (NDIM == 2)
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const double*,
                    const double*,
                    const double*,
                    const double*,
                    double*,
                    double*
#endif
#if (NDIM == 3)
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const double*,
                    const double*,
                    const double*,
                    const double*,
                    const double*,
                    const double*,
                    double*,
                    double*,
                    double*
#endif
                    );

void C_TO_F_CWISE_INTERP_2ND_FC(
#if (NDIM == 2)
    double*,
    double*,
    const int&,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&
#endif
#if (NDIM == 3)
    double*,
    double*,
    double*,
    const int&,
    const double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&
#endif
    );

void F_TO_C_DIV_FC(double*,
                   const int&,
                   const double&,
#if (NDIM == 2)
                   const double*,
                   const double*,
                   const int&,
                   const int&,
                   const int&,
                   const int&,
                   const int&,
#endif
#if (NDIM == 3)
                   const double*,
                   const double*,
                   const double*,
                   const int&,
                   const int&,
                   const int&,
                   const int&,
                   const int&,
                   const int&,
                   const int&,
#endif
                   const double*);

void F_TO_C_DIV_ADD_FC(double*,
                       const int&,
                       const double&,
#if (NDIM == 2)
                       const double*,
                       const double*,
                       const int&,
                       const double&,
                       const double*,
                       const int&,
                       const int&,
                       const int&,
                       const int&,
                       const int&,
#endif
#if (NDIM == 3)
                       const double*,
                       const double*,
                       const double*,
                       const int&,
                       const double&,
                       const double*,
                       const int&,
                       const int&,
                       const int&,
                       const int&,
                       const int&,
                       const int&,
                       const int&,
#endif
                       const double*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int GADVECTG = 1;

// Timers.
static boost::shared_ptr<Timer> t_apply_convective_operator;
static boost::shared_ptr<Timer> t_apply;
static boost::shared_ptr<Timer> t_initialize_operator_state;
static boost::shared_ptr<Timer> t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSCollocatedCenteredConvectiveOperator::INSCollocatedCenteredConvectiveOperator(
    const std::string& object_name,
    boost::shared_ptr<Database> input_db,
    const ConvectiveDifferencingType difference_form,
    const std::vector<RobinBcCoefStrategy*>& /*bc_coefs*/)
    : ConvectiveOperator(object_name, difference_form), d_ghostfill_alg(NULL), d_ghostfill_scheds(),
      d_bdry_extrap_type("CONSTANT"), d_hierarchy(NULL), d_coarsest_ln(-1), d_finest_ln(-1), d_U_var(NULL),
      d_U_scratch_idx(-1), d_u_extrap_var(NULL), d_u_flux_var(NULL), d_u_extrap_idx(-1), d_u_flux_idx(-1)
{
    if (d_difference_form != ADVECTIVE && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC)
    {
        TBOX_ERROR(
            "INSCollocatedCenteredConvectiveOperator::INSCollocatedCenteredConvectiveOperator("
            "):\n"
            << "  unsupported differencing form: " << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
            << " \n"
            << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
    }

    auto var_db = VariableDatabase::getDatabase();
    auto context = var_db->getContext("INSCollocatedCenteredConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSCollocatedCenteredConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var)
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
    else
    {
        d_U_var = boost::make_shared<CellVariable<double> >(DIM, U_var_name, NDIM);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector(DIM, GADVECTG));
    }
    TBOX_ASSERT(d_U_scratch_idx >= 0);
    const std::string u_extrap_var_name = "INSCollocatedCenteredConvectiveOperator::u_extrap";
    d_u_extrap_var = var_db->getVariable(u_extrap_var_name);
    if (d_u_extrap_var)
    {
        d_u_extrap_idx = var_db->mapVariableAndContextToIndex(d_u_extrap_var, context);
    }
    else
    {
        d_u_extrap_var = boost::make_shared<FaceVariable<double> >(DIM, u_extrap_var_name, NDIM);
        d_u_extrap_idx = var_db->registerVariableAndContext(d_u_extrap_var, context, IntVector::getZero(DIM));
    }
    TBOX_ASSERT(d_u_extrap_idx >= 0);
    const std::string u_flux_var_name = "INSCollocatedCenteredConvectiveOperator::u_flux";
    d_u_flux_var = var_db->getVariable(u_flux_var_name);
    if (d_u_flux_var)
    {
        d_u_flux_idx = var_db->mapVariableAndContextToIndex(d_u_flux_var, context);
    }
    else
    {
        d_u_flux_var = boost::make_shared<FaceVariable<double> >(DIM, u_flux_var_name, NDIM);
        d_u_flux_idx = var_db->registerVariableAndContext(d_u_flux_var, context, IntVector::getZero(DIM));
    }
    TBOX_ASSERT(d_u_flux_idx >= 0);

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator = TimerManager::getManager()->getTimer(
                      "IBAMR::INSCollocatedCenteredConvectiveOperator::applyConvectiveOperator()");
                  t_apply =
                      TimerManager::getManager()->getTimer("IBAMR::INSCollocatedCenteredConvectiveOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::INSCollocatedCenteredConvectiveOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::INSCollocatedCenteredConvectiveOperator::deallocateOperatorState()"););
    return;
}

INSCollocatedCenteredConvectiveOperator::~INSCollocatedCenteredConvectiveOperator()
{
    deallocateOperatorState();
    return;
}

void INSCollocatedCenteredConvectiveOperator::applyConvectiveOperator(const int U_idx, const int N_idx)
{
    IBAMR_TIMER_START(t_apply_convective_operator);
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSCollocatedCenteredConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }

    // Setup communications algorithm.
    auto grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
    auto refine_alg = boost::make_shared<RefineAlgorithm>();
    auto refine_op = grid_geom->lookupRefineOperator(d_U_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg->registerRefine(d_U_scratch_idx, U_idx, d_U_scratch_idx, refine_op);

    // Extrapolate from cell centers to cell faces.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        refine_alg->resetSchedule(d_ghostfill_scheds[ln]);
        d_ghostfill_scheds[ln]->fillData(d_solution_time);
        d_ghostfill_alg->resetSchedule(d_ghostfill_scheds[ln]);
        auto level = d_hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;

            const Box& patch_box = patch->getBox();
            const IntVector& patch_lower = patch_box.lower();
            const IntVector& patch_upper = patch_box.upper();

            auto U_data = BOOST_CAST<CellData<double> >(patch->getPatchData(d_U_scratch_idx));
            const IntVector& U_data_gcw = U_data->getGhostCellWidth();
            TBOX_ASSERT(U_data_gcw.min() == U_data_gcw.max());
            auto u_ADV_data = BOOST_CAST<FaceData<double> >(patch->getPatchData(d_u_idx));
            const IntVector& u_ADV_data_gcw = u_ADV_data->getGhostCellWidth();
            TBOX_ASSERT(u_ADV_data_gcw.min() == u_ADV_data_gcw.max());
            auto u_extrap_data = BOOST_CAST<FaceData<double> >(patch->getPatchData(d_u_extrap_idx));
            const IntVector& u_extrap_data_gcw = u_extrap_data->getGhostCellWidth();
            TBOX_ASSERT(u_extrap_data_gcw.min() == u_extrap_data_gcw.max());

            // Interpolate from cell centers to cell faces.
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                C_TO_F_CWISE_INTERP_2ND_FC(
#if (NDIM == 2)
                    u_extrap_data->getPointer(0, axis), u_extrap_data->getPointer(1, axis), u_extrap_data_gcw.min(),
                    U_data->getPointer(axis), U_data_gcw.min(), patch_lower(0), patch_upper(0), patch_lower(1),
                    patch_upper(1)
#endif
#if (NDIM == 3)
                        u_extrap_data->getPointer(0, axis),
                    u_extrap_data->getPointer(1, axis), u_extrap_data->getPointer(2, axis), u_extrap_data_gcw.min(),
                    U_data->getPointer(axis), U_data_gcw.min(), patch_lower(0), patch_upper(0), patch_lower(1),
                    patch_upper(1), patch_lower(2), patch_upper(2)
#endif
                        );
            }

            // If we are using conservative or skew-symmetric differencing,
            // compute the advective fluxes.  These need to be synchronized on
            // the patch hierarchy.
            if (d_difference_form == CONSERVATIVE || d_difference_form == SKEW_SYMMETRIC)
            {
                auto u_flux_data = BOOST_CAST<FaceData<double> >(patch->getPatchData(d_u_flux_idx));
                const IntVector& u_flux_data_gcw = u_flux_data->getGhostCellWidth();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    static const double dt = 1.0;
                    ADVECT_FLUX_FC(dt,
#if (NDIM == 2)
                                   patch_lower(0), patch_upper(0), patch_lower(1), patch_upper(1), u_ADV_data_gcw(0),
                                   u_ADV_data_gcw(1), u_extrap_data_gcw(0), u_extrap_data_gcw(1), u_flux_data_gcw(0),
                                   u_flux_data_gcw(1), u_ADV_data->getPointer(0), u_ADV_data->getPointer(1),
                                   u_extrap_data->getPointer(0, axis), u_extrap_data->getPointer(1, axis),
                                   u_flux_data->getPointer(0, axis), u_flux_data->getPointer(1, axis)
#endif
#if (NDIM == 3)
                                                                         patch_lower(0),
                                   patch_upper(0), patch_lower(1), patch_upper(1), patch_lower(2), patch_upper(2),
                                   u_ADV_data_gcw(0), u_ADV_data_gcw(1), u_ADV_data_gcw(2), u_extrap_data_gcw(0),
                                   u_extrap_data_gcw(1), u_extrap_data_gcw(2), u_flux_data_gcw(0), u_flux_data_gcw(1),
                                   u_flux_data_gcw(2), u_ADV_data->getPointer(0), u_ADV_data->getPointer(1),
                                   u_ADV_data->getPointer(2), u_extrap_data->getPointer(0, axis),
                                   u_extrap_data->getPointer(1, axis), u_extrap_data->getPointer(2, axis),
                                   u_flux_data->getPointer(0, axis), u_flux_data->getPointer(1, axis),
                                   u_flux_data->getPointer(2, axis)
#endif
                                       );
                }
            }
        }
    }

    // Synchronize data on the patch hierarchy.
    for (int ln = d_finest_ln; ln > d_coarsest_ln; --ln)
    {
        d_coarsen_scheds[ln]->coarsenData();
    }

    // Difference values on the patches.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;

            const Box& patch_box = patch->getBox();
            const IntVector& patch_lower = patch_box.lower();
            const IntVector& patch_upper = patch_box.upper();

            auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
            const double* const dx = pgeom->getDx();

            auto N_data = BOOST_CAST<CellData<double> >(patch->getPatchData(N_idx));
            const IntVector& N_data_gcw = N_data->getGhostCellWidth();

            if (d_difference_form == ADVECTIVE || d_difference_form == SKEW_SYMMETRIC)
            {
                auto u_ADV_data = BOOST_CAST<FaceData<double> >(patch->getPatchData(d_u_idx));
                const IntVector& u_ADV_data_gcw = u_ADV_data->getGhostCellWidth();
                auto u_extrap_data = BOOST_CAST<FaceData<double> >(patch->getPatchData(d_u_extrap_idx));
                const IntVector& u_extrap_data_gcw = u_extrap_data->getGhostCellWidth();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    ADVECT_DERIVATIVE_FC(
                        dx,
#if (NDIM == 2)
                        patch_lower(0), patch_upper(0), patch_lower(1), patch_upper(1), u_ADV_data_gcw(0),
                        u_ADV_data_gcw(1), u_extrap_data_gcw(0), u_extrap_data_gcw(1), u_ADV_data->getPointer(0),
                        u_ADV_data->getPointer(1), u_extrap_data->getPointer(0, axis),
                        u_extrap_data->getPointer(1, axis), N_data_gcw(0), N_data_gcw(1),
#endif
#if (NDIM == 3)
                        patch_lower(0), patch_upper(0), patch_lower(1), patch_upper(1), patch_lower(2), patch_upper(2),
                        u_ADV_data_gcw(0), u_ADV_data_gcw(1), u_ADV_data_gcw(2), u_extrap_data_gcw(0),
                        u_extrap_data_gcw(1), u_extrap_data_gcw(2), u_ADV_data->getPointer(0),
                        u_ADV_data->getPointer(1), u_ADV_data->getPointer(2), u_extrap_data->getPointer(0, axis),
                        u_extrap_data->getPointer(1, axis), u_extrap_data->getPointer(2, axis), N_data_gcw(0),
                        N_data_gcw(1), N_data_gcw(2),
#endif
                        N_data->getPointer(axis));
                }
            }

            if (d_difference_form == CONSERVATIVE)
            {
                auto u_flux_data = BOOST_CAST<FaceData<double> >(patch->getPatchData(d_u_flux_idx));
                const IntVector& u_flux_data_gcw = u_flux_data->getGhostCellWidth();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    static const double alpha = 1.0;
                    F_TO_C_DIV_FC(N_data->getPointer(axis), N_data_gcw.min(), alpha,
#if (NDIM == 2)
                                  u_flux_data->getPointer(0, axis), u_flux_data->getPointer(1, axis),
                                  u_flux_data_gcw.min(), patch_lower(0), patch_upper(0), patch_lower(1), patch_upper(1),
#endif
#if (NDIM == 3)
                                  u_flux_data->getPointer(0, axis), u_flux_data->getPointer(1, axis),
                                  u_flux_data->getPointer(2, axis), u_flux_data_gcw.min(), patch_lower(0),
                                  patch_upper(0), patch_lower(1), patch_upper(1), patch_lower(2), patch_upper(2),
#endif
                                  dx);
                }
            }

            if (d_difference_form == SKEW_SYMMETRIC)
            {
                auto u_flux_data = BOOST_CAST<FaceData<double> >(patch->getPatchData(d_u_flux_idx));
                const IntVector& u_flux_data_gcw = u_flux_data->getGhostCellWidth();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    static const double alpha = 0.5;
                    static const double beta = 0.5;
                    F_TO_C_DIV_ADD_FC(N_data->getPointer(axis), N_data_gcw.min(), alpha,
#if (NDIM == 2)
                                      u_flux_data->getPointer(0, axis), u_flux_data->getPointer(1, axis),
                                      u_flux_data_gcw.min(), beta, N_data->getPointer(axis), N_data_gcw.min(),
                                      patch_lower(0), patch_upper(0), patch_lower(1), patch_upper(1),
#endif
#if (NDIM == 3)
                                      u_flux_data->getPointer(0, axis), u_flux_data->getPointer(1, axis),
                                      u_flux_data->getPointer(2, axis), u_flux_data_gcw.min(), beta,
                                      N_data->getPointer(axis), N_data_gcw.min(), patch_lower(0), patch_upper(0),
                                      patch_lower(1), patch_upper(1), patch_lower(2), patch_upper(2),
#endif
                                      dx);
                }
            }
        }
    }
    IBAMR_TIMER_STOP(t_apply_convective_operator);
    return;
}

void INSCollocatedCenteredConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<double>& in,
                                                                      const SAMRAIVectorReal<double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    if (d_is_initialized) deallocateOperatorState();

    // Get the hierarchy configuration.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
    auto grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());

    // Setup the coarsen algorithm, operator, and schedules.
    auto coarsen_op = grid_geom->lookupCoarsenOperator(d_u_flux_var, "CONSERVATIVE_COARSEN");
    d_coarsen_alg = boost::make_shared<CoarsenAlgorithm>(DIM);
    if (d_difference_form == ADVECTIVE || d_difference_form == SKEW_SYMMETRIC)
        d_coarsen_alg->registerCoarsen(d_u_extrap_idx, d_u_extrap_idx, coarsen_op);
    if (d_difference_form == CONSERVATIVE || d_difference_form == SKEW_SYMMETRIC)
        d_coarsen_alg->registerCoarsen(d_u_flux_idx, d_u_flux_idx, coarsen_op);
    d_coarsen_scheds.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        auto coarser_level = d_hierarchy->getPatchLevel(ln - 1);
        d_coarsen_scheds[ln] = d_coarsen_alg->createSchedule(coarser_level, level);
    }

    // Setup the refine algorithm, operator, patch strategy, and schedules.
    auto refine_op = grid_geom->lookupRefineOperator(d_U_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ghostfill_alg = boost::make_shared<RefineAlgorithm>();
    d_ghostfill_alg->registerRefine(d_U_scratch_idx, in.getComponentDescriptorIndex(0), d_U_scratch_idx, refine_op);
    d_ghostfill_strategy = boost::make_shared<CartExtrapPhysBdryOp>(d_U_scratch_idx, d_bdry_extrap_type);
    d_ghostfill_scheds.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        d_ghostfill_scheds[ln] = d_ghostfill_alg->createSchedule(level, ln - 1, d_hierarchy, d_ghostfill_strategy);
    }

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_U_scratch_idx))
        {
            level->allocatePatchData(d_U_scratch_idx);
            level->allocatePatchData(d_u_extrap_idx);
            if (d_difference_form == CONSERVATIVE || d_difference_form == SKEW_SYMMETRIC)
                level->allocatePatchData(d_u_flux_idx);
        }
    }
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
}

void INSCollocatedCenteredConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_U_scratch_idx))
        {
            level->deallocatePatchData(d_U_scratch_idx);
        }
        if (level->checkAllocated(d_u_extrap_idx))
        {
            level->deallocatePatchData(d_u_extrap_idx);
        }
        if (level->checkAllocated(d_u_flux_idx))
        {
            level->deallocatePatchData(d_u_flux_idx);
        }
    }

    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_ghostfill_alg.reset();
    d_ghostfill_strategy.reset();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_ghostfill_scheds[ln].reset();
    }
    d_ghostfill_scheds.clear();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
