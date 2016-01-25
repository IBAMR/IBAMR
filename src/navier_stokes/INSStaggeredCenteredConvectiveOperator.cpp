// Filename: INSStaggeredCenteredConvectiveOperator.cpp
// Created on 30 Oct 2008 by Boyce Griffith
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

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSStaggeredCenteredConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE_FC                                                                      \
    IBAMR_FC_FUNC_(navier_stokes_staggered_adv_derivative2d, NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE2D)
#define NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE_FC                                                                      \
    IBAMR_FC_FUNC_(navier_stokes_staggered_div_derivative2d, NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE2D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE_FC                                                                 \
    IBAMR_FC_FUNC_(navier_stokes_staggered_skew_sym_derivative2d, NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE2D)
#endif

#if (NDIM == 3)
#define NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE_FC                                                                      \
    IBAMR_FC_FUNC_(navier_stokes_staggered_adv_derivative3d, NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE3D)
#define NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE_FC                                                                      \
    IBAMR_FC_FUNC_(navier_stokes_staggered_div_derivative3d, NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE3D)
#define NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE_FC                                                                 \
    IBAMR_FC_FUNC_(navier_stokes_staggered_skew_sym_derivative3d, NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE3D)
#endif

extern "C" {
void NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE_FC(const double*,
#if (NDIM == 2)
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const double*,
                                               const double*,
                                               const int&,
                                               const int&,
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
                                               const double*,
                                               const double*,
                                               const double*,
                                               const int&,
                                               const int&,
                                               const int&,
                                               double*,
                                               double*,
                                               double*
#endif
                                               );

void NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE_FC(const double*,
#if (NDIM == 2)
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const int&,
                                               const double*,
                                               const double*,
                                               const int&,
                                               const int&,
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
                                               const double*,
                                               const double*,
                                               const double*,
                                               const int&,
                                               const int&,
                                               const int&,
                                               double*,
                                               double*,
                                               double*
#endif
                                               );

void NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE_FC(const double*,
#if (NDIM == 2)
                                                    const int&,
                                                    const int&,
                                                    const int&,
                                                    const int&,
                                                    const int&,
                                                    const int&,
                                                    const double*,
                                                    const double*,
                                                    const int&,
                                                    const int&,
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
                                                    const double*,
                                                    const double*,
                                                    const double*,
                                                    const int&,
                                                    const int&,
                                                    const int&,
                                                    double*,
                                                    double*,
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
// Number of ghosts cells used for each variable quantity.
static const int GADVECTG = 1;

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredCenteredConvectiveOperator::INSStaggeredCenteredConvectiveOperator(
    const std::string& object_name,
    Pointer<Database> input_db,
    const ConvectiveDifferencingType difference_form,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    : ConvectiveOperator(object_name, difference_form),
      d_bc_coefs(bc_coefs),
      d_bdry_extrap_type("CONSTANT"),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_U_var(NULL),
      d_U_scratch_idx(-1)
{
    if (d_difference_form != ADVECTIVE && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC)
    {
        TBOX_ERROR(
            "INSStaggeredCenteredConvectiveOperator::INSStaggeredCenteredConvectiveOperator():"
            "\n"
            << "  unsupported differencing form: "
            << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
            << " \n"
            << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSStaggeredCenteredConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSStaggeredCenteredConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var)
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
    else
    {
        d_U_var = new SideVariable<NDIM, double>(U_var_name);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector<NDIM>(GADVECTG));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator = TimerManager::getManager()->getTimer(
                      "IBAMR::INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator()");
                  t_apply =
                      TimerManager::getManager()->getTimer("IBAMR::INSStaggeredCenteredConvectiveOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::INSStaggeredCenteredConvectiveOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::INSStaggeredCenteredConvectiveOperator::deallocateOperatorState()"););
    return;
} // INSStaggeredCenteredConvectiveOperator

INSStaggeredCenteredConvectiveOperator::~INSStaggeredCenteredConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~INSStaggeredCenteredConvectiveOperator

void
INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator(const int U_idx, const int N_idx)
{
    IBAMR_TIMER_START(t_apply_convective_operator);
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }
    TBOX_ASSERT(U_idx == d_u_idx);
#endif

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_U_scratch_idx))
        {
            level->allocatePatchData(d_U_scratch_idx);
        }
    }

    // Fill ghost cell values for all components.
    static const bool homogeneous_bc = false;
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(1);
    transaction_comps[0] = InterpolationTransactionComponent(d_U_scratch_idx,
                                                             U_idx,
                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                             false,
                                                             "CONSERVATIVE_COARSEN",
                                                             d_bdry_extrap_type,
                                                             false,
                                                             d_bc_coefs);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_bc_coefs, NULL, d_U_scratch_idx, -1, homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_bc_coefs, NULL);
    //  d_bc_helper->enforceDivergenceFreeConditionAtBoundary(d_U_scratch_idx);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the convective derivative.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            Pointer<SideData<NDIM, double> > N_data = patch->getPatchData(N_idx);
            Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(d_U_scratch_idx);

            const IntVector<NDIM>& N_data_gcw = N_data->getGhostCellWidth();
            const IntVector<NDIM>& U_data_gcw = U_data->getGhostCellWidth();

            switch (d_difference_form)
            {
            case CONSERVATIVE:
                NAVIER_STOKES_STAGGERED_DIV_DERIVATIVE_FC(dx,
#if (NDIM == 2)
                                                          patch_lower(0),
                                                          patch_upper(0),
                                                          patch_lower(1),
                                                          patch_upper(1),
                                                          U_data_gcw(0),
                                                          U_data_gcw(1),
                                                          U_data->getPointer(0),
                                                          U_data->getPointer(1),
                                                          N_data_gcw(0),
                                                          N_data_gcw(1),
                                                          N_data->getPointer(0),
                                                          N_data->getPointer(1)
#endif
#if (NDIM == 3)
                                                              patch_lower(0),
                                                          patch_upper(0),
                                                          patch_lower(1),
                                                          patch_upper(1),
                                                          patch_lower(2),
                                                          patch_upper(2),
                                                          U_data_gcw(0),
                                                          U_data_gcw(1),
                                                          U_data_gcw(2),
                                                          U_data->getPointer(0),
                                                          U_data->getPointer(1),
                                                          U_data->getPointer(2),
                                                          N_data_gcw(0),
                                                          N_data_gcw(1),
                                                          N_data_gcw(2),
                                                          N_data->getPointer(0),
                                                          N_data->getPointer(1),
                                                          N_data->getPointer(2)
#endif
                                                              );
                break;
            case ADVECTIVE:
                NAVIER_STOKES_STAGGERED_ADV_DERIVATIVE_FC(dx,
#if (NDIM == 2)
                                                          patch_lower(0),
                                                          patch_upper(0),
                                                          patch_lower(1),
                                                          patch_upper(1),
                                                          U_data_gcw(0),
                                                          U_data_gcw(1),
                                                          U_data->getPointer(0),
                                                          U_data->getPointer(1),
                                                          N_data_gcw(0),
                                                          N_data_gcw(1),
                                                          N_data->getPointer(0),
                                                          N_data->getPointer(1)
#endif
#if (NDIM == 3)
                                                              patch_lower(0),
                                                          patch_upper(0),
                                                          patch_lower(1),
                                                          patch_upper(1),
                                                          patch_lower(2),
                                                          patch_upper(2),
                                                          U_data_gcw(0),
                                                          U_data_gcw(1),
                                                          U_data_gcw(2),
                                                          U_data->getPointer(0),
                                                          U_data->getPointer(1),
                                                          U_data->getPointer(2),
                                                          N_data_gcw(0),
                                                          N_data_gcw(1),
                                                          N_data_gcw(2),
                                                          N_data->getPointer(0),
                                                          N_data->getPointer(1),
                                                          N_data->getPointer(2)
#endif
                                                              );
                break;
            case SKEW_SYMMETRIC:
                NAVIER_STOKES_STAGGERED_SKEW_SYM_DERIVATIVE_FC(dx,
#if (NDIM == 2)
                                                               patch_lower(0),
                                                               patch_upper(0),
                                                               patch_lower(1),
                                                               patch_upper(1),
                                                               U_data_gcw(0),
                                                               U_data_gcw(1),
                                                               U_data->getPointer(0),
                                                               U_data->getPointer(1),
                                                               N_data_gcw(0),
                                                               N_data_gcw(1),
                                                               N_data->getPointer(0),
                                                               N_data->getPointer(1)
#endif
#if (NDIM == 3)
                                                                   patch_lower(0),
                                                               patch_upper(0),
                                                               patch_lower(1),
                                                               patch_upper(1),
                                                               patch_lower(2),
                                                               patch_upper(2),
                                                               U_data_gcw(0),
                                                               U_data_gcw(1),
                                                               U_data_gcw(2),
                                                               U_data->getPointer(0),
                                                               U_data->getPointer(1),
                                                               U_data->getPointer(2),
                                                               N_data_gcw(0),
                                                               N_data_gcw(1),
                                                               N_data_gcw(2),
                                                               N_data->getPointer(0),
                                                               N_data->getPointer(1),
                                                               N_data->getPointer(2)
#endif
                                                                   );
                break;
            default:
                TBOX_ERROR("INSStaggeredCenteredConvectiveOperator::applyConvectiveOperator():\n"
                           << "  unsupported differencing form: "
                           << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
                           << " \n"
                           << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
            }
        }
    }

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_U_scratch_idx))
        {
            level->deallocatePatchData(d_U_scratch_idx);
        }
    }

    IBAMR_TIMER_STOP(t_apply_convective_operator);
    return;
} // applyConvectiveOperator

void
INSStaggeredCenteredConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                                const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    if (d_is_initialized) deallocateOperatorState();

    // Get the hierarchy configuration.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#else
    NULL_USE(out);
#endif

    // Setup the interpolation transaction information.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    d_transaction_comps.resize(1);
    d_transaction_comps[0] = InterpolationTransactionComponent(d_U_scratch_idx,
                                                               in.getComponentDescriptorIndex(0),
                                                               "CONSERVATIVE_LINEAR_REFINE",
                                                               false,
                                                               "CONSERVATIVE_COARSEN",
                                                               d_bdry_extrap_type,
                                                               false,
                                                               d_bc_coefs);

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy);

    // Initialize the BC helper.
    d_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    d_bc_helper->cacheBcCoefData(d_bc_coefs, d_solution_time, d_hierarchy);

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
INSStaggeredCenteredConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_hier_bdry_fill.setNull();
    d_bc_helper.setNull();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
