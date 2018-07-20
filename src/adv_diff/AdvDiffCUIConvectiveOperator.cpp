// Filename: AdvDiffCUIConvectiveOperator.cpp
// Created on 14 May 2018 by Nishant Nangia
//
// Copyright (c) 2002-2018, Nishant Nangia
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

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/AdvDiffCUIConvectiveOperator.h"
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartExtrapPhysBdryOp.h"
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
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux2d, ADVECT_FLUX2D)
#define F_TO_C_DIV_FC IBAMR_FC_FUNC_(ftocdiv2d, FTOCDIV2D)
#define F_TO_C_DIV_ADD_FC IBAMR_FC_FUNC_(ftocdivadd2d, FTOCDIVADD2D)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate2d, CUI_EXTRAPOLATE2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux3d, ADVECT_FLUX3D)
#define F_TO_C_DIV_FC IBAMR_FC_FUNC_(ftocdiv3d, FTOCDIV3D)
#define F_TO_C_DIV_ADD_FC IBAMR_FC_FUNC_(ftocdivadd3d, FTOCDIVADD3D)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate3d, CUI_EXTRAPOLATE3D)
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

void CUI_EXTRAPOLATE_FC(
#if (NDIM == 2)
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
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
    const double*,
    double*,
    double*,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const int&,
    const double*,
    const double*,
    const double*,
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
// NOTE: The number of ghost cells required by the advection scheme
// These values were chosen to work with CUI (the cubic interpolation
// upwind method of Waterson and Deconinck).
static const int GADVECTG = 2;

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffCUIConvectiveOperator::AdvDiffCUIConvectiveOperator(const std::string& object_name,
                                                           Pointer<CellVariable<NDIM, double> > Q_var,
                                                           Pointer<Database> input_db,
                                                           const ConvectiveDifferencingType difference_form,
                                                           std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs)
    : ConvectiveOperator(object_name, difference_form),
      d_ghostfill_alg(nullptr),
      d_ghostfill_scheds(),
      d_bc_coefs(std::move(bc_coefs)),
      d_outflow_bdry_extrap_type("CONSTANT"),
      d_hierarchy(nullptr),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_Q_var(Q_var),
      d_Q_data_depth(0),
      d_Q_scratch_idx(-1),
      d_q_extrap_var(nullptr),
      d_q_flux_var(nullptr),
      d_q_extrap_idx(-1),
      d_q_flux_idx(-1)
{
    if (d_difference_form != ADVECTIVE && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC)
    {
        TBOX_ERROR("AdvDiffCUIConvectiveOperator::AdvDiffCUIConvectiveOperator():\n"
                   << "  unsupported differencing form: "
                   << enum_to_string<ConvectiveDifferencingType>(d_difference_form)
                   << " \n"
                   << "  valid choices are: ADVECTIVE, CONSERVATIVE, SKEW_SYMMETRIC\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("outflow_bdry_extrap_type"))
            d_outflow_bdry_extrap_type = input_db->getString("outflow_bdry_extrap_type");
        if (input_db->keyExists("bdry_extrap_type"))
        {
            TBOX_ERROR("AdvDiffCUIConvectiveOperator::AdvDiffCUIConvectiveOperator():\n"
                       << "  input database key ``bdry_extrap_type'' has been changed to "
                          "``outflow_bdry_extrap_type''\n");
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
    d_Q_scratch_idx = var_db->registerVariableAndContext(d_Q_var, context, GADVECTG);
    Pointer<CellDataFactory<NDIM, double> > Q_pdat_fac = d_Q_var->getPatchDataFactory();
    d_Q_data_depth = Q_pdat_fac->getDefaultDepth();
    const std::string q_extrap_var_name = d_object_name + "::q_extrap";
    d_q_extrap_var = var_db->getVariable(q_extrap_var_name);
    if (d_q_extrap_var)
    {
        d_q_extrap_idx = var_db->mapVariableAndContextToIndex(d_q_extrap_var, context);
    }
    else
    {
        d_q_extrap_var = new FaceVariable<NDIM, double>(q_extrap_var_name, d_Q_data_depth);
        d_q_extrap_idx = var_db->registerVariableAndContext(d_q_extrap_var, context, IntVector<NDIM>(0));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_q_extrap_idx >= 0);
#endif
    const std::string q_flux_var_name = d_object_name + "::q_flux";
    d_q_flux_var = var_db->getVariable(q_flux_var_name);
    if (d_q_flux_var)
    {
        d_q_flux_idx = var_db->mapVariableAndContextToIndex(d_q_flux_var, context);
    }
    else
    {
        d_q_flux_var = new FaceVariable<NDIM, double>(q_flux_var_name, d_Q_data_depth);
        d_q_flux_idx = var_db->registerVariableAndContext(d_q_flux_var, context, IntVector<NDIM>(0));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_q_flux_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_apply_convective_operator =
            TimerManager::getManager()->getTimer("IBAMR::AdvDiffCUIConvectiveOperator::applyConvectiveOperator()");
        t_apply = TimerManager::getManager()->getTimer("IBAMR::AdvDiffCUIConvectiveOperator::apply()");
        t_initialize_operator_state =
            TimerManager::getManager()->getTimer("IBAMR::AdvDiffCUIConvectiveOperator::initializeOperatorState()");
        t_deallocate_operator_state =
            TimerManager::getManager()->getTimer("IBAMR::AdvDiffCUIConvectiveOperator::deallocateOperatorState()"););
    return;
} // AdvDiffCUIConvectiveOperator

AdvDiffCUIConvectiveOperator::~AdvDiffCUIConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~AdvDiffCUIConvectiveOperator

void
AdvDiffCUIConvectiveOperator::applyConvectiveOperator(const int Q_idx, const int N_idx)
{
    IBAMR_TIMER_START(t_apply_convective_operator);
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("AdvDiffCUIConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }
#endif

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_Q_scratch_idx);
        level->allocatePatchData(d_q_extrap_idx);
        if (d_difference_form == CONSERVATIVE || d_difference_form == SKEW_SYMMETRIC)
            level->allocatePatchData(d_q_flux_idx);
    }

    // Setup communications algorithm.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(d_Q_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg->registerRefine(d_Q_scratch_idx, Q_idx, d_Q_scratch_idx, refine_op);

    // Extrapolate from cell centers to cell faces.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        refine_alg->resetSchedule(d_ghostfill_scheds[ln]);
        d_ghostfill_scheds[ln]->fillData(d_solution_time);
        d_ghostfill_alg->resetSchedule(d_ghostfill_scheds[ln]);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(d_Q_scratch_idx);
            const IntVector<NDIM>& Q_data_gcw = Q_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(Q_data_gcw.min() == Q_data_gcw.max());
#endif
            Pointer<FaceData<NDIM, double> > u_ADV_data = patch->getPatchData(d_u_idx);
            const IntVector<NDIM>& u_ADV_data_gcw = u_ADV_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(u_ADV_data_gcw.min() == u_ADV_data_gcw.max());
#endif
            Pointer<FaceData<NDIM, double> > q_extrap_data = patch->getPatchData(d_q_extrap_idx);
            const IntVector<NDIM>& q_extrap_data_gcw = q_extrap_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(q_extrap_data_gcw.min() == q_extrap_data_gcw.max());
#endif
            CellData<NDIM, double>& Q0_data = *Q_data;
            CellData<NDIM, double> Q1_data(patch_box, 1, Q_data_gcw);
#if (NDIM == 3)
            CellData<NDIM, double> Q2_data(patch_box, 1, Q_data_gcw);
#endif

            // Enforce physical boundary conditions at inflow boundaries.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q_data,
                u_ADV_data,
                patch,
                d_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_outflow_bdry_extrap_type != "NONE",
                d_homogeneous_bc);

            // Extrapolate from cell centers to cell faces.
            for (unsigned int d = 0; d < d_Q_data_depth; ++d)
            {
                CUI_EXTRAPOLATE_FC(
#if (NDIM == 2)
                    patch_lower(0),
                    patch_upper(0),
                    patch_lower(1),
                    patch_upper(1),
                    Q_data_gcw(0),
                    Q_data_gcw(1),
                    Q0_data.getPointer(d),
                    Q1_data.getPointer(),
                    u_ADV_data_gcw(0),
                    u_ADV_data_gcw(1),
                    q_extrap_data_gcw(0),
                    q_extrap_data_gcw(1),
                    u_ADV_data->getPointer(0),
                    u_ADV_data->getPointer(1),
                    q_extrap_data->getPointer(0, d),
                    q_extrap_data->getPointer(1, d)
#endif
#if (NDIM == 3)
                        patch_lower(0),
                    patch_upper(0),
                    patch_lower(1),
                    patch_upper(1),
                    patch_lower(2),
                    patch_upper(2),
                    Q_data_gcw(0),
                    Q_data_gcw(1),
                    Q_data_gcw(2),
                    Q0_data.getPointer(d),
                    Q1_data.getPointer(),
                    Q2_data.getPointer(),
                    u_ADV_data_gcw(0),
                    u_ADV_data_gcw(1),
                    u_ADV_data_gcw(2),
                    q_extrap_data_gcw(0),
                    q_extrap_data_gcw(1),
                    q_extrap_data_gcw(2),
                    u_ADV_data->getPointer(0),
                    u_ADV_data->getPointer(1),
                    u_ADV_data->getPointer(2),
                    q_extrap_data->getPointer(0, d),
                    q_extrap_data->getPointer(1, d),
                    q_extrap_data->getPointer(2, d)
#endif
                        );
            }

            // If we are using conservative or skew-symmetric differencing,
            // compute the advective fluxes.  These need to be synchronized on
            // the patch hierarchy.
            if (d_difference_form == CONSERVATIVE || d_difference_form == SKEW_SYMMETRIC)
            {
                Pointer<FaceData<NDIM, double> > u_ADV_data = patch->getPatchData(d_u_idx);
                const IntVector<NDIM>& u_ADV_data_gcw = u_ADV_data->getGhostCellWidth();
                Pointer<FaceData<NDIM, double> > q_flux_data = patch->getPatchData(d_q_flux_idx);
                const IntVector<NDIM>& q_flux_data_gcw = q_flux_data->getGhostCellWidth();
                for (unsigned int d = 0; d < d_Q_data_depth; ++d)
                {
                    static const double dt = 1.0;
                    ADVECT_FLUX_FC(dt,
#if (NDIM == 2)
                                   patch_lower(0),
                                   patch_upper(0),
                                   patch_lower(1),
                                   patch_upper(1),
                                   u_ADV_data_gcw(0),
                                   u_ADV_data_gcw(1),
                                   q_extrap_data_gcw(0),
                                   q_extrap_data_gcw(1),
                                   q_flux_data_gcw(0),
                                   q_flux_data_gcw(1),
                                   u_ADV_data->getPointer(0),
                                   u_ADV_data->getPointer(1),
                                   q_extrap_data->getPointer(0, d),
                                   q_extrap_data->getPointer(1, d),
                                   q_flux_data->getPointer(0, d),
                                   q_flux_data->getPointer(1, d)
#endif
#if (NDIM == 3)
                                       patch_lower(0),
                                   patch_upper(0),
                                   patch_lower(1),
                                   patch_upper(1),
                                   patch_lower(2),
                                   patch_upper(2),
                                   u_ADV_data_gcw(0),
                                   u_ADV_data_gcw(1),
                                   u_ADV_data_gcw(2),
                                   q_extrap_data_gcw(0),
                                   q_extrap_data_gcw(1),
                                   q_extrap_data_gcw(2),
                                   q_flux_data_gcw(0),
                                   q_flux_data_gcw(1),
                                   q_flux_data_gcw(2),
                                   u_ADV_data->getPointer(0),
                                   u_ADV_data->getPointer(1),
                                   u_ADV_data->getPointer(2),
                                   q_extrap_data->getPointer(0, d),
                                   q_extrap_data->getPointer(1, d),
                                   q_extrap_data->getPointer(2, d),
                                   q_flux_data->getPointer(0, d),
                                   q_flux_data->getPointer(1, d),
                                   q_flux_data->getPointer(2, d)
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > N_data = patch->getPatchData(N_idx);
            const IntVector<NDIM>& N_data_gcw = N_data->getGhostCellWidth();

            if (d_difference_form == ADVECTIVE || d_difference_form == SKEW_SYMMETRIC)
            {
                Pointer<FaceData<NDIM, double> > u_ADV_data = patch->getPatchData(d_u_idx);
                const IntVector<NDIM>& u_ADV_data_gcw = u_ADV_data->getGhostCellWidth();
                Pointer<FaceData<NDIM, double> > q_extrap_data = patch->getPatchData(d_q_extrap_idx);
                const IntVector<NDIM>& q_extrap_data_gcw = q_extrap_data->getGhostCellWidth();
                for (unsigned int d = 0; d < d_Q_data_depth; ++d)
                {
                    ADVECT_DERIVATIVE_FC(dx,
#if (NDIM == 2)
                                         patch_lower(0),
                                         patch_upper(0),
                                         patch_lower(1),
                                         patch_upper(1),
                                         u_ADV_data_gcw(0),
                                         u_ADV_data_gcw(1),
                                         q_extrap_data_gcw(0),
                                         q_extrap_data_gcw(1),
                                         u_ADV_data->getPointer(0),
                                         u_ADV_data->getPointer(1),
                                         q_extrap_data->getPointer(0, d),
                                         q_extrap_data->getPointer(1, d),
                                         N_data_gcw(0),
                                         N_data_gcw(1),
#endif
#if (NDIM == 3)
                                         patch_lower(0),
                                         patch_upper(0),
                                         patch_lower(1),
                                         patch_upper(1),
                                         patch_lower(2),
                                         patch_upper(2),
                                         u_ADV_data_gcw(0),
                                         u_ADV_data_gcw(1),
                                         u_ADV_data_gcw(2),
                                         q_extrap_data_gcw(0),
                                         q_extrap_data_gcw(1),
                                         q_extrap_data_gcw(2),
                                         u_ADV_data->getPointer(0),
                                         u_ADV_data->getPointer(1),
                                         u_ADV_data->getPointer(2),
                                         q_extrap_data->getPointer(0, d),
                                         q_extrap_data->getPointer(1, d),
                                         q_extrap_data->getPointer(2, d),
                                         N_data_gcw(0),
                                         N_data_gcw(1),
                                         N_data_gcw(2),
#endif
                                         N_data->getPointer(d));
                }
            }

            if (d_difference_form == CONSERVATIVE)
            {
                Pointer<FaceData<NDIM, double> > q_flux_data = patch->getPatchData(d_q_flux_idx);
                const IntVector<NDIM>& q_flux_data_gcw = q_flux_data->getGhostCellWidth();
                for (unsigned int d = 0; d < d_Q_data_depth; ++d)
                {
                    static const double alpha = 1.0;
                    F_TO_C_DIV_FC(N_data->getPointer(d),
                                  N_data_gcw.min(),
                                  alpha,
#if (NDIM == 2)
                                  q_flux_data->getPointer(0, d),
                                  q_flux_data->getPointer(1, d),
                                  q_flux_data_gcw.min(),
                                  patch_lower(0),
                                  patch_upper(0),
                                  patch_lower(1),
                                  patch_upper(1),
#endif
#if (NDIM == 3)
                                  q_flux_data->getPointer(0, d),
                                  q_flux_data->getPointer(1, d),
                                  q_flux_data->getPointer(2, d),
                                  q_flux_data_gcw.min(),
                                  patch_lower(0),
                                  patch_upper(0),
                                  patch_lower(1),
                                  patch_upper(1),
                                  patch_lower(2),
                                  patch_upper(2),
#endif
                                  dx);
                }
            }

            if (d_difference_form == SKEW_SYMMETRIC)
            {
                Pointer<FaceData<NDIM, double> > q_flux_data = patch->getPatchData(d_q_flux_idx);
                const IntVector<NDIM>& q_flux_data_gcw = q_flux_data->getGhostCellWidth();
                for (unsigned int d = 0; d < d_Q_data_depth; ++d)
                {
                    static const double alpha = 0.5;
                    static const double beta = 0.5;
                    F_TO_C_DIV_ADD_FC(N_data->getPointer(d),
                                      N_data_gcw.min(),
                                      alpha,
#if (NDIM == 2)
                                      q_flux_data->getPointer(0, d),
                                      q_flux_data->getPointer(1, d),
                                      q_flux_data_gcw.min(),
                                      beta,
                                      N_data->getPointer(d),
                                      N_data_gcw.min(),
                                      patch_lower(0),
                                      patch_upper(0),
                                      patch_lower(1),
                                      patch_upper(1),
#endif
#if (NDIM == 3)
                                      q_flux_data->getPointer(0, d),
                                      q_flux_data->getPointer(1, d),
                                      q_flux_data->getPointer(2, d),
                                      q_flux_data_gcw.min(),
                                      beta,
                                      N_data->getPointer(d),
                                      N_data_gcw.min(),
                                      patch_lower(0),
                                      patch_upper(0),
                                      patch_lower(1),
                                      patch_upper(1),
                                      patch_lower(2),
                                      patch_upper(2),
#endif
                                      dx);
                }
            }
        }
    }

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_Q_scratch_idx);
        level->deallocatePatchData(d_q_extrap_idx);
        if (d_difference_form == CONSERVATIVE || d_difference_form == SKEW_SYMMETRIC)
            level->deallocatePatchData(d_q_flux_idx);
    }

    IBAMR_TIMER_STOP(t_apply_convective_operator);
    return;
} // applyConvectiveOperator

void
AdvDiffCUIConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
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
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Setup the coarsen algorithm, operator, and schedules.
    Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_q_flux_var, "CONSERVATIVE_COARSEN");
    d_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    if (d_difference_form == ADVECTIVE || d_difference_form == SKEW_SYMMETRIC)
        d_coarsen_alg->registerCoarsen(d_q_extrap_idx, d_q_extrap_idx, coarsen_op);
    if (d_difference_form == CONSERVATIVE || d_difference_form == SKEW_SYMMETRIC)
        d_coarsen_alg->registerCoarsen(d_q_flux_idx, d_q_flux_idx, coarsen_op);
    d_coarsen_scheds.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln - 1);
        d_coarsen_scheds[ln] = d_coarsen_alg->createSchedule(coarser_level, level);
    }

    // Setup the refine algorithm, operator, patch strategy, and schedules.
    Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(d_Q_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ghostfill_alg = new RefineAlgorithm<NDIM>();
    d_ghostfill_alg->registerRefine(d_Q_scratch_idx, in.getComponentDescriptorIndex(0), d_Q_scratch_idx, refine_op);
    if (d_outflow_bdry_extrap_type != "NONE")
        d_ghostfill_strategy = new CartExtrapPhysBdryOp(d_Q_scratch_idx, d_outflow_bdry_extrap_type);
    d_ghostfill_scheds.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_ghostfill_scheds[ln] = d_ghostfill_alg->createSchedule(level, ln - 1, d_hierarchy, d_ghostfill_strategy);
    }

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
AdvDiffCUIConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

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
} // deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
