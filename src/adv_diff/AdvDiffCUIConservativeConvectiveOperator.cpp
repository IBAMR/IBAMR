// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/AdvDiffCUIConservativeConvectiveOperator.h"
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartExtrapPhysBdryOp.h"

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
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux2d, ADVECT_FLUX2D)
#define F_TO_C_DIV_FC IBAMR_FC_FUNC_(ftocdiv2d, FTOCDIV2D)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate2d, CUI_EXTRAPOLATE2D)
#endif

#if (NDIM == 3)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux3d, ADVECT_FLUX3D)
#define F_TO_C_DIV_FC IBAMR_FC_FUNC_(ftocdiv3d, FTOCDIV3D)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate3d, CUI_EXTRAPOLATE3D)
#endif

extern "C"
{
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
} // namespace
/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffCUIConservativeConvectiveOperator::AdvDiffCUIConservativeConvectiveOperator(
    std::string object_name,
    Pointer<CellVariable<NDIM, double> > Q_var,
    Pointer<CellVariable<NDIM, double> > Q1_var,
    Pointer<Database> input_db,
    const ConvectiveDifferencingType difference_form,
    std::vector<RobinBcCoefStrategy<NDIM>*> Q_bc_coefs,
    std::vector<RobinBcCoefStrategy<NDIM>*> Q1_bc_coefs)
    : AdvDiffCUIConvectiveOperator(std::move(object_name), Q_var, input_db, difference_form, Q_bc_coefs),
      d_Q1_var(Q1_var),
      d_Q1_bc_coefs(std::move(Q1_bc_coefs))
{
    if (d_difference_form != CONSERVATIVE)
    {
        TBOX_ERROR("AdvDiffCUIConservativeConvectiveOperator::AdvDiffCUIConservativeConvectiveOperator():\n"
                   << "  unsupported differencing form: "
                   << enum_to_string<ConvectiveDifferencingType>(d_difference_form) << " \n"
                   << "  valid choice is: CONSERVATIVE\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("outflow_bdry_extrap_type"))
            d_outflow_bdry_extrap_type = input_db->getString("outflow_bdry_extrap_type");
        if (input_db->keyExists("bdry_extrap_type"))
        {
            TBOX_ERROR("AdvDiffCUIConservativeConvectiveOperator::AdvDiffCUIConservativeConvectiveOperator():\n"
                       << "  input database key ``bdry_extrap_type'' has been changed to "
                          "``outflow_bdry_extrap_type''\n");
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
    d_Q1_scratch_idx = var_db->registerVariableAndContext(d_Q1_var, context, GADVECTG);

    Pointer<CellDataFactory<NDIM, double> > Q1_pdat_fac = d_Q1_var->getPatchDataFactory();
    d_Q1_data_depth = Q1_pdat_fac->getDefaultDepth();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q1_data_depth == d_Q_data_depth);
#endif

    const std::string q1_extrap_var_name = d_object_name + "::q1_extrap";
    d_q1_extrap_var = var_db->getVariable(q1_extrap_var_name);
    if (d_q1_extrap_var)
    {
        d_q1_extrap_idx = var_db->mapVariableAndContextToIndex(d_q1_extrap_var, context);
    }
    else
    {
        d_q1_extrap_var = new FaceVariable<NDIM, double>(q1_extrap_var_name, d_Q1_data_depth);
        d_q1_extrap_idx = var_db->registerVariableAndContext(d_q1_extrap_var, context, IntVector<NDIM>(0));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_q1_extrap_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator = TimerManager::getManager()->getTimer(
                      "IBAMR::AdvDiffCUIConservativeConvectiveOperator::applyConvectiveOperator()");
                  t_apply =
                      TimerManager::getManager()->getTimer("IBAMR::AdvDiffCUIConservativeConvectiveOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::AdvDiffCUIConservativeConvectiveOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::AdvDiffCUIConservativeConvectiveOperator::deallocateOperatorState()"););
    return;
} // AdvDiffCUIConservativeConvectiveOperator

AdvDiffCUIConservativeConvectiveOperator::~AdvDiffCUIConservativeConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~AdvDiffCUIConservativeConvectiveOperator

void
AdvDiffCUIConservativeConvectiveOperator::applyConvectiveOperator(const int Q_idx, const int N_idx)
{
    TBOX_ERROR(
        "AdvDiffCUIConservativeConvectiveOperator::applyConvectiveOperator can handle only conservative form of "
        "convective operator with three arguments.\n");
    return;
}

void
AdvDiffCUIConservativeConvectiveOperator::applyConvectiveOperator(const int Q_idx, const int Q1_idx, const int N_idx)
{
    IBAMR_TIMER_START(t_apply_convective_operator);
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("AdvDiffCUIConservativeConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to "
                      "applyConvectiveOperator\n");
    }
#endif

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_Q_scratch_idx);
        level->allocatePatchData(d_q_extrap_idx);
        level->allocatePatchData(d_q_flux_idx);
        level->allocatePatchData(d_Q1_scratch_idx);
        level->allocatePatchData(d_q1_extrap_idx);
    }

    // Setup communications algorithm.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(d_Q_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg->registerRefine(d_Q_scratch_idx, Q_idx, d_Q_scratch_idx, refine_op);

    Pointer<RefineAlgorithm<NDIM> > refine_alg1 = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op1 = grid_geom->lookupRefineOperator(d_Q1_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg1->registerRefine(d_Q1_scratch_idx, Q1_idx, d_Q1_scratch_idx, refine_op1);

    // Extrapolate from cell centers to cell faces.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        refine_alg->resetSchedule(d_ghostfill_scheds[ln]);
        d_ghostfill_scheds[ln]->fillData(d_solution_time);
        d_ghostfill_alg->resetSchedule(d_ghostfill_scheds[ln]);

        refine_alg1->resetSchedule(d_ghostfill_scheds1[ln]);
        d_ghostfill_scheds1[ln]->fillData(d_solution_time);
        d_ghostfill_alg1->resetSchedule(d_ghostfill_scheds1[ln]);

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

            Pointer<CellData<NDIM, double> > Q1_data = patch->getPatchData(d_Q1_scratch_idx);
            const IntVector<NDIM>& Q1_data_gcw = Q1_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(Q1_data_gcw.min() == Q1_data_gcw.max());
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

            Pointer<FaceData<NDIM, double> > q1_extrap_data = patch->getPatchData(d_q1_extrap_idx);
            const IntVector<NDIM>& q1_extrap_data_gcw = q1_extrap_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(q1_extrap_data_gcw.min() == q1_extrap_data_gcw.max());
#endif

            CellData<NDIM, double>& Q00_data = *Q_data;
            CellData<NDIM, double> Q01_data(patch_box, 1, Q_data_gcw);
#if (NDIM == 3)
            CellData<NDIM, double> Q02_data(patch_box, 1, Q_data_gcw);
#endif

            CellData<NDIM, double>& Q10_data = *Q1_data;
            CellData<NDIM, double> Q11_data(patch_box, 1, Q1_data_gcw);
#if (NDIM == 3)
            CellData<NDIM, double> Q12_data(patch_box, 1, Q1_data_gcw);
#endif

            // Enforce physical boundary conditions at inflow boundaries for Q1 variable.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q_data,
                u_ADV_data,
                patch,
                d_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_outflow_bdry_extrap_type != "NONE",
                d_homogeneous_bc);

            // Enforce physical boundary conditions at inflow boundaries for Q2 variable.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q1_data,
                u_ADV_data,
                patch,
                d_Q1_bc_coefs,
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
                    Q00_data.getPointer(d),
                    Q01_data.getPointer(),
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
                    Q00_data.getPointer(d),
                    Q01_data.getPointer(),
                    Q02_data.getPointer(),
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

            for (unsigned int d = 0; d < d_Q1_data_depth; ++d)
            {
                CUI_EXTRAPOLATE_FC(
#if (NDIM == 2)
                    patch_lower(0),
                    patch_upper(0),
                    patch_lower(1),
                    patch_upper(1),
                    Q1_data_gcw(0),
                    Q1_data_gcw(1),
                    Q10_data.getPointer(d),
                    Q11_data.getPointer(),
                    u_ADV_data_gcw(0),
                    u_ADV_data_gcw(1),
                    q1_extrap_data_gcw(0),
                    q1_extrap_data_gcw(1),
                    u_ADV_data->getPointer(0),
                    u_ADV_data->getPointer(1),
                    q1_extrap_data->getPointer(0, d),
                    q1_extrap_data->getPointer(1, d)
#endif
#if (NDIM == 3)
                        patch_lower(0),
                    patch_upper(0),
                    patch_lower(1),
                    patch_upper(1),
                    patch_lower(2),
                    patch_upper(2),
                    Q1_data_gcw(0),
                    Q1_data_gcw(1),
                    Q1_data_gcw(2),
                    Q10_data.getPointer(d),
                    Q11_data.getPointer(),
                    Q12_data.getPointer(),
                    u_ADV_data_gcw(0),
                    u_ADV_data_gcw(1),
                    u_ADV_data_gcw(2),
                    q1_extrap_data_gcw(0),
                    q1_extrap_data_gcw(1),
                    q1_extrap_data_gcw(2),
                    u_ADV_data->getPointer(0),
                    u_ADV_data->getPointer(1),
                    u_ADV_data->getPointer(2),
                    q1_extrap_data->getPointer(0, d),
                    q1_extrap_data->getPointer(1, d),
                    q1_extrap_data->getPointer(2, d)
#endif
                );
            }
        }
    }

    HierarchyFaceDataOpsReal<NDIM, double> hier_fc_data_ops(d_hierarchy, d_coarsest_ln, d_finest_ln);
    // q_extrap_idx = q_extrap_idx*q1_extrap_idx
    hier_fc_data_ops.multiply(d_q_extrap_idx, d_q1_extrap_idx, d_q_extrap_idx);

    // Extrapolate from cell centers to cell faces.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            // Since we are using conservative differencing,
            // compute the advective fluxes.  These need to be synchronized on
            // the patch hierarchy.
            Pointer<FaceData<NDIM, double> > q_extrap_data = patch->getPatchData(d_q_extrap_idx);
            const IntVector<NDIM>& q_extrap_data_gcw = q_extrap_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(q_extrap_data_gcw.min() == q_extrap_data_gcw.max());
#endif

            Pointer<FaceData<NDIM, double> > u_ADV_data = patch->getPatchData(d_u_idx);
            const IntVector<NDIM>& u_ADV_data_gcw = u_ADV_data->getGhostCellWidth();
            Pointer<FaceData<NDIM, double> > q_flux_data = patch->getPatchData(d_q_flux_idx);
            const IntVector<NDIM>& q_flux_data_gcw = q_flux_data->getGhostCellWidth();
            for (unsigned int d = 0; d < d_Q1_data_depth; ++d)
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

            Pointer<FaceData<NDIM, double> > q_flux_data = patch->getPatchData(d_q_flux_idx);
            const IntVector<NDIM>& q_flux_data_gcw = q_flux_data->getGhostCellWidth();
            for (unsigned int d = 0; d < d_Q1_data_depth; ++d)
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
    }

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_Q_scratch_idx);
        level->deallocatePatchData(d_q_extrap_idx);
        level->deallocatePatchData(d_q_flux_idx);
        level->deallocatePatchData(d_Q1_scratch_idx);
        level->deallocatePatchData(d_q1_extrap_idx);
    }

    IBAMR_TIMER_STOP(t_apply_convective_operator);
    return;
} // applyConvectiveOperator

void
AdvDiffCUIConservativeConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                                  const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    if (d_is_initialized) deallocateOperatorState();

    AdvDiffCUIConvectiveOperator::initializeOperatorState(in, out);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Setup the refine algorithm, operator, patch strategy, and schedules.
    d_ghostfill_alg1 = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op1 = grid_geom->lookupRefineOperator(d_Q1_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ghostfill_alg1->registerRefine(d_Q1_scratch_idx, in.getComponentDescriptorIndex(0), d_Q1_scratch_idx, refine_op1);
    if (d_outflow_bdry_extrap_type != "NONE")
        d_ghostfill_strategy1 = new CartExtrapPhysBdryOp(d_Q1_scratch_idx, d_outflow_bdry_extrap_type);

    d_ghostfill_scheds1.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_ghostfill_scheds1[ln] = d_ghostfill_alg1->createSchedule(level, ln - 1, d_hierarchy, d_ghostfill_strategy1);
    }

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
AdvDiffCUIConservativeConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    AdvDiffCUIConvectiveOperator::deallocateOperatorState();

    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_ghostfill_alg1.setNull();
    d_ghostfill_strategy1.setNull();
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_ghostfill_scheds1[ln].setNull();
    }
    d_ghostfill_scheds1.clear();

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
