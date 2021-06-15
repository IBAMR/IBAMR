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

#include "ibamr/AdvDiffConservativeCUIConvectiveOperator.h"
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/ConvectiveOperator.h"
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

extern "C"
{
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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffConservativeCUIConvectiveOperator::AdvDiffConservativeCUIConvectiveOperator(
    std::string object_name,
    Pointer<CellVariable<NDIM, double> > Q1_var,
    Pointer<Database> input_db,
    const ConvectiveDifferencingType difference_form,
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs)
    : ConvectiveOperator(std::move(object_name), difference_form), d_bc_coefs(std::move(bc_coefs)), d_Q1_var(Q1_var)
{
    if (d_difference_form != CONSERVATIVE)
    {
        TBOX_ERROR("AdvDiffConservativeCUIConvectiveOperator::AdvDiffConservativeCUIConvectiveOperator():\n"
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
            TBOX_ERROR("AdvDiffConservativeCUIConvectiveOperator::AdvDiffConservativeCUIConvectiveOperator():\n"
                       << "  input database key ``bdry_extrap_type'' has been changed to "
                          "``outflow_bdry_extrap_type''\n");
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
    d_Q1_scratch_idx = var_db->registerVariableAndContext(d_Q1_var, context, GADVECTG);

    Pointer<CellDataFactory<NDIM, double> > Q1_pdat_fac = d_Q1_var->getPatchDataFactory();
    d_Q1_data_depth = Q1_pdat_fac->getDefaultDepth();

    const std::string q_extrap_var_name = d_object_name + "::q_extrap";
    d_q_extrap_var = var_db->getVariable(q_extrap_var_name);
    if (d_q_extrap_var)
    {
        d_q_extrap_idx = var_db->mapVariableAndContextToIndex(d_q_extrap_var, context);
    }
    else
    {
        d_q_extrap_var = new FaceVariable<NDIM, double>(q_extrap_var_name, d_Q1_data_depth);
        d_q_extrap_idx = var_db->registerVariableAndContext(d_q_extrap_var, context, IntVector<NDIM>(0));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_q_extrap_idx >= 0);
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

    const std::string q2_extrap_var_name = d_object_name + "::q2_extrap";
    d_q2_extrap_var = var_db->getVariable(q2_extrap_var_name);
    if (d_q2_extrap_var)
    {
        d_q2_extrap_idx = var_db->mapVariableAndContextToIndex(d_q2_extrap_var, context);
    }
    else
    {
        d_q2_extrap_var = new FaceVariable<NDIM, double>(q2_extrap_var_name, d_Q1_data_depth);
        d_q2_extrap_idx = var_db->registerVariableAndContext(d_q2_extrap_var, context, IntVector<NDIM>(0));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_q2_extrap_idx >= 0);
#endif

    const std::string q_flux_var_name = d_object_name + "::q_flux";
    d_q_flux_var = var_db->getVariable(q_flux_var_name);
    if (d_q_flux_var)
    {
        d_q_flux_idx = var_db->mapVariableAndContextToIndex(d_q_flux_var, context);
    }
    else
    {
        d_q_flux_var = new FaceVariable<NDIM, double>(q_flux_var_name, d_Q1_data_depth);
        d_q_flux_idx = var_db->registerVariableAndContext(d_q_flux_var, context, IntVector<NDIM>(0));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_q_flux_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator = TimerManager::getManager()->getTimer(
                      "IBAMR::AdvDiffConservativeCUIConvectiveOperator::applyConvectiveOperator()");
                  t_apply =
                      TimerManager::getManager()->getTimer("IBAMR::AdvDiffConservativeCUIConvectiveOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::AdvDiffConservativeCUIConvectiveOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::AdvDiffConservativeCUIConvectiveOperator::deallocateOperatorState()"););
    return;
} // AdvDiffConservativeCUIConvectiveOperator

AdvDiffConservativeCUIConvectiveOperator::~AdvDiffConservativeCUIConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~AdvDiffConservativeCUIConvectiveOperator

void
AdvDiffConservativeCUIConvectiveOperator::applyConvectiveOperator(const int Q_idx, const int N_idx)
{
    return;
}

void
AdvDiffConservativeCUIConvectiveOperator::applyConvectiveOperator(const int Q1_idx, const int Q2_idx, const int N_idx)
{
    IBAMR_TIMER_START(t_apply_convective_operator);
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("AdvDiffConservativeCUIConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to "
                      "applyConvectiveOperator\n");
    }
#endif
    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_Q1_scratch_idx);
        level->allocatePatchData(d_Q2_scratch_idx);
        level->allocatePatchData(d_q_extrap_idx);
        level->allocatePatchData(d_q1_extrap_idx);
        level->allocatePatchData(d_q2_extrap_idx);
        if (d_difference_form == CONSERVATIVE) level->allocatePatchData(d_q_flux_idx);
    }

    // Setup communications algorithm.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op1 = grid_geom->lookupRefineOperator(d_Q1_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg->registerRefine(d_Q1_scratch_idx, Q1_idx, d_Q1_scratch_idx, refine_op1);

    Pointer<RefineAlgorithm<NDIM> > refine_alg1 = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op2 = grid_geom->lookupRefineOperator(d_Q2_var, "CONSERVATIVE_LINEAR_REFINE");
    refine_alg1->registerRefine(d_Q2_scratch_idx, Q2_idx, d_Q2_scratch_idx, refine_op2);

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

            Pointer<CellData<NDIM, double> > Q1_data = patch->getPatchData(d_Q1_scratch_idx);
            const IntVector<NDIM>& Q1_data_gcw = Q1_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(Q1_data_gcw.min() == Q1_data_gcw.max());
#endif

            Pointer<CellData<NDIM, double> > Q2_data = patch->getPatchData(d_Q2_scratch_idx);
            const IntVector<NDIM>& Q2_data_gcw = Q2_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(Q2_data_gcw.min() == Q2_data_gcw.max());
#endif

            Pointer<FaceData<NDIM, double> > u_ADV_data = patch->getPatchData(d_u_idx);
            const IntVector<NDIM>& u_ADV_data_gcw = u_ADV_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(u_ADV_data_gcw.min() == u_ADV_data_gcw.max());
#endif

            Pointer<FaceData<NDIM, double> > q1_extrap_data = patch->getPatchData(d_q1_extrap_idx);
            const IntVector<NDIM>& q1_extrap_data_gcw = q1_extrap_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(q1_extrap_data_gcw.min() == q1_extrap_data_gcw.max());
#endif

            Pointer<FaceData<NDIM, double> > q2_extrap_data = patch->getPatchData(d_q2_extrap_idx);
            const IntVector<NDIM>& q2_extrap_data_gcw = q2_extrap_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(q2_extrap_data_gcw.min() == q2_extrap_data_gcw.max());
#endif

            CellData<NDIM, double>& Q10_data = *Q1_data;
            CellData<NDIM, double> Q11_data(patch_box, 1, Q1_data_gcw);
#if (NDIM == 3)
            CellData<NDIM, double> Q12_data(patch_box, 1, Q1_data_gcw);
#endif

            CellData<NDIM, double>& Q20_data = *Q2_data;
            CellData<NDIM, double> Q21_data(patch_box, 1, Q2_data_gcw);
#if (NDIM == 3)
            CellData<NDIM, double> Q22_data(patch_box, 1, Q2_data_gcw);
#endif

            // Enforce physical boundary conditions at inflow boundaries.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q1_data,
                u_ADV_data,
                patch,
                d_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_outflow_bdry_extrap_type != "NONE",
                d_homogeneous_bc);

            // Enforce physical boundary conditions at inflow boundaries.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q2_data,
                u_ADV_data,
                patch,
                d_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_outflow_bdry_extrap_type != "NONE",
                d_homogeneous_bc);

            // Extrapolate from cell centers to cell faces.
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

            for (unsigned int d = 0; d < d_Q2_data_depth; ++d)
            {
                CUI_EXTRAPOLATE_FC(
#if (NDIM == 2)
                    patch_lower(0),
                    patch_upper(0),
                    patch_lower(1),
                    patch_upper(1),
                    Q2_data_gcw(0),
                    Q2_data_gcw(1),
                    Q20_data.getPointer(d),
                    Q21_data.getPointer(),
                    u_ADV_data_gcw(0),
                    u_ADV_data_gcw(1),
                    q2_extrap_data_gcw(0),
                    q2_extrap_data_gcw(1),
                    u_ADV_data->getPointer(0),
                    u_ADV_data->getPointer(1),
                    q2_extrap_data->getPointer(0, d),
                    q2_extrap_data->getPointer(1, d)
#endif
#if (NDIM == 3)
                        patch_lower(0),
                    patch_upper(0),
                    patch_lower(1),
                    patch_upper(1),
                    patch_lower(2),
                    patch_upper(2),
                    Q2_data_gcw(0),
                    Q2_data_gcw(1),
                    Q2_data_gcw(2),
                    Q20_data.getPointer(d),
                    Q21_data.getPointer(),
                    Q22_data.getPointer(),
                    u_ADV_data_gcw(0),
                    u_ADV_data_gcw(1),
                    u_ADV_data_gcw(2),
                    q2_extrap_data_gcw(0),
                    q2_extrap_data_gcw(1),
                    q2_extrap_data_gcw(2),
                    u_ADV_data->getPointer(0),
                    u_ADV_data->getPointer(1),
                    u_ADV_data->getPointer(2),
                    q2_extrap_data->getPointer(0, d),
                    q2_extrap_data->getPointer(1, d),
                    q2_extrap_data->getPointer(2, d)
#endif
                );
            }
        }
    }

    HierarchyFaceDataOpsReal<NDIM, double> hier_fc_data_ops(d_hierarchy, d_coarsest_ln, d_finest_ln);
    hier_fc_data_ops.multiply(d_q_extrap_idx, d_q1_extrap_idx, d_q2_extrap_idx);

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

            // If we are using conservative or skew-symmetric differencing,
            // compute the advective fluxes.  These need to be synchronized on
            // the patch hierarchy.
            if (d_difference_form == CONSERVATIVE)
            {
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

            if (d_difference_form == CONSERVATIVE)
            {
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
    }

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_Q1_scratch_idx);
        level->deallocatePatchData(d_Q2_scratch_idx);
        level->deallocatePatchData(d_q_extrap_idx);
        level->deallocatePatchData(d_q1_extrap_idx);
        level->deallocatePatchData(d_q2_extrap_idx);
        if (d_difference_form == CONSERVATIVE) level->deallocatePatchData(d_q_flux_idx);
    }

    IBAMR_TIMER_STOP(t_apply_convective_operator);
    return;
} // applyConvectiveOperator

void
AdvDiffConservativeCUIConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
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
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
    d_Q2_scratch_idx = var_db->registerVariableAndContext(d_Q2_var, context, GADVECTG);
    Pointer<CellDataFactory<NDIM, double> > Q2_pdat_fac = d_Q2_var->getPatchDataFactory();
    d_Q2_data_depth = Q2_pdat_fac->getDefaultDepth();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q1_data_depth == d_Q2_data_depth);
#endif

    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Setup the coarsen algorithm, operator, and schedules.
    Pointer<CoarsenOperator<NDIM> > coarsen_op = grid_geom->lookupCoarsenOperator(d_q_flux_var, "CONSERVATIVE_COARSEN");
    d_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    if (d_difference_form == CONSERVATIVE) d_coarsen_alg->registerCoarsen(d_q_flux_idx, d_q_flux_idx, coarsen_op);
    d_coarsen_scheds.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln + 1; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln - 1);
        d_coarsen_scheds[ln] = d_coarsen_alg->createSchedule(coarser_level, level);
    }

    // Setup the refine algorithm, operator, patch strategy, and schedules.
    d_ghostfill_alg = new RefineAlgorithm<NDIM>();

    Pointer<RefineOperator<NDIM> > refine_op1 = grid_geom->lookupRefineOperator(d_Q1_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ghostfill_alg->registerRefine(d_Q1_scratch_idx, in.getComponentDescriptorIndex(0), d_Q1_scratch_idx, refine_op1);
    if (d_outflow_bdry_extrap_type != "NONE")
        d_ghostfill_strategy = new CartExtrapPhysBdryOp(d_Q1_scratch_idx, d_outflow_bdry_extrap_type);

    d_ghostfill_scheds.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_ghostfill_scheds[ln] = d_ghostfill_alg->createSchedule(level, ln - 1, d_hierarchy, d_ghostfill_strategy);
    }

    d_ghostfill_alg1 = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op2 = grid_geom->lookupRefineOperator(d_Q2_var, "CONSERVATIVE_LINEAR_REFINE");
    d_ghostfill_alg1->registerRefine(d_Q2_scratch_idx, in.getComponentDescriptorIndex(0), d_Q2_scratch_idx, refine_op2);
    if (d_outflow_bdry_extrap_type != "NONE")
        d_ghostfill_strategy1 = new CartExtrapPhysBdryOp(d_Q2_scratch_idx, d_outflow_bdry_extrap_type);

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
AdvDiffConservativeCUIConvectiveOperator::deallocateOperatorState()
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

void
AdvDiffConservativeCUIConvectiveOperator::setHeavysideIndex(int H_idx)
{
    d_H_idx = H_idx;
    return;
} // setMassDensity

void
AdvDiffConservativeCUIConvectiveOperator::setHeavisideBoundaryConditions(RobinBcCoefStrategy<NDIM>* H_bc_coef)
{
    //    for (unsigned int d = 0; d < NDIM; ++d) d_H_bc_coefs[d] = H_bc_coef;
    return;
} // setMassDensityBoundaryConditions
void
AdvDiffConservativeCUIConvectiveOperator::setHeavisideVariable(Pointer<Variable<NDIM> > H_var)
{
    d_Q2_var = H_var;
    return;
} // setMassDensityVariable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
