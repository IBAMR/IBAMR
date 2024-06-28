// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/CellConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartExtrapPhysBdryOp.h"

#include "SAMRAIVectorReal.h"
#include "tbox/TimerManager.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux2d, ADVECT_FLUX2D)
#define F_TO_C_DIV_FC IBAMR_FC_FUNC(ftocdiv2d, FTOCDIV2D)
#define F_TO_C_DIV_ADD_FC IBAMR_FC_FUNC(ftocdivadd2d, FTOCDIVADD2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux3d, ADVECT_FLUX3D)
#define F_TO_C_DIV_FC IBAMR_FC_FUNC(ftocdiv3d, FTOCDIV3D)
#define F_TO_C_DIV_ADD_FC IBAMR_FC_FUNC(ftocdivadd3d, FTOCDIVADD3D)
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
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CellConvectiveOperator::CellConvectiveOperator(std::string object_name,
                                               Pointer<CellVariable<NDIM, double> > Q_cell_var,
                                               const int Q_min_ghost_cell_width,
                                               Pointer<Database> input_db,
                                               const ConvectiveDifferencingType difference_form,
                                               std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs)

    : ConvectiveOperator(std::move(object_name), difference_form),
      d_bc_coefs(std::move(bc_coefs)),
      d_Q_cell_var(Q_cell_var)
{
    if (input_db)
    {
        if (input_db->keyExists("outflow_bdry_extrap_type"))
            d_outflow_bdry_extrap_type = input_db->getString("outflow_bdry_extrap_type");
        if (input_db->keyExists("bdry_extrap_type"))
        {
            TBOX_ERROR(
                "AdvDiffCenteredCellConvectiveOperator::"
                "AdvDiffCenteredCellConvectiveOperator():\n"
                << "  input database key ``bdry_extrap_type'' has been changed to "
                   "``outflow_bdry_extrap_type''\n");
        }
    }

    if (d_Q_cell_var)
    {
        Pointer<CellDataFactory<NDIM, double> > Q_pdat_fac = d_Q_cell_var->getPatchDataFactory();
        const int Q_depth = Q_pdat_fac->getDefaultDepth();

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
        d_Q_scratch_idx = var_db->registerVariableAndContext(d_Q_cell_var, context, Q_min_ghost_cell_width);
        d_Q_ghost_idx = var_db->registerClonedPatchDataIndex(d_Q_cell_var, d_Q_scratch_idx);

        const std::string q_flux_var_name = d_object_name + "::q_flux";
        d_q_flux_var = var_db->getVariable(q_flux_var_name);
        if (d_q_flux_var)
        {
            d_q_flux_idx = var_db->mapVariableAndContextToIndex(d_q_flux_var, context);
        }
        else
        {
            d_q_flux_var = new FaceVariable<NDIM, double>(q_flux_var_name, Q_depth);
            d_q_flux_idx = var_db->registerVariableAndContext(d_q_flux_var, context, IntVector<NDIM>(0));
        }
        TBOX_ASSERT(d_q_flux_idx >= 0);

        const std::string q_interp_var_name = d_object_name + "::q_interp";
        d_q_interp_var = var_db->getVariable(q_interp_var_name);
        if (d_q_interp_var)
        {
            d_q_interp_idx = var_db->mapVariableAndContextToIndex(d_q_interp_var, context);
        }
        else
        {
            d_q_interp_var = new FaceVariable<NDIM, double>(q_interp_var_name, Q_depth);
            d_q_interp_idx = var_db->registerVariableAndContext(d_q_interp_var, context, IntVector<NDIM>(0));
        }
        TBOX_ASSERT(d_q_interp_idx >= 0);
    }

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator =
                      TimerManager::getManager()->getTimer("IBAMR::CellConvectiveOperator::applyConvectiveOperator()");
                  t_initialize_operator_state =
                      TimerManager::getManager()->getTimer("IBAMR::CellConvectiveOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::CellConvectiveOperator::deallocateOperatorState()"););
    return;
} // CellConvectiveOperator

void
CellConvectiveOperator::interpolateToFaceOnHierarchy(int q_interp_idx, int Q_cell_idx, int u_idx, bool synch_cf_bdry)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("CellConvectiveOperator::interpolateToFaceOnHierarchy():\n"
                   << "  operator must be initialized\n");
    }

    // Setup communications algorithms.
    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(d_Q_ghost_idx, Q_cell_idx, d_Q_scratch_idx, d_Q_cell_refine_op);
    CoarsenAlgorithm<NDIM> coarsen_alg;
    coarsen_alg.registerCoarsen(q_interp_idx, q_interp_idx, d_q_interp_coarsen_op);

    // Interpolate data from cell centers to cell faces.
    auto coarsest_ln = 0, finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        refine_alg.resetSchedule(d_Q_cell_refine_scheds[ln]);
        d_Q_cell_refine_scheds[ln]->fillData(d_solution_time);
        d_Q_cell_refine_alg->resetSchedule(d_Q_cell_refine_scheds[ln]);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            auto patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > Q_ghost_data = patch->getPatchData(d_Q_ghost_idx);
            Pointer<FaceData<NDIM, double> > q_interp_data = patch->getPatchData(q_interp_idx);
            Pointer<FaceData<NDIM, double> > u_data = patch->getPatchData(u_idx);

            // Enforce physical boundary conditions at inflow boundaries.
            //
            // TODO: Consider making this more modular.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q_ghost_data,
                u_data,
                patch,
                d_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_outflow_bdry_extrap_type != "NONE",
                d_homogeneous_bc);

            // Evaluate the implementation-specific interpolation routine.
            interpolateToFaceOnPatch(*q_interp_data, *Q_ghost_data, *u_data, *patch);
        }
    }

    // Synchronize data as requested.
    if (synch_cf_bdry)
    {
        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            coarsen_alg.resetSchedule(d_q_interp_coarsen_scheds[ln]);
            d_q_interp_coarsen_scheds[ln]->coarsenData();
            d_q_interp_coarsen_alg->resetSchedule(d_q_interp_coarsen_scheds[ln]);
        }
    }
    return;
} // interpolateToFaceOnHierarchy

void
CellConvectiveOperator::evaluateAdvectiveFluxOnHierarchy(int q_flux_idx, int Q_cell_idx, int u_idx, bool synch_cf_bdry)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("CellConvectiveOperator::evaluateAdvectiveFluxOnHierarchy():\n"
                   << "  operator must be initialized\n");
    }

    // Setup communications algorithms.
    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(d_Q_ghost_idx, Q_cell_idx, d_Q_scratch_idx, d_Q_cell_refine_op);
    CoarsenAlgorithm<NDIM> coarsen_alg;
    coarsen_alg.registerCoarsen(q_flux_idx, q_flux_idx, d_q_flux_coarsen_op);

    // Interpolate data from cell centers to cell faces.
    auto coarsest_ln = 0, finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        refine_alg.resetSchedule(d_Q_cell_refine_scheds[ln]);
        d_Q_cell_refine_scheds[ln]->fillData(d_solution_time);
        d_Q_cell_refine_alg->resetSchedule(d_Q_cell_refine_scheds[ln]);
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            auto patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > Q_ghost_data = patch->getPatchData(d_Q_ghost_idx);
            Pointer<FaceData<NDIM, double> > q_flux_data = patch->getPatchData(q_flux_idx);
            Pointer<FaceData<NDIM, double> > u_data = patch->getPatchData(u_idx);

            // Enforce physical boundary conditions at inflow boundaries.
            //
            // TODO: Consider making this more modular.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q_ghost_data,
                u_data,
                patch,
                d_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_outflow_bdry_extrap_type != "NONE",
                d_homogeneous_bc);

            // Evaluate the implementation-specific interpolation routine.
            evaluateAdvectiveFluxOnPatch(*q_flux_data, *Q_ghost_data, *u_data, *patch);
        }
    }

    // Synchronize data as requested.
    if (synch_cf_bdry)
    {
        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            coarsen_alg.resetSchedule(d_q_flux_coarsen_scheds[ln]);
            d_q_flux_coarsen_scheds[ln]->coarsenData();
            d_q_flux_coarsen_alg->resetSchedule(d_q_flux_coarsen_scheds[ln]);
        }
    }
    return;
} // evaluateAdvectiveFluxOnHierarchy

void
CellConvectiveOperator::computeAdvectiveDerivativeOnHierarchy(int N_cell_idx,
                                                              int q_interp_idx,
                                                              int u_idx,
                                                              bool synch_cf_bdry)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("CellConvectiveOperator::computeAdvectiveDerivativeOnHierarchy():\n"
                   << "  operator must be initialized\n");
    }

    // Setup communications algorithms.
    CoarsenAlgorithm<NDIM> coarsen_alg;
    coarsen_alg.registerCoarsen(q_interp_idx, q_interp_idx, d_q_interp_coarsen_op);

    // Difference values on the patches.
    auto coarsest_ln = 0, finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        // Synchronize data as requested.
        if (synch_cf_bdry && ln > coarsest_ln)
        {
            coarsen_alg.resetSchedule(d_q_interp_coarsen_scheds[ln]);
            d_q_interp_coarsen_scheds[ln]->coarsenData();
            d_q_interp_coarsen_alg->resetSchedule(d_q_interp_coarsen_scheds[ln]);
        }

        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            auto patch = level->getPatch(p());

            const auto& patch_box = patch->getBox();
            const auto& patch_lower = patch_box.lower();
            const auto& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > N_cell_data = patch->getPatchData(N_cell_idx);
            const auto& N_cell_data_gcw = N_cell_data->getGhostCellWidth();
            Pointer<FaceData<NDIM, double> > u_data = patch->getPatchData(u_idx);
            const auto& u_data_gcw = u_data->getGhostCellWidth();
            Pointer<FaceData<NDIM, double> > q_interp_data = patch->getPatchData(q_interp_idx);
            const auto& q_interp_data_gcw = q_interp_data->getGhostCellWidth();

            for (int d = 0; d < N_cell_data->getDepth(); ++d)
            {
                ADVECT_DERIVATIVE_FC(dx,
#if (NDIM == 2)
                                     patch_lower(0),
                                     patch_upper(0),
                                     patch_lower(1),
                                     patch_upper(1),
                                     u_data_gcw(0),
                                     u_data_gcw(1),
                                     q_interp_data_gcw(0),
                                     q_interp_data_gcw(1),
                                     u_data->getPointer(0),
                                     u_data->getPointer(1),
                                     q_interp_data->getPointer(0, d),
                                     q_interp_data->getPointer(1, d),
                                     N_cell_data_gcw(0),
                                     N_cell_data_gcw(1),
#endif
#if (NDIM == 3)
                                     patch_lower(0),
                                     patch_upper(0),
                                     patch_lower(1),
                                     patch_upper(1),
                                     patch_lower(2),
                                     patch_upper(2),
                                     u_data_gcw(0),
                                     u_data_gcw(1),
                                     u_data_gcw(2),
                                     q_interp_data_gcw(0),
                                     q_interp_data_gcw(1),
                                     q_interp_data_gcw(2),
                                     u_data->getPointer(0),
                                     u_data->getPointer(1),
                                     u_data->getPointer(2),
                                     q_interp_data->getPointer(0, d),
                                     q_interp_data->getPointer(1, d),
                                     q_interp_data->getPointer(2, d),
                                     N_cell_data_gcw(0),
                                     N_cell_data_gcw(1),
                                     N_cell_data_gcw(2),
#endif
                                     N_cell_data->getPointer(d));
            }
        }
    }
    return;
} // computeAdvectiveDerivativeOnHierarchy

void
CellConvectiveOperator::computeConservativeDerivativeOnHierarchy(int N_cell_idx, int q_flux_idx, bool synch_cf_bdry)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("CellConvectiveOperator::computeConservativeDerivativeOnHierarchy():\n"
                   << "  operator must be initialized\n");
    }

    // Setup communications algorithms.
    CoarsenAlgorithm<NDIM> coarsen_alg;
    coarsen_alg.registerCoarsen(q_flux_idx, q_flux_idx, d_q_flux_coarsen_op);

    // Difference values on the patches.
    auto coarsest_ln = 0, finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        // Synchronize data as requested.
        if (synch_cf_bdry && ln > coarsest_ln)
        {
            coarsen_alg.resetSchedule(d_q_flux_coarsen_scheds[ln]);
            d_q_flux_coarsen_scheds[ln]->coarsenData();
            d_q_flux_coarsen_alg->resetSchedule(d_q_flux_coarsen_scheds[ln]);
        }

        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            auto patch = level->getPatch(p());

            const auto& patch_box = patch->getBox();
            const auto& patch_lower = patch_box.lower();
            const auto& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > N_cell_data = patch->getPatchData(N_cell_idx);
            const auto& N_cell_data_gcw = N_cell_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(N_cell_data_gcw.min() == N_cell_data_gcw.max());
#endif
            Pointer<FaceData<NDIM, double> > q_flux_data = patch->getPatchData(q_flux_idx);
            const auto& q_flux_data_gcw = q_flux_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(q_flux_data_gcw.min() == q_flux_data_gcw.max());
#endif
            for (int d = 0; d < N_cell_data->getDepth(); ++d)
            {
                static const double alpha = 1.0;
                F_TO_C_DIV_FC(N_cell_data->getPointer(d),
                              N_cell_data_gcw.min(),
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
    return;
} // computeConservativeDerivativeOnHierarchy

void
CellConvectiveOperator::computeSkewSymmetricDerivativeOnHierarchy(int N_cell_idx,
                                                                  int q_interp_idx,
                                                                  int q_flux_idx,
                                                                  int u_idx,
                                                                  bool synch_cf_bdry)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("CellConvectiveOperator::computeConservativeDerivativeOnHierarchy():\n"
                   << "  operator must be initialized\n");
    }

    // Setup communications algorithms.
    CoarsenAlgorithm<NDIM> flux_coarsen_alg, interp_coarsen_alg;
    flux_coarsen_alg.registerCoarsen(q_flux_idx, q_flux_idx, d_q_flux_coarsen_op);
    interp_coarsen_alg.registerCoarsen(q_interp_idx, q_interp_idx, d_q_interp_coarsen_op);

    // Difference values on the patches.
    auto coarsest_ln = 0, finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        // Synchronize data as requested.
        if (synch_cf_bdry && ln > coarsest_ln)
        {
            flux_coarsen_alg.resetSchedule(d_q_flux_coarsen_scheds[ln]);
            interp_coarsen_alg.resetSchedule(d_q_interp_coarsen_scheds[ln]);
            d_q_flux_coarsen_scheds[ln]->coarsenData();
            d_q_interp_coarsen_scheds[ln]->coarsenData();
            d_q_flux_coarsen_alg->resetSchedule(d_q_flux_coarsen_scheds[ln]);
            d_q_interp_coarsen_alg->resetSchedule(d_q_interp_coarsen_scheds[ln]);
        }

        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            auto patch = level->getPatch(p());

            const auto& patch_box = patch->getBox();
            const auto& patch_lower = patch_box.lower();
            const auto& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > N_cell_data = patch->getPatchData(N_cell_idx);
            const auto& N_cell_data_gcw = N_cell_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(N_cell_data_gcw.min() == N_cell_data_gcw.max());
#endif
            Pointer<FaceData<NDIM, double> > q_flux_data = patch->getPatchData(q_flux_idx);
            const auto& q_flux_data_gcw = q_flux_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(q_flux_data_gcw.min() == q_flux_data_gcw.max());
#endif
            Pointer<FaceData<NDIM, double> > q_interp_data = patch->getPatchData(q_interp_idx);
            const auto& q_interp_data_gcw = q_interp_data->getGhostCellWidth();
#if !defined(NDEBUG)
            TBOX_ASSERT(q_interp_data_gcw.min() == q_interp_data_gcw.max());
#endif
            Pointer<FaceData<NDIM, double> > u_data = patch->getPatchData(u_idx);
            const auto& u_data_gcw = u_data->getGhostCellWidth();
            for (int d = 0; d < N_cell_data->getDepth(); ++d)
            {
                ADVECT_DERIVATIVE_FC(dx,
#if (NDIM == 2)
                                     patch_lower(0),
                                     patch_upper(0),
                                     patch_lower(1),
                                     patch_upper(1),
                                     u_data_gcw(0),
                                     u_data_gcw(1),
                                     q_interp_data_gcw(0),
                                     q_interp_data_gcw(1),
                                     u_data->getPointer(0),
                                     u_data->getPointer(1),
                                     q_interp_data->getPointer(0, d),
                                     q_interp_data->getPointer(1, d),
                                     N_cell_data_gcw(0),
                                     N_cell_data_gcw(1),
#endif
#if (NDIM == 3)
                                     patch_lower(0),
                                     patch_upper(0),
                                     patch_lower(1),
                                     patch_upper(1),
                                     patch_lower(2),
                                     patch_upper(2),
                                     u_data_gcw(0),
                                     u_data_gcw(1),
                                     u_data_gcw(2),
                                     q_interp_data_gcw(0),
                                     q_interp_data_gcw(1),
                                     q_interp_data_gcw(2),
                                     u_data->getPointer(0),
                                     u_data->getPointer(1),
                                     u_data->getPointer(2),
                                     q_interp_data->getPointer(0, d),
                                     q_interp_data->getPointer(1, d),
                                     q_interp_data->getPointer(2, d),
                                     N_cell_data_gcw(0),
                                     N_cell_data_gcw(1),
                                     N_cell_data_gcw(2),
#endif
                                     N_cell_data->getPointer(d));

                static const double alpha = 0.5;
                F_TO_C_DIV_ADD_FC(N_cell_data->getPointer(d),
                                  N_cell_data_gcw.min(),
                                  alpha,
#if (NDIM == 2)
                                  q_flux_data->getPointer(0, d),
                                  q_flux_data->getPointer(1, d),
                                  q_flux_data_gcw.min(),
                                  alpha,
                                  N_cell_data->getPointer(d),
                                  N_cell_data_gcw.min(),
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
                                  alpha,
                                  N_cell_data->getPointer(d),
                                  N_cell_data_gcw.min(),
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
    return;
} // computeSkewSymmetricDerivativeOnHierarchy

void
CellConvectiveOperator::evaluateAdvectiveFluxOnPatch(FaceData<NDIM, double>& q_flux_data,
                                                     const CellData<NDIM, double>& Q_cell_data,
                                                     const FaceData<NDIM, double>& u_data,
                                                     const Patch<NDIM>& patch)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("CellConvectiveOperator::evaluateAdvectiveFluxOnPatch():\n"
                   << "  operator must be initialized\n");
    }

    // Setup temporary data.
    const auto& patch_box = patch.getBox();
    const auto& patch_lower = patch_box.lower();
    const auto& patch_upper = patch_box.upper();

    const auto& u_data_gcw = u_data.getGhostCellWidth();
    const IntVector<NDIM>& q_flux_data_gcw = q_flux_data.getGhostCellWidth();
    FaceData<NDIM, double> q_interp_data(patch_box, q_flux_data.getDepth(), q_flux_data_gcw);
    const IntVector<NDIM>& q_interp_data_gcw = q_interp_data.getGhostCellWidth();

    // Interpolate from cell centers to cell faces.
    interpolateToFaceOnPatch(q_interp_data, Q_cell_data, u_data, patch);

    // Evaluate the advective flux.
    for (int d = 0; d < Q_cell_data.getDepth(); ++d)
    {
        static const double dt = 1.0;
        ADVECT_FLUX_FC(dt,
#if (NDIM == 2)
                       patch_lower(0),
                       patch_upper(0),
                       patch_lower(1),
                       patch_upper(1),
                       u_data_gcw(0),
                       u_data_gcw(1),
                       q_interp_data_gcw(0),
                       q_interp_data_gcw(1),
                       q_flux_data_gcw(0),
                       q_flux_data_gcw(1),
                       u_data.getPointer(0),
                       u_data.getPointer(1),
                       q_interp_data.getPointer(0, d),
                       q_interp_data.getPointer(1, d),
                       q_flux_data.getPointer(0, d),
                       q_flux_data.getPointer(1, d)
#endif
#if (NDIM == 3)
                           patch_lower(0),
                       patch_upper(0),
                       patch_lower(1),
                       patch_upper(1),
                       patch_lower(2),
                       patch_upper(2),
                       u_data_gcw(0),
                       u_data_gcw(1),
                       u_data_gcw(2),
                       q_interp_data_gcw(0),
                       q_interp_data_gcw(1),
                       q_interp_data_gcw(2),
                       q_flux_data_gcw(0),
                       q_flux_data_gcw(1),
                       q_flux_data_gcw(2),
                       u_data.getPointer(0),
                       u_data.getPointer(1),
                       u_data.getPointer(2),
                       q_interp_data.getPointer(0, d),
                       q_interp_data.getPointer(1, d),
                       q_interp_data.getPointer(2, d),
                       q_flux_data.getPointer(0, d),
                       q_flux_data.getPointer(1, d),
                       q_flux_data.getPointer(2, d)
#endif
        );
    }
    return;
} // evaluateAdvectiveFluxOnPatch

void
CellConvectiveOperator::applyConvectiveOperator(int Q_idx, int N_idx)
{
    if (!d_is_initialized)
    {
        TBOX_ERROR("CellConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized\n");
    }

    switch (d_difference_form)
    {
    case ADVECTIVE:
        interpolateToFaceOnHierarchy(d_q_interp_idx, Q_idx, d_u_idx, /*synch_cf_bdry*/ false);
        computeAdvectiveDerivativeOnHierarchy(N_idx, d_q_interp_idx, d_u_idx, /*synch_cf_bdry*/ true);
        break;
    case CONSERVATIVE:
        evaluateAdvectiveFluxOnHierarchy(d_q_flux_idx, Q_idx, d_u_idx, /*synch_cf_bdry*/ false);
        computeConservativeDerivativeOnHierarchy(N_idx, d_q_flux_idx, /*synch_cf_bdry*/ true);
        break;
    case SKEW_SYMMETRIC:
        interpolateToFaceOnHierarchy(d_q_interp_idx, Q_idx, d_u_idx, /*synch_cf_bdry*/ false);
        evaluateAdvectiveFluxOnHierarchy(d_q_flux_idx, Q_idx, d_u_idx, /*synch_cf_bdry*/ false);
        computeSkewSymmetricDerivativeOnHierarchy(N_idx, d_q_interp_idx, d_q_flux_idx, d_u_idx, /*synch_cf_bdry*/ true);
        break;
    default:
        TBOX_ERROR("CellConvectiveOperator::applyConvectiveOperator(): unsupported ConvectiveDifferencingType "
                   << enum_to_string(d_difference_form) << "\n");
    }
} // applyConvectiveOperator

void
CellConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    if (d_is_initialized) deallocateOperatorState();

    // Get the hierarchy configuration.
    d_hierarchy = in.getPatchHierarchy();
    const int coarsest_ln = 0, finest_ln = d_hierarchy->getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(in.getPatchHierarchy() == out.getPatchHierarchy());
    TBOX_ASSERT(in.getCoarsestLevelNumber() == out.getCoarsestLevelNumber());
    TBOX_ASSERT(in.getFinestLevelNumber() == out.getFinestLevelNumber());
#else
    NULL_USE(out);
#endif
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Setup the communication operators, algorithms, and schedules.
    //
    // TODO: Make the choice of refine and coarsen operators more flexible.
    d_Q_cell_refine_op = grid_geom->lookupRefineOperator(d_Q_cell_var, "CONSERVATIVE_LINEAR_REFINE");
    d_Q_cell_refine_alg = new RefineAlgorithm<NDIM>();
    d_Q_cell_refine_alg->registerRefine(
        d_Q_ghost_idx, in.getComponentDescriptorIndex(0), d_Q_scratch_idx, d_Q_cell_refine_op);
    if (d_outflow_bdry_extrap_type != "NONE")
    {
        d_Q_cell_refine_bdry_op = new CartExtrapPhysBdryOp(
            d_Q_scratch_idx, d_outflow_bdry_extrap_type); // TODO: check if this should be ghost or scratch!
    }

    d_q_flux_coarsen_op = grid_geom->lookupCoarsenOperator(d_q_flux_var, "CONSERVATIVE_COARSEN");
    d_q_flux_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_q_flux_coarsen_alg->registerCoarsen(d_q_flux_idx, d_q_flux_idx, d_q_flux_coarsen_op);

    d_q_interp_coarsen_op = grid_geom->lookupCoarsenOperator(d_q_interp_var, "CONSERVATIVE_COARSEN");
    d_q_interp_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_q_interp_coarsen_alg->registerCoarsen(d_q_interp_idx, d_q_interp_idx, d_q_interp_coarsen_op);

    d_Q_cell_refine_scheds.resize(finest_ln + 1);
    d_q_interp_coarsen_scheds.resize(finest_ln + 1);
    d_q_flux_coarsen_scheds.resize(finest_ln + 1);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > coarser_level = (ln > coarsest_ln) ? d_hierarchy->getPatchLevel(ln - 1) : nullptr;
        level->allocatePatchData(d_Q_ghost_idx);
        level->allocatePatchData(d_Q_scratch_idx);
        level->allocatePatchData(d_q_flux_idx);
        level->allocatePatchData(d_q_interp_idx);
        d_Q_cell_refine_scheds[ln] =
            d_Q_cell_refine_alg->createSchedule(level, ln - 1, d_hierarchy, d_Q_cell_refine_bdry_op);
        if (ln > coarsest_ln)
        {
            d_q_flux_coarsen_scheds[ln] = d_q_flux_coarsen_alg->createSchedule(coarser_level, level);
            d_q_interp_coarsen_scheds[ln] = d_q_interp_coarsen_alg->createSchedule(coarser_level, level);
        }
    }

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
CellConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate the refine algorithm, operator, patch strategy, and schedules.
    d_Q_cell_refine_op.setNull();
    d_Q_cell_refine_alg.setNull();
    d_Q_cell_refine_bdry_op.setNull();
    d_Q_cell_refine_scheds.clear();

    d_q_flux_coarsen_op.setNull();
    d_q_flux_coarsen_alg.setNull();
    d_q_flux_coarsen_scheds.clear();

    d_q_interp_coarsen_op.setNull();
    d_q_interp_coarsen_alg.setNull();
    d_q_interp_coarsen_scheds.clear();

    // Deallocate scratch data.
    const int coarsest_ln = 0, finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_Q_ghost_idx);
        level->deallocatePatchData(d_Q_scratch_idx);
        level->deallocatePatchData(d_q_flux_idx);
        level->deallocatePatchData(d_q_interp_idx);
    }

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
