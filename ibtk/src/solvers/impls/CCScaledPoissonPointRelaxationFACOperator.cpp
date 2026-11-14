// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#include <ibtk/CCScaledPoissonPETScLevelSolver.h>
#include <ibtk/CCScaledPoissonPointRelaxationFACOperator.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/CoarseFineBoundaryRefinePatchStrategy.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

#include <ArrayData.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellDataFactory.h>
#include <CellVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchDescriptor.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <cmath>
#include <string>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// FORTRAN INTERFACE /////////////////////////////

#if (NDIM == 2)
#define SCALED_GS_FC IBTK_FC_FUNC_(scaled_gs_2d, SCALED_GS_2D)
#elif (NDIM == 3)
#define SCALED_GS_FC IBTK_FC_FUNC_(scaled_gs_3d, SCALED_GS_3D)
#endif

extern "C"
{
    void SCALED_GS_FC(double* U,
                      const int& U_gcw,
                      const double* F,
                      const int& F_gcw,
                      const double* C,
                      const int& C_gcw,
                      const double& C_constant,
                      const int& C_is_variable,
                      const double* S,
                      const int& S_gcw,
                      const double& S_constant,
                      const int& S_is_variable,
                      const double* D0,
                      const double* D1,
#if (NDIM == 3)
                      const double* D2,
#endif
                      const int& D_gcw,
                      const double& D_constant,
                      const int& D_is_variable,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1,
#if (NDIM == 3)
                      const int& ilower2,
                      const int& iupper2,
#endif
                      const double* dx,
                      const int& color);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static Timer* t_smooth_error;
static Timer* t_compute_residual;

int
scalar_ghost_width(const IntVector<NDIM>& gcw)
{
#if !defined(NDEBUG)
    if (gcw.min() != gcw.max())
    {
        TBOX_ERROR("CCScaledPoissonPointRelaxationFACOperator:\n"
                   << "  Fortran relaxation requires uniform ghost widths in all directions.\n");
    }
#endif
    return gcw.max();
} // scalar_ghost_width
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCScaledPoissonPointRelaxationFACOperator::CCScaledPoissonPointRelaxationFACOperator(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : CCPoissonPointRelaxationFACOperator(object_name, input_db, default_options_prefix)
{
    // An ordinary unscaled coarse solver is not a suitable default.
    // Preserve an explicitly requested backend; otherwise use coarse relaxation.
    if (!input_db || !input_db->keyExists("coarse_solver_type"))
    {
        setCoarseSolverType("PATCH_GAUSS_SEIDEL");
    }

    if (!input_db || !input_db->keyExists("coarse_solver_max_iterations"))
    {
        d_coarse_solver_max_iterations = 30;
    }

    IBTK_DO_ONCE(t_smooth_error = TimerManager::getManager()->getTimer(
                     "IBTK::CCScaledPoissonPointRelaxationFACOperator::smoothError()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::CCScaledPoissonPointRelaxationFACOperator::computeResidual()"););

    return;
} // CCScaledPoissonPointRelaxationFACOperator

void
CCScaledPoissonPointRelaxationFACOperator::setScalePatchDataIndex(const int scale_idx)
{
#if !defined(NDEBUG)
    if (scale_idx < 0)
    {
        TBOX_ERROR(d_object_name << "::setScalePatchDataIndex():\n"
                                 << "  scale_idx must be a valid patch data index.\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellDataFactory<NDIM, double>> factory = var_db->getPatchDescriptor()->getPatchDataFactory(scale_idx);

    if (!factory || factory->getDefaultDepth() != 1)
    {
        TBOX_ERROR(d_object_name << "::setScalePatchDataIndex():\n"
                                 << "  scale must refer to scalar cell-centered double data.\n");
    }
#endif
    d_scale_idx = scale_idx;

    return;
} // setScalePatchDataIndex

void
CCScaledPoissonPointRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                       const SAMRAIVectorReal<NDIM, double>& residual,
                                                       const int level_num,
                                                       const int num_sweeps,
                                                       const bool /*performing_pre_sweeps*/,
                                                       const bool /*performing_post_sweeps*/)
{
    if (num_sweeps <= 0) return;

    IBTK_TIMER_START(t_smooth_error);

    const int error_idx = error.getComponentDescriptorIndex(0);
    const int residual_idx = residual.getComponentDescriptorIndex(0);

    const int C_is_variable = d_poisson_spec.cIsVariable() ? 1 : 0;
    const int S_is_variable = d_scale_idx != invalid_index ? 1 : 0;
    const int D_is_variable = d_poisson_spec.dIsConstant() ? 0 : 1;

    const double C_constant = (C_is_variable || d_poisson_spec.cIsZero()) ? 0.0 : d_poisson_spec.getCConstant();
    const double D_constant = D_is_variable ? 0.0 : d_poisson_spec.getDConstant();

    const std::string& smoother_type = level_num == d_coarsest_ln ? d_coarse_solver_type : d_smoother_type;

    const bool red_black = smoother_type == "RED_BLACK_GAUSS_SEIDEL";
    const bool update_local_data = red_black || smoother_type == "PROCESSOR_GAUSS_SEIDEL";
    const int num_passes = red_black ? 2 * num_sweeps : num_sweeps;

    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);

    // Cache incoming coarse-fine ghost values using the inherited overlap maps.
    if (level_num > d_coarsest_ln && num_passes > 1)
    {
        int patch_counter = 0;

        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> error_data = patch->getPatchData(error_idx);
            Pointer<CellData<NDIM, double>> scratch_data = patch->getPatchData(d_scratch_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
            scratch_data->getArrayData().copy(
                error_data->getArrayData(), d_patch_bc_box_overlap[level_num][patch_counter], IntVector<NDIM>(0));
        }
    }

    for (int pass = 0; pass < num_passes; ++pass)
    {
        if (level_num > d_coarsest_ln)
        {
            if (pass > 0)
            {
                // Copy the coarse-fine interface ghost cell values which are
                // cached in the scratch data into the error data.
                int patch_counter = 0;
                for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double>> error_data = patch->getPatchData(error_idx);
                    Pointer<CellData<NDIM, double>> scratch_data = patch->getPatchData(d_scratch_idx);
#if !defined(NDEBUG)
                    const Box<NDIM>& ghost_box = error_data->getGhostBox();
                    TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
                    TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    error_data->getArrayData().copy(scratch_data->getArrayData(),
                                                    d_patch_bc_box_overlap[level_num][patch_counter],
                                                    IntVector<NDIM>(0));
                }

                // Fill the non-coarse-fine interface ghost cell values.
                xeqScheduleGhostFillNoCoarse(error_idx, level_num);
            }

            d_cf_bdry_op->setPatchDataIndex(error_idx);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();

            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                d_cf_bdry_op->computeNormalExtension(*patch, ratio, d_gcw);
            }
        }
        else if (pass > 0)
        {
            xeqScheduleGhostFillNoCoarse(error_idx, level_num);
        }

        const int color = red_black ? pass % 2 : -1;

        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> error_data = patch->getPatchData(error_idx);
            Pointer<CellData<NDIM, double>> residual_data = patch->getPatchData(residual_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == residual_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(error_data->getDepth() == residual_data->getDepth());
#endif

            // Copy updated values from neighboring local patches.
            if (update_local_data)
            {
                const std::map<int, Box<NDIM>> neighbor_overlap = d_patch_neighbor_overlap[level_num][patch_counter];
                for (const auto& pair : neighbor_overlap)
                {
                    const int src_patch_num = pair.first;
                    const Box<NDIM>& overlap = pair.second;
                    Pointer<Patch<NDIM>> src_patch = level->getPatch(src_patch_num);
                    Pointer<CellData<NDIM, double>> src_error_data = error.getComponentPatchData(0, *src_patch);
                    error_data->getArrayData().copy(src_error_data->getArrayData(), overlap, IntVector<NDIM>(0));
                }
            }

            Pointer<CellData<NDIM, double>> C_data;
            Pointer<CellData<NDIM, double>> S_data;
            Pointer<SideData<NDIM, double>> D_data;

            if (C_is_variable) C_data = patch->getPatchData(d_poisson_spec.getCPatchDataId());
            if (S_is_variable) S_data = patch->getPatchData(d_scale_idx);
            if (D_is_variable) D_data = patch->getPatchData(d_poisson_spec.getDPatchDataId());

            if ((C_is_variable && (!C_data || C_data->getDepth() != 1)) ||
                (S_is_variable && (!S_data || S_data->getDepth() != 1)) ||
                (D_is_variable && (!D_data || D_data->getDepth() != 1)))
            {
                TBOX_ERROR(d_object_name << "::smoothError():\n"
                                         << "  scalar coefficient data must be allocated on every patch.\n");
            }

            const int U_gcw = scalar_ghost_width(error_data->getGhostCellWidth());
            const int F_gcw = scalar_ghost_width(residual_data->getGhostCellWidth());
            const int C_gcw = C_data ? scalar_ghost_width(C_data->getGhostCellWidth()) : U_gcw;
            const int S_gcw = S_data ? scalar_ghost_width(S_data->getGhostCellWidth()) : U_gcw;
            const int D_gcw = D_data ? scalar_ghost_width(D_data->getGhostCellWidth()) : 0;

            // Constant-coefficient branches do not dereference the dummy arrays.
            const double* D0 = D_data ? D_data->getPointer(0, 0) : &D_constant;
            const double* D1 = D_data ? D_data->getPointer(1, 0) : &D_constant;
#if (NDIM == 3)
            const double* D2 = D_data ? D_data->getPointer(2, 0) : &D_constant;
#endif

            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* dx = pgeom->getDx();

            for (int depth = 0; depth < error_data->getDepth(); ++depth)
            {
                double* U = error_data->getPointer(depth);
                const double* F = residual_data->getPointer(depth);
                const double* C = C_data ? C_data->getPointer(0) : U;
                const double* S = S_data ? S_data->getPointer(0) : U;

                SCALED_GS_FC(U,
                             U_gcw,
                             F,
                             F_gcw,
                             C,
                             C_gcw,
                             C_constant,
                             C_is_variable,
                             S,
                             S_gcw,
                             d_scale_constant,
                             S_is_variable,
                             D0,
                             D1,
#if (NDIM == 3)
                             D2,
#endif
                             D_gcw,
                             D_constant,
                             D_is_variable,
                             patch_box.lower(0),
                             patch_box.upper(0),
                             patch_box.lower(1),
                             patch_box.upper(1),
#if (NDIM == 3)
                             patch_box.lower(2),
                             patch_box.upper(2),
#endif
                             dx,
                             color);
            }
        }
    }

    IBTK_TIMER_STOP(t_smooth_error);

    return;
} // smoothError

void
CCScaledPoissonPointRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                           const SAMRAIVectorReal<NDIM, double>& solution,
                                                           const SAMRAIVectorReal<NDIM, double>& rhs,
                                                           const int coarsest_level_num,
                                                           const int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    const int sol_idx = solution.getComponentDescriptorIndex(0);
    const int rhs_idx = rhs.getComponentDescriptorIndex(0);
    const int res_idx = residual.getComponentDescriptorIndex(0);

    Pointer<CellVariable<NDIM, double>> sol_var = solution.getComponentVariable(0);

    // Reuse the parent's residual ghost-fill infrastructure.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    Pointer<VariableFillPattern<NDIM>> fill_pattern = nullptr;
    if (d_poisson_spec.dIsConstant())
    {
        fill_pattern = new CellNoCornersFillPattern(1, /*overwrite_interior*/ false);
    }

    InterpolationTransactionComponent transaction_comp(sol_idx,
                                                       d_data_refine_type,
                                                       d_use_cf_interpolation,
                                                       d_data_coarsen_type,
                                                       d_bdry_extrap_type,
                                                       d_use_consistent_type_2_bdry,
                                                       d_bc_coefs,
                                                       fill_pattern);

    if (d_level_bdry_fill_ops[finest_level_num])
    {
        d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponent(transaction_comp);
    }
    else
    {
        d_level_bdry_fill_ops[finest_level_num] = new HierarchyGhostCellInterpolation();
        d_level_bdry_fill_ops[finest_level_num]->initializeOperatorState(
            transaction_comp, d_hierarchy, coarsest_level_num, finest_level_num);
    }

    d_level_bdry_fill_ops[finest_level_num]->setHomogeneousBc(true);
    d_level_bdry_fill_ops[finest_level_num]->fillData(d_solution_time);

    InterpolationTransactionComponent default_transaction_comp(d_solution->getComponentDescriptorIndex(0),
                                                               d_data_refine_type,
                                                               d_use_cf_interpolation,
                                                               d_data_coarsen_type,
                                                               d_bdry_extrap_type,
                                                               d_use_consistent_type_2_bdry,
                                                               d_bc_coefs,
                                                               fill_pattern);

    d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponent(default_transaction_comp);

    // Reuse the parent's hierarchy math operations cache.
    if (!d_level_math_ops[finest_level_num])
    {
        d_level_math_ops[finest_level_num] =
            new HierarchyMathOps(d_object_name + "::hier_math_ops_" + std::to_string(finest_level_num),
                                 d_hierarchy,
                                 coarsest_level_num,
                                 finest_level_num);
    }

    for (unsigned int depth = 0; depth < d_bc_coefs.size(); ++depth)
    {
        d_level_math_ops[finest_level_num]->scaled_laplace(d_scratch_idx,
                                                           sol_var,
                                                           d_poisson_spec,
                                                           d_scale_idx,
                                                           sol_idx,
                                                           sol_var,
                                                           nullptr,
                                                           d_solution_time,
                                                           0.0,
                                                           invalid_index,
                                                           nullptr,
                                                           depth,
                                                           depth,
                                                           0,
                                                           0);
    }

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_cc_data_ops.subtract(res_idx, rhs_idx, d_scratch_idx);

    IBTK_TIMER_STOP(t_compute_residual);

    return;
} // computeResidual

/////////////////////////////// PROTECTED ////////////////////////////////////

void
CCScaledPoissonPointRelaxationFACOperator::configureCoarseSolverSpecialized()
{
    if (!d_coarse_solver) return;

    Pointer<CCScaledPoissonPETScLevelSolver> scaled_solver = d_coarse_solver;

    if (!scaled_solver)
    {
        TBOX_ERROR(
            d_object_name << "::configureCoarseSolverSpecialized():\n"
                          << "  the selected coarse solver is not of the type CCScaledPoissonPETScLevelSolver.\n");
    }

    if (d_scale_idx != invalid_index) scaled_solver->setScalePatchDataIndex(d_scale_idx);

    return;
} // configureCoarseSolverSpecialized

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

#undef SCALED_GS_FC
