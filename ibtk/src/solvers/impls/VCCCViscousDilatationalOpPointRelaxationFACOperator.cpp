// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
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

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartCellDoubleCubicCoarsen.h"
#include "ibtk/CartCellDoubleQuadraticCFInterpolation.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/VCCCViscousDilatationalOpPointRelaxationFACOperator.h"
#include "ibtk/ibtk_utilities.h"

#include "ArrayData.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CoarsenOperator.h"
#include "HierarchyCellDataOpsReal.h"
#include "Patch.h"
#include "PatchDescriptor.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "ProcessorMapping.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define VC_CC_GS_SMOOTH_FC IBTK_FC_FUNC(vcccgssmooth2d, VCGSSMOOTH2D)
#define VC_CC_RB_GS_SMOOTH_FC IBTK_FC_FUNC(vcccrbgssmooth2d, VCRBGSSMOOTH2D)
#endif
#if (NDIM == 3)
#define VC_CC_GS_SMOOTH_FC IBTK_FC_FUNC(vcccgssmooth3d, VCGSSMOOTH3D)
#define VC_CC_RB_GS_SMOOTH_FC IBTK_FC_FUNC(vcccrbgssmooth3d, VCRBGSSMOOTH3D)
#endif

// Function interfaces
extern "C"
{
    void VC_CC_GS_SMOOTH_FC(double* U0,
                            double* U1,
#if (NDIM == 3)
                            double* U2,
#endif
                            const int& U_gcw,
                            const double* F0,
                            const double* F1,
#if (NDIM == 3)
                            const double* F2,
#endif
                            const int& F_gcw,
                            const double* C0,
                            const double* C1,
#if (NDIM == 3)
                            const double* C2,
#endif
                            const int& C_gcw,
                            const double* mu,
                            const int& mu_gcw,
                            const double* lambda,
                            const int& lambda_gcw,
                            const int& ilower0,
                            const int& iupper0,
                            const int& ilower1,
                            const int& iupper1,
#if (NDIM == 3)
                            const int& ilower2,
                            const int& iupper2,
#endif
                            const double* dx,
                            const int& is_C_const,
                            const double& C_const,
                            const int& use_harmonic_interp);

    void VC_CC_RB_GS_SMOOTH_FC(double* U0,
                               double* U1,
#if (NDIM == 3)
                               double* U2,
#endif
                               const int& U_gcw,
                               const double* F0,
                               const double* F1,
#if (NDIM == 3)
                               const double* F2,
#endif
                               const int& F_gcw,
                               const double* C0,
                               const double* C1,
#if (NDIM == 3)
                               const double* C2,
#endif
                               const int& C_gcw,
                               const double* mu,
                               const int& mu_gcw,
                               const double* lambda,
                               const int& lambda_gcw,
                               const int& ilower0,
                               const int& iupper0,
                               const int& ilower1,
                               const int& iupper1,
#if (NDIM == 3)
                               const int& ilower2,
                               const int& iupper2,
#endif
                               const double* dx,
                               const int& is_C_const,
                               const double& C_const,
                               const int& use_harmonic_interp,
                               const int& red_or_black);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_restrict_residual;
static Timer* t_smooth_error;
static Timer* t_solve_coarsest_level;
static Timer* t_compute_residual;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "LINEAR_REFINE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries; used only to evaluate
// composite grid residuals.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells; used only to evaluate composite grid residuals.
static const bool CONSISTENT_TYPE_2_BDRY = false;

enum SmootherType
{
    PATCH_GAUSS_SEIDEL,
    PROCESSOR_GAUSS_SEIDEL,
    RED_BLACK_GAUSS_SEIDEL,
    UNKNOWN = -1
};

inline SmootherType
get_smoother_type(const std::string& smoother_type_string)
{
    if (smoother_type_string == "PATCH_GAUSS_SEIDEL") return PATCH_GAUSS_SEIDEL;
    if (smoother_type_string == "PROCESSOR_GAUSS_SEIDEL") return PROCESSOR_GAUSS_SEIDEL;
    if (smoother_type_string == "RED_BLACK_GAUSS_SEIDEL")
        return RED_BLACK_GAUSS_SEIDEL;
    else
        return UNKNOWN;
} // get_smoother_type

inline bool
use_red_black_ordering(SmootherType smoother_type)
{
    if (smoother_type == RED_BLACK_GAUSS_SEIDEL)
    {
        return true;
    }
    else
    {
        return false;
    }
} // use_red_black_ordering

inline bool
do_local_data_update(SmootherType smoother_type)
{
    if (smoother_type == PROCESSOR_GAUSS_SEIDEL || smoother_type == RED_BLACK_GAUSS_SEIDEL)
    {
        return true;
    }
    else
    {
        return false;
    }
} // do_local_data_update
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

VCCCViscousDilatationalOpPointRelaxationFACOperator::VCCCViscousDilatationalOpPointRelaxationFACOperator(
    std::string object_name,
    const Pointer<Database> input_db,
    std::string default_options_prefix)
    : CCPoissonPointRelaxationFACOperator(std::move(object_name), input_db, std::move(default_options_prefix))
{
    // Setup Timers.
    IBTK_DO_ONCE(t_restrict_residual = TimerManager::getManager()->getTimer(
                     "IBTK::VCCCViscousDilatationalOpPointRelaxationFACOperator::restrictResidual()");
                 t_smooth_error = TimerManager::getManager()->getTimer(
                     "IBTK::VCCCViscousDilatationalOpPointRelaxationFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "IBTK::VCCCViscousDilatationalOpPointRelaxationFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::VCCCViscousDilatationalOpPointRelaxationFACOperator::computeResidual()"););

    // Set a default interpolation type.
    d_D_interp_type = IBTK::VC_HARMONIC_INTERP;

    return;
} // VCCCViscousDilatationalOpPointRelaxationFACOperator

VCCCViscousDilatationalOpPointRelaxationFACOperator::~VCCCViscousDilatationalOpPointRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~VCCCViscousDilatationalOpPointRelaxationFACOperator

void
VCCCViscousDilatationalOpPointRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                                 const SAMRAIVectorReal<NDIM, double>& residual,
                                                                 int level_num,
                                                                 int num_sweeps,
                                                                 bool /*performing_pre_sweeps*/,
                                                                 bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBTK_TIMER_START(t_smooth_error);

    const auto& vc_op_spec = static_cast<const VCViscousDilatationalOpSpec&>(*d_problem_spec);

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int error_idx = error.getComponentDescriptorIndex(0);
    const int scratch_idx = d_scratch_idx;

    // Determine the smoother type.
    const std::string& smoother_type_string = (level_num == d_coarsest_ln ? d_coarse_solver_type : d_smoother_type);
    const SmootherType smoother_type = get_smoother_type(smoother_type_string);
#if !defined(NDEBUG)
    TBOX_ASSERT(smoother_type != UNKNOWN);
#endif
    const bool red_black_ordering = use_red_black_ordering(smoother_type);
    const bool update_local_data = do_local_data_update(smoother_type);
    const bool use_harmonic_interp = (d_D_interp_type == VC_HARMONIC_INTERP);

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
            Pointer<CellData<NDIM, double> > scratch_data = patch->getPatchData(scratch_idx);
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

    // Smooth the error by the specified number of sweeps.
    if (red_black_ordering) num_sweeps *= 2;
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        // Re-fill ghost cell data as needed.
        if (level_num > d_coarsest_ln)
        {
            if (isweep > 0)
            {
                // Copy the coarse-fine interface ghost cell values which are
                // cached in the scratch data into the error data.
                int patch_counter = 0;
                for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
                    Pointer<CellData<NDIM, double> > scratch_data = patch->getPatchData(scratch_idx);
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

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_cf_bdry_op->setPatchDataIndex(error_idx);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const IntVector<NDIM>& ghost_width_to_fill = d_gcw;
                d_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }
        else if (isweep > 0)
        {
            xeqScheduleGhostFillNoCoarse(error_idx, level_num);
        }

        // Smooth the error on the patches.
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
            Pointer<CellData<NDIM, double> > residual_data = residual.getComponentPatchData(0, *patch);
#if !defined(NDEBUG)
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == residual_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(error_data->getDepth() == residual_data->getDepth());
#endif
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            // Copy updated values from neighboring local patches.
            if (update_local_data)
            {
                const std::map<int, Box<NDIM> > neighbor_overlap = d_patch_neighbor_overlap[level_num][patch_counter];
                for (const auto& pair : neighbor_overlap)
                {
                    const int src_patch_num = pair.first;
                    const Box<NDIM>& overlap = pair.second;
                    Pointer<Patch<NDIM> > src_patch = level->getPatch(src_patch_num);
                    Pointer<CellData<NDIM, double> > src_error_data = error.getComponentPatchData(0, *src_patch);
                    error_data->getArrayData().copy(src_error_data->getArrayData(), overlap, IntVector<NDIM>(0));
                }
            }

            // Smooth the error using Gauss-Seidel.
            const bool C_is_const = vc_op_spec.d_C_is_const;
            const double C_const = vc_op_spec.d_C_const;
            const int C_idx = vc_op_spec.d_C_idx;
            int C_data_depth = 0;
            Pointer<CellData<NDIM, double> > C_data = nullptr;
            double* C0 = nullptr;
            double* C1 = nullptr;
            double* C2 = nullptr;
            int C_gcw = 0;
            if (!C_is_const)
            {
                C_data = patch->getPatchData(C_idx);
                C_data_depth = C_data->getDepth();
                C_gcw = (C_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
                TBOX_ASSERT(C_data_depth == 1 || C_data_depth == NDIM);
#endif
                if (C_data_depth == 1)
                {
                    C0 = C_data->getPointer(0);
                    C1 = C0;
#if (NDIM == 3)
                    C2 = C0;
#endif
                }
                else
                {
                    C0 = C_data->getPointer(0);
                    C1 = C_data->getPointer(1);
#if (NDIM == 3)
                    C2 = C_data->getPointer(2);
#endif
                }
            }

            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(vc_op_spec.d_D_idx);
            Pointer<CellData<NDIM, double> > lambda_data = patch->getPatchData(vc_op_spec.d_L_idx);

#if !defined(NDEBUG)
            TBOX_ASSERT(!mu_data.isNull());
            TBOX_ASSERT(!lambda_data.isNull());
            if (!C_is_const) TBOX_ASSERT(!C_data.isNull());
#endif

            double* const U0 = error_data->getPointer(0);
            double* const U1 = error_data->getPointer(1);
#if (NDIM == 3)
            double* const U2 = error_data->getPointer(2);
#endif
            const int U_gcw = (error_data->getGhostCellWidth()).max();

            const double* const F0 = residual_data->getPointer(0);
            const double* const F1 = residual_data->getPointer(1);
#if (NDIM == 3)
            const double* const F2 = residual_data->getPointer(2);
#endif
            const int F_gcw = (residual_data->getGhostCellWidth()).max();

            const double* const mu = mu_data->getPointer();
            const int mu_gcw = (mu_data->getGhostCellWidth()).max();

            const double* const lambda = lambda_data->getPointer();
            const int lambda_gcw = (lambda_data->getGhostCellWidth()).max();

            if (red_black_ordering)
            {
                int red_or_black = isweep % 2; // "red" = 0, "black" = 1
                VC_CC_RB_GS_SMOOTH_FC(U0,
                                      U1,
#if (NDIM == 3)
                                      U2,
#endif
                                      U_gcw,
                                      F0,
                                      F1,
#if (NDIM == 3)
                                      F2,
#endif
                                      F_gcw,
                                      C0,
                                      C1,
#if (NDIM == 3)
                                      C2,
#endif
                                      C_gcw,
                                      mu,
                                      mu_gcw,
                                      lambda,
                                      lambda_gcw,
                                      patch_box.lower(0),
                                      patch_box.upper(0),
                                      patch_box.lower(1),
                                      patch_box.upper(1),
#if (NDIM == 3)
                                      patch_box.lower(2),
                                      patch_box.upper(2),
#endif
                                      dx,
                                      C_is_const,
                                      C_const,
                                      use_harmonic_interp,
                                      red_or_black);
            }
            else
            {
                VC_CC_GS_SMOOTH_FC(U0,
                                   U1,
#if (NDIM == 3)
                                   U2,
#endif
                                   U_gcw,
                                   F0,
                                   F1,
#if (NDIM == 3)
                                   F2,
#endif
                                   F_gcw,
                                   C0,
                                   C1,
#if (NDIM == 3)
                                   C2,
#endif
                                   C_gcw,
                                   mu,
                                   mu_gcw,
                                   lambda,
                                   lambda_gcw,
                                   patch_box.lower(0),
                                   patch_box.upper(0),
                                   patch_box.lower(1),
                                   patch_box.upper(1),
#if (NDIM == 3)
                                   patch_box.lower(2),
                                   patch_box.upper(2),
#endif
                                   dx,
                                   C_is_const,
                                   C_const,
                                   use_harmonic_interp);
            }
        }
    }
    IBTK_TIMER_STOP(t_smooth_error);
    return;
} // smoothError

void
VCCCViscousDilatationalOpPointRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                                     const SAMRAIVectorReal<NDIM, double>& solution,
                                                                     const SAMRAIVectorReal<NDIM, double>& rhs,
                                                                     int coarsest_level_num,
                                                                     int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    auto& vc_op_spec = static_cast<const VCViscousDilatationalOpSpec&>(*d_problem_spec);

    const int res_idx = residual.getComponentDescriptorIndex(0);
    const int sol_idx = solution.getComponentDescriptorIndex(0);
    const int rhs_idx = rhs.getComponentDescriptorIndex(0);

    const Pointer<CellVariable<NDIM, double> > res_var = residual.getComponentVariable(0);
    const Pointer<CellVariable<NDIM, double> > sol_var = solution.getComponentVariable(0);

    // Fill ghost-cell values.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    Pointer<CellNoCornersFillPattern> fill_pattern = nullptr;
    InterpolationTransactionComponent transaction_comp(sol_idx,
                                                       DATA_REFINE_TYPE,
                                                       USE_CF_INTERPOLATION,
                                                       DATA_COARSEN_TYPE,
                                                       BDRY_EXTRAP_TYPE,
                                                       CONSISTENT_TYPE_2_BDRY,
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
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_bc_coefs,
                                                               fill_pattern);
    d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponent(default_transaction_comp);

    // Compute the residual, r = f - A*u.
    if (!d_level_math_ops[finest_level_num])
    {
        d_level_math_ops[finest_level_num] =
            new HierarchyMathOps(d_object_name + "::hier_math_ops_" + std::to_string(finest_level_num),
                                 d_hierarchy,
                                 coarsest_level_num,
                                 finest_level_num);
    }

    // A_U = (C*I + L1(mu))*U
    double alpha = 1.0;
    double beta = 1.0;
    if (vc_op_spec.d_C_is_const)
    {
        beta = vc_op_spec.d_C_const;
    }
    d_level_math_ops[finest_level_num]->vc_laplace(res_idx,
                                                   res_var,
                                                   alpha,
                                                   beta,
                                                   vc_op_spec.d_D_idx,
                                                   Pointer<CellVariable<NDIM, double> >(nullptr),
                                                   sol_idx,
                                                   sol_var,
                                                   nullptr,
                                                   d_solution_time,
                                                   d_D_interp_type,
                                                   vc_op_spec.d_C_idx,
                                                   Pointer<CellVariable<NDIM, double> >(nullptr));

    // A_U += L2(lambda)*U
    if (vc_op_spec.d_L_idx > 0)
    {
        d_level_math_ops[finest_level_num]->vc_dilatational(res_idx,
                                                            res_var,
                                                            1.0,
                                                            vc_op_spec.d_L_idx,
                                                            Pointer<CellVariable<NDIM, double> >(nullptr),
                                                            sol_idx,
                                                            sol_var,
                                                            nullptr,
                                                            d_solution_time,
                                                            1.0,
                                                            res_idx,
                                                            res_var);
    }

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_cc_data_ops.axpy(res_idx, -1.0, res_idx, rhs_idx, false);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
} // computeResidual

void
VCCCViscousDilatationalOpPointRelaxationFACOperator::setDPatchDataInterpolationType(
    const IBTK::VCInterpType D_interp_type)
{
    d_D_interp_type = D_interp_type;
    return;
} // setDPatchDataInterpolationType

Pointer<PoissonSolver>
VCCCViscousDilatationalOpPointRelaxationFACOperator::getCoarseSolver()
{
    return d_coarse_solver;
} // getCoarseSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

void
VCCCViscousDilatationalOpPointRelaxationFACOperator::initializeOperatorStateSpecialized(
    const SAMRAIVectorReal<NDIM, double>& solution,
    const SAMRAIVectorReal<NDIM, double>& rhs,
    const int coarsest_reset_ln,
    const int finest_reset_ln)
{
    // Set the problem specification object for the coarse solver
    if (coarsest_reset_ln == d_coarsest_ln && d_coarse_solver)
    {
        d_coarse_solver->setProblemSpecification(d_problem_spec);
    }

    CCPoissonPointRelaxationFACOperator::initializeOperatorStateSpecialized(
        solution, rhs, coarsest_reset_ln, finest_reset_ln);

    // Set stencil fill pattern to nullptr if the base class
    // has set it to non-null.
    d_op_stencil_fill_pattern = nullptr;

    return;
} // initializeOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
