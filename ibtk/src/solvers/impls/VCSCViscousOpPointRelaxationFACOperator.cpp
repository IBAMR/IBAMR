// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
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

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/SCPoissonPointRelaxationFACOperator.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "ibtk/VCSCViscousOpPointRelaxationFACOperator.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/ibtk_utilities.h"

#include "ArrayData.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "EdgeData.h"
#include "EdgeVariable.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "MultiblockDataTranslator.h"
#include "NodeData.h"
#include "NodeVariable.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "VariableFillPattern.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"

#include <map>
#include <memory>

#include "ibtk/namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define VC_GS_SMOOTH_FC IBTK_FC_FUNC(vcgssmooth2d, VCGSSMOOTH2D)
#define VC_GS_SMOOTH_MASK_FC IBTK_FC_FUNC(vcgssmoothmask2d, VCGSSMOOTHMASK2D)
#define VC_RB_GS_SMOOTH_FC IBTK_FC_FUNC(vcrbgssmooth2d, VCRBGSSMOOTH2D)
#define VC_RB_GS_SMOOTH_MASK_FC IBTK_FC_FUNC(vcrbgssmoothmask2d, VCRBGSSMOOTHMASK2D)
#endif
#if (NDIM == 3)
#define VC_GS_SMOOTH_FC IBTK_FC_FUNC(vcgssmooth3d, VCGSSMOOTH3D)
#define VC_GS_SMOOTH_MASK_FC IBTK_FC_FUNC(vcgssmoothmask3d, VCGSSMOOTHMASK3D)
#define VC_RB_GS_SMOOTH_FC IBTK_FC_FUNC(vcrbgssmooth3d, VCRBGSSMOOTH3D)
#define VC_RB_GS_SMOOTH_MASK_FC IBTK_FC_FUNC(vcrbgssmoothmask3d, VCRBGSSMOOTHMASK3D)
#endif

// Function interfaces
extern "C"
{
    void VC_GS_SMOOTH_FC(double* U0,
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
#if (NDIM == 2)
                         const double* mu,
#endif
#if (NDIM == 3)
                         const double* mu0,
                         const double* mu1,
                         const double* mu2,
#endif
                         const int& mu_gcw,
                         const double& alpha,
                         const double& beta,
                         const int& ilower0,
                         const int& iupper0,
                         const int& ilower1,
                         const int& iupper1,
#if (NDIM == 3)
                         const int& ilower2,
                         const int& iupper2,
#endif
                         const double* dx,
                         const int& var_c,
                         const int& use_harmonic_interp);

    void VC_GS_SMOOTH_MASK_FC(double* U0,
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
                              const int* mask0,
                              const int* mask1,
#if (NDIM == 3)
                              const int* mask2,
#endif
                              const int& mask_gcw,
                              const double* C0,
                              const double* C1,
#if (NDIM == 3)
                              const double* C2,
#endif
                              const int& C_gcw,
#if (NDIM == 2)
                              const double* mu,
#endif
#if (NDIM == 3)
                              const double* mu0,
                              const double* mu1,
                              const double* mu2,
#endif
                              const int& mu_gcw,
                              const double& alpha,
                              const double& beta,
                              const int& ilower0,
                              const int& iupper0,
                              const int& ilower1,
                              const int& iupper1,
#if (NDIM == 3)
                              const int& ilower2,
                              const int& iupper2,
#endif
                              const double* dx,
                              const int& var_c,
                              const int& use_harmonic_interp);

    void VC_RB_GS_SMOOTH_FC(double* U0,
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
#if (NDIM == 2)
                            const double* mu,
#endif
#if (NDIM == 3)
                            const double* mu0,
                            const double* mu1,
                            const double* mu2,
#endif
                            const int& mu_gcw,
                            const double& alpha,
                            const double& beta,
                            const int& ilower0,
                            const int& iupper0,
                            const int& ilower1,
                            const int& iupper1,
#if (NDIM == 3)
                            const int& ilower2,
                            const int& iupper2,
#endif
                            const double* dx,
                            const int& var_c,
                            const int& use_harmonic_interp,
                            const int& red_or_black);

    void VC_RB_GS_SMOOTH_MASK_FC(double* U0,
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
                                 const int* mask0,
                                 const int* mask1,
#if (NDIM == 3)
                                 const int* mask2,
#endif
                                 const int& mask_gcw,
                                 const double* C0,
                                 const double* C1,
#if (NDIM == 3)
                                 const double* C2,
#endif
                                 const int& C_gcw,
#if (NDIM == 2)
                                 const double* mu,
#endif
#if (NDIM == 3)
                                 const double* mu0,
                                 const double* mu1,
                                 const double* mu2,
#endif
                                 const int& mu_gcw,
                                 const double& alpha,
                                 const double& beta,
                                 const int& ilower0,
                                 const int& iupper0,
                                 const int& ilower1,
                                 const int& iupper1,
#if (NDIM == 3)
                                 const int& ilower2,
                                 const int& iupper2,
#endif
                                 const double* dx,
                                 const int& var_c,
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
static const std::string DATA_REFINE_TYPE = "NONE";
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

VCSCViscousOpPointRelaxationFACOperator::VCSCViscousOpPointRelaxationFACOperator(std::string object_name,
                                                                                 const Pointer<Database> input_db,
                                                                                 std::string default_options_prefix)
    : SCPoissonPointRelaxationFACOperator(std::move(object_name), input_db, std::move(default_options_prefix))
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_restrict_residual =
            TimerManager::getManager()->getTimer("IBTK::VCSCViscousOpPointRelaxationFACOperator::restrictResidual()");
        t_smooth_error =
            TimerManager::getManager()->getTimer("IBTK::VCSCViscousOpPointRelaxationFACOperator::smoothError()");
        t_solve_coarsest_level =
            TimerManager::getManager()->getTimer("IBTK::VCSCViscousOpPointRelaxationFACOperator::solveCoarsestLevel()");
        t_compute_residual =
            TimerManager::getManager()->getTimer("IBTK::VCSCViscousOpPointRelaxationFACOperator::computeResidual()"););

    // Set a default interpolation type.
    d_D_interp_type = IBTK::VC_HARMONIC_INTERP;

    return;
} // VCSCViscousOpPointRelaxationFACOperator

VCSCViscousOpPointRelaxationFACOperator::~VCSCViscousOpPointRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~VCSCViscousOpPointRelaxationFACOperator

void
VCSCViscousOpPointRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                     const SAMRAIVectorReal<NDIM, double>& residual,
                                                     int level_num,
                                                     int num_sweeps,
                                                     bool /*performing_pre_sweeps*/,
                                                     bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBTK_TIMER_START(t_smooth_error);

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

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double> > scratch_data = patch->getPatchData(scratch_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                scratch_data->getArrayData(axis).copy(error_data->getArrayData(axis),
                                                      d_patch_bc_box_overlap[level_num][patch_counter][axis],
                                                      IntVector<NDIM>(0));
            }
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
                    Pointer<SideData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM, double> > scratch_data = patch->getPatchData(scratch_idx);
#if !defined(NDEBUG)
                    const Box<NDIM>& ghost_box = error_data->getGhostBox();
                    TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
                    TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        error_data->getArrayData(axis).copy(scratch_data->getArrayData(axis),
                                                            d_patch_bc_box_overlap[level_num][patch_counter][axis],
                                                            IntVector<NDIM>(0));
                    }
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
            Pointer<SideData<NDIM, double> > error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double> > residual_data = residual.getComponentPatchData(0, *patch);
#if (NDIM == 2)
            Pointer<NodeData<NDIM, double> > mu_data = patch->getPatchData(d_poisson_spec.getDPatchDataId());
#endif
#if (NDIM == 3)
            Pointer<EdgeData<NDIM, double> > mu_data = patch->getPatchData(d_poisson_spec.getDPatchDataId());
#endif
#if !defined(NDEBUG)
            const Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == residual_data->getGhostBox());
            TBOX_ASSERT(ghost_box == mu_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(mu_data->getGhostCellWidth() >= d_gcw);
            TBOX_ASSERT(error_data->getDepth() == residual_data->getDepth());
            TBOX_ASSERT(error_data->getDepth() == mu_data->getDepth());
#endif
            Pointer<SideData<NDIM, double> > C_data = nullptr;
            if (d_poisson_spec.cIsVariable())
            {
                C_data = patch->getPatchData(d_poisson_spec.getCPatchDataId());
            }

            Pointer<SideData<NDIM, int> > mask_data = patch->getPatchData(d_mask_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            // Copy updated values from neighboring local patches.
            if (update_local_data)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const std::map<int, Box<NDIM> > neighbor_overlap =
                        d_patch_neighbor_overlap[level_num][patch_counter][axis];
                    for (const auto& pair : neighbor_overlap)
                    {
                        const int src_patch_num = pair.first;
                        const Box<NDIM>& overlap = pair.second;
                        Pointer<Patch<NDIM> > src_patch = level->getPatch(src_patch_num);
                        Pointer<SideData<NDIM, double> > src_error_data = error.getComponentPatchData(0, *src_patch);
                        error_data->getArrayData(axis).copy(
                            src_error_data->getArrayData(axis), overlap, IntVector<NDIM>(0));
                    }
                }
            }

            // Enforce any Dirichlet boundary conditions.
            const bool patch_has_dirichlet_bdry = d_bc_helper->patchTouchesDirichletBoundary(patch);
            if (patch_has_dirichlet_bdry)
            {
                d_bc_helper->copyDataAtDirichletBoundaries(error_data, residual_data, patch);
            }

            // What type of averaging to use for the patch smoothers
            const bool use_harmonic_interp = (d_D_interp_type == VC_HARMONIC_INTERP);

            // Smooth the error using Gauss-Seidel.
            double alpha = 1.0;
            double beta = 0.0;
            const int C_is_varying = d_poisson_spec.cIsVariable();
            if (d_poisson_spec.cIsConstant())
            {
                beta = d_poisson_spec.getCConstant();
            }
            else if (d_poisson_spec.cIsVariable())
            {
                beta = 1.0;
            }

            for (int depth = 0; depth < error_data->getDepth(); ++depth)
            {
                double* const U0 = error_data->getPointer(0, depth);
                double* const U1 = error_data->getPointer(1, depth);
#if (NDIM == 3)
                double* const U2 = error_data->getPointer(2, depth);
#endif
                const int U_ghosts = (error_data->getGhostCellWidth()).max();

                const double* const F0 = residual_data->getPointer(0, depth);
                const double* const F1 = residual_data->getPointer(1, depth);
#if (NDIM == 3)
                const double* const F2 = residual_data->getPointer(2, depth);
#endif
                const int F_ghosts = (residual_data->getGhostCellWidth()).max();

                const int* const mask0 = mask_data->getPointer(0, depth);
                const int* const mask1 = mask_data->getPointer(1, depth);
#if (NDIM == 3)
                const int* const mask2 = mask_data->getPointer(2, depth);
#endif
                const int mask_ghosts = (mask_data->getGhostCellWidth()).max();

#if (NDIM == 2)
                const double* const mu = mu_data->getPointer(depth);
#endif
#if (NDIM == 3)
                const double* const mu0 = mu_data->getPointer(0, depth);
                const double* const mu1 = mu_data->getPointer(1, depth);
                const double* const mu2 = mu_data->getPointer(2, depth);
#endif
                const int mu_ghosts = (mu_data->getGhostCellWidth()).max();

                const double* C0 = nullptr;
                const double* C1 = nullptr;
#if (NDIM == 3)
                const double* C2 = nullptr;
#endif
                int C_ghosts = 0;
                if (d_poisson_spec.cIsVariable())
                {
                    C0 = C_data->getPointer(0, depth);
                    C1 = C_data->getPointer(1, depth);
#if (NDIM == 3)
                    C2 = C_data->getPointer(2, depth);
#endif
                    C_ghosts = (C_data->getGhostCellWidth()).max();
                }

                if (patch_has_dirichlet_bdry)
                {
                    if (red_black_ordering)
                    {
                        int red_or_black = isweep % 2; // "red" = 0, "black" = 1
                        VC_RB_GS_SMOOTH_MASK_FC(U0,
                                                U1,
#if (NDIM == 3)
                                                U2,
#endif
                                                U_ghosts,
                                                F0,
                                                F1,
#if (NDIM == 3)
                                                F2,
#endif
                                                F_ghosts,
                                                mask0,
                                                mask1,
#if (NDIM == 3)
                                                mask2,
#endif
                                                mask_ghosts,
                                                C0,
                                                C1,
#if (NDIM == 3)
                                                C2,
#endif
                                                C_ghosts,
#if (NDIM == 2)
                                                mu,
#endif
#if (NDIM == 3)
                                                mu0,
                                                mu1,
                                                mu2,
#endif
                                                mu_ghosts,
                                                alpha,
                                                beta,
                                                patch_box.lower(0),
                                                patch_box.upper(0),
                                                patch_box.lower(1),
                                                patch_box.upper(1),
#if (NDIM == 3)
                                                patch_box.lower(2),
                                                patch_box.upper(2),
#endif
                                                dx,
                                                C_is_varying,
                                                use_harmonic_interp,
                                                red_or_black);
                    }
                    else
                    {
                        VC_GS_SMOOTH_MASK_FC(U0,
                                             U1,
#if (NDIM == 3)
                                             U2,
#endif
                                             U_ghosts,
                                             F0,
                                             F1,
#if (NDIM == 3)
                                             F2,
#endif
                                             F_ghosts,
                                             mask0,
                                             mask1,
#if (NDIM == 3)
                                             mask2,
#endif
                                             mask_ghosts,
                                             C0,
                                             C1,
#if (NDIM == 3)
                                             C2,
#endif
                                             C_ghosts,
#if (NDIM == 2)
                                             mu,
#endif
#if (NDIM == 3)
                                             mu0,
                                             mu1,
                                             mu2,
#endif
                                             mu_ghosts,
                                             alpha,
                                             beta,
                                             patch_box.lower(0),
                                             patch_box.upper(0),
                                             patch_box.lower(1),
                                             patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2),
                                             patch_box.upper(2),
#endif
                                             dx,
                                             C_is_varying,
                                             use_harmonic_interp);
                    }
                }
                else
                {
                    if (red_black_ordering)
                    {
                        int red_or_black = isweep % 2; // "red" = 0, "black" = 1
                        VC_RB_GS_SMOOTH_FC(U0,
                                           U1,
#if (NDIM == 3)
                                           U2,
#endif
                                           U_ghosts,
                                           F0,
                                           F1,
#if (NDIM == 3)
                                           F2,
#endif
                                           F_ghosts,
                                           C0,
                                           C1,
#if (NDIM == 3)
                                           C2,
#endif
                                           C_ghosts,
#if (NDIM == 2)
                                           mu,
#endif
#if (NDIM == 3)
                                           mu0,
                                           mu1,
                                           mu2,
#endif
                                           mu_ghosts,
                                           alpha,
                                           beta,
                                           patch_box.lower(0),
                                           patch_box.upper(0),
                                           patch_box.lower(1),
                                           patch_box.upper(1),
#if (NDIM == 3)
                                           patch_box.lower(2),
                                           patch_box.upper(2),
#endif
                                           dx,
                                           C_is_varying,
                                           use_harmonic_interp,
                                           red_or_black);
                    }
                    else
                    {
                        VC_GS_SMOOTH_FC(U0,
                                        U1,
#if (NDIM == 3)
                                        U2,
#endif
                                        U_ghosts,
                                        F0,
                                        F1,
#if (NDIM == 3)
                                        F2,
#endif
                                        F_ghosts,
                                        C0,
                                        C1,
#if (NDIM == 3)
                                        C2,
#endif
                                        C_ghosts,
#if (NDIM == 2)
                                        mu,
#endif
#if (NDIM == 3)
                                        mu0,
                                        mu1,
                                        mu2,
#endif
                                        mu_ghosts,
                                        alpha,
                                        beta,
                                        patch_box.lower(0),
                                        patch_box.upper(0),
                                        patch_box.lower(1),
                                        patch_box.upper(1),
#if (NDIM == 3)
                                        patch_box.lower(2),
                                        patch_box.upper(2),
#endif
                                        dx,
                                        C_is_varying,
                                        use_harmonic_interp);
                    }
                }
            }
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleDataSynch(error_idx, level_num);
    IBTK_TIMER_STOP(t_smooth_error);
    return;
} // smoothError

void
VCSCViscousOpPointRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                         const SAMRAIVectorReal<NDIM, double>& solution,
                                                         const SAMRAIVectorReal<NDIM, double>& rhs,
                                                         int coarsest_level_num,
                                                         int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    const int res_idx = residual.getComponentDescriptorIndex(0);
    const int sol_idx = solution.getComponentDescriptorIndex(0);
    const int rhs_idx = rhs.getComponentDescriptorIndex(0);

    const Pointer<SideVariable<NDIM, double> > res_var = residual.getComponentVariable(0);
    const Pointer<SideVariable<NDIM, double> > sol_var = solution.getComponentVariable(0);

    // Fill ghost-cell values.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    Pointer<SideNoCornersFillPattern> fill_pattern = nullptr;
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

    double alpha = 1.0;
    double beta = 1.0;
    if (d_poisson_spec.cIsZero() || d_poisson_spec.cIsConstant())
    {
        beta = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();
    }

    d_level_math_ops[finest_level_num]->vc_laplace(res_idx,
                                                   res_var,
                                                   alpha,
                                                   beta,
                                                   d_poisson_spec.getDPatchDataId(),
#if (NDIM == 2)
                                                   Pointer<NodeVariable<NDIM, double> >(nullptr),
#endif
#if (NDIM == 3)
                                                   Pointer<EdgeVariable<NDIM, double> >(nullptr),
#endif
                                                   sol_idx,
                                                   sol_var,
                                                   Pointer<HierarchyGhostCellInterpolation>(nullptr),
                                                   d_solution_time,
                                                   d_D_interp_type,
                                                   d_poisson_spec.cIsVariable() ? d_poisson_spec.getCPatchDataId() :
                                                                                  -1);

    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_sc_data_ops.axpy(res_idx, -1.0, res_idx, rhs_idx, false);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
} // computeResidual

void
VCSCViscousOpPointRelaxationFACOperator::restrictResidual(const SAMRAIVectorReal<NDIM, double>& src,
                                                          SAMRAIVectorReal<NDIM, double>& dst,
                                                          int dst_ln)
{
    IBTK_TIMER_START(t_restrict_residual);

    const int src_idx = src.getComponentDescriptorIndex(0);
    const int dst_idx = dst.getComponentDescriptorIndex(0);

    // Copy src into scratch and rescale
    d_level_data_ops[dst_ln + 1]->copyData(d_scratch_idx, src_idx, /*interior_only*/ false);
    d_level_data_ops[dst_ln + 1]->scale(
        d_scratch_idx, d_A_scale[dst_ln] / d_A_scale[dst_ln + 1], d_scratch_idx, /*interior_only*/ false);

    if (src_idx != dst_idx)
    {
        d_level_data_ops[dst_ln]->copyData(dst_idx, src_idx, /*interior_only*/ false);
    }
    xeqScheduleRestriction(dst_idx, d_scratch_idx, dst_ln);

    IBTK_TIMER_STOP(t_restrict_residual);
    return;
} // restrictResidual

void
VCSCViscousOpPointRelaxationFACOperator::setDPatchDataInterpolationType(const IBTK::VCInterpType D_interp_type)
{
    d_D_interp_type = D_interp_type;
    return;
} // setDPatchDataInterpolationType

void
VCSCViscousOpPointRelaxationFACOperator::setOperatorScaling(const Array<double> A_scale)
{
    d_A_scale = A_scale;
    return;
} // setOperatorScaling

Pointer<PoissonSolver>
VCSCViscousOpPointRelaxationFACOperator::getCoarseSolver()
{
    return d_coarse_solver;
} // getCoarseSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

void
VCSCViscousOpPointRelaxationFACOperator::initializeOperatorStateSpecialized(
    const SAMRAIVectorReal<NDIM, double>& solution,
    const SAMRAIVectorReal<NDIM, double>& rhs,
    const int coarsest_reset_ln,
    const int finest_reset_ln)
{
    SCPoissonPointRelaxationFACOperator::initializeOperatorStateSpecialized(
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
