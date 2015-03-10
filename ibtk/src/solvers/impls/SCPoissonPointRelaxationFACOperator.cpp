// Filename: SCPoissonPointRelaxationFACOperator.cpp
// Created on 13 Nov 2008 by Boyce Griffith
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
#include <algorithm>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"
#include "IBTK_config.h"
#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideDataFactory.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/xfer/VariableFillPattern.h"
#include "boost/array.hpp"
#include "ibtk/CartSideDoubleCubicCoarsen.h"
#include "ibtk/CartSideDoubleQuadraticCFInterpolation.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/SCPoissonPointRelaxationFACOperator.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define GS_SMOOTH_FC IBTK_FC_FUNC(gssmooth2d, GSSMOOTH2D)
#define GS_SMOOTH_MASK_FC IBTK_FC_FUNC(gssmoothmask2d, GSSMOOTHMASK2D)
#define RB_GS_SMOOTH_FC IBTK_FC_FUNC(rbgssmooth2d, RBGSSMOOTH2D)
#define RB_GS_SMOOTH_MASK_FC IBTK_FC_FUNC(rbgssmoothmask2d, RBGSSMOOTHMASK2D)
#endif
#if (NDIM == 3)
#define GS_SMOOTH_FC IBTK_FC_FUNC(gssmooth3d, GSSMOOTH3D)
#define GS_SMOOTH_MASK_FC IBTK_FC_FUNC(gssmoothmask3d, GSSMOOTHMASK3D)
#define RB_GS_SMOOTH_FC IBTK_FC_FUNC(rbgssmooth3d, RBGSSMOOTH3D)
#define RB_GS_SMOOTH_MASK_FC IBTK_FC_FUNC(rbgssmoothmask3d, RBGSSMOOTHMASK3D)
#endif

// Function interfaces
extern "C" {
void GS_SMOOTH_FC(double* U,
                  const int& U_gcw,
                  const double& alpha,
                  const double& beta,
                  const double* F,
                  const int& F_gcw,
                  const int& ilower0,
                  const int& iupper0,
                  const int& ilower1,
                  const int& iupper1,
#if (NDIM == 3)
                  const int& ilower2,
                  const int& iupper2,
#endif
                  const double* dx);

void GS_SMOOTH_MASK_FC(double* U,
                       const int& U_gcw,
                       const double& alpha,
                       const double& beta,
                       const double* F,
                       const int& F_gcw,
                       const int* mask,
                       const int& mask_gcw,
                       const int& ilower0,
                       const int& iupper0,
                       const int& ilower1,
                       const int& iupper1,
#if (NDIM == 3)
                       const int& ilower2,
                       const int& iupper2,
#endif
                       const double* dx);

void RB_GS_SMOOTH_FC(double* U,
                     const int& U_gcw,
                     const double& alpha,
                     const double& beta,
                     const double* F,
                     const int& F_gcw,
                     const int& ilower0,
                     const int& iupper0,
                     const int& ilower1,
                     const int& iupper1,
#if (NDIM == 3)
                     const int& ilower2,
                     const int& iupper2,
#endif
                     const double* dx,
                     const int& red_or_black);

void RB_GS_SMOOTH_MASK_FC(double* U,
                          const int& U_gcw,
                          const double& alpha,
                          const double& beta,
                          const double* F,
                          const int& F_gcw,
                          const int* mask,
                          const int& mask_gcw,
                          const int& ilower0,
                          const int& iupper0,
                          const int& ilower1,
                          const int& iupper1,
#if (NDIM == 3)
                          const int& ilower2,
                          const int& iupper2,
#endif
                          const double* dx,
                          const int& red_or_black);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static boost::shared_ptr<Timer> t_smooth_error;
static boost::shared_ptr<Timer> t_solve_coarsest_level;
static boost::shared_ptr<Timer> t_compute_residual;

// Default data depth.
static const int DEFAULT_DATA_DEPTH = 1;

// Number of ghosts cells used for each variable quantity.
static const int SIDEG = 1;

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

inline SmootherType get_smoother_type(const std::string& smoother_type_string)
{
    if (smoother_type_string == "PATCH_GAUSS_SEIDEL") return PATCH_GAUSS_SEIDEL;
    if (smoother_type_string == "PROCESSOR_GAUSS_SEIDEL") return PROCESSOR_GAUSS_SEIDEL;
    if (smoother_type_string == "RED_BLACK_GAUSS_SEIDEL")
        return RED_BLACK_GAUSS_SEIDEL;
    else
        return UNKNOWN;
}

inline bool use_red_black_ordering(SmootherType smoother_type)
{
    if (smoother_type == RED_BLACK_GAUSS_SEIDEL)
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline bool do_local_data_update(SmootherType smoother_type)
{
    if (smoother_type == PROCESSOR_GAUSS_SEIDEL || smoother_type == RED_BLACK_GAUSS_SEIDEL)
    {
        return true;
    }
    else
    {
        return false;
    }
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SCPoissonPointRelaxationFACOperator::SCPoissonPointRelaxationFACOperator(const std::string& object_name,
                                                                         const boost::shared_ptr<Database> input_db,
                                                                         const std::string& default_options_prefix)
    : PoissonFACPreconditionerStrategy(
          object_name,
          boost::make_shared<SideVariable<double> >(DIM, object_name + "::side_scratch", DEFAULT_DATA_DEPTH),
          SIDEG,
          input_db,
          default_options_prefix),
      d_coarse_solver(NULL), d_coarse_solver_db(), d_patch_bc_box_overlap(), d_patch_neighbor_overlap()
{
    // Set some default values.
    d_smoother_type = "PATCH_GAUSS_SEIDEL";
    d_prolongation_method = "CONSTANT_REFINE";
    d_restriction_method = "CONSERVATIVE_COARSEN";
    d_coarse_solver_type = SCPoissonSolverManager::HYPRE_LEVEL_SOLVER;
    d_coarse_solver_rel_residual_tol = 1.0e-5;
    d_coarse_solver_abs_residual_tol = 1.0e-50;
    d_coarse_solver_max_iterations = 1;
    d_coarse_solver_db = boost::make_shared<MemoryDatabase>(object_name + "::coarse_solver_db");
    d_coarse_solver_db->putString("solver_type", "Split");
    d_coarse_solver_db->putString("split_solver_type", "PFMG");

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("smoother_type")) d_smoother_type = input_db->getString("smoother_type");
        if (input_db->keyExists("prolongation_method"))
            d_prolongation_method = input_db->getString("prolongation_method");
        if (input_db->keyExists("restriction_method")) d_restriction_method = input_db->getString("restriction_method");
        if (input_db->keyExists("coarse_solver_type")) d_coarse_solver_type = input_db->getString("coarse_solver_type");
        if (input_db->keyExists("coarse_solver_rel_residual_tol"))
            d_coarse_solver_rel_residual_tol = input_db->getDouble("coarse_solver_rel_residual_tol");
        if (input_db->keyExists("coarse_solver_abs_residual_tol"))
            d_coarse_solver_abs_residual_tol = input_db->getDouble("coarse_solver_abs_residual_tol");
        if (input_db->keyExists("coarse_solver_max_iterations"))
            d_coarse_solver_max_iterations = input_db->getInteger("coarse_solver_max_iterations");
        if (input_db->isDatabase("coarse_solver_db"))
        {
            d_coarse_solver_db = input_db->getDatabase("coarse_solver_db");
        }
        if (input_db->isDatabase("bottom_solver"))
        {
            tbox::pout << "WARNING: ``bottom_solver'' input entry is no longer used by class "
                          "SCPoissonPointRelaxationFACOperator.\n"
                       << "         use ``coarse_solver_db'' input entry instead.\n";
        }
        if (input_db->isDatabase("hypre_solver"))
        {
            tbox::pout << "WARNING: ``hypre_solver'' input entry is no longer used by class "
                          "SCPoissonPointRelaxationFACOperator.\n"
                       << "         use ``coarse_solver_db'' input entry instead.\n";
        }
        if (input_db->isDatabase("petsc_solver"))
        {
            tbox::pout << "WARNING: ``petsc_solver'' input entry is no longer used by class "
                          "SCPoissonPointRelaxationFACOperator.\n"
                       << "         use ``coarse_solver_db'' input entry instead.\n";
        }
    }

    // Configure the coarse level solver.
    setCoarseSolverType(d_coarse_solver_type);

    // Construct a variable to store any needed masking data.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    auto mask_var =
        boost::make_shared<SideVariable<int> >(DIM, object_name + "::mask");
    if (var_db->checkVariableExists(mask_var->getName()))
    {
        mask_var = var_db->getVariable(mask_var->getName());
        d_mask_idx = var_db->mapVariableAndContextToIndex(mask_var, d_context);
        var_db->removePatchDataIndex(d_mask_idx);
    }
    IntVector no_ghosts = IntVector::getZero(DIM);
    d_mask_idx = var_db->registerVariableAndContext(mask_var, d_context, no_ghosts);

    // Setup Timers.
    IBTK_DO_ONCE(t_smooth_error =
                     TimerManager::getManager()->getTimer("IBTK::SCPoissonPointRelaxationFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "IBTK::SCPoissonPointRelaxationFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::SCPoissonPointRelaxationFACOperator::computeResidual()"););
    return;
}

SCPoissonPointRelaxationFACOperator::~SCPoissonPointRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}

void SCPoissonPointRelaxationFACOperator::setSmootherType(const std::string& smoother_type)
{
    d_smoother_type = smoother_type;
    return;
}

void SCPoissonPointRelaxationFACOperator::setCoarseSolverType(const std::string& coarse_solver_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setCoarseSolverType():\n"
                                 << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (d_coarse_solver_type != coarse_solver_type) d_coarse_solver.reset();
    d_coarse_solver_type = coarse_solver_type;
    if (get_smoother_type(d_coarse_solver_type) == UNKNOWN && !d_coarse_solver)
    {
        d_coarse_solver = SCPoissonSolverManager::getManager()->allocateSolver(d_coarse_solver_type,
                                                                               d_object_name + "::coarse_solver",
                                                                               d_coarse_solver_db,
                                                                               d_coarse_solver_default_options_prefix);
    }
    return;
}

void SCPoissonPointRelaxationFACOperator::smoothError(SAMRAIVectorReal<double>& error,
                                                      const SAMRAIVectorReal<double>& residual,
                                                      int level_num,
                                                      int num_sweeps,
                                                      bool /*performing_pre_sweeps*/,
                                                      bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBTK_TIMER_START(t_smooth_error);

    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    const int mpi_rank = comm.getRank();

    auto level =d_hierarchy->getPatchLevel(level_num);
    const int error_idx = error.getComponentDescriptorIndex(0);
    const int scratch_idx = d_scratch_idx;

    // Determine the smoother type.
    const std::string& smoother_type_string = (level_num == d_coarsest_ln ? d_coarse_solver_type : d_smoother_type);
    const SmootherType smoother_type = get_smoother_type(smoother_type_string);
    TBOX_ASSERT(smoother_type != UNKNOWN);
    const bool red_black_ordering = use_red_black_ordering(smoother_type);
    const bool update_local_data = do_local_data_update(smoother_type);

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
        {
            auto patch =*p;
            boost::shared_ptr<SideData<double> > error_data = error.getComponentPatchData(0, *patch);
            boost::shared_ptr<SideData<double> > scratch_data = patch->getPatchData(scratch_idx);
            const Box& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                scratch_data->getArrayData(axis).copy(error_data->getArrayData(axis),
                                                      d_patch_bc_box_overlap[level_num][patch_counter][axis],
                                                      IntVector::getZero(DIM));
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
                for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
                {
                    auto patch =*p;
                    boost::shared_ptr<SideData<double> > error_data = error.getComponentPatchData(0, *patch);
                    boost::shared_ptr<SideData<double> > scratch_data = patch->getPatchData(scratch_idx);
                    const Box& ghost_box = error_data->getGhostBox();
                    TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
                    TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        error_data->getArrayData(axis).copy(scratch_data->getArrayData(axis),
                                                            d_patch_bc_box_overlap[level_num][patch_counter][axis],
                                                            IntVector::getZero(DIM));
                    }
                }

                // Fill the non-coarse-fine interface ghost cell values.
                xeqScheduleGhostFillNoCoarse(error_idx, level_num);
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_cf_bdry_op->setPatchDataIndex(error_idx);
            const IntVector& ratio = level->getRatioToCoarserLevel();
            for (auto p = level->begin(); p != level->end(); ++p)
            {
                auto patch =*p;
                const IntVector& ghost_width_to_fill = d_gcw;
                d_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }
        else if (isweep > 0)
        {
            xeqScheduleGhostFillNoCoarse(error_idx, level_num);
        }

        // Smooth the error on the patches.
        int patch_counter = 0;
        for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
        {
            auto patch =*p;
            boost::shared_ptr<SideData<double> > error_data = error.getComponentPatchData(0, *patch);
            boost::shared_ptr<SideData<double> > residual_data = residual.getComponentPatchData(0, *patch);
            const Box& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == residual_data->getGhostBox());
            TBOX_ASSERT(error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(error_data->getDepth() == residual_data->getDepth());
            boost::shared_ptr<SideData<int> > mask_data = patch->getPatchData(d_mask_idx);
            const Box& patch_box = patch->getBox();
            const auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
            const double* const dx = pgeom->getDx();

            // Copy updated values from neighboring local patches.
            if (update_local_data)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const std::map<int, Box> neighbor_overlap =
                        d_patch_neighbor_overlap[level_num][patch_counter][axis];
                    for (auto cit = neighbor_overlap.begin();
                         cit != neighbor_overlap.end();
                         ++cit)
                    {
                        const GlobalId src_patch_id(LocalId(cit->first), mpi_rank);
                        const Box& overlap = cit->second;
                        auto src_patch = level->getPatch(src_patch_id);
                        boost::shared_ptr<SideData<double> > src_error_data =
                            error.getComponentPatchData(0, *src_patch);
                        error_data->getArrayData(axis)
                            .copy(src_error_data->getArrayData(axis), overlap, IntVector::getZero(DIM));
                    }
                }
            }

            // Enforce any Dirichlet boundary conditions.
            const bool patch_has_dirichlet_bdry = d_bc_helper->patchTouchesDirichletBoundary(patch);
            if (patch_has_dirichlet_bdry)
            {
                d_bc_helper->copyDataAtDirichletBoundaries(error_data, residual_data, patch);
            }

            // Smooth the error using Gauss-Seidel.
            const double& alpha = d_poisson_spec.getDConstant();
            const double beta = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box side_patch_box = SideGeometry::toSideBox(patch_box, axis);
                for (int depth = 0; depth < error_data->getDepth(); ++depth)
                {
                    double* const U = error_data->getPointer(axis, depth);
                    const int U_ghosts = (error_data->getGhostCellWidth()).max();
                    const double* const F = residual_data->getPointer(axis, depth);
                    const int F_ghosts = (residual_data->getGhostCellWidth()).max();
                    const int* const mask = mask_data->getPointer(axis, depth);
                    const int mask_ghosts = (mask_data->getGhostCellWidth()).max();
                    if (patch_has_dirichlet_bdry && d_bc_helper->patchTouchesDirichletBoundaryAxis(patch, axis))
                    {
                        if (red_black_ordering)
                        {
                            int red_or_black = isweep % 2; // "red" = 0, "black" = 1
                            RB_GS_SMOOTH_MASK_FC(U,
                                                 U_ghosts,
                                                 alpha,
                                                 beta,
                                                 F,
                                                 F_ghosts,
                                                 mask,
                                                 mask_ghosts,
                                                 side_patch_box.lower(0),
                                                 side_patch_box.upper(0),
                                                 side_patch_box.lower(1),
                                                 side_patch_box.upper(1),
#if (NDIM == 3)
                                                 side_patch_box.lower(2),
                                                 side_patch_box.upper(2),
#endif
                                                 dx,
                                                 red_or_black);
                        }
                        else
                        {
                            GS_SMOOTH_MASK_FC(U,
                                              U_ghosts,
                                              alpha,
                                              beta,
                                              F,
                                              F_ghosts,
                                              mask,
                                              mask_ghosts,
                                              side_patch_box.lower(0),
                                              side_patch_box.upper(0),
                                              side_patch_box.lower(1),
                                              side_patch_box.upper(1),
#if (NDIM == 3)
                                              side_patch_box.lower(2),
                                              side_patch_box.upper(2),
#endif
                                              dx);
                        }
                    }
                    else
                    {
                        if (red_black_ordering)
                        {
                            int red_or_black = isweep % 2; // "red" = 0, "black" = 1
                            RB_GS_SMOOTH_FC(U,
                                            U_ghosts,
                                            alpha,
                                            beta,
                                            F,
                                            F_ghosts,
                                            side_patch_box.lower(0),
                                            side_patch_box.upper(0),
                                            side_patch_box.lower(1),
                                            side_patch_box.upper(1),
#if (NDIM == 3)
                                            side_patch_box.lower(2),
                                            side_patch_box.upper(2),
#endif
                                            dx,
                                            red_or_black);
                        }
                        else
                        {
                            GS_SMOOTH_FC(U,
                                         U_ghosts,
                                         alpha,
                                         beta,
                                         F,
                                         F_ghosts,
                                         side_patch_box.lower(0),
                                         side_patch_box.upper(0),
                                         side_patch_box.lower(1),
                                         side_patch_box.upper(1),
#if (NDIM == 3)
                                         side_patch_box.lower(2),
                                         side_patch_box.upper(2),
#endif
                                         dx);
                        }
                    }
                }
            }
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleDataSynch(error_idx, level_num);
    IBTK_TIMER_STOP(t_smooth_error);
    return;
}

bool SCPoissonPointRelaxationFACOperator::solveCoarsestLevel(SAMRAIVectorReal<double>& error,
                                                             const SAMRAIVectorReal<double>& residual,
                                                             int coarsest_ln)
{
    IBTK_TIMER_START(t_solve_coarsest_level);

    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
    if (d_coarse_solver)
    {
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setMaxIterations(d_coarse_solver_max_iterations);
        d_coarse_solver->setAbsoluteTolerance(d_coarse_solver_abs_residual_tol);
        d_coarse_solver->setRelativeTolerance(d_coarse_solver_rel_residual_tol);
        auto p_coarse_solver = dynamic_cast<LinearSolver*>(d_coarse_solver.getPointer());
        if (p_coarse_solver) p_coarse_solver->setInitialGuessNonzero(true);
        d_coarse_solver->solveSystem(*getLevelSAMRAIVectorReal(error, d_coarsest_ln),
                                     *getLevelSAMRAIVectorReal(residual, d_coarsest_ln));
        xeqScheduleDataSynch(error.getComponentDescriptorIndex(0), coarsest_ln);
    }
    else
    {
        TBOX_ASSERT(get_smoother_type(d_coarse_solver_type) != UNKNOWN);
        smoothError(error, residual, coarsest_ln, d_coarse_solver_max_iterations, false, false);
    }
    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
}

void SCPoissonPointRelaxationFACOperator::computeResidual(SAMRAIVectorReal<double>& residual,
                                                          const SAMRAIVectorReal<double>& solution,
                                                          const SAMRAIVectorReal<double>& rhs,
                                                          int coarsest_level_num,
                                                          int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    const int res_idx = residual.getComponentDescriptorIndex(0);
    const int sol_idx = solution.getComponentDescriptorIndex(0);
    const int rhs_idx = rhs.getComponentDescriptorIndex(0);

    auto res_var = BOOST_CAST<SideVariable<double> >(residual.getComponentVariable(0));
    auto sol_var = BOOST_CAST<SideVariable<double> >(solution.getComponentVariable(0));
    auto rhs_var = BOOST_CAST<SideVariable<double> >(rhs.getComponentVariable(0));

    // Fill ghost-cell values.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    auto fill_pattern = boost::make_shared<SideNoCornersFillPattern>(SIDEG, false, false, true);
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
        d_level_bdry_fill_ops[finest_level_num] = boost::make_shared<HierarchyGhostCellInterpolation>();
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
        std::ostringstream stream;
        stream << d_object_name << "::hier_math_ops_" << finest_level_num;
        d_level_math_ops[finest_level_num] =
            boost::make_shared<HierarchyMathOps>(stream.str(), d_hierarchy, coarsest_level_num, finest_level_num);
    }
    boost::shared_ptr<HierarchyGhostCellInterpolation> no_fill_op;
    d_level_math_ops[finest_level_num]->laplace(
        res_idx, res_var, d_poisson_spec, sol_idx, sol_var, no_fill_op, d_solution_time);
    HierarchySideDataOpsReal<double> hier_sc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    hier_sc_data_ops.axpy(res_idx, -1.0, res_idx, rhs_idx, false);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void SCPoissonPointRelaxationFACOperator::initializeOperatorStateSpecialized(const SAMRAIVectorReal<double>& solution,
                                                                             const SAMRAIVectorReal<double>& rhs,
                                                                             const int coarsest_reset_ln,
                                                                             const int finest_reset_ln)
{
    // Setup solution and rhs vectors.
    auto solution_var = BOOST_CAST<SideVariable<double> >(solution.getComponentVariable(0));
    auto rhs_var = BOOST_CAST<SideVariable<double> >(rhs.getComponentVariable(0));
    TBOX_ASSERT(solution_var);
    TBOX_ASSERT(rhs_var);
    if (solution_var->getDepth() != rhs_var->getDepth())
    {
        TBOX_ERROR("SCPoissonPointRelaxationFACOperator::initializeOperatorState()\n"
                   << "  solution and rhs vectors must have the same data depths\n"
                   << "  solution data depth = " << solution_var->getDepth() << "\n"
                   << "  rhs      data depth = " << rhs_var->getDepth() << std::endl);
    }

    VariableDatabase* var_db = VariableDatabase::getDatabase();
    boost::shared_ptr<Variable> scratch_var;
    var_db->mapIndexToVariable(d_scratch_idx, scratch_var);
    auto scratch_sc_var = BOOST_CAST<SideVariable<double> >(scratch_var);
    TBOX_ASSERT(scratch_sc_var);
    const int depth = solution_var->getDepth();
    if (scratch_sc_var->getDepth() != depth)
    {
        var_db->removePatchDataIndex(d_scratch_idx);
        const IntVector ghosts = d_gcw;
        scratch_var = boost::make_shared<SideVariable<double> >(scratch_var->getDim(), scratch_var->getName(), depth);
        d_scratch_idx = var_db->registerVariableAndContext(scratch_var, d_context, ghosts);
    }

    // Setup cached BC data.
    d_bc_helper = boost::make_shared<StaggeredPhysicalBoundaryHelper>();
    d_bc_helper->cacheBcCoefData(d_bc_coefs, d_solution_time, d_hierarchy);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        auto level =d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_mask_idx)) level->allocatePatchData(d_mask_idx);
    }
    d_bc_helper->setupMaskingFunction(d_mask_idx);

    // Initialize the coarse level solvers when needed.
    if (coarsest_reset_ln == d_coarsest_ln && d_coarse_solver)
    {
        // Note that since the coarse level solver is solving for the error, it
        // must always employ homogeneous boundary conditions.
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setPoissonSpecifications(d_poisson_spec);
        d_coarse_solver->setPhysicalBcCoefs(d_bc_coefs);
        d_coarse_solver->setHomogeneousBc(true);
        d_coarse_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, d_coarsest_ln),
                                               *getLevelSAMRAIVectorReal(*d_rhs, d_coarsest_ln));
    }

    // Setup specialized transfer operators.
    auto geometry = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());
    IBTK_DO_ONCE(geometry->addSpatialCoarsenOperator(boost::make_shared<CartSideDoubleCubicCoarsen>()););

    // Setup coarse-fine interface and physical boundary operators.
    d_cf_bdry_op = boost::make_shared<CartSideDoubleQuadraticCFInterpolation>();
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);
    d_bc_op = boost::make_shared<CartSideRobinPhysBdryOp>(d_scratch_idx, d_bc_coefs, false);

    // Setup fill pattern spec objects.
    if (d_poisson_spec.dIsConstant())
    {
        d_op_stencil_fill_pattern = boost::make_shared<SideNoCornersFillPattern>(SIDEG, true, false, false);
    }
    else
    {
        d_op_stencil_fill_pattern = NULL;
    }
    d_synch_fill_pattern = boost::make_shared<SideSynchCopyFillPattern>();

    // Get overlap information for setting patch boundary conditions.
    d_patch_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        auto level =d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (auto p = level->begin(); p != level->end(); ++p, ++patch_counter)
        {
            auto patch =*p;
            const Box& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box side_box = SideGeometry::toSideBox(patch_box, axis);
                const Box side_ghost_box = Box::grow(side_box, IntVector::getOne(DIM));
                d_patch_bc_box_overlap[ln][patch_counter][axis] = BoxContainer(side_ghost_box);
                d_patch_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    // Get overlap information for re-setting patch boundary conditions during
    // smoothing.
    d_patch_neighbor_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        auto level =d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_neighbor_overlap[ln].resize(num_local_patches);
        int patch_counter1 = 0;
        for (auto p1(level); p1; p1++, ++patch_counter1)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                d_patch_neighbor_overlap[ln][patch_counter1][axis].clear();
            }
            auto dst_patch = p1();
            const Box& dst_patch_box = dst_patch->getBox();
            const Box& dst_ghost_box = Box::grow(dst_patch_box, IntVector::getOne(DIM));
            int patch_counter2 = 0;
            for (auto p2(level); patch_counter2 < patch_counter1; p2++, ++patch_counter2)
            {
                auto src_patch = p2();
                const Box& src_patch_box = src_patch->getBox();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const Box overlap =
                        SideGeometry::toSideBox(dst_ghost_box, axis) * SideGeometry::toSideBox(src_patch_box, axis);
                    if (!overlap.empty())
                    {
                        d_patch_neighbor_overlap[ln][patch_counter1][axis].insert(
                            std::make_pair(src_patch->getLocalId().getValue(), overlap));
                    }
                }
            }
        }
    }
    return;
}

void SCPoissonPointRelaxationFACOperator::deallocateOperatorStateSpecialized(const int /*coarsest_reset_ln*/,
                                                                             const int /*finest_reset_ln*/)
{
    if (!d_is_initialized) return;

    if (!d_in_initialize_operator_state)
    {
        d_patch_bc_box_overlap.clear();
        d_patch_neighbor_overlap.clear();
        if (d_coarse_solver) d_coarse_solver->deallocateSolverState();
    }
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
