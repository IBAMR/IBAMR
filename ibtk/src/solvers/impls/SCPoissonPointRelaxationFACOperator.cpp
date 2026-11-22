// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleQuadraticCFInterpolation.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/CoarseFineBoundaryRefinePatchStrategy.h>
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/PoissonFACPreconditionerStrategy.h>
#include <ibtk/PoissonSolver.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/SCPoissonPointRelaxationFACOperator.h>
#include <ibtk/SCPoissonSolverManager.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/StaggeredPhysicalBoundaryHelper.h>
#include <ibtk/VCSCViscousPETScLevelSolver.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/solver_utilities.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/MemoryDatabase.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxList.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarseFineBoundary.h>
#include <CoarsenOperator.h>
#include <EdgeData.h>
#include <MultiblockDataTranslator.h>
#include <NodeData.h>
#include <Patch.h>
#include <PatchDescriptor.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <ProcessorMapping.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define SMOOTH_GS_CONST_DC_FC IBTK_FC_FUNC_(smooth_gs_const_dc_2d, SMOOTH_GS_CONST_DC_2D)
#define SMOOTH_GS_RB_CONST_DC_FC IBTK_FC_FUNC_(smooth_gs_rb_const_dc_2d, SMOOTH_GS_RB_CONST_DC_2D)
#define GS_SMOOTH_MASK_FC IBTK_FC_FUNC(gssmoothmask2d, GSSMOOTHMASK2D)
#define RB_GS_SMOOTH_MASK_FC IBTK_FC_FUNC(rbgssmoothmask2d, RBGSSMOOTHMASK2D)
#endif
#if (NDIM == 3)
#define SMOOTH_GS_CONST_DC_FC IBTK_FC_FUNC_(smooth_gs_const_dc_3d, SMOOTH_GS_CONST_DC_3D)
#define SMOOTH_GS_RB_CONST_DC_FC IBTK_FC_FUNC_(smooth_gs_rb_const_dc_3d, SMOOTH_GS_RB_CONST_DC_3D)
#define GS_SMOOTH_MASK_FC IBTK_FC_FUNC(gssmoothmask3d, GSSMOOTHMASK3D)
#define RB_GS_SMOOTH_MASK_FC IBTK_FC_FUNC(rbgssmoothmask3d, RBGSSMOOTHMASK3D)
#endif

// Function interfaces
extern "C"
{
    void SMOOTH_GS_CONST_DC_FC(double* U,
                               const int& U_gcw,
                               const double& D,
                               const double& C,
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

    void SMOOTH_GS_RB_CONST_DC_FC(double* U,
                                  const int& U_gcw,
                                  const double& D,
                                  const double& C,
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
static Timer* t_smooth_error;
static Timer* t_solve_coarsest_level;
static Timer* t_compute_residual;

// Default data depth.
static const int DEFAULT_DATA_DEPTH = 1;

// Number of ghosts cells used for each variable quantity.
static const int SIDEG = 1;

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

SCPoissonPointRelaxationFACOperator::SCPoissonPointRelaxationFACOperator(const std::string& object_name,
                                                                         const Pointer<Database> input_db,
                                                                         std::string default_options_prefix)
    : PoissonFACPreconditionerStrategy(
          object_name,
          new SideVariable<NDIM, double>(object_name + "::side_scratch", DEFAULT_DATA_DEPTH),
          SIDEG,
          input_db,
          std::move(default_options_prefix))
{
    // Set some default values.
    d_smoother_type = "PATCH_GAUSS_SEIDEL";
    d_prolongation_method = "CONSTANT_REFINE";
    d_restriction_method = "CONSERVATIVE_COARSEN";
    d_coarse_solver_type = SCPoissonSolverManager::HYPRE_LEVEL_SOLVER;
    d_coarse_solver_rel_residual_tol = 1.0e-5;
    d_coarse_solver_abs_residual_tol = 1.0e-50;
    d_coarse_solver_max_iterations = 1;
    d_coarse_solver_db = new MemoryDatabase(object_name + "::coarse_solver_db");
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
        if (input_db->keyExists("coarse_solver_prefix"))
            d_coarse_solver_default_options_prefix = input_db->getString("coarse_solver_prefix");
        if (input_db->isDatabase("coarse_solver_db"))
        {
            d_coarse_solver_db = input_db->getDatabase("coarse_solver_db");
        }
        if (input_db->isDatabase("bottom_solver"))
        {
            TBOX_WARNING(d_object_name << "::SCPoissonPointRelaxationFACOperator():\n"
                                       << "  ``bottom_solver'' input entry is no longer used\n"
                                       << "  use ``coarse_solver_db'' input entry instead" << std::endl);
        }
        if (input_db->isDatabase("hypre_solver"))
        {
            TBOX_WARNING(d_object_name << "::SCPoissonPointRelaxationFACOperator():\n"
                                       << "  ``hypre_solver'' input entry is no longer used\n"
                                       << "  use ``coarse_solver_db'' input entry instead" << std::endl);
        }
        if (input_db->isDatabase("petsc_solver"))
        {
            TBOX_WARNING(d_object_name << "::SCPoissonPointRelaxationFACOperator():\n"
                                       << "  ``petsc_solver'' input entry is no longer used\n"
                                       << "  use ``coarse_solver_db'' input entry instead" << std::endl);
        }
        if (input_db->isBool("reuse_unchanged_coarse_solver"))
        {
            d_reuse_unchanged_coarse_solver = input_db->getBool("reuse_unchanged_coarse_solver");
        }
        if (input_db->keyExists("coarse_solver_reuse_rel_tol"))
        {
            d_coarse_solver_reuse_rel_tol = input_db->getDouble("coarse_solver_reuse_rel_tol");
            if (!std::isfinite(d_coarse_solver_reuse_rel_tol) || d_coarse_solver_reuse_rel_tol < 0.0)
            {
                TBOX_ERROR(d_object_name << "::SCPoissonPointRelaxationFACOperator():\n"
                                         << "  coarse_solver_reuse_rel_tol must be finite and >= 0.0, got "
                                         << d_coarse_solver_reuse_rel_tol << "\n");
            }
        }
    }

    // Check the smoother type and configure the coarse level solver.
    if (get_smoother_type(d_smoother_type) == UNKNOWN)
    {
        TBOX_ERROR(d_object_name << "::SCPoissonPointRelaxationFACOperator():\n"
                                 << "  unsupported smoother type: " << d_smoother_type << "\n"
                                 << "  supported types: PATCH_GAUSS_SEIDEL, PROCESSOR_GAUSS_SEIDEL, "
                                    "RED_BLACK_GAUSS_SEIDEL"
                                 << std::endl);
    }
    if (get_smoother_type(d_coarse_solver_type) == UNKNOWN)
    {
        d_coarse_solver = SCPoissonSolverManager::getManager()->allocateSolver(d_coarse_solver_type,
                                                                               d_object_name + "::coarse_solver",
                                                                               d_coarse_solver_db,
                                                                               d_coarse_solver_default_options_prefix);
        set_fixed_iteration_ksp_defaults(d_coarse_solver.getPointer(), d_coarse_solver_db);
    }

    // Construct a variable to store any needed masking data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, int>> mask_var = new SideVariable<NDIM, int>(object_name + "::mask");
    if (var_db->checkVariableExists(mask_var->getName()))
    {
        mask_var = var_db->getVariable(mask_var->getName());
        d_mask_idx = var_db->mapVariableAndContextToIndex(mask_var, d_context);
        var_db->removePatchDataIndex(d_mask_idx);
    }
    IntVector<NDIM> no_ghosts = 0;
    d_mask_idx = var_db->registerVariableAndContext(mask_var, d_context, no_ghosts);

    // Setup Timers.
    IBTK_DO_ONCE(t_smooth_error =
                     TimerManager::getManager()->getTimer("IBTK::SCPoissonPointRelaxationFACOperator::smoothError()");
                 t_solve_coarsest_level = TimerManager::getManager()->getTimer(
                     "IBTK::SCPoissonPointRelaxationFACOperator::solveCoarsestLevel()");
                 t_compute_residual = TimerManager::getManager()->getTimer(
                     "IBTK::SCPoissonPointRelaxationFACOperator::computeResidual()"););
    return;
} // SCPoissonPointRelaxationFACOperator

namespace
{
inline bool
scalars_match_rel(const double a, const double b, const double rel_tol)
{
    if (!std::isfinite(a) || !std::isfinite(b))
    {
        return false;
    }
    const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= rel_tol * scale;
}

bool
bc_coefs_match_exact(const std::vector<double>& cur, const std::vector<double>& cached, double& max_rel_diff)
{
    max_rel_diff = 0.0;
    if (cur.size() != cached.size())
    {
        max_rel_diff = std::numeric_limits<double>::infinity();
        return false;
    }
    bool all_match = true;
    for (size_t i = 0; i < cur.size(); ++i)
    {
        const double a = cur[i];
        const double b = cached[i];
        if (!std::isfinite(a) || !std::isfinite(b))
        {
            max_rel_diff = std::numeric_limits<double>::infinity();
            all_match = false;
            continue;
        }
        if (a != b)
        {
            all_match = false;
            const double denom = std::max(std::abs(a), std::abs(b));
            const double rel = (denom > 0.0) ? (std::abs(a - b) / denom) : 1.0;
            if (rel > max_rel_diff)
            {
                max_rel_diff = rel;
            }
        }
    }
    return all_match;
}

template <typename EnumType>
EnumType
classify_coarse_solver(const PoissonSolver* solver,
                       const EnumType unsupported_val,
                       const EnumType vcsc_petsc_val,
                       IBTK::VCInterpType& mu_interp_out)
{
    mu_interp_out = IBTK::VC_HARMONIC_INTERP;
    if (!solver)
    {
        return unsupported_val;
    }
    const std::type_info& ti = typeid(*solver);
    if (ti == typeid(VCSCViscousPETScLevelSolver))
    {
        mu_interp_out = static_cast<const VCSCViscousPETScLevelSolver*>(solver)->getViscosityInterpolationType();
        return vcsc_petsc_val;
    }
    return unsupported_val;
}

bool
coeff_vectors_match_rel(const std::vector<double>& cur,
                        const std::vector<double>& cached,
                        const double rel_tol,
                        double& max_rel_diff)
{
    max_rel_diff = 0.0;
    if (cur.size() != cached.size())
    {
        max_rel_diff = std::numeric_limits<double>::infinity();
        return false;
    }
    bool all_match = true;
    for (size_t i = 0; i < cur.size(); ++i)
    {
        const double a = cur[i];
        const double b = cached[i];
        if (!std::isfinite(a) || !std::isfinite(b))
        {
            max_rel_diff = std::numeric_limits<double>::infinity();
            all_match = false;
            continue;
        }
        const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
        const double rel = std::abs(a - b) / scale;
        if (rel > max_rel_diff)
        {
            max_rel_diff = rel;
        }
        if (rel > rel_tol)
        {
            all_match = false;
        }
    }
    return all_match;
}

inline void
append_array_data(const SAMRAI::pdat::ArrayData<NDIM, double>& arr, std::vector<double>& out)
{
    const int n = arr.getOffset() * arr.getDepth();
    if (n > 0)
    {
        const double* ptr = arr.getPointer(0);
        out.insert(out.end(), ptr, ptr + n);
    }
}

bool
append_patch_data_values(const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>& patch,
                         const int data_id,
                         std::vector<double>& out)
{
    if (!patch || data_id < 0 || !patch->checkAllocated(data_id))
    {
        return false;
    }
    const SAMRAI::hier::PatchData<NDIM>* raw = patch->getPatchData(data_id).getPointer();
    if (const auto* sd = dynamic_cast<const SAMRAI::pdat::SideData<NDIM, double>*>(raw))
    {
        for (int d = 0; d < NDIM; ++d)
        {
            append_array_data(sd->getArrayData(d), out);
        }
        return true;
    }
    if (const auto* ed = dynamic_cast<const SAMRAI::pdat::EdgeData<NDIM, double>*>(raw))
    {
        for (int d = 0; d < NDIM; ++d)
        {
            append_array_data(ed->getArrayData(d), out);
        }
        return true;
    }
    if (const auto* nd = dynamic_cast<const SAMRAI::pdat::NodeData<NDIM, double>*>(raw))
    {
        append_array_data(nd->getArrayData(), out);
        return true;
    }
    return false;
}

bool
snapshot_l0_coeff_values(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>& l0_level,
                         const int data_id,
                         std::vector<double>& out)
{
    out.clear();
    if (!l0_level)
    {
        return false;
    }
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(l0_level); p; p++)
    {
        if (!append_patch_data_values(l0_level->getPatch(p()), data_id, out))
        {
            return false;
        }
    }
    return true;
}

void
capture_l0_boxes_and_procs(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>& l0_level,
                           std::vector<SAMRAI::hier::Box<NDIM>>& boxes_out,
                           std::vector<int>& procs_out)
{
    const SAMRAI::hier::BoxArray<NDIM>& boxes = l0_level->getBoxes();
    const SAMRAI::hier::ProcessorMapping& procs = l0_level->getProcessorMapping();
    const int n_boxes = boxes.size();
    boxes_out.resize(n_boxes);
    procs_out.resize(n_boxes);
    for (int i = 0; i < n_boxes; ++i)
    {
        boxes_out[i] = boxes[i];
        procs_out[i] = procs.getProcessorAssignment(i);
    }
}

bool
matches_l0_boxes_and_procs(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>& l0_level,
                           const std::vector<SAMRAI::hier::Box<NDIM>>& cached_boxes,
                           const std::vector<int>& cached_procs)
{
    const SAMRAI::hier::BoxArray<NDIM>& boxes = l0_level->getBoxes();
    const SAMRAI::hier::ProcessorMapping& procs = l0_level->getProcessorMapping();
    const int n_boxes = boxes.size();
    if (static_cast<int>(cached_boxes.size()) != n_boxes || static_cast<int>(cached_procs.size()) != n_boxes)
    {
        return false;
    }
    for (int i = 0; i < n_boxes; ++i)
    {
        if (boxes[i] != cached_boxes[i] || procs.getProcessorAssignment(i) != cached_procs[i])
        {
            return false;
        }
    }
    return true;
}

bool
snapshot_l0_robin_bc_matrix_coefs(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>& l0_level,
                                  const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                  const double solution_time,
                                  std::vector<double>& bc_vals_out)
{
    bc_vals_out.clear();
    if (!l0_level)
    {
        return false;
    }

    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(l0_level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = l0_level->getPatch(p());
        if (!patch)
        {
            return false;
        }
        const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>> physical_codim1_boxes =
            IBTK::PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();
        if (n_physical_codim1_boxes == 0)
        {
            continue;
        }

        if (static_cast<int>(bc_coefs.size()) < NDIM)
        {
            return false;
        }
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        const double* const patch_x_lower = pgeom->getXLower();
        const double* const patch_x_upper = pgeom->getXUpper();
        const SAMRAI::hier::IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
        SAMRAI::tbox::Array<SAMRAI::tbox::Array<bool>> touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            touches_regular_bdry[d].resizeArray(2);
            touches_periodic_bdry[d].resizeArray(2);
            for (int ul = 0; ul < 2; ++ul)
            {
                touches_regular_bdry[d][ul] = pgeom->getTouchesRegularBoundary(d, ul);
                touches_periodic_bdry[d][ul] = pgeom->getTouchesPeriodicBoundary(d, ul);
            }
        }

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const bc_coef = bc_coefs[axis];
            if (!bc_coef)
            {
                return false;
            }
            auto* const extended_bc_coef = dynamic_cast<IBTK::ExtendedRobinBcCoefStrategy*>(bc_coef);
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const SAMRAI::hier::BoundaryBox<NDIM> trimmed_bdry_box =
                    IBTK::PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                SAMRAI::hier::Box<NDIM> bc_coef_box =
                    IBTK::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                const bool shifted = (bdry_normal_axis != axis);
                if (shifted)
                {
                    bc_coef_box.upper()(axis) += 1;
                    std::array<double, NDIM> shifted_lower, shifted_upper;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        shifted_lower[d] = patch_x_lower[d];
                        shifted_upper[d] = patch_x_upper[d];
                    }
                    shifted_lower[axis] -= 0.5 * dx[axis];
                    shifted_upper[axis] -= 0.5 * dx[axis];
                    patch->setPatchGeometry(new SAMRAI::geom::CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                                           touches_regular_bdry,
                                                                                           touches_periodic_bdry,
                                                                                           dx,
                                                                                           shifted_lower.data(),
                                                                                           shifted_upper.data()));
                }

                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>> acoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM, double>(bc_coef_box, 1);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>> bcoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM, double>(bc_coef_box, 1);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>> gcoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM, double>(bc_coef_box, 1);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(true);
                }
                bc_coef->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, solution_time);
                if (shifted)
                {
                    patch->setPatchGeometry(pgeom);
                }

                for (SAMRAI::hier::Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const SAMRAI::hier::Index<NDIM>& idx = bc();
                    const double a = (*acoef_data)(idx, 0);
                    const double b = (*bcoef_data)(idx, 0);
                    if (!std::isfinite(a) || !std::isfinite(b) || (a == 0.0 && b == 0.0))
                    {
                        return false;
                    }
                    bc_vals_out.push_back(a);
                    bc_vals_out.push_back(b);
                }
            }
        }
    }
    return true;
}
} // namespace

SCPoissonPointRelaxationFACOperator::~SCPoissonPointRelaxationFACOperator()
{
    d_in_destructor = true;
    if (d_is_initialized)
    {
        deallocateOperatorState();
    }
    clearCachedCoarseState(/*deallocate_coarse_solver*/ true);
    if (d_reuse_unchanged_coarse_solver && (d_coarse_solver_reuse_count > 0 || d_coarse_solver_reinit_count > 0))
    {
        if (IBTK_MPI::getRank() == 0)
        {
            SAMRAI::tbox::pout << d_object_name << ": coarse solver reuse_count=" << d_coarse_solver_reuse_count
                               << " reinit_count=" << d_coarse_solver_reinit_count << "\n";
        }
    }
    return;
} // ~SCPoissonPointRelaxationFACOperator

void
SCPoissonPointRelaxationFACOperator::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    if (poisson_spec.cIsVariable())
    {
        TBOX_ERROR(d_object_name << "::setPoissonSpecifications():\n"
                                 << "  the patch data for C cannot be used in SCPoissonPointRelaxationFACOperator\n"
                                 << "  C must be zero or constant, and D must be constant\n");
    }
    if (poisson_spec.dIsVariable())
    {
        TBOX_ERROR(d_object_name << "::setPoissonSpecifications():\n"
                                 << "  the patch data for D cannot be used in SCPoissonPointRelaxationFACOperator\n"
                                 << "  C must be zero or constant, and D must be constant\n");
    }
    PoissonFACPreconditionerStrategy::setPoissonSpecifications(poisson_spec);
    return;
} // setPoissonSpecifications

long
SCPoissonPointRelaxationFACOperator::getCoarseSolverReuseCount() const
{
    return d_coarse_solver_reuse_count;
} // getCoarseSolverReuseCount

long
SCPoissonPointRelaxationFACOperator::getCoarseSolverReinitCount() const
{
    return d_coarse_solver_reinit_count;
} // getCoarseSolverReinitCount

void
SCPoissonPointRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                 SAMRAIVectorReal<NDIM, double>& residual,
                                                 int level_num,
                                                 int num_sweeps,
                                                 bool /*performing_pre_sweeps*/,
                                                 bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBTK_TIMER_START(t_smooth_error);

    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(level_num);
    const int error_idx = error.getComponentDescriptorIndex(0);
    const int scratch_idx = d_scratch_idx;

    // Determine the smoother type.
    const std::string& smoother_type_string = (level_num == d_coarsest_ln ? d_coarse_solver_type : d_smoother_type);
    const SmootherType smoother_type = get_smoother_type(smoother_type_string);
    const bool red_black_ordering = use_red_black_ordering(smoother_type);
    const bool update_local_data = do_local_data_update(smoother_type);

    // Count passes, not sweeps, for the caching test below: red-black ordering makes two passes per sweep.
    if (red_black_ordering) num_sweeps *= 2;
    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double>> scratch_data = patch->getPatchData(scratch_idx);
#if !defined(NDEBUG)
            TBOX_ASSERT(error_data->getGhostCellWidth() >= d_gcw);
            TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                scratch_data->getArrayData(axis).copy(error_data->getArrayData(axis),
                                                      d_patch_cf_bdry_ghost_boxes[level_num][patch_counter][axis],
                                                      IntVector<NDIM>(0));
            }
        }
    }

    // Smooth the error by the specified number of sweeps.
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
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<SideData<NDIM, double>> error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM, double>> scratch_data = patch->getPatchData(scratch_idx);
#if !defined(NDEBUG)
                    TBOX_ASSERT(error_data->getGhostCellWidth() >= d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        error_data->getArrayData(axis).copy(scratch_data->getArrayData(axis),
                                                            d_patch_cf_bdry_ghost_boxes[level_num][patch_counter][axis],
                                                            IntVector<NDIM>(0));
                    }
                }

                // Synchronize the copies of the values on shared faces, which differ after a patch sweep.
                xeqScheduleDataSynch(error_idx, level_num);

                // Fill the non-coarse-fine interface ghost cell values.
                xeqScheduleGhostFillNoCoarse(error_idx, level_num);
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_cf_bdry_op->setPatchDataIndex(error_idx);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const IntVector<NDIM>& ghost_width_to_fill = d_gcw;
                d_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }
        else if (isweep > 0)
        {
            xeqScheduleDataSynch(error_idx, level_num);
            xeqScheduleGhostFillNoCoarse(error_idx, level_num);
        }

        // Smooth the error on the patches.
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double>> residual_data = residual.getComponentPatchData(0, *patch);
#if !defined(NDEBUG)
            TBOX_ASSERT(error_data->getGhostCellWidth() >= d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() >= d_gcw);
            TBOX_ASSERT(error_data->getDepth() == residual_data->getDepth());
#endif
            Pointer<SideData<NDIM, int>> mask_data = patch->getPatchData(d_mask_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            // Copy updated values from neighboring local patches.
            if (update_local_data)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const std::map<int, Box<NDIM>>& neighbor_overlap =
                        d_patch_neighbor_overlap[level_num][patch_counter][axis];
                    for (const auto& pair : neighbor_overlap)
                    {
                        const int src_patch_num = pair.first;
                        const Box<NDIM>& overlap = pair.second;
                        Pointer<Patch<NDIM>> src_patch = level->getPatch(src_patch_num);
                        Pointer<SideData<NDIM, double>> src_error_data = error.getComponentPatchData(0, *src_patch);
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

            // Smooth the error using Gauss-Seidel.
            const double& D = d_poisson_spec.getDConstant();
            const double& C = d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_patch_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
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
                                                 D,
                                                 C,
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
                                              D,
                                              C,
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
                            SMOOTH_GS_RB_CONST_DC_FC(U,
                                                     U_ghosts,
                                                     D,
                                                     C,
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
                            SMOOTH_GS_CONST_DC_FC(U,
                                                  U_ghosts,
                                                  D,
                                                  C,
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
} // smoothError

bool
SCPoissonPointRelaxationFACOperator::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                        SAMRAIVectorReal<NDIM, double>& residual,
                                                        int coarsest_ln)
{
    IBTK_TIMER_START(t_solve_coarsest_level);

#if !defined(NDEBUG)
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif
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
    }
    else
    {
        smoothError(error, residual, coarsest_ln, d_coarse_solver_max_iterations, false, false);
    }
    IBTK_TIMER_STOP(t_solve_coarsest_level);
    return true;
} // solveCoarsestLevel

void
SCPoissonPointRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                     SAMRAIVectorReal<NDIM, double>& solution,
                                                     const SAMRAIVectorReal<NDIM, double>& rhs,
                                                     int coarsest_level_num,
                                                     int finest_level_num)
{
    IBTK_TIMER_START(t_compute_residual);

    const int res_idx = residual.getComponentDescriptorIndex(0);
    const int sol_idx = solution.getComponentDescriptorIndex(0);
    const int rhs_idx = rhs.getComponentDescriptorIndex(0);

    const Pointer<SideVariable<NDIM, double>> res_var = residual.getComponentVariable(0);
    const Pointer<SideVariable<NDIM, double>> sol_var = solution.getComponentVariable(0);
    const Pointer<SideVariable<NDIM, double>> rhs_var = rhs.getComponentVariable(0);

    // Fill ghost-cell values.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    Pointer<SideNoCornersFillPattern> fill_pattern =
        new SideNoCornersFillPattern(d_gcw.max(), /*overwrite_interior*/ false);
    InterpolationTransactionComponent transaction_comp(sol_idx,
                                                       d_data_refine_type,
                                                       d_use_cf_interpolation,
                                                       d_data_coarsen_type,
                                                       d_bdry_extrap_type,
                                                       /*consistent_type_2_bdry*/ false,
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
    d_level_bdry_fill_ops[finest_level_num]->setHomogeneousBc(ALWAYS_HOMOGENEOUS_BC);
    d_level_bdry_fill_ops[finest_level_num]->fillData(d_solution_time);
    InterpolationTransactionComponent default_transaction_comp(d_solution->getComponentDescriptorIndex(0),
                                                               d_data_refine_type,
                                                               d_use_cf_interpolation,
                                                               d_data_coarsen_type,
                                                               d_bdry_extrap_type,
                                                               /*consistent_type_2_bdry*/ false,
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

    // Negating the coefficients and passing f as the additive source gives f - A*u.
    PoissonSpecifications negated_spec(d_object_name + "::negated_spec");
    negatePoissonSpecifications(negated_spec);
    d_level_math_ops[finest_level_num]->laplace(
        res_idx, res_var, negated_spec, sol_idx, sol_var, nullptr, d_solution_time, 1.0, rhs_idx, rhs_var);

    IBTK_TIMER_STOP(t_compute_residual);
    return;
} // computeResidual

/////////////////////////////// PROTECTED ////////////////////////////////////

void
SCPoissonPointRelaxationFACOperator::initializeOperatorStateSpecialized(const SAMRAIVectorReal<NDIM, double>& solution,
                                                                        const SAMRAIVectorReal<NDIM, double>& rhs,
                                                                        const int coarsest_reset_ln,
                                                                        const int finest_reset_ln)
{
    // Setup solution and rhs vectors.
    Pointer<SideVariable<NDIM, double>> solution_var = solution.getComponentVariable(0);
    Pointer<SideVariable<NDIM, double>> rhs_var = rhs.getComponentVariable(0);

    Pointer<SideDataFactory<NDIM, double>> solution_pdat_fac = solution_var->getPatchDataFactory();
    Pointer<SideDataFactory<NDIM, double>> rhs_pdat_fac = rhs_var->getPatchDataFactory();

#if !defined(NDEBUG)
    TBOX_ASSERT(solution_var);
    TBOX_ASSERT(rhs_var);
    TBOX_ASSERT(solution_pdat_fac);
    TBOX_ASSERT(rhs_pdat_fac);
#endif

    if (solution_pdat_fac->getDefaultDepth() != rhs_pdat_fac->getDefaultDepth())
    {
        TBOX_ERROR(d_object_name << "::initializeOperatorState():\n"
                                 << "  solution and rhs vectors must have the same data depths\n"
                                 << "  solution data depth = " << solution_pdat_fac->getDefaultDepth() << "\n"
                                 << "  rhs      data depth = " << rhs_pdat_fac->getDefaultDepth() << std::endl);
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideDataFactory<NDIM, double>> scratch_pdat_fac =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_scratch_idx);
    scratch_pdat_fac->setDefaultDepth(solution_pdat_fac->getDefaultDepth());

    // Setup cached BC data.
    d_bc_helper = new StaggeredPhysicalBoundaryHelper();
    d_bc_helper->cacheBcCoefData(d_bc_coefs, d_solution_time, d_hierarchy);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_mask_idx)) level->allocatePatchData(d_mask_idx);
    }
    d_bc_helper->setupMaskingFunction(d_mask_idx, d_coarsest_ln, d_finest_ln);

    // Initialize the coarse level solvers when needed.
    if (coarsest_reset_ln == d_coarsest_ln && d_coarse_solver)
    {
        bool reused = false;
        bool cache_ready = false;
        const bool was_cached = d_has_cached_coarse_state;
        const double rel_tol = d_coarse_solver_reuse_rel_tol;
        double max_rel_dC = 0.0;
        double max_rel_dD = 0.0;
        double max_rel_dBC = 0.0;
        CachedCoarseSolverKind cur_solver_kind = CachedCoarseSolverKind::UNSUPPORTED;
        IBTK::VCInterpType cur_mu_interp_type = IBTK::VC_HARMONIC_INTERP;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l0_level;
        double cur_dt = 0.0;
        bool c_is_var = false;
        bool c_is_zero = true;
        double c_const = 0.0;
        bool d_is_var = false;
        double d_const = 0.0;
        int cur_x_depth = -1;
        SAMRAI::hier::IntVector<NDIM> cur_x_gcw(-1);
        int cur_C_id = IBTK::invalid_index;
        int cur_D_id = IBTK::invalid_index;
        std::vector<double> cur_C_vals;
        std::vector<double> cur_D_vals;
        std::vector<double> cur_bc_vals;

        if (d_reuse_unchanged_coarse_solver && d_coarsest_ln == 0 && d_hierarchy)
        {
            cur_solver_kind = classify_coarse_solver(d_coarse_solver.getPointer(),
                                                     CachedCoarseSolverKind::UNSUPPORTED,
                                                     CachedCoarseSolverKind::VCSC_PETSC,
                                                     cur_mu_interp_type);
            l0_level = d_hierarchy->getPatchLevel(0);
            const bool have_dt = std::isfinite(d_current_time) && std::isfinite(d_new_time);
            cur_dt = have_dt ? (d_new_time - d_current_time) : 0.0;
            c_is_var = d_poisson_spec.cIsVariable();
            c_is_zero = d_poisson_spec.cIsZero();
            c_const = (!c_is_var && !c_is_zero) ? d_poisson_spec.getCConstant() : 0.0;
            d_is_var = d_poisson_spec.dIsVariable();
            const bool d_is_const = !d_is_var && d_poisson_spec.dIsConstant();
            d_const = d_is_const ? d_poisson_spec.getDConstant() : 0.0;
            cur_x_depth = solution_pdat_fac->getDefaultDepth();
            cur_x_gcw = solution_pdat_fac->getGhostCellWidth();
            cur_C_id = c_is_var ? d_poisson_spec.getCPatchDataId() : IBTK::invalid_index;
            cur_D_id = d_is_var ? d_poisson_spec.getDPatchDataId() : IBTK::invalid_index;

            bool local_snapshot_ok = (cur_solver_kind != CachedCoarseSolverKind::UNSUPPORTED) &&
                                     static_cast<bool>(l0_level) && have_dt && (d_is_var || d_is_const);
            if (local_snapshot_ok && c_is_var)
            {
                local_snapshot_ok = snapshot_l0_coeff_values(l0_level, cur_C_id, cur_C_vals);
            }
            if (local_snapshot_ok && d_is_var)
            {
                local_snapshot_ok = snapshot_l0_coeff_values(l0_level, cur_D_id, cur_D_vals);
            }
            if (local_snapshot_ok)
            {
                local_snapshot_ok =
                    snapshot_l0_robin_bc_matrix_coefs(l0_level, d_bc_coefs, d_solution_time, cur_bc_vals);
            }

            int snapshot_fail = local_snapshot_ok ? 0 : 1;
            snapshot_fail = IBTK_MPI::maxReduction(snapshot_fail);
            cache_ready = (snapshot_fail == 0);

            if (cache_ready && was_cached)
            {
                bool local_match = (cur_solver_kind == d_cached_l0_solver_kind) &&
                                   (cur_solver_kind != CachedCoarseSolverKind::VCSC_PETSC ||
                                    cur_mu_interp_type == d_cached_l0_mu_interp_type) &&
                                   (l0_level == d_cached_l0_level) &&
                                   matches_l0_boxes_and_procs(l0_level, d_cached_l0_boxes, d_cached_l0_procs) &&
                                   scalars_match_rel(cur_dt, d_cached_l0_dt, rel_tol) &&
                                   (cur_x_depth == d_cached_l0_x_depth) && (cur_x_gcw == d_cached_l0_x_gcw) &&
                                   (c_is_var == d_cached_l0_C_is_var) && (c_is_zero == d_cached_l0_C_is_zero) &&
                                   (d_is_var == d_cached_l0_D_is_var);
                if (local_match && !c_is_var && !c_is_zero)
                {
                    local_match = scalars_match_rel(c_const, d_cached_l0_C_const, rel_tol);
                }
                if (local_match && !d_is_var)
                {
                    local_match = scalars_match_rel(d_const, d_cached_l0_D_const, rel_tol);
                }
                const bool c_vals_match =
                    (!c_is_var || !d_cached_l0_C_is_var) ?
                        true :
                        coeff_vectors_match_rel(cur_C_vals, d_cached_l0_C_vals, rel_tol, max_rel_dC);
                const bool d_vals_match =
                    (!d_is_var || !d_cached_l0_D_is_var) ?
                        true :
                        coeff_vectors_match_rel(cur_D_vals, d_cached_l0_D_vals, rel_tol, max_rel_dD);
                const bool bc_vals_match = bc_coefs_match_exact(cur_bc_vals, d_cached_l0_bc_vals, max_rel_dBC);
                if (local_match && c_is_var)
                {
                    local_match = c_vals_match;
                }
                if (local_match && d_is_var)
                {
                    local_match = d_vals_match;
                }
                if (local_match)
                {
                    local_match = bc_vals_match;
                }

                int local_mismatch = local_match ? 0 : 1;
                local_mismatch = IBTK_MPI::maxReduction(local_mismatch);
                reused = (local_mismatch == 0);
                if (!reused)
                {
                    max_rel_dC = IBTK_MPI::maxReduction(max_rel_dC);
                    max_rel_dD = IBTK_MPI::maxReduction(max_rel_dD);
                    max_rel_dBC = IBTK_MPI::maxReduction(max_rel_dBC);
                }
            }
        }

        if (reused)
        {
            d_coarse_solver->setSolutionTime(d_solution_time);
            d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
            d_coarse_solver->setPoissonSpecifications(d_poisson_spec);
            d_coarse_solver->setPhysicalBcCoefs(d_bc_coefs);
            ++d_coarse_solver_reuse_count;
            if (d_reuse_unchanged_coarse_solver && d_coarsest_ln == 0 && IBTK_MPI::getRank() == 0)
            {
                SAMRAI::tbox::plog << d_object_name
                                   << ": coarse solver reused (reuse_count=" << d_coarse_solver_reuse_count
                                   << ", reinit_count=" << d_coarse_solver_reinit_count << ")\n";
            }
        }
        else
        {
            clearCachedCoarseState(/*deallocate_coarse_solver*/ true);
            // Note that since the coarse level solver is solving for the error, it
            // must always employ homogeneous boundary conditions.
            d_coarse_solver->setSolutionTime(d_solution_time);
            d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
            d_coarse_solver->setPoissonSpecifications(d_poisson_spec);
            d_coarse_solver->setPhysicalBcCoefs(d_bc_coefs);
            d_coarse_solver->setHomogeneousBc(ALWAYS_HOMOGENEOUS_BC);

            Pointer<SAMRAIVectorReal<NDIM, double>> l0_solution = getLevelSAMRAIVectorReal(*d_solution, d_coarsest_ln);
            Pointer<SAMRAIVectorReal<NDIM, double>> l0_rhs = getLevelSAMRAIVectorReal(*d_rhs, d_coarsest_ln);
            if (cache_ready)
            {
                d_cached_l0_solution = l0_solution->cloneVector(l0_solution->getName() + "::cached");
                d_cached_l0_rhs = l0_rhs->cloneVector(l0_rhs->getName() + "::cached");
                d_coarse_solver->initializeSolverState(*d_cached_l0_solution, *d_cached_l0_rhs);
                d_has_cached_coarse_state = true;
                d_cached_l0_solver_kind = cur_solver_kind;
                d_cached_l0_mu_interp_type = cur_mu_interp_type;
                d_cached_l0_level = l0_level;
                capture_l0_boxes_and_procs(l0_level, d_cached_l0_boxes, d_cached_l0_procs);
                d_cached_l0_dt = cur_dt;
                d_cached_l0_C_is_var = c_is_var;
                d_cached_l0_C_is_zero = c_is_zero;
                d_cached_l0_C_const = c_const;
                d_cached_l0_D_is_var = d_is_var;
                d_cached_l0_D_const = d_const;
                d_cached_l0_x_depth = cur_x_depth;
                d_cached_l0_x_gcw = cur_x_gcw;
                d_cached_l0_C_vals.swap(cur_C_vals);
                d_cached_l0_D_vals.swap(cur_D_vals);
                d_cached_l0_bc_vals.swap(cur_bc_vals);
            }
            else
            {
                d_coarse_solver->initializeSolverState(*l0_solution, *l0_rhs);
            }

            if (d_reuse_unchanged_coarse_solver && d_coarsest_ln == 0)
            {
                ++d_coarse_solver_reinit_count;
                if (IBTK_MPI::getRank() == 0)
                {
                    SAMRAI::tbox::plog << d_object_name
                                       << ": coarse solver reinitialized (reuse_count=" << d_coarse_solver_reuse_count
                                       << ", reinit_count=" << d_coarse_solver_reinit_count;
                    if (was_cached)
                    {
                        SAMRAI::tbox::plog << ", cache_miss=1, max_rel_dC=" << max_rel_dC
                                           << ", max_rel_dD=" << max_rel_dD << ", max_rel_dBC=" << max_rel_dBC;
                    }
                    SAMRAI::tbox::plog << ")\n";
                }
            }
        }
    }

    // Setup specialized transfer operators.
    Pointer<CartesianGridGeometry<NDIM>> geometry = d_hierarchy->getGridGeometry();
    IBTK_DO_ONCE(geometry->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen()););

    // Setup coarse-fine interface and physical boundary operators.
    d_cf_bdry_op = new CartSideDoubleQuadraticCFInterpolation();
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);
    d_bc_op = new CartSideRobinPhysBdryOp(d_scratch_idx, d_bc_coefs, false);

    // Setup fill pattern spec objects.
    d_op_stencil_fill_pattern = new SideNoCornersFillPattern(d_gcw.max(), /*overwrite_interior*/ true);
    d_synch_fill_pattern = new SideSynchCopyFillPattern();

    // Find the ghost sides along the coarse-fine interface of each patch.
    d_patch_cf_bdry_ghost_boxes.resize(d_finest_ln + 1);
    for (int ln = std::max(coarsest_reset_ln, d_coarsest_ln + 1); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_cf_bdry_ghost_boxes[ln].assign(num_local_patches, std::array<BoxList<NDIM>, NDIM>());
        const CoarseFineBoundary<NDIM> cf_boundary(*d_hierarchy, ln, IntVector<NDIM>(1));
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Array<BoundaryBox<NDIM>>& cf_bdry_codim1_boxes =
                cf_boundary.getBoundaries(patch->getPatchNumber(), /*boundary type*/ 1);
            for (int k = 0; k < cf_bdry_codim1_boxes.size(); ++k)
            {
                const BoundaryBox<NDIM>& bdry_box = cf_bdry_codim1_boxes[k];

                // Remove from the boundary box the corner cells, which are not part of the interface. The location
                // index is 2 * (axis normal to the face) + side.
                const int bdry_axis = bdry_box.getLocationIndex() / 2;
                Box<NDIM> face_box = patch_box;
                face_box.grow(bdry_axis, 1);
                const Box<NDIM> cf_box = bdry_box.getBox() * face_box;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    // Only the sides outside the patch box are ghost sides: the sides on the patch boundary belong to
                    // the patch.
                    BoxList<NDIM> ghost_sides(SideGeometry<NDIM>::toSideBox(cf_box, axis));
                    ghost_sides.removeIntersections(SideGeometry<NDIM>::toSideBox(patch_box, axis));
                    for (BoxList<NDIM>::Iterator b(ghost_sides); b; b++)
                    {
                        d_patch_cf_bdry_ghost_boxes[ln][patch_counter][axis].addItem(b());
                    }
                }
            }
        }
    }

    // Get overlap information for re-setting patch boundary conditions during
    // smoothing.
    d_patch_neighbor_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_neighbor_overlap[ln].resize(num_local_patches);
        int patch_counter1 = 0;
        for (PatchLevel<NDIM>::Iterator p1(level); p1; p1++, ++patch_counter1)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                d_patch_neighbor_overlap[ln][patch_counter1][axis].clear();
            }
            Pointer<Patch<NDIM>> dst_patch = level->getPatch(p1());
            const Box<NDIM>& dst_patch_box = dst_patch->getBox();
            const Box<NDIM>& dst_ghost_box = Box<NDIM>::grow(dst_patch_box, 1);
            int patch_counter2 = 0;
            for (PatchLevel<NDIM>::Iterator p2(level); patch_counter2 < patch_counter1; p2++, ++patch_counter2)
            {
                Pointer<Patch<NDIM>> src_patch = level->getPatch(p2());
                const Box<NDIM>& src_patch_box = src_patch->getBox();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const Box<NDIM> overlap = SideGeometry<NDIM>::toSideBox(dst_ghost_box, axis) *
                                              SideGeometry<NDIM>::toSideBox(src_patch_box, axis);
                    if (!overlap.empty())
                    {
                        d_patch_neighbor_overlap[ln][patch_counter1][axis][p2()] = overlap;
                    }
                }
            }
        }
    }
    return;
} // initializeOperatorStateSpecialized

void
SCPoissonPointRelaxationFACOperator::deallocateOperatorStateSpecialized(const int /*coarsest_reset_ln*/,
                                                                        const int /*finest_reset_ln*/)
{
    if (!d_is_initialized) return;

    if (!d_in_initialize_operator_state)
    {
        d_patch_cf_bdry_ghost_boxes.clear();
        d_patch_neighbor_overlap.clear();
        if (d_coarse_solver && (d_in_destructor || !d_reuse_unchanged_coarse_solver || !d_has_cached_coarse_state))
        {
            d_coarse_solver->deallocateSolverState();
            clearCachedCoarseState(/*deallocate_coarse_solver*/ false);
        }
    }
    return;
} // deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SCPoissonPointRelaxationFACOperator::clearCachedCoarseState(const bool deallocate_coarse_solver)
{
    if (deallocate_coarse_solver && d_coarse_solver && d_has_cached_coarse_state)
    {
        d_coarse_solver->deallocateSolverState();
    }
    d_has_cached_coarse_state = false;
    if (d_cached_l0_solution)
    {
        free_vector_components(*d_cached_l0_solution);
        d_cached_l0_solution.setNull();
    }
    if (d_cached_l0_rhs)
    {
        free_vector_components(*d_cached_l0_rhs);
        d_cached_l0_rhs.setNull();
    }
    d_cached_l0_solver_kind = CachedCoarseSolverKind::UNSUPPORTED;
    d_cached_l0_mu_interp_type = IBTK::VC_HARMONIC_INTERP;
    d_cached_l0_level.setNull();
    d_cached_l0_boxes.clear();
    d_cached_l0_procs.clear();
    d_cached_l0_x_depth = -1;
    d_cached_l0_x_gcw = SAMRAI::hier::IntVector<NDIM>(-1);
    d_cached_l0_C_vals.clear();
    d_cached_l0_D_vals.clear();
    d_cached_l0_bc_vals.clear();
    return;
} // clearCachedCoarseState

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
