// Filename: CCPoissonHypreSolver.C
// Last modified: <04.Sep.2006 03:03:48 boyce@bigboy.nyconnect.com>
// Created on 20 Jan 2006 by Boyce Griffith (boyce@boyce.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "CCPoissonHypreSolver.h"

// HYPRE INCLUDES
#ifndef included_sstruct_ls
#define included_sstruct_ls
#include <sstruct_ls.h>
#endif

// IBTK INCLUDES
#ifndef included_IBTK_config
#include <IBTK_config.h>
#endif

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <Index.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <VariableDatabase.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

#define USE_HO_FLUXES 1

namespace
{
// hypre data.
static const int VAR = 0;
static const int NVARS = 1;
static const int STENCIL_SIZE = (NDIM == 2 ? 5 : 7);

// Ghost cell width.
static const int GHOST_CELL_WIDTH = 1;

// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_set_matrix_coefficients;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_solve_system;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_solve_system_hypre;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_solver_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_deallocate_solver_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCPoissonHypreSolver::CCPoissonHypreSolver(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_context(SAMRAI::hier::VariableDatabase<NDIM>::getDatabase()->getContext(object_name+"::CONTEXT")),
      d_is_refined_idx(
          SAMRAI::hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
              new SAMRAI::pdat::CellVariable<NDIM,int>(object_name+"::is_refined"), d_context, GHOST_CELL_WIDTH)),
      d_scratch_idx(
          SAMRAI::hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
              new SAMRAI::pdat::CellVariable<NDIM,int>(object_name+"::scratch"), d_context, GHOST_CELL_WIDTH)),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_cf_boundary(),
      d_in_initialize_solver_state(false),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_poisson_spec(d_object_name+"::poisson_spec"),
      d_sol_depth(0),
      d_rhs_depth(0),
      d_grid(NULL),
      d_graph(NULL),
      d_stencil(NULL),
      d_matrix(NULL),
      d_rhs_vec(NULL),
      d_sol_vec(NULL),
      d_solver(NULL),
      d_precond(NULL),
      d_nparts(-1),
      d_plevels(NULL),
      d_rfactors(NULL),
      d_solver_type("SysPFMG"),
      d_precond_type("none"),
      d_max_iterations(10),
      d_rel_residual_tol(1.0e-6),
      d_initial_guess_nonzero(false),
      d_rel_change(0),
      d_num_pre_relax_steps(1),
      d_num_post_relax_steps(1),
      d_memory_use(0),
      d_relax_type(1),         /* SysPFMG:
                                * 0 = Jacobi
                                * 1 = weighted Jacobi (default)
                                * 2 = red-black Gauss-Seidel
                                * 3 = red-black Gauss-Seidel (nonsymmetric)
                                * FAC:
                                * 1 = Jacobi
                                * 2 = Gauss-Seidel */
      d_skip_relax(1),         /* skip levels in SysPFMG (0 or 1) */
      d_coarse_solver_type(2), /* 1 = PCG + SysPFMG
                                * 2 = SysPFMG */
      d_two_norm(1),
      d_current_its(-1),
      d_current_residual_norm(-1.0),
      d_synch_soln(false),
      d_urestriction_coarsen_operator(),
      d_urestriction_coarsen_algorithm(),
      d_urestriction_coarsen_schedules(),
      d_print_solver_info(false),
      d_enable_logging(false)
{
    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Get values from the input database.
    if (!input_db.isNull())
    {
        d_print_solver_info = input_db->getBoolWithDefault("print_solver_info", d_print_solver_info);
        d_enable_logging = input_db->getBoolWithDefault("enable_logging", d_enable_logging);
        d_synch_soln = input_db->getBoolWithDefault("synch_soln", d_synch_soln);
        d_solver_type = input_db->getStringWithDefault("solver_type", d_solver_type);
        d_precond_type = input_db->getStringWithDefault("precond_type", d_precond_type);
        d_max_iterations = input_db->getIntegerWithDefault("max_iterations", d_max_iterations);
        d_rel_residual_tol = input_db->getDoubleWithDefault("relative_residual_tol", d_rel_residual_tol);
        d_initial_guess_nonzero = input_db->getBoolWithDefault("initial_guess_nonzero", d_initial_guess_nonzero);
        d_rel_change = input_db->getIntegerWithDefault("rel_change", d_rel_change);

        if (d_solver_type == "SysPFMG" || d_solver_type == "FAC" ||
            d_precond_type == "SysPFMG" || d_precond_type == "FAC" ||
            d_precond_type == "Split-SMG" || d_precond_type == "Split-PFMG")
        {
            d_num_pre_relax_steps = input_db->getIntegerWithDefault("num_pre_relax_steps", d_num_pre_relax_steps);
            d_num_post_relax_steps = input_db->getIntegerWithDefault("num_post_relax_steps", d_num_post_relax_steps);
        }

        if (d_solver_type == "SysPFMG" || d_solver_type == "FAC" ||
            d_precond_type == "SysPFMG" || d_precond_type == "FAC")
        {
            d_relax_type = input_db->getIntegerWithDefault("relax_type", d_relax_type);
            d_skip_relax = input_db->getIntegerWithDefault("skip_relax", d_skip_relax);
        }

        if (d_solver_type == "SysPFMG" || d_precond_type == "SysPFMG")
        {
            d_skip_relax = input_db->getIntegerWithDefault("skip_relax", d_skip_relax);
        }

        if (d_solver_type == "FAC" || d_precond_type == "FAC")
        {
            d_coarse_solver_type = input_db->getIntegerWithDefault("coarse_solver_type", d_coarse_solver_type);
        }

        if (d_solver_type == "PCG")
        {
            d_two_norm = input_db->getIntegerWithDefault("two_norm", d_two_norm);
        }
    }

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_set_matrix_coefficients = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBTK::CCPoissonHypreSolver::setMatrixCoefficients()");
        t_solve_system = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBTK::CCPoissonHypreSolver::solveSystem()");
        t_solve_system_hypre = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBTK::CCPoissonHypreSolver::solveSystem()[hypre]");
        t_initialize_solver_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBTK::CCPoissonHypreSolver::initializeSolverState()");
        t_deallocate_solver_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBTK::CCPoissonHypreSolver::deallocateSolverState()");
        timers_need_init = false;
    }
    return;
}// CCPoissonHypreSolver

CCPoissonHypreSolver::~CCPoissonHypreSolver()
{
    // Deallocate all data.
    deallocateSolverState();
    return;
}// ~CCPoissonHypreSolver

void
CCPoissonHypreSolver::setPoissonSpecifications(
    const SAMRAI::solv::PoissonSpecifications* poisson_spec)
{
    if (poisson_spec != NULL)
    {
        d_poisson_spec = *poisson_spec;
    }
    else
    {
        d_poisson_spec.setCZero();
        d_poisson_spec.setDConstant(-1.0);
    }
    return;
}// setPoissonSpecifications

void
CCPoissonHypreSolver::setResetLevels(
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert((coarsest_ln == -1 && finest_ln == -1) ||
           (coarsest_ln >=  0 && finest_ln >= coarsest_ln));
#endif
    if (d_is_initialized)
    {
        d_coarsest_reset_ln = coarsest_ln;
        d_finest_reset_ln = finest_ln;
    }
    return;
}// setResetLevels

bool
CCPoissonHypreSolver::solveSystem(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b)
{
    t_solve_system->start();

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = d_hierarchy.isNull();

    if (d_hierarchy.isNull())
    {
        initializeSolverState(x,b);
    }

    // Solve the system using the "legacy" solve method.
    const bool converged = solveSystem(x.getComponentDescriptorIndex(0),
                                       b.getComponentDescriptorIndex(0));

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve)
    {
        deallocateSolverState();
    }

    t_solve_system->stop();
    return converged;
}// solveSystem

namespace
{
inline
int getLevelPart(const int ln, const int coarsest_ln, const int finest_ln)
{
    return ln-coarsest_ln;
}
}

void
CCPoissonHypreSolver::initializeSolverState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b)
{
    t_initialize_solver_state->start();

    // Rudimentary error checking.
#ifdef DEBUG_CHECK_ASSERTIONS
    if (x.getNumberOfComponents() != b.getNumberOfComponents())
    {
        TBOX_ERROR(d_object_name
                   << ": vectors must have the same number of components.\n");
    }

    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >& patch_hierarchy =
        x.getPatchHierarchy();

    if (patch_hierarchy != b.getPatchHierarchy())
    {
        TBOX_ERROR(d_object_name
                   << ": vectors must have the same hierarchy.\n");
    }

    const int coarsest_ln = x.getCoarsestLevelNumber();

    if (coarsest_ln < 0)
    {
        TBOX_ERROR(d_object_name
                   << ": coarsest level number must not be negative.\n");
    }

    if (coarsest_ln != b.getCoarsestLevelNumber())
    {
        TBOX_ERROR(d_object_name
                   << ": vectors must have same coarsest level number.\n");
    }

    const int finest_ln = x.getFinestLevelNumber();

    if (finest_ln != b.getFinestLevelNumber())
    {
        TBOX_ERROR(d_object_name
                   << ": vectors must have same finest level number.\n");
    }

    if (finest_ln < coarsest_ln)
    {
        TBOX_ERROR(d_object_name
                   << ": finest level number must be >= coarsest level number.\n");
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (patch_hierarchy->getPatchLevel(ln).isNull())
        {
            TBOX_ERROR(d_object_name
                       <<": hierarchy level " << ln << " does not exist\n");
        }
    }
#endif

    // Indicate that we are in the initializeSolverState() method and
    // cache the level range to be reset.
    //
    // NOTE: The range of levels to be reset is itself reset by
    // deallocateSolverState().
    d_in_initialize_solver_state = true;
    const int coarsest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1
         ? d_coarsest_reset_ln
         : x.getCoarsestLevelNumber());
    const int finest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1
         ? d_finest_reset_ln
         : x.getFinestLevelNumber());

    // Deallocate (some) of the solver state.
    if (d_is_initialized) deallocateSolverState();

    // Reset the hierarchy configuration.
    d_hierarchy   = x.getPatchHierarchy();
    d_coarsest_ln = x.getCoarsestLevelNumber();
    d_finest_ln   = x.getFinestLevelNumber();

    // Setup the hypre level refinement data.
    d_nparts = d_finest_ln-d_coarsest_ln+1;
    d_plevels  = new int[d_nparts];
    d_rfactors = new int[d_nparts][3];
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int part = getLevelPart(ln,d_coarsest_ln,d_finest_ln);
        d_plevels[part] = ln-d_coarsest_ln;  // level 0 is the coarsest level

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& ratio_to_coarser = level->getRatioToCoarserLevel();
        for (int d = 0; d < NDIM; ++d)
        {
            d_rfactors[part][d] = (ln == d_coarsest_ln ? 1 : ratio_to_coarser[d]);
        }
        for (int d = NDIM; d < 3; ++d)
        {
            d_rfactors[part][d] = 1;
        }
    }

    // Initialize the coarse-fine boundary description for the
    // hierarchy.
    d_cf_boundary.resize(d_finest_ln+1);
    const SAMRAI::hier::IntVector<NDIM> max_gcw = GHOST_CELL_WIDTH;
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        d_cf_boundary[ln].computeFromHierarchy(*d_hierarchy, ln, max_gcw);
    }

    // Allocate patch data.
    for (int ln = SAMRAI::tbox::Utilities::imax(d_coarsest_ln, coarsest_reset_ln);
         ln <= finest_reset_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_is_refined_idx);
    }

    for (int ln = d_coarsest_ln; ln <= finest_reset_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->setTime(0.0,d_is_refined_idx);
        level->allocatePatchData(d_scratch_idx);
    }

    // Determine the coarsened-fine cells (i.e., the cells on coarser
    // level of the patch hierarchy that contain refined cells on some
    // finer level of the hierarchy).
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::hier::BoxArray<NDIM> refined_region_boxes;

        if (ln < d_finest_ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > next_finer_level =
                d_hierarchy->getPatchLevel(ln+1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->
                                         getRatioToCoarserLevel());
        }

        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > is_refined_data =
                patch->getPatchData(d_is_refined_idx);

            is_refined_data->fillAll(0);

            if (ln < d_finest_ln)
            {
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const SAMRAI::hier::Box<NDIM>& refined_box = refined_region_boxes(i);
                    const SAMRAI::hier::Box<NDIM> intersection = patch_box * refined_box;
                    if (!intersection.empty())
                    {
                        is_refined_data->fillAll(1, intersection);
                    }
                }
            }
        }
    }

    // Fill ghost cell values for the is_refined data.
    //
    // NOTE: It would be possible to avoid filling ghost cell values
    // here by intersecting the refined_region_boxes above with the
    // "ghost box" instead of the "patch box".  If this approach is
    // followed, however, care must be taken in dealing with periodic
    // boundaries.  It seems less error prone to fill the ghost cell
    // values via a refinement algorithm and schedule.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var;

    var_db->mapIndexToVariable(d_is_refined_idx, var);
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > is_refined_fill_op =
        geometry->lookupRefineOperator(var, "CONSTANT_REFINE");

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > is_refined_fill_alg =
        new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    is_refined_fill_alg->
        registerRefine(d_is_refined_idx,
                       d_is_refined_idx,
                       d_scratch_idx,
                       is_refined_fill_op);

    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        is_refined_fill_alg->createSchedule(level)->fillData(0.0);
    }

    // Walk the coarse-fine interface to reset ghost cells values on
    // the "coarse side" of the coarse-fine interface.
    for (int ln = d_coarsest_ln+1; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > is_refined_data =
                patch->getPatchData(d_is_refined_idx);
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::Box<NDIM>& ghost_box = is_refined_data->getGhostBox();

            const int bdry_type = 1;
            const int patch_ln  = patch->getPatchLevelNumber();
            const int patch_num = patch->getPatchNumber();

            const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >& bdry_boxes =
                d_cf_boundary[patch_ln].getBoundaries(patch_num,bdry_type);
            const int num_boxes = bdry_boxes.getSize();

            for (int box_num = 0; box_num < num_boxes; ++box_num)
            {
                const int loc_index = bdry_boxes[box_num].getLocationIndex();
                const SAMRAI::hier::Box<NDIM> bdry_box = trimBoundaryBox(bdry_boxes[box_num].getBox(),patch_box,loc_index);
                const SAMRAI::hier::Box<NDIM> intersection = ghost_box * bdry_box;

                if (!intersection.empty())
                {
                    is_refined_data->fillAll(-1, intersection);
                }
            }
        }
    }

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_scratch_idx);
    }

    // Reset current solution data.
    d_current_its = -1;
    d_current_residual_norm = -1.0;

    // Allocate the hypre solver data and set the hypre matrix
    // coefficients.
    allocateHypreData();
    setMatrixCoefficients();

    // Initialize the transfer operators, algorithms, and schedules.
    if (d_synch_soln)
    {
        var = x.getComponentVariable(0);
        d_urestriction_coarsen_operator = geometry->
            lookupCoarsenOperator(var, "CONSERVATIVE_COARSEN");

        d_urestriction_coarsen_algorithm = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
        d_urestriction_coarsen_algorithm->
            registerCoarsen(x.getComponentDescriptorIndex(0),
                            x.getComponentDescriptorIndex(0),
                            d_urestriction_coarsen_operator);

        d_urestriction_coarsen_schedules.resize(d_finest_ln+1);
        for (int dst_ln = d_coarsest_ln; dst_ln < d_finest_ln; ++dst_ln)
        {
            d_urestriction_coarsen_schedules[dst_ln] =
                d_urestriction_coarsen_algorithm->
                createSchedule(d_hierarchy->getPatchLevel(dst_ln),
                               d_hierarchy->getPatchLevel(dst_ln+1));
        }
    }

    // Indicate that we are not in the initializeSolverState() method.
    d_in_initialize_solver_state = false;

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    t_initialize_solver_state->stop();
    return;
}// initializeSolverState

void
CCPoissonHypreSolver::deallocateSolverState()
{
    if (d_is_initialized && !d_in_initialize_solver_state &&
        (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
    {
        return;
    }

    t_deallocate_solver_state->start();

    if (d_is_initialized)
    {
        const int coarsest_reset_ln =
            (d_in_initialize_solver_state &&
             (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
            ? d_coarsest_reset_ln : d_coarsest_ln;
        const int finest_reset_ln =
            (d_in_initialize_solver_state &&
             (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
            ? d_finest_reset_ln : d_finest_ln;

        for (int ln = coarsest_reset_ln;
             ln <= SAMRAI::tbox::Utilities::imin(d_finest_ln, finest_reset_ln); ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            if (level->checkAllocated(d_is_refined_idx))
            {
                level->deallocatePatchData(d_is_refined_idx);
            }
        }

        if (!d_in_initialize_solver_state ||
            (d_coarsest_reset_ln == -1) || (d_finest_reset_ln == -1))
        {
            d_cf_boundary.resize(0);

            d_hierarchy.setNull();
            d_coarsest_ln = -1;
            d_finest_ln   = -1;
        }
    }

    // Deallocate the hypre data.
    deallocateHypreData();

    // Clear the "reset level" range.
    d_coarsest_reset_ln = -1;
    d_finest_reset_ln = -1;

    // Indicate that the solver is not initialized.
    d_is_initialized = false;

    t_deallocate_solver_state->stop();
    return;
}// deallocateSolverState

void
CCPoissonHypreSolver::enableLogging(
    bool enabled)
{
    /*!
     * \todo implement this.
     */
    TBOX_WARNING("CCPoissonHypreSolver: " <<
                 "enableLogging() not implemented.\n");
    d_enable_logging = enabled;
    return;
}// enableLogging

void
CCPoissonHypreSolver::printClassData(
    std::ostream& os) const
{
    /*!
     * \todo implement this.
     */
    TBOX_WARNING("CCPoissonHypreSolver: " <<
                 "printClassData() not implemented.\n");
    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CCPoissonHypreSolver::setupHypreSolver()
{
    // Get the MPI communicator.
#ifdef HAVE_MPI
    MPI_Comm communicator = SAMRAI::tbox::MPI::getCommunicator();
#else
    MPI_Comm communicator;
#endif

    // Setup the hypre solver and preconditioner.
    if (d_solver_type == "PCG"     ||
        d_solver_type == "GMRES"   ||
        d_solver_type == "BiCGSTAB")
    {
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructSysPFMGCreate(communicator, &d_precond);
            HYPRE_SStructSysPFMGSetMaxIter(d_precond, 1);
            HYPRE_SStructSysPFMGSetTol(d_precond, 0.0);
            HYPRE_SStructSysPFMGSetRelaxType(d_precond, d_relax_type);
            HYPRE_SStructSysPFMGSetNumPreRelax(d_precond, d_num_pre_relax_steps);
            HYPRE_SStructSysPFMGSetNumPostRelax(d_precond, d_num_post_relax_steps);
            HYPRE_SStructSysPFMGSetSkipRelax(d_precond, d_skip_relax);
        }
        else if (d_precond_type == "FAC")
        {
            HYPRE_SStructFACCreate(communicator, &d_precond);
            HYPRE_SStructFACSetMaxLevels(d_precond, d_nparts);
            HYPRE_SStructFACSetMaxIter(d_precond, 1);
            HYPRE_SStructFACSetTol(d_precond, 0.0);
            HYPRE_SStructFACSetZeroGuess(d_precond);
            HYPRE_SStructFACSetRelChange(d_precond, d_rel_change);
            HYPRE_SStructFACSetRelaxType(d_precond, d_relax_type);
            HYPRE_SStructFACSetNumPreRelax(d_precond, d_num_pre_relax_steps);
            HYPRE_SStructFACSetNumPostRelax(d_precond, d_num_post_relax_steps);
            HYPRE_SStructFACSetCoarseSolverType(d_precond, d_coarse_solver_type);
            HYPRE_SStructFACSetPLevels(d_precond, d_nparts, d_plevels);
            HYPRE_SStructFACSetPRefinements(d_precond, d_nparts, d_rfactors);
            if (d_initial_guess_nonzero)
            {
                HYPRE_SStructFACSetNonZeroGuess(d_precond);
            }
            else
            {
                HYPRE_SStructFACSetZeroGuess(d_precond);
            }
        }
        else if ((d_precond_type == "Split-SMG") ||
                 (d_precond_type == "Split-PFMG"))
        {
            HYPRE_SStructSplitCreate(communicator, &d_precond);
            HYPRE_SStructSplitSetMaxIter(d_precond, 1);
            HYPRE_SStructSplitSetTol(d_precond, 0.0);
            HYPRE_SStructSplitSetZeroGuess(d_precond);
            if (d_precond_type == "Split-SMG")
            {
                HYPRE_SStructSplitSetStructSolver(d_precond, HYPRE_SMG);
            }
            else if (d_precond_type == "Split-PFMG")
            {
                HYPRE_SStructSplitSetStructSolver(d_precond, HYPRE_PFMG);
            }
        }
    }
    else
    {
        if (d_precond_type != "none")
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  invalid preconditioner type: " << d_precond_type << endl
                       << "  no preconditioners are supported for solver type: " << d_solver_type << endl);
        }
    }

    if (d_solver_type == "SysPFMG")
    {
        HYPRE_SStructSysPFMGCreate(communicator, &d_solver);
        HYPRE_SStructSysPFMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructSysPFMGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructSysPFMGSetRelChange(d_solver, d_rel_change);
        HYPRE_SStructSysPFMGSetRelaxType(d_solver, d_relax_type);
        HYPRE_SStructSysPFMGSetNumPreRelax(d_solver, d_num_pre_relax_steps);
        HYPRE_SStructSysPFMGSetNumPostRelax(d_solver, d_num_post_relax_steps);
        HYPRE_SStructSysPFMGSetSkipRelax(d_solver, d_skip_relax);
        if (d_initial_guess_nonzero)
        {
            HYPRE_SStructSysPFMGSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_SStructSysPFMGSetZeroGuess(d_solver);
        }
        HYPRE_SStructSysPFMGSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "FAC")
    {
        HYPRE_SStructFACCreate(communicator, &d_solver);
        HYPRE_SStructFACSetMaxLevels(d_solver, d_nparts);
        HYPRE_SStructFACSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructFACSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructFACSetRelChange(d_solver, d_rel_change);
        HYPRE_SStructFACSetRelaxType(d_solver, d_relax_type);
        HYPRE_SStructFACSetNumPreRelax(d_solver, d_num_pre_relax_steps);
        HYPRE_SStructFACSetNumPostRelax(d_solver, d_num_post_relax_steps);
        HYPRE_SStructFACSetCoarseSolverType(d_solver, d_coarse_solver_type);
        HYPRE_SStructFACSetPLevels(d_solver, d_nparts, d_plevels);
        HYPRE_SStructFACSetPRefinements(d_solver, d_nparts, d_rfactors);
        if (d_initial_guess_nonzero)
        {
            HYPRE_SStructFACSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_SStructFACSetZeroGuess(d_solver);
        }
        HYPRE_SStructFACSetup2(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_SStructPCGCreate(communicator, &d_solver);
        HYPRE_SStructPCGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructPCGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructPCGSetTwoNorm(d_solver, d_two_norm);
        HYPRE_SStructPCGSetRelChange(d_solver, d_rel_change);
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructPCGSetPrecond(d_solver,
                                       HYPRE_SStructSysPFMGSolve,
                                       HYPRE_SStructSysPFMGSetup,
                                       d_precond);
        }
        else if (d_precond_type == "FAC")
        {
            HYPRE_SStructPCGSetPrecond(d_solver,
                                       HYPRE_SStructFACSolve3,
                                       HYPRE_SStructFACSetup2,
                                       d_precond);
        }
        else if ((d_precond_type == "Split-SMG") ||
                 (d_precond_type == "Split-PFMG"))
        {
            HYPRE_SStructPCGSetPrecond(d_solver,
                                       HYPRE_SStructSplitSolve,
                                       HYPRE_SStructSplitSetup,
                                       d_precond);
        }
        else if (d_precond_type == "diagonal_scaling")
        {
            d_precond = NULL;
            HYPRE_SStructPCGSetPrecond(d_solver,
                                       HYPRE_SStructDiagScale,
                                       HYPRE_SStructDiagScaleSetup,
                                       d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = NULL;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown preconditioner type: " << d_precond_type << endl);
        }
        HYPRE_SStructPCGSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_SStructGMRESCreate(communicator, &d_solver);
        HYPRE_SStructGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructGMRESSetTol(d_solver, d_rel_residual_tol);
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructGMRESSetPrecond(d_solver,
                                         HYPRE_SStructSysPFMGSolve,
                                         HYPRE_SStructSysPFMGSetup,
                                         d_precond);
        }
        else if (d_precond_type == "FAC")
        {
            HYPRE_SStructGMRESSetPrecond(d_solver,
                                         HYPRE_SStructFACSolve3,
                                         HYPRE_SStructFACSetup2,
                                         d_precond);
        }
        else if ((d_precond_type == "Split-SMG") ||
                 (d_precond_type == "Split-PFMG"))
        {
            HYPRE_SStructGMRESSetPrecond(d_solver,
                                         HYPRE_SStructSplitSolve,
                                         HYPRE_SStructSplitSetup,
                                         d_precond);
        }
        else if (d_precond_type == "diagonal_scaling")
        {
            d_precond = NULL;
            HYPRE_SStructGMRESSetPrecond(d_solver,
                                         HYPRE_SStructDiagScale,
                                         HYPRE_SStructDiagScaleSetup,
                                         d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = NULL;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown preconditioner type: " << d_precond_type << endl);
        }
        HYPRE_SStructGMRESSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_SStructBiCGSTABCreate(communicator, &d_solver);
        HYPRE_SStructBiCGSTABSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructBiCGSTABSetTol(d_solver, d_rel_residual_tol);
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructBiCGSTABSetPrecond(d_solver,
                                            HYPRE_SStructSysPFMGSolve,
                                            HYPRE_SStructSysPFMGSetup,
                                            d_precond);
        }
        else if (d_precond_type == "FAC")
        {
            HYPRE_SStructBiCGSTABSetPrecond(d_solver,
                                            HYPRE_SStructFACSolve3,
                                            HYPRE_SStructFACSetup2,
                                            d_precond);
        }
        else if ((d_precond_type == "Split-SMG") ||
                 (d_precond_type == "Split-PFMG"))
        {
            HYPRE_SStructBiCGSTABSetPrecond(d_solver,
                                            HYPRE_SStructSplitSolve,
                                            HYPRE_SStructSplitSetup,
                                            d_precond);
        }
        else if (d_precond_type == "diagonal_scaling")
        {
            d_precond = NULL;
            HYPRE_SStructBiCGSTABSetPrecond(d_solver,
                                            HYPRE_SStructDiagScale,
                                            HYPRE_SStructDiagScaleSetup,
                                            d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = NULL;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown preconditioner type: " << d_precond_type << endl);
        }
        HYPRE_SStructBiCGSTABSetup(d_solver,d_matrix, d_rhs_vec, d_sol_vec);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  unknown solver type: " << d_solver_type << endl);
    }

    return;
}// setupHypreSolver

namespace
{
inline
int
coarsen(const int index, const int ratio)
{
    return index < 0 ? (index+1)/ratio-1 : index/ratio;
}// coarsen

inline
SAMRAI::hier::Index<NDIM>
coarsen(const SAMRAI::hier::Index<NDIM>& index, const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    SAMRAI::hier::Index<NDIM> coarsened_index;
    for (int d = 0; d < NDIM; ++d)
    {
        coarsened_index(d) = coarsen(index(d), ratio(d));
    }
    return coarsened_index;
}// coarsen

inline
int
refine(const int index, const int ratio)
{
    return index*ratio;
}// refine

inline
SAMRAI::hier::Index<NDIM> refine(const SAMRAI::hier::Index<NDIM>& index, const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    SAMRAI::hier::Index<NDIM> refined_index;
    for (int d = 0; d < NDIM; ++d)
    {
        refined_index(d) = refine(index(d), ratio(d));
    }
    return refined_index;
}// refine

inline
SAMRAI::hier::Index<NDIM>
periodicIndexFix(const SAMRAI::hier::Index<NDIM>& index, const SAMRAI::hier::BoxArray<NDIM>& domain, const SAMRAI::hier::IntVector<NDIM>& periodic_shift)
{
    if (domain.contains(index))
    {
        return index;
    }

    SAMRAI::hier::IntVector<NDIM> shift;
#if (NDIM == 3)
    for (int k = -1; k <= 1; ++k)
    {
        shift(2) = k*periodic_shift(2);
#endif
        for (int j = -1; j <= 1; ++j)
        {
            shift(1) = j*periodic_shift(1);
            for (int i = -1; i <= 1; ++i)
            {
                shift(0) = i*periodic_shift(0);

                const SAMRAI::hier::Index<NDIM> shifted_index = index + shift;
                if (domain.contains(shifted_index))
                {
                    return shifted_index;
                }
            }
        }
#if (NDIM == 3)
    }
#endif
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ERROR("this statement should not be reached!\n");
#endif
    return index;
}// periodicIndexFix
}

void
CCPoissonHypreSolver::allocateHypreData()
{
    // Get the MPI communicator.
#ifdef HAVE_MPI
    MPI_Comm communicator = SAMRAI::tbox::MPI::getCommunicator();
#else
    MPI_Comm communicator;
#endif

    // Setup the hypre grid.
    static const int nvars = NVARS;
    HYPRE_SStructVariable vartypes[nvars] = {HYPRE_SSTRUCT_VARIABLE_CELL};
    static const int var = VAR;

    HYPRE_SStructGridCreate(communicator, NDIM, d_nparts, &d_grid);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int part = getLevelPart(ln,d_coarsest_ln,d_finest_ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const SAMRAI::hier::Box<NDIM>& patch_box = level->getPatch(p())->getBox();
            SAMRAI::hier::Index<NDIM> lower = patch_box.lower();
            SAMRAI::hier::Index<NDIM> upper = patch_box.upper();
            HYPRE_SStructGridSetExtents(d_grid, part, lower, upper);
        }

//    The following code only allows the COARSEST level to touch the
//    periodic boundary.
//
//    if (ln == d_coarsest_ln)
//      {
//          SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
//          const SAMRAI::hier::IntVector<NDIM>& ratio = level->getRatio();
//          const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift(ratio);
//
//          if (periodic_shift != SAMRAI::hier::IntVector<NDIM>(0))
//          {
//              int hypre_periodic_shift[3] = {0 , 0 , 0};
//              for (int d = 0; d < NDIM; ++d)
//              {
//                  hypre_periodic_shift[d] = periodic_shift(d);
//              }
//              HYPRE_SStructGridSetPeriodic(d_grid, part, hypre_periodic_shift);
//          }
//      }
//      else
//      {
//          for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
//          {
//              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
//              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
//              if (pgeom->getTouchesPeriodicBoundary())
//              {
//                  TBOX_ERROR(d_object_name << ": Invalid hierarchy configuration:\n"
//                             << "patch " << patch->getPatchNumber() << " on level " << ln << "\n"
//                             << "touches the periodic boundary, but level " << ln << "\n"
//                             << "is not the coarsest level in the hierarchy!\n");
//              }
//          }
//      }

        const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift(level->getRatio());

        if (periodic_shift != SAMRAI::hier::IntVector<NDIM>(0))
        {
            int hypre_periodic_shift[3] = {0 , 0 , 0};
            for (int d = 0; d < NDIM; ++d)
            {
                hypre_periodic_shift[d] = periodic_shift(d);
            }
            HYPRE_SStructGridSetPeriodic(d_grid, part, hypre_periodic_shift);
        }

        HYPRE_SStructGridSetVariables(d_grid, part, nvars, vartypes);
    }

    HYPRE_SStructGridAssemble(d_grid);

    // Allocate graph and stencil data and set stencil offsets.
    HYPRE_SStructGraphCreate(communicator, d_grid, &d_graph);

    static const int stencil_size = STENCIL_SIZE;
#if (NDIM == 2)
    int stencil_offsets[stencil_size][2] = {
        {-1,0}, {0,-1}, {0,0}, {1,0}, {0,1}
    };
#endif
#if (NDIM == 3)
    int stencil_offsets[stencil_size][3] = {
        {-1,0,0}, {0,-1,0}, {0,0,-1}, {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}
    };
#endif

    HYPRE_SStructStencilCreate(NDIM, stencil_size, &d_stencil);
    for (int s = 0; s < stencil_size; ++s)
    {
        HYPRE_SStructStencilSetEntry(d_stencil, s, stencil_offsets[s], var);
    }

    // Setup the structured part of the graph (corresponding to
    // stencils that do not cross the coarse-fine interface).
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int part = getLevelPart(ln,d_coarsest_ln,d_finest_ln);
        HYPRE_SStructGraphSetStencil(d_graph, part, var, d_stencil);
    }

    // Setup the "unstructured" coarse-to-fine part of the graph
    // (corresponding to stencil entries that cross the coarse-fine
    // interface from the coarse grid to the fine grid).
    for (int ln = d_coarsest_ln; ln < d_finest_ln; ++ln)
    {
        const int f_part = getLevelPart(ln+1,d_coarsest_ln,d_finest_ln);
        const int c_part = getLevelPart(ln  ,d_coarsest_ln,d_finest_ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > finer_level = d_hierarchy->getPatchLevel(ln+1);
        const SAMRAI::hier::IntVector<NDIM>& ratio_to_finer = finer_level->getRatioToCoarserLevel();
        const int r_x = ratio_to_finer(0);
        const int r_y = ratio_to_finer(1);
#if (NDIM == 3)
        const int r_z = ratio_to_finer(2);
#endif
        const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift(finer_level->getRatio());
        SAMRAI::hier::BoxArray<NDIM> domain;
        grid_geometry->computePhysicalDomain(domain, finer_level->getRatio());

        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > is_refined_data =
                patch->getPatchData(d_is_refined_idx);
            for (SAMRAI::hier::Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const SAMRAI::hier::Index<NDIM>& i = b();
                if ((*is_refined_data)(i) == 0)
                {
                    SAMRAI::hier::Index<NDIM> f_idx;
                    SAMRAI::hier::Index<NDIM> c_idx = i;
                    SAMRAI::hier::Index<NDIM> refined_c_idx = refine(c_idx,ratio_to_finer);

                    // The x-faces.
                    SAMRAI::hier::Index<NDIM> ixl = i;  --ixl(0);
                    if ((*is_refined_data)(ixl) == 1)
                    {
                        f_idx(0) = refined_c_idx(0)-1;
                        for (int p_j = 0; p_j < r_y; ++p_j)
                        {
                            f_idx(1) = refined_c_idx(1) + p_j;
#if (NDIM == 3)
                            for (int p_k = 0; p_k < r_z; ++p_k)
                            {
                                f_idx(2) = refined_c_idx(2) + p_k;
#endif
                                f_idx = periodicIndexFix(f_idx,domain,periodic_shift);
                                HYPRE_SStructGraphAddEntries(d_graph,
                                                             c_part, c_idx, var,
                                                             f_part, f_idx, var);
#if (NDIM == 3)
                            }
#endif
                        }
                    }

                    SAMRAI::hier::Index<NDIM> ixu = i;  ++ixu(0);
                    if ((*is_refined_data)(ixu) == 1)
                    {
                        f_idx(0) = refined_c_idx(0)+r_x;
                        for (int p_j = 0; p_j < r_y; ++p_j)
                        {
                            f_idx(1) = refined_c_idx(1) + p_j;
#if (NDIM == 3)
                            for (int p_k = 0; p_k < r_z; ++p_k)
                            {
                                f_idx(2) = refined_c_idx(2) + p_k;
#endif
                                f_idx = periodicIndexFix(f_idx,domain,periodic_shift);
                                HYPRE_SStructGraphAddEntries(d_graph,
                                                             c_part, c_idx, var,
                                                             f_part, f_idx, var);
#if (NDIM == 3)
                            }
#endif
                        }
                    }

                    // The y-faces.
                    SAMRAI::hier::Index<NDIM> iyl = i;  --iyl(1);
                    if ((*is_refined_data)(iyl) == 1)
                    {
                        f_idx(1) = refined_c_idx(1)-1;
                        for (int p_i = 0; p_i < r_x; ++p_i)
                        {
                            f_idx(0) = refined_c_idx(0) + p_i;
#if (NDIM == 3)
                            for (int p_k = 0; p_k < r_z; ++p_k)
                            {
                                f_idx(2) = refined_c_idx(2) + p_k;
#endif
                                f_idx = periodicIndexFix(f_idx,domain,periodic_shift);
                                HYPRE_SStructGraphAddEntries(d_graph,
                                                             c_part, c_idx, var,
                                                             f_part, f_idx, var);
#if (NDIM == 3)
                            }
#endif
                        }
                    }

                    SAMRAI::hier::Index<NDIM> iyu = i;  ++iyu(1);
                    if ((*is_refined_data)(iyu) == 1)
                    {
                        f_idx(1) = refined_c_idx(1)+r_y;
                        for (int p_i = 0; p_i < r_x; ++p_i)
                        {
                            f_idx(0) = refined_c_idx(0) + p_i;
#if (NDIM == 3)
                            for (int p_k = 0; p_k < r_z; ++p_k)
                            {
                                f_idx(2) = refined_c_idx(2) + p_k;
#endif
                                f_idx = periodicIndexFix(f_idx,domain,periodic_shift);
                                HYPRE_SStructGraphAddEntries(d_graph,
                                                             c_part, c_idx, var,
                                                             f_part, f_idx, var);
#if (NDIM == 3)
                            }
#endif
                        }
                    }
#if (NDIM == 3)
                    // The z-faces.
                    SAMRAI::hier::Index<NDIM> izl = i;  --izl(2);
                    if ((*is_refined_data)(izl) == 1)
                    {
                        f_idx(2) = refined_c_idx(2)-1;
                        for (int p_i = 0; p_i < r_x; ++p_i)
                        {
                            f_idx(0) = refined_c_idx(0) + p_i;
                            for (int p_j = 0; p_j < r_y; ++p_j)
                            {
                                f_idx(1) = refined_c_idx(1) + p_j;

                                f_idx = periodicIndexFix(f_idx,domain,periodic_shift);
                                HYPRE_SStructGraphAddEntries(d_graph,
                                                             c_part, c_idx, var,
                                                             f_part, f_idx, var);
                            }
                        }
                    }

                    SAMRAI::hier::Index<NDIM> izu = i;  ++izu(2);
                    if ((*is_refined_data)(izu) == 1)
                    {
                        f_idx(2) = refined_c_idx(2)+r_z;
                        for (int p_i = 0; p_i < r_x; ++p_i)
                        {
                            f_idx(0) = refined_c_idx(0) + p_i;
                            for (int p_j = 0; p_j < r_y; ++p_j)
                            {
                                f_idx(1) = refined_c_idx(1) + p_j;

                                f_idx = periodicIndexFix(f_idx,domain,periodic_shift);
                                HYPRE_SStructGraphAddEntries(d_graph,
                                                             c_part, c_idx, var,
                                                             f_part, f_idx, var);
                            }
                        }
                    }
#endif
                }
            }
        }
    }

    // Setup the "unstructured" fine-to-coarse part of the graph
    // (corresponding to stencil entries that cross the coarse-fine
    // interface from the fine grid to the coarse grid).
    for (int ln = d_coarsest_ln+1; ln <= d_finest_ln; ++ln)
    {
        const int f_part = getLevelPart(ln  ,d_coarsest_ln,d_finest_ln);
        const int c_part = getLevelPart(ln-1,d_coarsest_ln,d_finest_ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level = d_hierarchy->getPatchLevel(ln-1);
        const SAMRAI::hier::IntVector<NDIM>& ratio_to_coarser = level->getRatioToCoarserLevel();
        const int r_x = ratio_to_coarser(0);
        const int r_y = ratio_to_coarser(1);
#if (NDIM == 3)
        const int r_z = ratio_to_coarser(2);
#endif
        const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift(coarser_level->getRatio());
        SAMRAI::hier::BoxArray<NDIM> domain;
        grid_geometry->computePhysicalDomain(domain, coarser_level->getRatio());

        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

            const int bdry_type = 1;
            const int patch_ln  = patch->getPatchLevelNumber();
            const int patch_num = patch->getPatchNumber();

            const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >& bdry_boxes =
                d_cf_boundary[patch_ln].getBoundaries(patch_num,bdry_type);
            const int num_boxes = bdry_boxes.getSize();

            for (int box_num = 0; box_num < num_boxes; ++box_num)
            {
                const int loc_index = bdry_boxes[box_num].getLocationIndex();
                const SAMRAI::hier::Box<NDIM> bdry_box = trimBoundaryBox(bdry_boxes[box_num].getBox(),patch_box,loc_index);
                const SAMRAI::hier::Index<NDIM>& blower = bdry_box.lower();
                const SAMRAI::hier::Index<NDIM>& bupper = bdry_box.upper();

                SAMRAI::hier::Index<NDIM> p_idx;

                switch (loc_index)
                {
                    case 0:   // x_lo
                    case 1:   // x_hi
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(blower(0) == bupper(0));
#endif
                        p_idx(0) = blower(0) + (loc_index == 0 ? 1 : -1);
                        for (int c_j = coarsen(blower(1),r_y);
                             c_j <= coarsen(bupper(1),r_y); ++c_j)
                        {
                            for (int p_j = 0; p_j < r_y; ++p_j)
                            {
                                p_idx(1) = r_y*c_j + p_j;
#if (NDIM == 3)
                                for (int c_k = coarsen(blower(2),r_z);
                                     c_k <= coarsen(bupper(2),r_z); ++c_k)
                                {
                                    for (int p_k = 0; p_k < r_z; ++p_k)
                                    {
                                        p_idx(2) = r_z*c_k + p_k;
#endif
                                        SAMRAI::hier::Index<NDIM> fine_c_idx = p_idx;
                                        fine_c_idx(0) += (loc_index == 0 ? -1 : 1);
                                        SAMRAI::hier::Index<NDIM> c_idx = coarsen(fine_c_idx,ratio_to_coarser);

                                        c_idx = periodicIndexFix(c_idx,domain,periodic_shift);
                                        HYPRE_SStructGraphAddEntries(d_graph,
                                                                     f_part, p_idx, var,
                                                                     c_part, c_idx, var);
#if (NDIM == 3)
                                    }
                                }
#endif
                            }
                        }
                        break;
                    case 2:   // y_lo
                    case 3:   // y_hi
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(blower(1) == bupper(1));
#endif
                        p_idx(1) = blower(1) + (loc_index == 2 ? 1 : -1);
                        for (int c_i = coarsen(blower(0),r_x);
                             c_i <= coarsen(bupper(0),r_x); ++c_i)
                        {
                            for (int p_i = 0; p_i < r_x; ++p_i)
                            {
                                p_idx(0) = r_x*c_i + p_i;
#if (NDIM == 3)
                                for (int c_k = coarsen(blower(2),r_z);
                                     c_k <= coarsen(bupper(2),r_z); ++c_k)
                                {
                                    for (int p_k = 0; p_k < r_z; ++p_k)
                                    {
                                        p_idx(2) = r_z*c_k + p_k;
#endif
                                        SAMRAI::hier::Index<NDIM> fine_c_idx = p_idx;
                                        fine_c_idx(1) += (loc_index == 2 ? -1 : 1);
                                        SAMRAI::hier::Index<NDIM> c_idx = coarsen(fine_c_idx,ratio_to_coarser);

                                        c_idx = periodicIndexFix(c_idx,domain,periodic_shift);
                                        HYPRE_SStructGraphAddEntries(d_graph,
                                                                     f_part, p_idx, var,
                                                                     c_part, c_idx, var);
#if (NDIM == 3)
                                    }
                                }
#endif
                            }
                        }
                        break;
#if (NDIM == 3)
                    case 4:   // z_lo
                    case 5:   // z_hi
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(blower(2) == bupper(2));
#endif
                        p_idx(2) = blower(2) + (loc_index == 4 ? 1 : -1);
                        for (int c_i = coarsen(blower(0),r_x);
                             c_i <= coarsen(bupper(0),r_x); ++c_i)
                        {
                            for (int p_i = 0; p_i < r_x; ++p_i)
                            {
                                p_idx(0) = r_x*c_i + p_i;
                                for (int c_j = coarsen(blower(1),r_y);
                                     c_j <= coarsen(bupper(1),r_y); ++c_j)
                                {
                                    for (int p_j = 0; p_j < r_y; ++p_j)
                                    {
                                        p_idx(1) = r_y*c_j + p_j;
                                        SAMRAI::hier::Index<NDIM> fine_c_idx = p_idx;
                                        fine_c_idx(2) += (loc_index == 4 ? -1 : 1);
                                        SAMRAI::hier::Index<NDIM> c_idx = coarsen(fine_c_idx,ratio_to_coarser);

                                        c_idx = periodicIndexFix(c_idx,domain,periodic_shift);
                                        HYPRE_SStructGraphAddEntries(d_graph,
                                                                     f_part, p_idx, var,
                                                                     c_part, c_idx, var);
                                    }
                                }
                            }
                        }
                        break;
#endif
                    default:
                        TBOX_ERROR("this statement should not be reached!\n");
                }
            }
        }
    }

    HYPRE_SStructGraphAssemble(d_graph);

    // Initialize the hypre matrix.
    HYPRE_SStructMatrixCreate(communicator, d_graph, &d_matrix);
    HYPRE_SStructMatrixInitialize(d_matrix);

    // Setup the hypre vectors
    HYPRE_SStructVectorCreate(communicator, d_grid, &d_sol_vec);
    HYPRE_SStructVectorInitialize(d_sol_vec);

    HYPRE_SStructVectorCreate(communicator, d_grid, &d_rhs_vec);
    HYPRE_SStructVectorInitialize(d_rhs_vec);

    return;
}// allocateHypreData

namespace
{
SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > getCData(
    const SAMRAI::solv::PoissonSpecifications& poisson_spec,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch)
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > C_data;
    if (poisson_spec.cIsZero() || poisson_spec.cIsConstant())
    {
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        C_data = new SAMRAI::pdat::CellData<NDIM,double>(patch_box, 1, 0);
        if (poisson_spec.cIsZero())
        {
            C_data->fillAll(0.0);
        }
        else
        {
            C_data->fillAll(poisson_spec.getCConstant());
        }
    }
    else
    {
        C_data = patch->getPatchData(poisson_spec.getCPatchDataId());
        if (C_data.isNull())
        {
            TBOX_ERROR("CCPoissonHypreSolver: Invalid patch descriptor index "
                       <<  poisson_spec.getCPatchDataId()
                       <<" for the C parameter.  It is not\n"
                       <<"cell-centered double data.");
        }
    }
    return C_data;
}// getCData

SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > getDData(
    const SAMRAI::solv::PoissonSpecifications& poisson_spec,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch)
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > D_data;
    if (poisson_spec.dIsConstant())
    {
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        D_data = new SAMRAI::pdat::SideData<NDIM,double>(patch_box, 1, 0);
        D_data->fillAll(poisson_spec.getDConstant());
    }
    else
    {
        D_data = patch->getPatchData(poisson_spec.getDPatchDataId());
        if (D_data.isNull())
        {
            TBOX_ERROR("CCPoissonHypreSolver: Invalid patch descriptor index "
                       <<  poisson_spec.getDPatchDataId()
                       << " for diffusion coefficient.  It is not\n"
                       << "side-centered double data.");
        }
    }
    return D_data;
}// getDData
}

void
CCPoissonHypreSolver::setMatrixCoefficients()
{
    t_set_matrix_coefficients->start();

    if (!d_in_initialize_solver_state)
    {
        TBOX_ERROR(d_object_name << ": Solver must be in initializeSolverState to make call to setMatrixCoefficients().\n");
    }

    // There is only one variable (scalar Poisson-type problem).
    const int var = VAR;

    // The size of the structured stencil is the same on all levels of
    // the composite grid.
    static const int stencil_size = STENCIL_SIZE;

    // Setup the structured hypre matrix entries (corresponding to
    // stencils that do not cross the coarse-fine interface).
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int part = getLevelPart(ln,d_coarsest_ln,d_finest_ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
#if (NDIM == 2)
            const double cell_vol = dx[0]*dx[1];
            const double xfactor  = dx[1]/dx[0];
            const double yfactor  = dx[0]/dx[1];
#endif
#if (NDIM == 3)
            const double cell_vol = dx[0]*dx[1]*dx[2];
            const double xfactor  = dx[1]*dx[2]/dx[0];
            const double yfactor  = dx[0]*dx[2]/dx[1];
            const double zfactor  = dx[0]*dx[1]/dx[2];
#endif
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > C_data = getCData(d_poisson_spec,patch);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > D_data = getDData(d_poisson_spec,patch);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > is_refined_data =
                patch->getPatchData(d_is_refined_idx);

            // Set the "structured" hypre matrix entries.
            for (SAMRAI::hier::Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                // Reset the stencil indices and the matrix entries.
                int stencil_indices[stencil_size];
                double mat_entries[stencil_size];
                for (int i = 0; i < stencil_size; ++i)
                {
                    stencil_indices[i] = i;
                    mat_entries[i] = 0.0;
                }

                SAMRAI::hier::Index<NDIM> i = b();
                if ((*is_refined_data)(i) != 0)
                {
                    mat_entries[NDIM] = 1.0;
                }
                else
                {
                    mat_entries[NDIM] = (*C_data)(i)*cell_vol;

                    static const int x_axis = 0;
                    static const int y_axis = 1;
                    static const int z_axis = 2;

                    (void) x_axis;
                    (void) y_axis;
                    (void) z_axis;

                    // The x-faces.
                    SAMRAI::hier::Index<NDIM> ixl = i;  --ixl(x_axis);
                    if ((*is_refined_data)(ixl) == 0)
                    {
                        const SAMRAI::pdat::SideIndex<NDIM> ixlower(
                            i, SAMRAI::pdat::SideIndex<NDIM>::X, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                        const double Dxl = (*D_data)(ixlower)*xfactor;
                        mat_entries[0]      = +Dxl;
                        mat_entries[NDIM]  += -Dxl;
                    }

                    SAMRAI::hier::Index<NDIM> ixu = i;  ++ixu(x_axis);
                    if ((*is_refined_data)(ixu) == 0)
                    {
                        const SAMRAI::pdat::SideIndex<NDIM> ixupper(
                            i, SAMRAI::pdat::SideIndex<NDIM>::X, SAMRAI::pdat::SideIndex<NDIM>::Upper);
                        const double Dxu = (*D_data)(ixupper)*xfactor;
                        mat_entries[NDIM+1] = +Dxu;
                        mat_entries[NDIM]  += -Dxu;
                    }

                    // The y-faces.
                    SAMRAI::hier::Index<NDIM> iyl = i;  --iyl(y_axis);
                    if ((*is_refined_data)(iyl) == 0)
                    {
                        const SAMRAI::pdat::SideIndex<NDIM> iylower(
                            i, SAMRAI::pdat::SideIndex<NDIM>::Y, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                        const double Dyl = (*D_data)(iylower)*yfactor;
                        mat_entries[1]      = +Dyl;
                        mat_entries[NDIM]  += -Dyl;
                    }

                    SAMRAI::hier::Index<NDIM> iyu = i;  ++iyu(y_axis);
                    if ((*is_refined_data)(iyu) == 0)
                    {
                        const SAMRAI::pdat::SideIndex<NDIM> iyupper(
                            i, SAMRAI::pdat::SideIndex<NDIM>::Y, SAMRAI::pdat::SideIndex<NDIM>::Upper);
                        const double Dyu = (*D_data)(iyupper)*yfactor;
                        mat_entries[NDIM+2] = +Dyu;
                        mat_entries[NDIM]  += -Dyu;
                    }
#if (NDIM == 3)
                    // The z-faces.
                    SAMRAI::hier::Index<NDIM> izl = i;  --izl(z_axis);
                    if ((*is_refined_data)(izl) == 0)
                    {
                        const SAMRAI::pdat::SideIndex<NDIM> izlower(
                            i, SAMRAI::pdat::SideIndex<NDIM>::Z, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                        const double Dzl = (*D_data)(izlower)*zfactor;
                        mat_entries[2]      = +Dzl;
                        mat_entries[NDIM]  += -Dzl;
                    }

                    SAMRAI::hier::Index<NDIM> izu = i;  ++izu(z_axis);
                    if ((*is_refined_data)(izu) == 0)
                    {
                        const SAMRAI::pdat::SideIndex<NDIM> izupper(
                            i, SAMRAI::pdat::SideIndex<NDIM>::Z, SAMRAI::pdat::SideIndex<NDIM>::Upper);
                        const double Dzu = (*D_data)(izupper)*zfactor;
                        mat_entries[NDIM+3] = +Dzu;
                        mat_entries[NDIM]  += -Dzu;
                    }
#endif
                }

                HYPRE_SStructMatrixSetValues(
                    d_matrix, part, i, var, stencil_size, stencil_indices, mat_entries);
            }
        }
    }

    // Setup the "unstructured" coarse-to-fine hypre matrix entries
    // (corresponding to stencil entries that cross the coarse-fine
    // interface from the coarse grid to the fine grid).
    for (int ln = d_coarsest_ln; ln < d_finest_ln; ++ln)
    {
        const int part = getLevelPart(ln,d_coarsest_ln,d_finest_ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > finer_level = d_hierarchy->getPatchLevel(ln+1);
        const SAMRAI::hier::IntVector<NDIM>& ratio_to_finer = finer_level->getRatioToCoarserLevel();
        const int r_x = ratio_to_finer(0);
        const int r_y = ratio_to_finer(1);
#if (NDIM == 3)
        const int r_z = ratio_to_finer(2);
#endif
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
#if (NDIM == 2)
            const double xfactor = (dx[1]/dx[0])*(static_cast<double>(r_x)/static_cast<double>(r_y));
            const double yfactor = (dx[0]/dx[1])*(static_cast<double>(r_y)/static_cast<double>(r_x));
#endif
#if (NDIM == 3)
            const double xfactor = (dx[1]*dx[2]/dx[0])*(static_cast<double>(r_x)/static_cast<double>(r_y*r_z));
            const double yfactor = (dx[0]*dx[2]/dx[1])*(static_cast<double>(r_y)/static_cast<double>(r_x*r_z));
            const double zfactor = (dx[0]*dx[1]/dx[2])*(static_cast<double>(r_z)/static_cast<double>(r_x*r_y));
#endif
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > D_data = getDData(d_poisson_spec,patch);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > is_refined_data =
                patch->getPatchData(d_is_refined_idx);

            for (SAMRAI::hier::Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const SAMRAI::hier::Index<NDIM>& idx = b();
                if ((*is_refined_data)(idx) == 0)
                {
                    SAMRAI::hier::Index<NDIM> c_idx = idx;
                    const int f_offset = stencil_size;

                    // The x-faces.
                    SAMRAI::hier::Index<NDIM> ixl = idx;  --ixl(0);
                    if ((*is_refined_data)(ixl) == 1)
                    {
                        for (int p_j = 0; p_j < r_y; ++p_j)
                        {
                            const int j = (p_j%2 == 0 ? +1 : -1);
#if (NDIM == 3)
                            for (int p_k = 0; p_k < r_z; ++p_k)
                            {
                                const int k = (p_k%2 == 0 ? +1 : -1);
#endif
                                const SAMRAI::pdat::SideIndex<NDIM> ixlower(
                                    c_idx, SAMRAI::pdat::SideIndex<NDIM>::X, SAMRAI::pdat::SideIndex<NDIM>::Lower);

                                double D_c =  -2.0/static_cast<double>(r_x+1);
                                double D_ppp = 2.0/static_cast<double>(r_x+1);
                                double D_pqp = 0.0;
                                double D_ppq = 0.0;
#if (USE_HO_FLUXES == 1)
                                D_ppp += static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_x));
                                D_pqp -= static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_x));
#if (NDIM == 3)
                                D_ppp += static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_x));
                                D_ppq -= static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_x));
#endif
#endif
                                const double Dx = +(*D_data)(ixlower)*xfactor;
                                D_c   *= Dx;
                                D_ppp *= Dx;
                                D_pqp *= Dx;
                                D_ppq *= Dx;

                                int entry_c = NDIM;
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_c, &D_c);
#if (NDIM == 2)
                                int entry_ppp = f_offset+p_j;
                                int entry_pqp = f_offset+p_j+j;
#endif
#if (NDIM == 3)
                                int entry_ppp = f_offset+p_k+r_z*p_j;
                                int entry_pqp = f_offset+p_k+r_z*(p_j+j);
                                int entry_ppq = f_offset+(p_k+k)+r_z*p_j;
#endif
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppp, &D_ppp);

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_pqp, &D_pqp);
#if (NDIM == 3)
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppq, &D_ppq);
#endif
#if (NDIM == 3)
                            }
#endif
                        }
                    }

                    SAMRAI::hier::Index<NDIM> ixu = idx;  ++ixu(0);
                    if ((*is_refined_data)(ixu) == 1)
                    {
                        for (int p_j = 0; p_j < r_y; ++p_j)
                        {
                            const int j = (p_j%2 == 0 ? +1 : -1);
#if (NDIM == 3)
                            for (int p_k = 0; p_k < r_z; ++p_k)
                            {
                                const int k = (p_k%2 == 0 ? +1 : -1);
#endif
                                const SAMRAI::pdat::SideIndex<NDIM> ixupper(
                                    c_idx, SAMRAI::pdat::SideIndex<NDIM>::X, SAMRAI::pdat::SideIndex<NDIM>::Upper);

                                double D_c =  -2.0/static_cast<double>(r_x+1);
                                double D_ppp = 2.0/static_cast<double>(r_x+1);
                                double D_pqp = 0.0;
                                double D_ppq = 0.0;
#if (USE_HO_FLUXES == 1)
                                D_ppp += static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_x));
                                D_pqp -= static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_x));
#if (NDIM == 3)
                                D_ppp += static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_x));
                                D_ppq -= static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_x));
#endif
#endif
                                const double Dx = +(*D_data)(ixupper)*xfactor;
                                D_c   *= Dx;
                                D_ppp *= Dx;
                                D_pqp *= Dx;
                                D_ppq *= Dx;

                                int entry_c = NDIM;
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_c, &D_c);
#if (NDIM == 2)
                                int entry_ppp = f_offset+p_j;
                                int entry_pqp = f_offset+p_j+j;
#endif
#if (NDIM == 3)
                                int entry_ppp = f_offset+p_k+r_z*p_j;
                                int entry_pqp = f_offset+p_k+r_z*(p_j+j);
                                int entry_ppq = f_offset+(p_k+k)+r_z*p_j;
#endif
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppp, &D_ppp);

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_pqp, &D_pqp);
#if (NDIM == 3)
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppq, &D_ppq);
#endif
#if (NDIM == 3)
                            }
#endif
                        }
                    }

                    // The y-faces.
                    SAMRAI::hier::Index<NDIM> iyl = idx;  --iyl(1);
                    if ((*is_refined_data)(iyl) == 1)
                    {
                        for (int p_i = 0; p_i < r_x; ++p_i)
                        {
                            const int i = (p_i%2 == 0 ? +1 : -1);
#if (NDIM == 3)
                            for (int p_k = 0; p_k < r_z; ++p_k)
                            {
                                const int k = (p_k%2 == 0 ? +1 : -1);
#endif
                                const SAMRAI::pdat::SideIndex<NDIM> iylower(
                                    c_idx, SAMRAI::pdat::SideIndex<NDIM>::Y, SAMRAI::pdat::SideIndex<NDIM>::Lower);

                                double D_c =  -2.0/static_cast<double>(r_y+1);
                                double D_ppp = 2.0/static_cast<double>(r_y+1);
                                double D_qpp = 0.0;
                                double D_ppq = 0.0;
#if (USE_HO_FLUXES == 1)
                                D_ppp += static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_y));
                                D_qpp -= static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_y));
#if (NDIM == 3)
                                D_ppp += static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_y));
                                D_ppq -= static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_y));
#endif
#endif
                                const double Dy = +(*D_data)(iylower)*yfactor;
                                D_c   *= Dy;
                                D_ppp *= Dy;
                                D_qpp *= Dy;
                                D_ppq *= Dy;

                                int entry_c = NDIM;
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_c, &D_c);
#if (NDIM == 2)
                                int entry_ppp = f_offset+p_i;
                                int entry_qpp = f_offset+p_i+i;
#endif
#if (NDIM == 3)
                                int entry_ppp = f_offset+p_k+r_z*p_i;
                                int entry_qpp = f_offset+p_k+r_z*(p_i+i);
                                int entry_ppq = f_offset+(p_k+k)+r_z*p_i;
#endif
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppp, &D_ppp);

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_qpp, &D_qpp);
#if (NDIM == 3)
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppq, &D_ppq);
#endif
#if (NDIM == 3)
                            }
#endif
                        }
                    }

                    SAMRAI::hier::Index<NDIM> iyu = idx;  ++iyu(1);
                    if ((*is_refined_data)(iyu) == 1)
                    {
                        for (int p_i = 0; p_i < r_x; ++p_i)
                        {
                            const int i = (p_i%2 == 0 ? +1 : -1);
#if (NDIM == 3)
                            for (int p_k = 0; p_k < r_z; ++p_k)
                            {
                                const int k = (p_k%2 == 0 ? +1 : -1);
#endif
                                const SAMRAI::pdat::SideIndex<NDIM> iyupper(
                                    c_idx, SAMRAI::pdat::SideIndex<NDIM>::Y, SAMRAI::pdat::SideIndex<NDIM>::Upper);

                                double D_c =  -2.0/static_cast<double>(r_y+1);
                                double D_ppp = 2.0/static_cast<double>(r_y+1);
                                double D_qpp = 0.0;
                                double D_ppq = 0.0;
#if (USE_HO_FLUXES == 1)
                                D_ppp += static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_y));
                                D_qpp -= static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_y));
#if (NDIM == 3)
                                D_ppp += static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_y));
                                D_ppq -= static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_y));
#endif
#endif
                                const double Dy = +(*D_data)(iyupper)*yfactor;
                                D_c   *= Dy;
                                D_ppp *= Dy;
                                D_qpp *= Dy;
                                D_ppq *= Dy;

                                int entry_c = NDIM;
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_c, &D_c);
#if (NDIM == 2)
                                int entry_ppp = f_offset+p_i;
                                int entry_qpp = f_offset+p_i+i;
#endif
#if (NDIM == 3)
                                int entry_ppp = f_offset+p_k+r_z*p_i;
                                int entry_qpp = f_offset+p_k+r_z*(p_i+i);
                                int entry_ppq = f_offset+(p_k+k)+r_z*p_i;
#endif
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppp, &D_ppp);

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_qpp, &D_qpp);
#if (NDIM == 3)
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppq, &D_ppq);
#endif
#if (NDIM == 3)
                            }
#endif
                        }
                    }
#if (NDIM == 3)
                    // The z-faces.
                    SAMRAI::hier::Index<NDIM> izl = idx;  --izl(2);
                    if ((*is_refined_data)(izl) == 1)
                    {
                        for (int p_i = 0; p_i < r_x; ++p_i)
                        {
                            const int i = (p_i%2 == 0 ? +1 : -1);
                            for (int p_j = 0; p_j < r_y; ++p_j)
                            {
                                const int j = (p_j%2 == 0 ? +1 : -1);

                                const SAMRAI::pdat::SideIndex<NDIM> izlower(
                                    c_idx, SAMRAI::pdat::SideIndex<NDIM>::Z, SAMRAI::pdat::SideIndex<NDIM>::Lower);

                                double D_c =  -2.0/static_cast<double>(r_z+1);
                                double D_ppp = 2.0/static_cast<double>(r_z+1);
                                double D_qpp = 0.0;
                                double D_pqp = 0.0;
#if (USE_HO_FLUXES == 1)
                                D_ppp += static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_z));
                                D_qpp -= static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_z));

                                D_ppp += static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_z));
                                D_pqp -= static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_z));
#endif
                                const double Dz = +(*D_data)(izlower)*zfactor;
                                D_c   *= Dz;
                                D_ppp *= Dz;
                                D_qpp *= Dz;
                                D_pqp *= Dz;

                                int entry_c = NDIM;
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_c, &D_c);

                                int entry_ppp = f_offset+p_j+r_y*p_i;
                                int entry_qpp = f_offset+p_j+r_y*(p_i+i);
                                int entry_pqp = f_offset+(p_j+j)+r_y*p_i;

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppp, &D_ppp);

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_qpp, &D_qpp);

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_pqp, &D_pqp);
                            }
                        }
                    }

                    SAMRAI::hier::Index<NDIM> izu = idx;  ++izu(2);
                    if ((*is_refined_data)(izu) == 1)
                    {
                        for (int p_i = 0; p_i < r_x; ++p_i)
                        {
                            const int i = (p_i%2 == 0 ? +1 : -1);
                            for (int p_j = 0; p_j < r_y; ++p_j)
                            {
                                const int j = (p_j%2 == 0 ? +1 : -1);

                                const SAMRAI::pdat::SideIndex<NDIM> izupper(
                                    c_idx, SAMRAI::pdat::SideIndex<NDIM>::Z, SAMRAI::pdat::SideIndex<NDIM>::Upper);

                                double D_c =  -2.0/static_cast<double>(r_z+1);
                                double D_ppp = 2.0/static_cast<double>(r_z+1);
                                double D_qpp = 0.0;
                                double D_pqp = 0.0;
#if (USE_HO_FLUXES == 1)
                                D_ppp += static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_z));
                                D_qpp -= static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_z));

                                D_ppp += static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_z));
                                D_pqp -= static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_z));
#endif
                                const double Dz = +(*D_data)(izupper)*zfactor;
                                D_c   *= Dz;
                                D_ppp *= Dz;
                                D_qpp *= Dz;
                                D_pqp *= Dz;

                                int entry_c = NDIM;
                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_c, &D_c);

                                int entry_ppp = f_offset+p_j+r_y*p_i;
                                int entry_qpp = f_offset+p_j+r_y*(p_i+i);
                                int entry_pqp = f_offset+(p_j+j)+r_y*p_i;

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_ppp, &D_ppp);

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_qpp, &D_qpp);

                                HYPRE_SStructMatrixAddToValues(d_matrix,
                                                               part, c_idx, var, 1,
                                                               &entry_pqp, &D_pqp);
                            }
                        }
                    }
#endif
                }
            }
        }
    }

    // Setup the "unstructured" fine-to-coarse hypre matrix entries
    // (corresponding to stencil entries that cross the coarse-fine
    // interface from the fine grid to the coarse grid).
    for (int ln = d_coarsest_ln+1; ln <= d_finest_ln; ++ln)
    {
        const int part = getLevelPart(ln,d_coarsest_ln,d_finest_ln);
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& ratio_to_coarser = level->getRatioToCoarserLevel();
        const int r_x = ratio_to_coarser(0);
        const int r_y = ratio_to_coarser(1);
#if (NDIM == 3)
        const int r_z = ratio_to_coarser(2);
#endif
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
#if (NDIM == 2)
            const double xfactor = dx[1]/dx[0];
            const double yfactor = dx[0]/dx[1];
#endif
#if (NDIM == 3)
            const double xfactor = dx[1]*dx[2]/dx[0];
            const double yfactor = dx[0]*dx[2]/dx[1];
            const double zfactor = dx[0]*dx[1]/dx[2];
#endif
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > D_data = getDData(d_poisson_spec,patch);
            const int c_offset = stencil_size;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > entries_data =
                new SAMRAI::pdat::CellData<NDIM,int>(patch_box, 1, 0);
            entries_data->fillAll(c_offset);

            const int bdry_type = 1;
            const int patch_ln  = patch->getPatchLevelNumber();
            const int patch_num = patch->getPatchNumber();

            const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >& bdry_boxes =
                d_cf_boundary[patch_ln].getBoundaries(patch_num,bdry_type);
            const int num_boxes = bdry_boxes.getSize();

            for (int box_num = 0; box_num < num_boxes; ++box_num)
            {
                const int loc_index = bdry_boxes[box_num].getLocationIndex();
                const SAMRAI::hier::Box<NDIM> bdry_box = trimBoundaryBox(bdry_boxes[box_num].getBox(),patch_box,loc_index);
                const SAMRAI::hier::Index<NDIM>& blower = bdry_box.lower();
                const SAMRAI::hier::Index<NDIM>& bupper = bdry_box.upper();

                SAMRAI::hier::Index<NDIM> p_idx;

                switch (loc_index)
                {
                    case 0:   // x_lo
                    case 1:   // x_hi
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(blower(0) == bupper(0));
#endif
                        p_idx(0) = blower(0) + (loc_index == 0 ? 1 : -1);
                        for (int c_j = coarsen(blower(1),r_y);
                             c_j <= coarsen(bupper(1),r_y); ++c_j)
                        {
                            for (int p_j = 0; p_j < r_y; ++p_j)
                            {
                                p_idx(1) = r_y*c_j + p_j;
                                const int j = (p_j%2 == 0 ? +1 : -1);
#if (NDIM == 3)
                                for (int c_k = coarsen(blower(2),r_z);
                                     c_k <= coarsen(bupper(2),r_z); ++c_k)
                                {
                                    for (int p_k = 0; p_k < r_z; ++p_k)
                                    {
                                        p_idx(2) = r_z*c_k + p_k;
                                        const int k = (p_k%2 == 0 ? +1 : -1);
#endif
                                        const SAMRAI::pdat::SideIndex<NDIM> ixlower(
                                            p_idx, SAMRAI::pdat::SideIndex<NDIM>::X, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                                        const SAMRAI::pdat::SideIndex<NDIM> ixupper(
                                            p_idx, SAMRAI::pdat::SideIndex<NDIM>::X, SAMRAI::pdat::SideIndex<NDIM>::Upper);

                                        double D_c =  -2.0/static_cast<double>(r_x+1);
                                        double D_ppp = 2.0/static_cast<double>(r_x+1);
                                        double D_ppq = 0.0;
                                        double D_pqp = 0.0;
#if (USE_HO_FLUXES == 1)
                                        D_ppp += static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_x));
                                        D_pqp -= static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_x));
#if (NDIM == 3)
                                        D_ppp += static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_x));
                                        D_ppq -= static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_x));
#endif
#endif
                                        const double Dx = -(*D_data)(loc_index == 0 ? ixlower : ixupper)*xfactor;
                                        D_c   *= Dx;
                                        D_ppp *= Dx;
                                        D_pqp *= Dx;
                                        D_ppq *= Dx;

                                        int entry_c = (*entries_data)(p_idx)++;
                                        HYPRE_SStructMatrixSetValues(d_matrix,
                                                                     part, p_idx, var, 1,
                                                                     &entry_c, &D_c);

                                        int entry_ppp = NDIM;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_ppp, &D_ppp);

                                        int entry_pqp = (j>0 ? NDIM+1 : 0) + 1;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_pqp, &D_pqp);
#if (NDIM == 3)
                                        int entry_ppq = (k>0 ? NDIM+1 : 0) + 2;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_ppq, &D_ppq);
#endif
#if (NDIM == 3)
                                    }
                                }
#endif
                            }
                        }
                        break;
                    case 2:   // y_lo
                    case 3:   // y_hi
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(blower(1) == bupper(1));
#endif
                        p_idx(1) = blower(1) + (loc_index == 2 ? 1 : -1);
                        for (int c_i = coarsen(blower(0),r_x);
                             c_i <= coarsen(bupper(0),r_x); ++c_i)
                        {
                            for (int p_i = 0; p_i < r_x; ++p_i)
                            {
                                p_idx(0) = r_x*c_i + p_i;
                                const int i = (p_i%2 == 0 ? +1 : -1);
#if (NDIM == 3)
                                for (int c_k = coarsen(blower(2),r_z);
                                     c_k <= coarsen(bupper(2),r_z); ++c_k)
                                {
                                    for (int p_k = 0; p_k < r_z; ++p_k)
                                    {
                                        p_idx(2) = r_z*c_k + p_k;
                                        const int k = (p_k%2 == 0 ? +1 : -1);
#endif
                                        const SAMRAI::pdat::SideIndex<NDIM> iylower(
                                            p_idx, SAMRAI::pdat::SideIndex<NDIM>::Y, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                                        const SAMRAI::pdat::SideIndex<NDIM> iyupper(
                                            p_idx, SAMRAI::pdat::SideIndex<NDIM>::Y, SAMRAI::pdat::SideIndex<NDIM>::Upper);

                                        double D_c =  -2.0/static_cast<double>(r_y+1);
                                        double D_ppp = 2.0/static_cast<double>(r_y+1);
                                        double D_qpp = 0.0;
                                        double D_ppq = 0.0;
#if (USE_HO_FLUXES == 1)
                                        D_ppp += static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_y));
                                        D_qpp -= static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_y));
#if (NDIM == 3)
                                        D_ppp += static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_y));
                                        D_ppq -= static_cast<double>(2*p_k+1-r_z)/static_cast<double>(k*(1+r_y));
#endif
#endif
                                        const double Dy = -(*D_data)(loc_index == 2 ? iylower : iyupper)*yfactor;
                                        D_c   *= Dy;
                                        D_ppp *= Dy;
                                        D_qpp *= Dy;
                                        D_ppq *= Dy;

                                        int entry_c = (*entries_data)(p_idx)++;
                                        HYPRE_SStructMatrixSetValues(d_matrix,
                                                                     part, p_idx, var, 1,
                                                                     &entry_c, &D_c);

                                        int entry_ppp = NDIM;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_ppp, &D_ppp);

                                        int entry_qpp = (i>0 ? NDIM+1 : 0) + 0;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_qpp, &D_qpp);
#if (NDIM == 3)
                                        int entry_ppq = (k>0 ? NDIM+1 : 0) + 2;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_ppq, &D_ppq);
#endif
#if (NDIM == 3)
                                    }
                                }
#endif
                            }
                        }
                        break;
#if (NDIM == 3)
                    case 4:   // z_lo
                    case 5:   // z_hi
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(blower(2) == bupper(2));
#endif
                        p_idx(2) = blower(2) + (loc_index == 4 ? 1 : -1);
                        for (int c_i = coarsen(blower(0),r_x);
                             c_i <= coarsen(bupper(0),r_x); ++c_i)
                        {
                            for (int p_i = 0; p_i < r_x; ++p_i)
                            {
                                p_idx(0) = r_x*c_i + p_i;
                                const int i = (p_i%2 == 0 ? +1 : -1);
                                for (int c_j = coarsen(blower(1),r_y);
                                     c_j <= coarsen(bupper(1),r_y); ++c_j)
                                {
                                    for (int p_j = 0; p_j < r_y; ++p_j)
                                    {
                                        p_idx(1) = r_y*c_j + p_j;
                                        const int j = (p_j%2 == 0 ? +1 : -1);
                                        const SAMRAI::pdat::SideIndex<NDIM> izlower(
                                            p_idx, SAMRAI::pdat::SideIndex<NDIM>::Y, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                                        const SAMRAI::pdat::SideIndex<NDIM> izupper(
                                            p_idx, SAMRAI::pdat::SideIndex<NDIM>::Y, SAMRAI::pdat::SideIndex<NDIM>::Upper);

                                        double D_c =  -2.0/static_cast<double>(r_z+1);
                                        double D_ppp = 2.0/static_cast<double>(r_z+1);
                                        double D_qpp = 0.0;
                                        double D_pqp = 0.0;
#if (USE_HO_FLUXES == 1)
                                        D_ppp += static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_z));
                                        D_qpp -= static_cast<double>(2*p_i+1-r_x)/static_cast<double>(i*(1+r_z));

                                        D_ppp += static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_z));
                                        D_pqp -= static_cast<double>(2*p_j+1-r_y)/static_cast<double>(j*(1+r_z));
#endif
                                        const double Dz = -(*D_data)(loc_index == 4 ? izlower : izupper)*zfactor;
                                        D_c   *= Dz;
                                        D_ppp *= Dz;
                                        D_qpp *= Dz;
                                        D_pqp *= Dz;

                                        int entry_c = (*entries_data)(p_idx)++;
                                        HYPRE_SStructMatrixSetValues(d_matrix,
                                                                     part, p_idx, var, 1,
                                                                     &entry_c, &D_c);

                                        int entry_ppp = NDIM;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_ppp, &D_ppp);

                                        int entry_qpp = (i>0 ? NDIM+1 : 0) + 0;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_qpp, &D_qpp);

                                        int entry_pqp = (j>0 ? NDIM+1 : 0) + 1;
                                        HYPRE_SStructMatrixAddToValues(d_matrix,
                                                                       part, p_idx, var, 1,
                                                                       &entry_pqp, &D_pqp);
                                    }
                                }
                            }
                        }
                        break;
#endif
                    default:
                        TBOX_ERROR("this statement should not be reached!\n");
                }
            }
        }
    }
    HYPRE_SStructMatrixAssemble(d_matrix);

    // Now that the matrix is assembled, setup the hypre solvers.
    setupHypreSolver();

    t_set_matrix_coefficients->stop();
    return;
}// setMatrixCoefficients

int
CCPoissonHypreSolver::solveSystem(
    const int u,
    const int f)
{
    // Copy solution and right-hand-side to hypre structures.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int part = getLevelPart(ln,d_coarsest_ln,d_finest_ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
#if (NDIM == 2)
            const double scale = dx[0]*dx[1];
#endif
#if (NDIM == 3)
            const double scale = dx[0]*dx[1]*dx[2];
#endif
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > u_data = patch->getPatchData(u);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_data = patch->getPatchData(f);

            copyToHypre(d_sol_vec,        *u_data, d_sol_depth, patch_box, part);
            copyToHypre(d_rhs_vec, scale, *f_data, d_rhs_depth, patch_box, part);
        }
    }

    // Zero out the right-hand-side values in the coarsened-fine
    // regions.
    hypre_ZeroAMRVectorData(d_rhs_vec, d_plevels, d_rfactors);
    hypre_ZeroAMRVectorData(d_sol_vec, d_plevels, d_rfactors);

    // Assemble the hypre vectors.
    HYPRE_SStructVectorAssemble(d_sol_vec);
    HYPRE_SStructVectorAssemble(d_rhs_vec);

    if (d_print_solver_info)
    {
        HYPRE_SStructMatrixPrint("A0.out"  ,d_matrix ,1);
        HYPRE_SStructVectorPrint("sol0.out",d_sol_vec,1);
        HYPRE_SStructVectorPrint("rhs0.out",d_rhs_vec,1);
    }

    // Solve the system.
    t_solve_system_hypre->start();

    if (d_solver_type == "SysPFMG")
    {
        HYPRE_SStructSysPFMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructSysPFMGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructSysPFMGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructSysPFMGGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "FAC")
    {
        HYPRE_SStructFACSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructFACSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructFACSolve3(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructFACGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructFACGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_SStructPCGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructPCGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructPCGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructPCGGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructPCGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_SStructGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructGMRESGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_SStructBiCGSTABSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructBiCGSTABSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructBiCGSTABSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructBiCGSTABGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }

    t_solve_system_hypre->stop();

    if (d_print_solver_info)
    {
        HYPRE_SStructMatrixPrint("A.out"  ,d_matrix ,1);
        HYPRE_SStructVectorPrint("sol.out",d_sol_vec,1);
        HYPRE_SStructVectorPrint("rhs.out",d_rhs_vec,1);
    }

    // Pull the solution vector out of the hypre structures.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int part = getLevelPart(ln,d_coarsest_ln,d_finest_ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > u_data = patch->getPatchData(u);
            copyFromHypre(*u_data, d_sol_depth, d_sol_vec, patch_box, part);
        }
    }

    // Synchronize the solution across the patch hierarchy.
    if (d_synch_soln)
    {
        SAMRAI::xfer::CoarsenAlgorithm<NDIM> coarsener;
        coarsener.registerCoarsen(u, u, d_urestriction_coarsen_operator);
        for (int dst_ln = d_finest_ln-1; dst_ln >= d_coarsest_ln; --dst_ln)
        {
            coarsener.resetSchedule(d_urestriction_coarsen_schedules[dst_ln]);
            d_urestriction_coarsen_schedules[dst_ln]->coarsenData();
            d_urestriction_coarsen_algorithm->resetSchedule(d_urestriction_coarsen_schedules[dst_ln]);
        }
    }

    return d_current_residual_norm <= d_rel_residual_tol;
}// solveSystem

void
CCPoissonHypreSolver::deallocateHypreData()
{
    if (d_grid   ) HYPRE_SStructGridDestroy(d_grid);
    if (d_graph  ) HYPRE_SStructGraphDestroy(d_graph);
    if (d_stencil) HYPRE_SStructStencilDestroy(d_stencil);
    if (d_matrix ) HYPRE_SStructMatrixDestroy(d_matrix);
    if (d_sol_vec) HYPRE_SStructVectorDestroy(d_sol_vec);
    if (d_rhs_vec) HYPRE_SStructVectorDestroy(d_rhs_vec);

    d_grid    = NULL;
    d_graph   = NULL;
    d_stencil = NULL;
    d_matrix  = NULL;
    d_sol_vec = NULL;
    d_rhs_vec = NULL;

    if (d_solver)
    {
        if (d_solver_type == "SysPFMG")
        {
            HYPRE_SStructSysPFMGDestroy(d_solver);
        }
        else if (d_solver_type == "FAC")
        {
            HYPRE_SStructFACDestroy2(d_solver);
        }
        else if (d_solver_type == "PCG")
        {
            HYPRE_SStructPCGDestroy(d_solver);
        }
        else if (d_solver_type == "GMRES")
        {
            HYPRE_SStructGMRESDestroy(d_solver);
        }
        else if (d_solver_type == "BiCGSTAB")
        {
            HYPRE_SStructBiCGSTABDestroy(d_solver);
        }
    }

    if (d_precond)
    {
        if (d_solver_type == "PCG"     ||
            d_solver_type == "GMRES"   ||
            d_solver_type == "BiCGSTAB")
        {
            if (d_precond_type == "SysPFMG")
            {
                HYPRE_SStructSysPFMGDestroy(d_precond);
            }
            else if (d_precond_type == "FAC")
            {
                HYPRE_SStructFACDestroy2(d_precond);
            }
            else if  ((d_precond_type == "Split-SMG") ||
                      (d_precond_type == "Split-PFMG"))
            {
                HYPRE_SStructSplitDestroy(d_precond);
            }
        }
    }

    d_solver  = NULL;
    d_precond = NULL;

    if (d_plevels  != NULL) delete[] d_plevels;
    if (d_rfactors != NULL) delete[] d_rfactors;

    d_nparts   = 0;
    d_plevels  = NULL;
    d_rfactors = NULL;

    return;
}// deallocateHypreData

void
CCPoissonHypreSolver::copyToHypre(
    HYPRE_SStructVector vector,
    SAMRAI::pdat::CellData<NDIM,double>& src,
    const int depth,
    const SAMRAI::hier::Box<NDIM>& box,
    const int part)
{
    static const int var = VAR;

    const SAMRAI::hier::Index<NDIM>& box_lower = box.lower();
    const SAMRAI::hier::Index<NDIM>& box_upper = box.upper();
    SAMRAI::hier::Index<NDIM> ilower = box_lower;
    SAMRAI::hier::Index<NDIM> iupper = box_upper;

#if (NDIM == 3)
    for (int k = box_lower(2); k <= box_upper(2); ++k)
    {
        ilower(2) = k;
        iupper(2) = k;
#endif
        for (int j = box_lower(1); j <= box_upper(1); ++j)
        {
            ilower(1) = j;
            iupper(1) = j;

            HYPRE_SStructVectorSetBoxValues(vector, part, ilower, iupper, var, &src(ilower,depth));
        }
#if (NDIM == 3)
    }
#endif
    return;
}// copyToHypre

void
CCPoissonHypreSolver::copyToHypre(
    HYPRE_SStructVector vector,
    const double scale,
    SAMRAI::pdat::CellData<NDIM,double>& src,
    const int depth,
    const SAMRAI::hier::Box<NDIM>& box,
    const int part)
{
    static const int var = VAR;

    const SAMRAI::hier::Index<NDIM>& box_lower = box.lower();
    const SAMRAI::hier::Index<NDIM>& box_upper = box.upper();
    SAMRAI::hier::Index<NDIM> index;

#if (NDIM == 3)
    for (int k = box_lower(2); k <= box_upper(2); ++k)
    {
        index(2) = k;
#endif
        for (int j = box_lower(1); j <= box_upper(1); ++j)
        {
            index(1) = j;
            for (int i = box_lower(0); i <= box_upper(0); ++i)
            {
                index(0) = i;
                double value = scale*src(index,depth);
                HYPRE_SStructVectorSetValues(vector, part, index, var, &value);
            }
        }
#if (NDIM == 3)
    }
#endif
    return;
}// copyToHypre

void
CCPoissonHypreSolver::copyFromHypre(
    SAMRAI::pdat::CellData<NDIM,double>& dst,
    const int depth,
    HYPRE_SStructVector vector,
    const SAMRAI::hier::Box<NDIM>& box,
    const int part)
{
    static const int var = VAR;

    const SAMRAI::hier::Index<NDIM>& box_lower = box.lower();
    const SAMRAI::hier::Index<NDIM>& box_upper = box.upper();
    SAMRAI::hier::Index<NDIM> ilower = box_lower;
    SAMRAI::hier::Index<NDIM> iupper = box_upper;

#if (NDIM == 3)
    for (int k = box_lower(2); k <= box_upper(2); ++k)
    {
        ilower(2) = k;
        iupper(2) = k;
#endif
        for (int j = box_lower(1); j <= box_upper(1); ++j)
        {
            ilower(1) = j;
            iupper(1) = j;

            HYPRE_SStructVectorGetBoxValues(vector, part, ilower, iupper, var, &dst(ilower,depth));
        }
#if (NDIM == 3)
    }
#endif
    return;
}// copyFromHypre

void
CCPoissonHypreSolver::copyFromHypre(
    const double scale,
    SAMRAI::pdat::CellData<NDIM,double>& dst,
    const int depth,
    HYPRE_SStructVector vector,
    const SAMRAI::hier::Box<NDIM>& box,
    const int part)
{
    static const int var = VAR;

    const SAMRAI::hier::Index<NDIM>& box_lower = box.lower();
    const SAMRAI::hier::Index<NDIM>& box_upper = box.upper();
    SAMRAI::hier::Index<NDIM> index;

#if (NDIM == 3)
    for (int k = box_lower(2); k <= box_upper(2); ++k)
    {
        index(2) = k;
#endif
        for (int j = box_lower(1); j <= box_upper(1); ++j)
        {
            index(1) = j;
            for (int i = box_lower(0); i <= box_upper(0); ++i)
            {
                index(0) = i;
                double& value = dst(index,depth);
                HYPRE_SStructVectorGetValues(vector, part, index, var, &value);
                value *= scale;
            }
        }
#if (NDIM == 3)
    }
#endif
    return;
}// copyFromHypre

SAMRAI::hier::Box<NDIM>
CCPoissonHypreSolver::trimBoundaryBox(
    const SAMRAI::hier::Box<NDIM>& bdry_box,
    const SAMRAI::hier::Box<NDIM>& patch_box,
    const int loc_index)
{
    SAMRAI::hier::Box<NDIM> trimmed_box = bdry_box;

    switch (loc_index)
    {
        case 0:   // x_lo
        case 1:   // x_hi
            trimmed_box.lower(1) = SAMRAI::tbox::Utilities::imax(trimmed_box.lower(1),patch_box.lower(1));
            trimmed_box.upper(1) = SAMRAI::tbox::Utilities::imin(trimmed_box.upper(1),patch_box.upper(1));
#if (NDIM == 3)
            trimmed_box.lower(2) = SAMRAI::tbox::Utilities::imax(trimmed_box.lower(2),patch_box.lower(2));
            trimmed_box.upper(2) = SAMRAI::tbox::Utilities::imin(trimmed_box.upper(2),patch_box.upper(2));
#endif
            break;
        case 2:   // y_lo
        case 3:   // y_hi
            trimmed_box.lower(0) = SAMRAI::tbox::Utilities::imax(trimmed_box.lower(0),patch_box.lower(0));
            trimmed_box.upper(0) = SAMRAI::tbox::Utilities::imin(trimmed_box.upper(0),patch_box.upper(0));
#if (NDIM == 3)
            trimmed_box.lower(2) = SAMRAI::tbox::Utilities::imax(trimmed_box.lower(2),patch_box.lower(2));
            trimmed_box.upper(2) = SAMRAI::tbox::Utilities::imin(trimmed_box.upper(2),patch_box.upper(2));
#endif
            break;
#if (NDIM == 3)
        case 4:   // z_lo
        case 5:   // z_hi
            trimmed_box.lower(0) = SAMRAI::tbox::Utilities::imax(trimmed_box.lower(0),patch_box.lower(0));
            trimmed_box.upper(0) = SAMRAI::tbox::Utilities::imin(trimmed_box.upper(0),patch_box.upper(0));

            trimmed_box.lower(1) = SAMRAI::tbox::Utilities::imax(trimmed_box.lower(1),patch_box.lower(1));
            trimmed_box.upper(1) = SAMRAI::tbox::Utilities::imin(trimmed_box.upper(1),patch_box.upper(1));

            break;
#endif
        default:
            TBOX_ERROR("this statement should not be reached!\n");
    }
    return trimmed_box;
}// trimBoundaryBox

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBTK::CCPoissonHypreSolver>;

//////////////////////////////////////////////////////////////////////////////
