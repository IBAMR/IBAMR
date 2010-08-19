// Filename: CCHierarchyProjector.C
// Created on 18 Feb 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "CCHierarchyProjector.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScSAMRAIVectorReal.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_project_hierarchy;
static Pointer<Timer> t_initialize_level_data;
static Pointer<Timer> t_reset_hierarchy_configuration;
static Pointer<Timer> t_put_to_database;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of CCHierarchyProjector restart file data.
static const int CC_HIERARCHY_PROJECTOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCHierarchyProjector::CCHierarchyProjector(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_do_log(false),
      d_hierarchy(hierarchy),
      d_grid_geom(d_hierarchy->getGridGeometry()),
      d_hier_cc_data_ops(NULL),
      d_hier_math_ops(NULL),
      d_is_managing_hier_math_ops(false),
      d_wgt_var(NULL),
      d_wgt_idx(-1),
      d_volume(0.0),
      d_context(NULL),
      d_F_var(NULL),
      d_F_idx(-1),
      d_null_space_idxs(int(pow(2,NDIM)),-1),
      d_petsc_null_space_vecs(int(pow(2,NDIM)),static_cast<Vec>(PETSC_NULL)),
      d_W_var(NULL),
      d_W_idx(-1),
      d_sol_vec(NULL),
      d_rhs_vec(NULL),
      d_null_space_vecs(int(pow(2,NDIM)),Pointer<SAMRAIVectorReal<NDIM,double> >(NULL)),
      d_max_iterations(50),
      d_abs_residual_tol(1.0e-12),
      d_rel_residual_tol(1.0e-8),
      d_initial_guess_nonzero(true),
      d_poisson_solver(NULL),
      d_cc_div_grad_op(NULL),
      d_cc_div_grad_hypre_solver(NULL),
      d_sol_var(NULL),
      d_rhs_var(NULL),
      d_sol_idx(-1),
      d_rhs_idx(-1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
#endif

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Initialize Variables and contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    d_context = var_db->getContext(d_object_name+"::CONTEXT");

    d_sol_var = new CellVariable<NDIM,double>(d_object_name+"::sol",1);
    d_rhs_var = new CellVariable<NDIM,double>(d_object_name+"::rhs",1);
    d_F_var   = new CellVariable<NDIM,double>(d_object_name+"::F",1);
    d_W_var   = new CellVariable<NDIM,double>(d_object_name+"::W",NDIM);

    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> no_ghosts = 0;

    d_sol_idx = var_db->registerVariableAndContext( d_sol_var, d_context, cell_ghosts);
    d_rhs_idx = var_db->registerVariableAndContext( d_rhs_var, d_context, cell_ghosts);
    d_F_idx   = var_db->registerVariableAndContext(   d_F_var, d_context, cell_ghosts);
    d_W_idx   = var_db->registerVariableAndContext(   d_W_var, d_context, cell_ghosts);

    d_null_space_var = d_F_var;
    for (int k = 0; k < int(pow(2,NDIM)); ++k)
    {
        d_null_space_idxs[k] = var_db->registerClonedPatchDataIndex(d_null_space_var, d_F_idx);
    }

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var =
        new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy);

    // Initialize the hypre preconditioner.
    Pointer<Database> hypre_db;
    if (input_db->isDatabase("hypre_solver"))
    {
        hypre_db = input_db->getDatabase("hypre_solver");
    }
    d_cc_div_grad_hypre_solver = new CCDivGradHypreLevelSolver(d_object_name+"::hypre D*G solver", hypre_db);

    // Initialize the Poisson solver.
    d_cc_div_grad_op = new CCDivGradOperator(d_object_name+"::-D*G Operator");
    d_poisson_solver = new PETScKrylovLinearSolver(d_object_name+"::PETSc Krylov solver", "cc_proj_");
    d_poisson_solver->setMaxIterations(d_max_iterations);
    d_poisson_solver->setAbsoluteTolerance(d_abs_residual_tol);
    d_poisson_solver->setRelativeTolerance(d_rel_residual_tol);
    d_poisson_solver->setInitialGuessNonzero(d_initial_guess_nonzero);
    d_poisson_solver->setOperator(d_cc_div_grad_op);
    d_poisson_solver->setPreconditioner(d_cc_div_grad_hypre_solver);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_project_hierarchy = TimerManager::getManager()->
            getTimer("IBAMR::CCHierarchyProjector::projectHierarchy");
        t_initialize_level_data = TimerManager::getManager()->
            getTimer("IBAMR::CCHierarchyProjector::initializeLevelData()");
        t_reset_hierarchy_configuration = TimerManager::getManager()->
            getTimer("IBAMR::CCHierarchyProjector::resetHierarchyConfiguration()");
        t_put_to_database = TimerManager::getManager()->
            getTimer("IBAMR::CCHierarchyProjector::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// CCHierarchyProjector

CCHierarchyProjector::~CCHierarchyProjector()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    // Deallocate null space object.
    int ierr;
    KSP petsc_ksp = dynamic_cast<PETScKrylovLinearSolver*>(d_poisson_solver.getPointer())->getPETScKSP();
    MatNullSpace petsc_nullsp;
    ierr = KSPGetNullSpace(petsc_ksp, &petsc_nullsp); IBTK_CHKERRQ(ierr);
    if (petsc_nullsp != PETSC_NULL)
    {
        ierr = MatNullSpaceDestroy(petsc_nullsp); IBTK_CHKERRQ(ierr);
    }
    for (int k = 0; k < int(pow(2,NDIM)); ++k)
    {
        if (d_petsc_null_space_vecs[k] != PETSC_NULL)
        {
            ierr = VecDestroy(d_petsc_null_space_vecs[k]); IBTK_CHKERRQ(ierr);
        }
    }

    // Deallocate solver.
    d_poisson_solver->deallocateSolverState();
    d_poisson_solver.setNull();
    return;
}// ~CCHierarchyProjector

const std::string&
CCHierarchyProjector::getName() const
{
    return d_object_name;
}// getName

///
///  The following routines:
///
///      getHierarchyMathOps(),
///      setHierarchyMathOps(),
///      isManagingHierarchyMathOps()
///
///  allow for the sharing of a single HierarchyMathOps object between multiple
///  HierarchyIntegrator objects.
///

Pointer<HierarchyMathOps>
CCHierarchyProjector::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
CCHierarchyProjector::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops,
    const bool manage_ops)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hier_math_ops.isNull());
#endif
    d_hier_math_ops = hier_math_ops;
    d_is_managing_hier_math_ops = manage_ops;
    return;
}// setHierarchyMathOps

bool
CCHierarchyProjector::isManagingHierarchyMathOps() const
{
    return d_is_managing_hier_math_ops;
}// isManagingHierarchyMathOps

///
///  The following routines:
///
///      getPoissonSolver()
///
///  allow other objects to access the Poisson solver and related data used by
///  this integrator.
///

Pointer<KrylovLinearSolver>
CCHierarchyProjector::getPoissonSolver() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_poisson_solver.isNull());
#endif
    return d_poisson_solver;
}// getPoissonSolver

///
///  The following routines:
///
///      projectHierarchy()
///
///  provide the projection functionality.
///

void
CCHierarchyProjector::projectHierarchy(
    const double rho,
    const double dt,
    const double time,
    const int U_idx,
    const Pointer<CellVariable<NDIM,double> >& U_var,
    const int Phi_idx,
    const Pointer<CellVariable<NDIM,double> >& Phi_var,
    const int Grad_Phi_idx,
    const Pointer<CellVariable<NDIM,double> >& Grad_Phi_var,
    const int W_idx,
    const Pointer<CellVariable<NDIM,double> >& W_var,
    const int Q_idx,
    const Pointer<CellVariable<NDIM,double> >& Q_var)
{
    t_project_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_F_idx);
    scratch_idxs.setFlag(d_W_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, time);
    }

    // Fill the intermediate velocity data.
    d_hier_cc_data_ops->copyData(d_W_idx, W_idx);
    d_W_hier_bdry_fill_op->fillData(time);

    // Setup the linear operator.
    d_cc_div_grad_op->setHierarchyMathOps(d_hier_math_ops);

    // Compute F = (rho/dt)*(Q - div W).
    d_hier_math_ops->div(
        d_F_idx, d_F_var,   // dst
        -rho/dt,            // alpha
        d_W_idx, d_W_var,   // src1
        d_no_fill_op,       // src1_bdry_fill
        time,               // src1_bdry_fill_time
        +rho/dt,            // beta
        Q_idx, Q_var     ); // src2

    // Solve -div grad Phi = F = (rho/dt)*(Q - div W).
    SAMRAIVectorReal<NDIM,double> sol_vec(
        d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(Phi_var, Phi_idx, d_wgt_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM,double> rhs_vec(
        d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_F_var, d_F_idx, d_wgt_idx, d_hier_cc_data_ops);

    d_poisson_solver->solveSystem(sol_vec,rhs_vec);

    if (d_do_log) plog << "CCHierarchyProjector::projectHierarchy(): linear solve number of iterations = " << d_poisson_solver->getNumIterations() << "\n";
    if (d_do_log) plog << "CCHierarchyProjector::projectHierarchy(): linear solve residual norm        = " << d_poisson_solver->getResidualNorm()  << "\n";

    if (d_poisson_solver->getNumIterations() == d_poisson_solver->getMaxIterations())
    {
        pout << d_object_name << "::projectHierarchy():"
             <<"  WARNING: linear solver iterations == max iterations\n";
    }

    // Setup the interpolation transaction information.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_transaction_comp(Phi_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    d_Phi_hier_bdry_fill_op->resetTransactionComponent(Phi_transaction_comp);

    // Fill the physical boundary conditions for Phi.
    d_Phi_hier_bdry_fill_op->setHomogeneousBc(false);
    d_Phi_hier_bdry_fill_op->fillData(time);

    // Set U = W - (dt/rho)*grad Phi.
    d_hier_math_ops->grad(
        Grad_Phi_idx, Grad_Phi_var, // dst
        1.0,                        // alpha
        Phi_idx, Phi_var,           // src
        d_no_fill_op,               // src_bdry_fill
        time              );        // src_bdry_fill_time
    d_hier_cc_data_ops->axpy(U_idx, -dt/rho, Grad_Phi_idx, W_idx);

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }

    // Reset the transaction components.
    InterpolationTransactionComponent d_Phi_transaction_comp(d_F_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    d_Phi_hier_bdry_fill_op->resetTransactionComponent(d_Phi_transaction_comp);

    InterpolationTransactionComponent d_W_transaction_comp(d_W_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    d_W_hier_bdry_fill_op->resetTransactionComponent(d_W_transaction_comp);

    t_project_hierarchy->stop();
    return;
}// projectHierarchy

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration()
///
///  are concrete implementations of functions declared in the
///  mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
CCHierarchyProjector::initializeLevelData(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    // Setup the 2^NDIM dimensional null space of the cell-centered exact
    // projection operator (the checkerboard modes).
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    if (allocate_data)
    {
        for (int k = 0; k < int(pow(2,NDIM)); ++k)
        {
            level->allocatePatchData(d_null_space_idxs[k]);
        }
    }
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> coarsened_box = Box<NDIM>::coarsen(patch_box,2);
        IntVector<NDIM> chkbrd_mode_id;
#if (NDIM > 2)
        for (chkbrd_mode_id(2) = 0; chkbrd_mode_id(2) < 2; ++chkbrd_mode_id(2))
        {
#endif
            for (chkbrd_mode_id(1) = 0; chkbrd_mode_id(1) < 2; ++chkbrd_mode_id(1))
            {
                for (chkbrd_mode_id(0) = 0; chkbrd_mode_id(0) < 2; ++chkbrd_mode_id(0))
                {
                    Pointer<CellData<NDIM,double> > null_space_data =
                        patch->getPatchData(d_null_space_idxs[chkbrd_mode_id(0) + 2*(chkbrd_mode_id(1) + 2*(NDIM > 2 ? chkbrd_mode_id(2) : 0))]);
                    null_space_data->fillAll(0.0);
                    for (CellData<NDIM,double>::Iterator it(coarsened_box); it; it++)
                    {
                        const Index<NDIM>& i = (*it);
                        (*null_space_data)(i*2+chkbrd_mode_id) = 1.0;
                    }
                }
            }
#if (NDIM > 2)
        }
#endif
    }

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
CCHierarchyProjector::resetHierarchyConfiguration(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level <= finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    // Reset the Hierarchy math operations for the new configuration.
    if (d_is_managing_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }
    else if (d_hier_math_ops.isNull())
    {
        d_hier_math_ops = new HierarchyMathOps(d_object_name+"::HierarchyMathOps", hierarchy);
        d_is_managing_hier_math_ops = true;
    }

    // Get the cell weights data.
    d_wgt_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Get the volume of the physical domain.
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Reset the solution and rhs vectors.
    d_sol_vec = new SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::sol_vec", d_hierarchy, 0, finest_hier_level);
    d_sol_vec->addComponent(d_sol_var,d_sol_idx,d_wgt_idx,d_hier_cc_data_ops);

    d_rhs_vec = new SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::rhs_vec", d_hierarchy, 0, finest_hier_level);
    d_rhs_vec->addComponent(d_rhs_var,d_rhs_idx,d_wgt_idx,d_hier_cc_data_ops);

    // (Re)-initialize the Poisson solver.
    d_cc_div_grad_op->setHierarchyMathOps(d_hier_math_ops);
    d_poisson_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);

    // Setup the nullspace object associated with the Poisson solver.
    int ierr;
    for (int k = 0; k < int(pow(2,NDIM)); ++k)
    {
        std::ostringstream stream;
        stream << k;
        d_null_space_vecs[k] = new SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::null_space_vec_"+stream.str(), d_hierarchy, 0, finest_hier_level);
        d_null_space_vecs[k]->addComponent(d_null_space_var,d_null_space_idxs[k],d_wgt_idx,d_hier_cc_data_ops);

        if (d_petsc_null_space_vecs[k] != PETSC_NULL)
        {
            ierr = VecDestroy(d_petsc_null_space_vecs[k]); IBTK_CHKERRQ(ierr);
        }
        d_petsc_null_space_vecs[k] = PETScSAMRAIVectorReal<double>::createPETScVector(d_null_space_vecs[k], PETSC_COMM_WORLD);

        double v_dot_v;
        ierr = VecDot(d_petsc_null_space_vecs[k], d_petsc_null_space_vecs[k], &v_dot_v); IBTK_CHKERRQ(ierr);
        ierr = VecScale(d_petsc_null_space_vecs[k], 1.0/v_dot_v); IBTK_CHKERRQ(ierr);
    }

    KSP petsc_ksp = dynamic_cast<PETScKrylovLinearSolver*>(d_poisson_solver.getPointer())->getPETScKSP();
    MatNullSpace petsc_nullsp;
    ierr = KSPGetNullSpace(petsc_ksp, &petsc_nullsp); IBTK_CHKERRQ(ierr);
    if (petsc_nullsp != PETSC_NULL)
    {
        ierr = MatNullSpaceDestroy(petsc_nullsp); IBTK_CHKERRQ(ierr);
    }
    static const PetscTruth has_cnst = PETSC_FALSE;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, has_cnst, int(pow(2,NDIM)), &d_petsc_null_space_vecs[0], &petsc_nullsp); IBTK_CHKERRQ(ierr);
    ierr = KSPSetNullSpace(petsc_ksp, petsc_nullsp); IBTK_CHKERRQ(ierr);

    // Initialize the interpolation operators.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;

    InterpolationTransactionComponent Phi_transaction_comp(d_F_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    d_Phi_hier_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_Phi_hier_bdry_fill_op->initializeOperatorState(Phi_transaction_comp, d_hierarchy);

    InterpolationTransactionComponent W_transaction_comp(d_W_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY);
    d_W_hier_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_W_hier_bdry_fill_op->initializeOperatorState(W_transaction_comp, d_hierarchy);

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  Serializable abstract base class.
///

void
CCHierarchyProjector::putToDatabase(
    Pointer<Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("CC_HIERARCHY_PROJECTOR_VERSION", CC_HIERARCHY_PROJECTOR_VERSION);

    t_put_to_database->stop();
    return;
}// putToDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CCHierarchyProjector::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    (void) is_from_restart;
    d_max_iterations = db->getIntegerWithDefault("max_iterations", d_max_iterations);
    d_abs_residual_tol = db->getDoubleWithDefault("abs_residual_tol", d_abs_residual_tol);
    d_rel_residual_tol = db->getDoubleWithDefault("rel_residual_tol", d_rel_residual_tol);
    d_initial_guess_nonzero = db->getBoolWithDefault("initial_guess_nonzero", d_initial_guess_nonzero);
    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);
    return;
}// getFromInput

void
CCHierarchyProjector::getFromRestart()
{
    Pointer<Database> restart_db =
        RestartManager::getManager()->getRootDatabase();

    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  \n"
                   << "Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("CC_HIERARCHY_PROJECTOR_VERSION");
    if (ver != CC_HIERARCHY_PROJECTOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::CCHierarchyProjector>;

//////////////////////////////////////////////////////////////////////////////
