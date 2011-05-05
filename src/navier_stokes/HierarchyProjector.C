// Filename: HierarchyProjector.C
// Created on 30 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "HierarchyProjector.h"

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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/PETScKrylovLinearSolver.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_project_hierarchy;
static Timer* t_initialize_level_data;
static Timer* t_reset_hierarchy_configuration;
static Timer* t_put_to_database;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int FACEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of HierarchyProjector restart file data.
static const int HIERARCHY_PROJECTOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

HierarchyProjector::HierarchyProjector(
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
      d_hier_fc_data_ops(NULL),
      d_hier_sc_data_ops(NULL),
      d_hier_math_ops(NULL),
      d_is_managing_hier_math_ops(false),
      d_wgt_var(NULL),
      d_wgt_idx(-1),
      d_volume(0.0),
      d_context(NULL),
      d_F_var(NULL),
      d_F_idx(-1),
      d_P_var(NULL),
      d_P_idx(-1),
      d_w_fc_var(NULL),
      d_w_fc_idx(-1),
      d_w_sc_var(NULL),
      d_w_sc_idx(-1),
      d_sol_vec(NULL),
      d_rhs_vec(NULL),
      d_max_iterations(50),
      d_abs_residual_tol(1.0e-12),
      d_rel_residual_tol(1.0e-8),
      d_initial_guess_nonzero(true),
      d_poisson_spec(d_object_name+"::Poisson spec"),
      d_u_bc_coefs(NDIM,static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_default_u_bc_coefs(),
      d_P_bc_coef(NULL),
      d_default_P_bc_coef(),
      d_Phi_bc_coef(NULL),
      d_poisson_solver(NULL),
      d_laplace_op(NULL),
      d_poisson_fac_op(NULL),
      d_poisson_fac_pc(NULL),
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

    d_sol_var  = new CellVariable<NDIM,double>(d_object_name+"::sol",1);
    d_rhs_var  = new CellVariable<NDIM,double>(d_object_name+"::rhs",1);
    d_F_var    = new CellVariable<NDIM,double>(d_object_name+"::F",1);
    d_P_var    = new CellVariable<NDIM,double>(d_object_name+"::P",1);
    d_w_fc_var = new FaceVariable<NDIM,double>(d_object_name+"::w_fc",1);
    d_w_sc_var = new SideVariable<NDIM,double>(d_object_name+"::w_sc",1);

    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> face_ghosts = FACEG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> no_ghosts = 0;

    d_sol_idx  = var_db->registerVariableAndContext( d_sol_var, d_context, cell_ghosts);
    d_rhs_idx  = var_db->registerVariableAndContext( d_rhs_var, d_context, cell_ghosts);
    d_F_idx    = var_db->registerVariableAndContext(   d_F_var, d_context, cell_ghosts);
    d_P_idx    = var_db->registerVariableAndContext(   d_P_var, d_context, cell_ghosts);
    d_w_fc_idx = var_db->registerVariableAndContext(d_w_fc_var, d_context, face_ghosts);
    d_w_sc_idx = var_db->registerVariableAndContext(d_w_sc_var, d_context, side_ghosts);

    // Initialize communications algorithms.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    Pointer<RefineOperator<NDIM> > refine_operator;

    refine_operator = grid_geom->lookupRefineOperator(d_w_fc_var, "CONSERVATIVE_LINEAR_REFINE");
    d_fc_velocity_ralg = new RefineAlgorithm<NDIM>();
    d_fc_velocity_ralg->registerRefine(d_w_fc_idx,  // destination
                                       d_w_fc_idx,  // source
                                       d_w_fc_idx,  // temporary work space
                                       refine_operator);
    d_fc_velocity_rstrategy = new CartExtrapPhysBdryOp(d_w_fc_idx, "LINEAR");

    refine_operator = grid_geom->lookupRefineOperator(d_w_sc_var, "CONSERVATIVE_LINEAR_REFINE");
    d_sc_velocity_ralg = new RefineAlgorithm<NDIM>();
    d_sc_velocity_ralg->registerRefine(d_w_sc_idx,  // destination
                                       d_w_sc_idx,  // source
                                       d_w_sc_idx,  // temporary work space
                                       refine_operator);
    d_sc_velocity_rstrategy = new CartExtrapPhysBdryOp(d_w_sc_idx, "LINEAR");

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();

    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    Pointer<FaceVariable<NDIM,double> > fc_var = new FaceVariable<NDIM,double>("fc_var");
    Pointer<SideVariable<NDIM,double> > sc_var = new SideVariable<NDIM,double>("sc_var");

    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy, true);
    d_hier_fc_data_ops = hier_ops_manager->getOperationsDouble(fc_var, hierarchy, true);
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, hierarchy, true);

    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup default velocity boundary condition objects to specify homogeneous
    // Dirichlet boundary conditions for all components of the velocity.
    for (int d = 0; d < NDIM; ++d)
    {
        std::ostringstream stream;
        stream << d;
        LocationIndexRobinBcCoefs<NDIM>* u_bc_coef = new LocationIndexRobinBcCoefs<NDIM>(
            d_object_name + "::default_u_bc_coef_" + stream.str(),
            Pointer<Database>(NULL));

        for (int location_index = 0; location_index < 2*NDIM; ++location_index)
        {
            u_bc_coef->setBoundaryValue(location_index,0.0);
        }

        d_default_u_bc_coefs.push_back(u_bc_coef);
    }

    // Setup a default pressure boundary condition object to specify homogeneous
    // Neumann boundary conditions for the pressure.
    LocationIndexRobinBcCoefs<NDIM>* P_bc_coef = new LocationIndexRobinBcCoefs<NDIM>(
        d_object_name + "::default_P_bc_coef", Pointer<Database>(NULL));

    for (int location_index = 0; location_index < 2*NDIM; ++location_index)
    {
        P_bc_coef->setBoundarySlope(location_index,0.0);
    }

    d_default_P_bc_coef = P_bc_coef;

    // Initialize the boundary conditions objects.
    d_Phi_bc_coef = new INSProjectionBcCoef(d_P_idx,d_default_P_bc_coef,PRESSURE_UPDATE,d_w_fc_idx,d_default_u_bc_coefs,true);
    setPressurePhysicalBcCoef(d_default_P_bc_coef);
    setVelocityPhysicalBcCoefs(d_default_u_bc_coefs);

    // Get initialization data for the FAC ops and FAC preconditioners and
    // initialize them.
    Pointer<Database> fac_op_db, fac_pc_db;

    if (input_db->keyExists("FACOp"))
    {
        fac_op_db = input_db->getDatabase("FACOp");
    }
    if (input_db->keyExists("FACPreconditioner"))
    {
        fac_pc_db = input_db->getDatabase("FACPreconditioner");
    }

    d_poisson_fac_op = new CCPoissonFACOperator(d_object_name+"::FAC Op", fac_op_db);
    d_poisson_fac_op->setPoissonSpecifications(d_poisson_spec);
    d_poisson_fac_op->setPhysicalBcCoef(d_Phi_bc_coef);
    d_poisson_fac_pc = new IBTK::FACPreconditioner(d_object_name+"::FAC Preconditioner", *d_poisson_fac_op, fac_pc_db);

    // Initialize the Poisson solver.
    static const bool homogeneous_bc = false;
    d_laplace_op = new CCLaplaceOperator(
        d_object_name+"::Laplace Operator",
        d_poisson_spec, d_Phi_bc_coef, homogeneous_bc);

    d_poisson_solver = new PETScKrylovLinearSolver(d_object_name+"::PETSc Krylov solver", "proj_");
    d_poisson_solver->setMaxIterations(d_max_iterations);
    d_poisson_solver->setAbsoluteTolerance(d_abs_residual_tol);
    d_poisson_solver->setRelativeTolerance(d_rel_residual_tol);
    d_poisson_solver->setInitialGuessNonzero(d_initial_guess_nonzero);
    d_poisson_solver->setOperator(d_laplace_op);
    d_poisson_solver->setPreconditioner(d_poisson_fac_pc);

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_project_hierarchy             = TimerManager::getManager()->getTimer("IBAMR::HierarchyProjector::projectHierarchy");
        t_initialize_level_data         = TimerManager::getManager()->getTimer("IBAMR::HierarchyProjector::initializeLevelData()");
        t_reset_hierarchy_configuration = TimerManager::getManager()->getTimer("IBAMR::HierarchyProjector::resetHierarchyConfiguration()");
        t_put_to_database               = TimerManager::getManager()->getTimer("IBAMR::HierarchyProjector::putToDatabase()");
                  );
    return;
}// HierarchyProjector

HierarchyProjector::~HierarchyProjector()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    // Deallocate solver.
    d_poisson_solver->deallocateSolverState();
    d_poisson_solver.setNull();

    // Deallocate other components.
    delete d_Phi_bc_coef;
    delete d_default_P_bc_coef;
    for (int d = 0; d < NDIM; ++d)
    {
        delete d_default_u_bc_coefs[d];
    }
    d_default_u_bc_coefs.clear();
    return;
}// ~HierarchyProjector

const std::string&
HierarchyProjector::getName() const
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
HierarchyProjector::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
HierarchyProjector::setHierarchyMathOps(
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
HierarchyProjector::isManagingHierarchyMathOps() const
{
    return d_is_managing_hier_math_ops;
}// isManagingHierarchyMathOps

///
///  The following routines:
///
///      setVelocityPhysicalBcCoefs(),
///      getVelocityPhysicalBcCoefs(),
///      setPressurePhysicalBcCoef(),
///      getPressurePhysicalBcCoef(),
///      getPoissonSolver()
///
///  allow other objects to access the Poisson solver and related data used by
///  this integrator.
///

void
HierarchyProjector::setVelocityPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR(d_object_name << "::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condition objects must be provided." << std::endl);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < u_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(u_bc_coefs[l] != NULL);
    }
#endif
    d_u_bc_coefs = u_bc_coefs;
    d_Phi_bc_coef->setVelocityPhysicalBcCoefs(d_u_bc_coefs);
    return;
}// setVelocityPhysicalBcCoefs

const std::vector<RobinBcCoefStrategy<NDIM>*>&
HierarchyProjector::getVelocityPhysicalBcCoefs() const
{
    return d_u_bc_coefs;
}// getVelocityPhysicalBcCoefs

void
HierarchyProjector::setPressurePhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* const P_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(P_bc_coef != NULL);
#endif
    d_P_bc_coef = P_bc_coef;
    d_Phi_bc_coef->setPressurePhysicalBcCoef(d_P_bc_coef);
    return;
}// setPressurePhysicalBcCoef

RobinBcCoefStrategy<NDIM>*
HierarchyProjector::getPressurePhysicalBcCoef() const
{
    return d_P_bc_coef;
}// getPressurePhysicalBcCoef

Pointer<KrylovLinearSolver>
HierarchyProjector::getPoissonSolver() const
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
HierarchyProjector::projectHierarchy(
    const double rho,
    const double dt,
    const double time,
    const ProjectionMethodType& projection_type,
    const int u_idx,
    const Pointer<FaceVariable<NDIM,double> >& u_var,
    const int P_idx,
    const Pointer<CellVariable<NDIM,double> >& P_var,
    const int Phi_idx,
    const Pointer<CellVariable<NDIM,double> >& Phi_var,
    const int grad_Phi_idx,
    const Pointer<FaceVariable<NDIM,double> >& grad_Phi_var,
    const int w_idx,
    const Pointer<FaceVariable<NDIM,double> >& w_var,
    const int Q_idx,
    const Pointer<CellVariable<NDIM,double> >& Q_var)
{
    t_project_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_F_idx);
    scratch_idxs.setFlag(d_P_idx);
    scratch_idxs.setFlag(d_w_fc_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, time);
    }

    // Fill the pressure data if we are using a pressure-increment projection.
    if (projection_type == PRESSURE_INCREMENT)
    {
        d_hier_cc_data_ops->copyData(d_P_idx, P_idx);
        d_P_hier_bdry_fill_op->fillData(time);
    }

    // Fill the intermediate velocity data.
    d_hier_fc_data_ops->copyData(d_w_fc_idx, w_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_fc_velocity_rscheds[ln]->fillData(time);
    }

    // Setup the boundary coefficient specification object.
    //
    // NOTE: These boundary coefficients are also used by the linear operator
    // and by the FAC preconditioner objects associated with this class.
    d_Phi_bc_coef->setProblemCoefs(rho, dt);
    d_Phi_bc_coef->setCurrentPressurePatchDataIndex(d_P_idx);
    d_Phi_bc_coef->setPressurePhysicalBcCoef(d_P_bc_coef);
    d_Phi_bc_coef->setProjectionType(projection_type);
    d_Phi_bc_coef->setIntermediateVelocityPatchDataIndex(d_w_fc_idx);
    d_Phi_bc_coef->setVelocityPhysicalBcCoefs(d_u_bc_coefs);

    // Setup the linear operator.
    d_laplace_op->setTime(time);
    d_laplace_op->setHierarchyMathOps(d_hier_math_ops);

    // Setup the preconditioner.
    d_poisson_fac_op->setTime(time);

    // Compute F = (rho/dt)*(Q - div w).
    const bool w_cf_bdry_synch = true;
    d_hier_math_ops->div(
        d_F_idx, d_F_var,   // dst
        -rho/dt,            // alpha
        w_idx, w_var,       // src1
        d_no_fill_op,       // src1_bdry_fill
        time,               // src1_bdry_fill_time
        w_cf_bdry_synch,    // src1_cf_bdry_synch
        +rho/dt,            // beta
        Q_idx, Q_var     ); // src2

    // Solve -div grad Phi = F = (rho/dt)*(Q - div w).
    SAMRAIVectorReal<NDIM,double> sol_vec(
        d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(Phi_var, Phi_idx, d_wgt_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM,double> rhs_vec(
        d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_F_var, d_F_idx, d_wgt_idx, d_hier_cc_data_ops);

    d_poisson_solver->solveSystem(sol_vec,rhs_vec);

    if (d_do_log) plog << "HierarchyProjector::projectHierarchy(): linear solve number of iterations = " << d_poisson_solver->getNumIterations() << "\n";
    if (d_do_log) plog << "HierarchyProjector::projectHierarchy(): linear solve residual norm        = " << d_poisson_solver->getResidualNorm()  << "\n";

    if (d_poisson_solver->getNumIterations() == d_poisson_solver->getMaxIterations())
    {
        pout << d_object_name << "::projectHierarchy():"
             <<"  WARNING: linear solver iterations == max iterations\n";
    }

    // Setup the interpolation transaction information.
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(CELLG, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_transaction_comp(Phi_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_Phi_bc_coef, cc_fill_pattern);
    d_Phi_hier_bdry_fill_op->resetTransactionComponent(Phi_transaction_comp);

    // Fill the physical boundary conditions for Phi.
    d_Phi_bc_coef->setHomogeneousBc(false);
    d_Phi_hier_bdry_fill_op->setHomogeneousBc(false);
    d_Phi_hier_bdry_fill_op->fillData(time);
    d_Phi_bc_coef->setHomogeneousBc(true);

    // Set u = w - (dt/rho)*grad Phi.
    const bool grad_Phi_cf_bdry_synch = true;
    d_hier_math_ops->grad(
        grad_Phi_idx, grad_Phi_var, // dst
        grad_Phi_cf_bdry_synch,     // dst_cf_bdry_synch
        1.0,                        // alpha
        Phi_idx, Phi_var,           // src
        d_no_fill_op,               // src_bdry_fill
        time              );        // src_bdry_fill_time
    d_hier_fc_data_ops->axpy(u_idx, -dt/rho, grad_Phi_idx, w_idx);

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }

    // Reset the transaction components.
    InterpolationTransactionComponent d_P_transaction_comp(d_P_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL, cc_fill_pattern);
    d_P_hier_bdry_fill_op->resetTransactionComponent(d_P_transaction_comp);

    InterpolationTransactionComponent d_Phi_transaction_comp(d_F_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_Phi_bc_coef, cc_fill_pattern);
    d_Phi_hier_bdry_fill_op->resetTransactionComponent(d_Phi_transaction_comp);

    t_project_hierarchy->stop();
    return;
}// projectHierarchy

void
HierarchyProjector::projectHierarchy(
    const double rho,
    const double dt,
    const double time,
    const ProjectionMethodType& projection_type,
    const int u_idx,
    const Pointer<SideVariable<NDIM,double> >& u_var,
    const int P_idx,
    const Pointer<CellVariable<NDIM,double> >& P_var,
    const int Phi_idx,
    const Pointer<CellVariable<NDIM,double> >& Phi_var,
    const int grad_Phi_idx,
    const Pointer<SideVariable<NDIM,double> >& grad_Phi_var,
    const int w_idx,
    const Pointer<SideVariable<NDIM,double> >& w_var,
    const int Q_idx,
    const Pointer<CellVariable<NDIM,double> >& Q_var)
{
    t_project_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_F_idx);
    scratch_idxs.setFlag(d_P_idx);
    scratch_idxs.setFlag(d_w_sc_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, time);
    }

    // Fill the pressure data if we are using a pressure-increment projection.
    if (projection_type == PRESSURE_INCREMENT)
    {
        d_hier_cc_data_ops->copyData(d_P_idx, P_idx);
        d_P_hier_bdry_fill_op->fillData(time);
    }

    // Fill the intermediate velocity data.
    d_hier_sc_data_ops->copyData(d_w_sc_idx, w_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_sc_velocity_rscheds[ln]->fillData(time);
    }

    // Setup the boundary coefficient specification object.
    //
    // NOTE: These boundary coefficients are also used by the linear operator
    // and by the FAC preconditioner objects associated with this class.
    d_Phi_bc_coef->setProblemCoefs(rho, dt);
    d_Phi_bc_coef->setCurrentPressurePatchDataIndex(d_P_idx);
    d_Phi_bc_coef->setPressurePhysicalBcCoef(d_P_bc_coef);
    d_Phi_bc_coef->setProjectionType(projection_type);
    d_Phi_bc_coef->setIntermediateVelocityPatchDataIndex(d_w_sc_idx);
    d_Phi_bc_coef->setVelocityPhysicalBcCoefs(d_u_bc_coefs);

    // Setup the linear operator.
    d_laplace_op->setTime(time);
    d_laplace_op->setHierarchyMathOps(d_hier_math_ops);

    // Setup the preconditioner.
    d_poisson_fac_op->setTime(time);

    // Compute F = (rho/dt)*(Q - div w).
    const bool w_cf_bdry_synch = true;
    d_hier_math_ops->div(
        d_F_idx, d_F_var,   // dst
        -rho/dt,            // alpha
        w_idx, w_var,       // src1
        d_no_fill_op,       // src1_bdry_fill
        time,               // src1_bdry_fill_time
        w_cf_bdry_synch,    // src1_cf_bdry_synch
        +rho/dt,            // beta
        Q_idx, Q_var     ); // src2

    // Solve -div grad Phi = F = (rho/dt)*(Q - div w).
    SAMRAIVectorReal<NDIM,double> sol_vec(
        d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(Phi_var, Phi_idx, d_wgt_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM,double> rhs_vec(
        d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_F_var, d_F_idx, d_wgt_idx, d_hier_cc_data_ops);

    d_poisson_solver->solveSystem(sol_vec,rhs_vec);

    if (d_do_log) plog << "HierarchyProjector::projectHierarchy(): linear solve number of iterations = " << d_poisson_solver->getNumIterations() << "\n";
    if (d_do_log) plog << "HierarchyProjector::projectHierarchy(): linear solve residual norm        = " << d_poisson_solver->getResidualNorm()  << "\n";

    if (d_poisson_solver->getNumIterations() == d_poisson_solver->getMaxIterations())
    {
        pout << d_object_name << "::projectHierarchy():"
             <<"  WARNING: linear solver iterations == max iterations\n";
    }

    // Setup the interpolation transaction information.
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(CELLG, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_transaction_comp(Phi_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_Phi_bc_coef, cc_fill_pattern);
    d_Phi_hier_bdry_fill_op->resetTransactionComponent(Phi_transaction_comp);

    // Fill the physical boundary conditions for Phi.
    d_Phi_bc_coef->setHomogeneousBc(false);
    d_Phi_hier_bdry_fill_op->setHomogeneousBc(false);
    d_Phi_hier_bdry_fill_op->fillData(time);
    d_Phi_bc_coef->setHomogeneousBc(true);

    // Set u = w - (dt/rho)*grad Phi.
    const bool grad_Phi_cf_bdry_synch = true;
    d_hier_math_ops->grad(
        grad_Phi_idx, grad_Phi_var, // dst
        grad_Phi_cf_bdry_synch,     // dst_cf_bdry_synch
        1.0,                        // alpha
        Phi_idx, Phi_var,           // src
        d_no_fill_op,               // src_bdry_fill
        time              );        // src_bdry_fill_time
    d_hier_sc_data_ops->axpy(u_idx, -dt/rho, grad_Phi_idx, w_idx);

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }

    // Reset the transaction components.
    InterpolationTransactionComponent d_P_transaction_comp(d_P_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL, cc_fill_pattern);
    d_P_hier_bdry_fill_op->resetTransactionComponent(d_P_transaction_comp);

    InterpolationTransactionComponent d_Phi_transaction_comp(d_F_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_Phi_bc_coef, cc_fill_pattern);
    d_Phi_hier_bdry_fill_op->resetTransactionComponent(d_Phi_transaction_comp);

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
HierarchyProjector::initializeLevelData(
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

    // intentionally blank

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
HierarchyProjector::resetHierarchyConfiguration(
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

    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

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
    d_laplace_op->setHierarchyMathOps(d_hier_math_ops);
    d_poisson_fac_op->setResetLevels(coarsest_level, finest_level);
    d_poisson_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);

    // Initialize the interpolation operators.
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(CELLG, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;

    InterpolationTransactionComponent P_transaction_comp(d_P_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL, cc_fill_pattern);
    d_P_hier_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_P_hier_bdry_fill_op->initializeOperatorState(P_transaction_comp, d_hierarchy);

    InterpolationTransactionComponent Phi_transaction_comp(d_F_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_Phi_bc_coef, cc_fill_pattern);
    d_Phi_hier_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_Phi_hier_bdry_fill_op->initializeOperatorState(Phi_transaction_comp, d_hierarchy);

    // (Re)build refine communication schedules.  These are created for all
    // levels in the hierarchy.
    d_fc_velocity_rscheds.resize(finest_hier_level+1);
    d_sc_velocity_rscheds.resize(finest_hier_level+1);
    for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        d_fc_velocity_rscheds[ln] = d_fc_velocity_ralg->createSchedule(
            level, ln-1, hierarchy, d_fc_velocity_rstrategy);
        d_sc_velocity_rscheds[ln] = d_sc_velocity_ralg->createSchedule(
            level, ln-1, hierarchy, d_sc_velocity_rstrategy);
    }

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
HierarchyProjector::putToDatabase(
    Pointer<Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("HIERARCHY_PROJECTOR_VERSION", HIERARCHY_PROJECTOR_VERSION);

    t_put_to_database->stop();
    return;
}// putToDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

void
HierarchyProjector::getFromInput(
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
HierarchyProjector::getFromRestart()
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

    int ver = db->getInteger("HIERARCHY_PROJECTOR_VERSION");
    if (ver != HIERARCHY_PROJECTOR_VERSION)
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
template class Pointer<IBAMR::HierarchyProjector>;

//////////////////////////////////////////////////////////////////////////////
