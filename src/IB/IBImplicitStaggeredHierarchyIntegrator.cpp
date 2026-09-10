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

#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#include <ibamr/IBImplicitStrategy.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/ibamr_enums.h>
#include <ibamr/private/StaggeredStokesIBTimeSteppingUtilities.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBKernelTensorProduct.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/KrylovLinearSolver.h>
#include <ibtk/PETScNewtonKrylovSolver.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Database.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>
#include <tbox/Utilities.h>

#include <petscmat.h>

#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchyDataOpsReal.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>

#include <cmath>
#include <string>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
constexpr int IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;

} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitStaggeredHierarchyIntegrator::IBImplicitStaggeredHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<IBImplicitStrategy> ib_method_ops,
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    bool register_for_restart)
    : IBHierarchyIntegrator(object_name, input_db, ib_method_ops, ins_hier_integrator, register_for_restart),
      d_ib_implicit_ops(ib_method_ops)
{
    d_ib_implicit_ops->setUseFixedLEOperators(true);

    d_use_structure_predictor = false;

    Pointer<Database> stokes_ib_pc_db = nullptr;
    if (input_db)
    {
        if (input_db->keyExists("use_structure_predictor"))
        {
            d_use_structure_predictor = input_db->getBool("use_structure_predictor");
        }
        if (input_db->keyExists("jacobian_delta_fcn"))
        {
            d_jac_kernel = IBKernelTensorProduct(input_db->getString("jacobian_delta_fcn"));
        }
        if (input_db->isDatabase("stokes_ib_precond_db"))
        {
            stokes_ib_pc_db = input_db->getDatabase("stokes_ib_precond_db");
        }
    }

    if (stokes_ib_pc_db)
    {
        d_has_velocity_nullspace = stokes_ib_pc_db->getBoolWithDefault("has_velocity_nullspace", false);
        d_has_pressure_nullspace = stokes_ib_pc_db->getBoolWithDefault("has_pressure_nullspace", true);
    }

    if (d_use_structure_predictor)
    {
        pout << "WARNING: explicit predictor for the structural configuration "
                "appears to be nonlinearly unstable!\n";
    }

    d_stokes_op = new StaggeredStokesOperator(d_object_name + "::stokes_op", false);
    d_ib_op = new StaggeredStokesIBOperator(d_object_name + "::ib_op", false);
    d_ib_jac_op = new StaggeredStokesIBJacobianOperator(d_object_name + "::ib_jac_op");
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> ib_fac_op = new StaggeredStokesIBLevelRelaxationFACOperator(
        d_object_name + "::ib_fac_op", stokes_ib_pc_db, "stokes_ib_pc_");
    d_ib_jac_pc = new StaggeredStokesIBJacobianFACPreconditioner(
        d_object_name + "::ib_fac_pc", ib_fac_op, stokes_ib_pc_db, "stokes_ib_pc_");
    d_stokes_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    d_ib_solver = new IBTK::PETScNewtonKrylovSolver(d_object_name + "::ib_solver", input_db, "ib_");
    d_ib_solver->setOperator(d_ib_op);
    d_ib_solver->setJacobian(d_ib_jac_op);
    d_ib_solver->getLinearSolver()->setPreconditioner(d_ib_jac_pc);

    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    return;
} // IBImplicitStaggeredHierarchyIntegrator

IBImplicitStaggeredHierarchyIntegrator::~IBImplicitStaggeredHierarchyIntegrator()
{
    deallocateOperatorsAndSolvers();
    if (d_rhs_vec)
    {
        free_vector_components(*d_rhs_vec);
    }
    return;
} // ~IBImplicitStaggeredHierarchyIntegrator

void
IBImplicitStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                     const double new_time,
                                                                     const int num_cycles)
{
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
    case TRAPEZOIDAL_RULE:
    case MIDPOINT_RULE:
        break;
    default:
        TBOX_ERROR(
            "IBImplicitStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy()"
            ": time_stepping_type = "
            << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
            << "  only supported time_stepping_types are:\n"
            << "    " << enum_to_string<TimeSteppingType>(BACKWARD_EULER) << "\n"
            << "    " << enum_to_string<TimeSteppingType>(TRAPEZOIDAL_RULE) << "\n"
            << "    " << enum_to_string<TimeSteppingType>(MIDPOINT_RULE) << "\n");
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    Pointer<PatchLevel<NDIM>> finest_level = d_hierarchy->getPatchLevel(finest_ln);
    if (!finest_level->checkAllocated(d_u_dof_index_idx))
    {
        finest_level->allocatePatchData(d_u_dof_index_idx, current_time);
    }
    if (!finest_level->checkAllocated(d_p_dof_index_idx))
    {
        finest_level->allocatePatchData(d_p_dof_index_idx, current_time);
    }
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, finest_level);

    setupSolverVectors(current_time, coarsest_ln, finest_ln);

    if (d_time_stepping_type == TRAPEZOIDAL_RULE)
    {
        d_ib_method_ops->computeLagrangianForce(current_time);
        d_hier_velocity_data_ops->setToScalar(d_f_current_idx, 0.0, false);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_current_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_current_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), current_time);
    }

    if (d_use_structure_predictor)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << "::preprocessIntegrateHierarchy(): Lagrangian forward Euler predictor\n";
        }
        d_ib_method_ops->forwardEulerStep(current_time, new_time);
    }

    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
IBImplicitStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                      const double new_time,
                                                                      const bool skip_synchronize_new_state_data,
                                                                      const int num_cycles)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());
    d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    if (d_enable_logging)
    {
        plog << d_object_name << "::postprocessIntegrateHierarchy(): interpolating velocity to the Lagrangian mesh\n";
    }
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_u_phys_bdry_op->setHomogeneousBc(false);
    d_ib_method_ops->interpolateVelocity(d_u_idx,
                                         getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                         getGhostfillRefineSchedules(d_object_name + "::u"),
                                         new_time);

    deallocateOperatorsAndSolvers();

    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> finest_level = d_hierarchy->getPatchLevel(finest_ln);
    if (finest_level->checkAllocated(d_u_dof_index_idx))
    {
        finest_level->deallocatePatchData(d_u_dof_index_idx);
    }
    if (finest_level->checkAllocated(d_p_dof_index_idx))
    {
        finest_level->deallocatePatchData(d_p_dof_index_idx);
    }
    if (d_rhs_vec)
    {
        d_rhs_vec->deallocateVectorData();
    }
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
IBImplicitStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                                      Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_integrator_is_initialized)
    {
        return;
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_u_dof_index_var = new SideVariable<NDIM, int>(d_object_name + "::u_dof_index");
    d_p_dof_index_var = new CellVariable<NDIM, int>(d_object_name + "::p_dof_index");
    const IntVector<NDIM> ib_ghosts = d_ib_method_ops->getMinimumGhostCellWidth();
    const IntVector<NDIM> no_ghosts = 0;
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, getScratchContext(), ib_ghosts);
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, getScratchContext(), no_ghosts);

    d_ins_hier_integrator->registerBodyForceFunction(d_body_force_fcn);

    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
    d_ib_solver->setHierarchyMathOps(d_hier_math_ops);
    return;
} // initializeHierarchyIntegrator

int
IBImplicitStaggeredHierarchyIntegrator::getNumberOfCycles() const
{
    return d_ins_hier_integrator->getNumberOfCycles();
} // getNumberOfCycles

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBImplicitStaggeredHierarchyIntegrator::integrateHierarchySpecialized(const double current_time,
                                                                      const double new_time,
                                                                      const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchySpecialized(current_time, new_time, cycle_num);

    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator = d_ins_hier_integrator;
    TBOX_ASSERT(ins_hier_integrator);

    const StaggeredStokesIBTimeStepParameters schedule = get_staggered_stokes_ib_time_step_parameters(
        d_time_stepping_type, current_time, new_time, d_object_name + "::integrateHierarchySpecialized()");

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getCurrentContext());
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Solve the coupled velocity-pressure equations.
    ins_hier_integrator->skipCycle(current_time, new_time, cycle_num);

    TBOX_ASSERT(d_sol_vec && d_rhs_vec);

    ins_hier_integrator->setupSolverVectors(d_sol_vec, d_rhs_vec, current_time, new_time, cycle_num);
    if (d_time_stepping_type == TRAPEZOIDAL_RULE)
    {
        d_hier_velocity_data_ops->axpy(
            d_rhs_vec->getComponentDescriptorIndex(0), 0.5, d_f_current_idx, d_rhs_vec->getComponentDescriptorIndex(0));
    }
    reinitializeOperatorsAndSolvers(current_time, new_time);

    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    // Keep the hierarchy Jacobian analytic and matrix-free; FAC owns assembled level coupling.
    d_ib_solver->initializeSolverState(*d_sol_vec, *d_rhs_vec);
    if (!d_ib_solver->solveSystem(*d_sol_vec, *d_rhs_vec))
    {
        TBOX_ERROR(d_object_name << "::integrateHierarchySpecialized(): nonlinear Stokes-IB solve failed\n");
    }

    d_stokes_op->imposeSolBcs(*d_sol_vec);
    if (d_has_velocity_nullspace || d_has_pressure_nullspace)
    {
        ins_hier_integrator->removeNullSpace(d_sol_vec);
    }
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    ins_hier_integrator->resetSolverVectors(d_sol_vec, d_rhs_vec, current_time, new_time, cycle_num);
    deallocateOperatorsAndSolvers();

    // Interpolate the Eulerian velocity to the Lagrangian mesh.
    if (d_enable_logging)
    {
        plog << d_object_name
             << "::integrateHierarchySpecialized(): interpolating "
                "Eulerian velocity to the Lagrangian mesh\n";
    }
    if (d_time_stepping_type == MIDPOINT_RULE)
    {
        d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    }
    else
    {
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    }
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_u_phys_bdry_op->setHomogeneousBc(false);
    d_ib_method_ops->interpolateVelocity(d_u_idx,
                                         getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                         getGhostfillRefineSchedules(d_object_name + "::u"),
                                         schedule.evaluation_time);

    // Update the positions of the Lagrangian structure.
    advance_staggered_stokes_ib_strategy(*d_ib_implicit_ops,
                                         d_time_stepping_type,
                                         current_time,
                                         new_time,
                                         d_object_name + "::integrateHierarchySpecialized()");
    return;
} // integrateHierarchySpecialized

void
IBImplicitStaggeredHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM>> hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    IBHierarchyIntegrator::resetHierarchyConfigurationSpecialized(hierarchy, coarsest_level, finest_level);

    d_num_dofs_per_proc.clear();
    d_vectors_need_init = true;
    deallocateOperatorsAndSolvers();
    return;
} // resetHierarchyConfigurationSpecialized

void
IBImplicitStaggeredHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    IBHierarchyIntegrator::putToDatabaseSpecialized(db);
    db->putInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",
                   IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBImplicitStaggeredHierarchyIntegrator::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    const int ver = db->getInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

void
IBImplicitStaggeredHierarchyIntegrator::setupSolverVectors(const double current_time,
                                                           const int coarsest_ln,
                                                           const int finest_ln)
{
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator = d_ins_hier_integrator;
    TBOX_ASSERT(ins_hier_integrator);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> scratch_ctx = ins_hier_integrator->getScratchContext();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    Pointer<Variable<NDIM>> u_var = ins_hier_integrator->getVelocityVariable();
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, scratch_ctx);
    Pointer<Variable<NDIM>> p_var = ins_hier_integrator->getPressureVariable();
    const int p_scratch_idx = var_db->mapVariableAndContextToIndex(p_var, scratch_ctx);

    const bool reset_solver_vecs =
        d_vectors_need_init || !d_sol_vec || !d_rhs_vec || d_sol_vec->getCoarsestLevelNumber() != coarsest_ln ||
        d_sol_vec->getFinestLevelNumber() != finest_ln || d_sol_vec->getComponentDescriptorIndex(0) != u_scratch_idx ||
        d_sol_vec->getComponentDescriptorIndex(1) != p_scratch_idx;
    if (reset_solver_vecs)
    {
        if (d_rhs_vec)
        {
            free_vector_components(*d_rhs_vec);
        }
        d_sol_vec = new SAMRAIVectorReal<NDIM, double>(
            d_object_name + "::eulerian_sol_vec", d_hierarchy, coarsest_ln, finest_ln);
        d_sol_vec->addComponent(u_var, u_scratch_idx, wgt_sc_idx, d_hier_velocity_data_ops);
        d_sol_vec->addComponent(p_var, p_scratch_idx, wgt_cc_idx, d_hier_pressure_data_ops);

        d_rhs_vec = d_sol_vec->cloneVector(d_object_name + "::eulerian_rhs_vec");
        d_vectors_need_init = false;
    }

    d_rhs_vec->allocateVectorData(current_time);
    return;
} // setupSolverVectors

void
IBImplicitStaggeredHierarchyIntegrator::reinitializeOperatorsAndSolvers(const double current_time,
                                                                        const double new_time)
{
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator = d_ins_hier_integrator;
    TBOX_ASSERT(ins_hier_integrator);

    if (d_ib_solver->getIsInitialized() || d_ib_force_jac || d_ib_interp_op)
    {
        deallocateOperatorsAndSolvers();
    }

    TBOX_ASSERT(d_sol_vec && d_rhs_vec);

    const double dt = new_time - current_time;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = ins_hier_integrator->getCurrentContext();
    Pointer<Variable<NDIM>> u_var = ins_hier_integrator->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);

    const StokesSpecifications* stokes_specs = ins_hier_integrator->getStokesSpecifications();
    TBOX_ASSERT(stokes_specs);
    PoissonSpecifications U_problem_coefs(d_object_name + "::U_problem_coefs");
    // Match the INS right-hand side, whose viscous rule is independent of the IB rule.
    double mass_factor = 1.0;
    double implicit_fraction = 1.0;
    switch (ins_hier_integrator->getViscousTimeSteppingType())
    {
    case BACKWARD_EULER:
        break;
    case FORWARD_EULER:
        implicit_fraction = 0.0;
        break;
    case TRAPEZOIDAL_RULE:
        implicit_fraction = 0.5;
        break;
    case BDF2:
        if (getIntegratorStep() > 0)
        {
            const double omega = dt / d_dt_previous[0];
            mass_factor = (1.0 + 2.0 * omega) / (1.0 + omega);
        }
        break;
    default:
        TBOX_ERROR(d_object_name << "::reinitializeOperatorsAndSolvers(): unsupported viscous time stepping type\n");
    }
    U_problem_coefs.setCConstant(mass_factor * stokes_specs->getRho() / dt +
                                 implicit_fraction * stokes_specs->getLambda());
    U_problem_coefs.setDConstant(-implicit_fraction * stokes_specs->getMu());

    d_ib_implicit_ops->updateFixedLEOperators();

    d_stokes_bc_helper->cacheBcCoefData(ins_hier_integrator->getVelocityBoundaryConditions(), new_time, d_hierarchy);

    d_stokes_op->setVelocityPoissonSpecifications(U_problem_coefs);
    d_stokes_op->setPhysicalBcCoefs(ins_hier_integrator->getVelocityBoundaryConditions(),
                                    ins_hier_integrator->getPressureBoundaryConditions());
    d_stokes_op->setPhysicalBoundaryHelper(d_stokes_bc_helper);
    d_stokes_op->setTimeInterval(current_time, new_time);
    d_stokes_op->setSolutionTime(new_time);

    StaggeredStokesIBOperator::Context ctx;
    ctx.ib_implicit_ops = d_ib_implicit_ops;
    ctx.stokes_op = d_stokes_op;
    ctx.u_phys_bdry_op = d_u_phys_bdry_op;
    ctx.hier_velocity_data_ops = d_hier_velocity_data_ops;
    ctx.u_synch_scheds = getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN");
    ctx.u_ghost_fill_scheds = getGhostfillRefineSchedules(d_object_name + "::u");
    ctx.f_prolongation_scheds = getProlongRefineSchedules(d_object_name + "::f");
    ctx.u_idx = d_u_idx;
    ctx.f_idx = d_f_idx;
    ctx.u_current_idx = u_current_idx;
    ctx.u_dof_index_idx = d_u_dof_index_idx;
    ctx.p_dof_index_idx = d_p_dof_index_idx;
    ctx.use_fixed_le_operators = true;
    ctx.time_stepping_type = d_time_stepping_type;

    d_ib_op->setOperatorContext(ctx);
    d_ib_jac_op->setOperatorContext(ctx);
    const StaggeredStokesIBTimeStepParameters schedule = get_staggered_stokes_ib_time_step_parameters(
        d_time_stepping_type, current_time, new_time, d_object_name + "::reinitializeOperatorsAndSolvers()");
    d_ib_implicit_ops->constructLagrangianForceJacobian(d_ib_force_jac, MATAIJ, schedule.evaluation_time);

    d_ib_implicit_ops->constructInterpOp(
        d_ib_interp_op, d_jac_kernel, d_num_dofs_per_proc, d_u_dof_index_idx, schedule.evaluation_time);

    d_ib_jac_pc->setVelocityPoissonSpecifications(U_problem_coefs);
    d_ib_jac_pc->setPhysicalBcCoefs(ins_hier_integrator->getIntermediateVelocityBoundaryConditions(),
                                    ins_hier_integrator->getProjectionBoundaryConditions());
    d_ib_jac_pc->setPhysicalBoundaryHelper(d_stokes_bc_helper);
    d_ib_jac_pc->setComponentsHaveNullSpace(d_has_velocity_nullspace, d_has_pressure_nullspace);
    d_ib_jac_pc->setIBTimeSteppingType(d_time_stepping_type);
    d_ib_jac_pc->setIBForceJacobian(d_ib_force_jac);
    d_ib_jac_pc->setIBInterpOp(d_ib_interp_op);
    d_ib_jac_pc->setIBImplicitStrategy(d_ib_implicit_ops);
    d_ib_jac_pc->setTimeInterval(current_time, new_time);
    d_ib_jac_pc->setSolutionTime(new_time);
    d_ib_jac_pc->setHomogeneousBc(true);

    d_ib_solver->setHomogeneousBc(false);
    d_ib_solver->setSolutionTime(new_time);
    d_ib_solver->setTimeInterval(current_time, new_time);
    d_ib_op->setHomogeneousBc(false);
    d_ib_op->setSolutionTime(new_time);
    d_ib_op->setTimeInterval(current_time, new_time);
    d_ib_jac_op->setHomogeneousBc(true);
    d_ib_jac_op->setSolutionTime(new_time);
    d_ib_jac_op->setTimeInterval(current_time, new_time);
    Pointer<IBTK::KrylovLinearSolver> linear_solver = d_ib_solver->getLinearSolver();
    linear_solver->setInitialGuessNonzero(false);

    return;
} // reinitializeOperatorsAndSolvers

void
IBImplicitStaggeredHierarchyIntegrator::deallocateOperatorsAndSolvers()
{
    if (d_ib_solver && d_ib_solver->getIsInitialized())
    {
        d_ib_solver->deallocateSolverState();
    }
    if (d_ib_force_jac)
    {
        const int ierr = MatDestroy(&d_ib_force_jac);
        IBTK_CHKERRQ(ierr);
        d_ib_force_jac = nullptr;
    }
    if (d_ib_interp_op)
    {
        const int ierr = MatDestroy(&d_ib_interp_op);
        IBTK_CHKERRQ(ierr);
        d_ib_interp_op = nullptr;
    }
    return;
} // deallocateOperatorsAndSolvers

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
