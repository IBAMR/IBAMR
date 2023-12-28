// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibamr/IBImplicitHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/IBTK_MPI.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBImplicitHierarchyIntegrator restart file data.
static const int IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitHierarchyIntegrator::IBImplicitHierarchyIntegrator(const std::string& object_name,
                                                             Pointer<Database> input_db,
                                                             Pointer<IBImplicitStrategy> ib_implicit_ops,
                                                             Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                                                             bool register_for_restart)
    : IBHierarchyIntegrator(object_name, input_db, ib_implicit_ops, ins_hier_integrator, register_for_restart),
      d_ib_implicit_ops(ib_implicit_ops)
{
    // Set default configuration options.
    d_use_structure_predictor = false;
    d_use_fixed_LE_operators = false;
    d_solve_for_position = false;
    d_eta = std::numeric_limits<double>::quiet_NaN();

    // Set options from input.
    if (input_db)
    {
        if (input_db->keyExists("use_structure_predictor"))
            d_use_structure_predictor = input_db->getBool("use_structure_predictor");
        if (input_db->keyExists("use_fixed_LE_operators"))
            d_use_fixed_LE_operators = input_db->getBool("use_fixed_LE_operators");
        if (input_db->keyExists("solve_for_position")) d_solve_for_position = input_db->getBool("solve_for_position");
        if (input_db->keyExists("eta")) d_eta = input_db->getDouble("eta");
        if (input_db->keyExists("f_evals_target_min"))
            d_snes_f_evals_target_min = input_db->getInteger("f_evals_target_min");
        if (input_db->keyExists("f_evals_target_max"))
            d_snes_f_evals_target_max = input_db->getInteger("f_evals_target_max");
    }
    d_ib_implicit_ops->setUseFixedLEOperators(d_use_fixed_LE_operators);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
} // IBImplicitHierarchyIntegrator

void
IBImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
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
        TBOX_ERROR("IBImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(): time_stepping_type = "
                   << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                   << "  only supported time_stepping_types are:\n"
                   << "    " << enum_to_string<TimeSteppingType>(BACKWARD_EULER) << "\n"
                   << "    " << enum_to_string<TimeSteppingType>(TRAPEZOIDAL_RULE) << "\n"
                   << "    " << enum_to_string<TimeSteppingType>(MIDPOINT_RULE) << "\n");
    }

    // Initialize IB data.
    d_ib_implicit_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize the fluid solver.
    d_ins_hier_integrator->setSkipEnforceNumCycles();

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case TRAPEZOIDAL_RULE:
        if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(current_time);
        if (d_enable_logging)
            plog << d_object_name
                 << "::preprocessIntegrateHierarchy(): spreading Lagrangian force "
                    "to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), current_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_hier_velocity_data_ops->copyData(d_f_current_idx, d_f_idx);
        break;
    case BACKWARD_EULER:
    case MIDPOINT_RULE:
        // intentionally blank
        break;
    default:
        TBOX_ERROR(
            d_object_name << "::preprocessIntegrateHierarchy():\n"
                          << "  unsupported time stepping type: "
                          << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                          << "  supported time stepping types are: BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute an initial prediction of the updated positions of the Lagrangian
    // structure.
    //
    // NOTE: The velocity should already have been interpolated to the
    // curvilinear mesh and should not need to be re-interpolated.
    if (d_use_structure_predictor)
    {
        if (d_enable_logging)
            plog << d_object_name << "::preprocessIntegrateHierarchy(): performing Lagrangian forward Euler step\n";
        d_ib_implicit_ops->forwardEulerStep(current_time, new_time);
        if (!d_solve_for_position) d_ib_method_ops->computeLagrangianForce(new_time);
    }

    // Zero out the function evaluation counter.
    d_snes_f_evals = 0;

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
IBImplicitHierarchyIntegrator::integrateHierarchy(const double current_time, const double new_time, const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    d_current_time = current_time;
    d_new_time = new_time;
    d_cycle_num = cycle_num;

    if (d_use_fixed_LE_operators) d_ib_implicit_ops->updateFixedLEOperators();

    // Setup Lagrangian vectors used in solving the implicit IB equations.
    PetscErrorCode ierr;
    Vec X, F, R, R_work;
    d_ib_implicit_ops->createSolverVecs(&X, &F, &R);
    d_ib_implicit_ops->setupSolverVecs(X, F, R);
    ierr = VecDuplicate(R, &R_work);
    IBTK_CHKERRQ(ierr);

    // Solve the implicit IB equations.
    d_ins_cycle_num = 0;

    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetFunction(snes, R_work, IBFunction, this);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(snes, "ib_");
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes);
    IBTK_CHKERRQ(ierr);

    // We either solve for position or force.
    Vec& Y = d_solve_for_position ? X : F;

    ierr = SNESSolve(snes, R, Y);
    IBTK_CHKERRQ(ierr);
    SNESConvergedReason converged_reason;
    ierr = SNESGetConvergedReason(snes, &converged_reason);
    IBTK_CHKERRQ(ierr);
    const bool snes_diverged = converged_reason <= 0;
    int f_evals;
    ierr = SNESGetNumberFunctionEvals(snes, &f_evals);
    d_snes_f_evals += f_evals;
    IBTK_CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);
    IBTK_CHKERRQ(ierr);

    // Check to see if we need to re-do this time step.
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    const auto u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                d_ins_hier_integrator->getNewContext());
    const double dt = new_time - current_time;
    const auto cfl = INSHierarchyIntegrator::compute_CFL_number(u_new_idx, dt, d_hierarchy);
    const auto cfl_max = d_ins_hier_integrator->getMaximumCFLNumber();
    const bool exceeded_cfl_number = cfl > d_cfl_max_tol * cfl_max;

    d_redo_time_step = snes_diverged || exceeded_cfl_number;

    // Ensure that the state variables are consistent with the solution.
    if (!d_redo_time_step) updateSolution(Y, nullptr);

    // Deallocate temporary data.
    ierr = VecDestroy(&X);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&F);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&R);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&R_work);
    IBTK_CHKERRQ(ierr);

    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
} // integrateHierarchy

void
IBImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                             const double new_time,
                                                             const bool skip_synchronize_new_state_data,
                                                             const int num_cycles)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    if (d_enable_logging)
        plog << d_object_name
             << "::postprocessIntegrateHierarchy(): interpolating Eulerian "
                "velocity to the Lagrangian mesh\n";
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_u_phys_bdry_op->setHomogeneousBc(false);
    d_ib_implicit_ops->interpolateVelocity(d_u_idx,
                                           getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                           getGhostfillRefineSchedules(d_object_name + "::u"),
                                           new_time);

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    // postprocess the objects this class manages...
    d_ib_implicit_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, ins_num_cycles);

    // ... and postprocess ourself.
    HierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Determine the CFL number.
    const auto cfl = INSHierarchyIntegrator::compute_CFL_number(u_new_idx, dt, d_hierarchy);
    d_regrid_fluid_cfl_estimate += cfl;

    // Not all IBStrategy objects implement this so make it optional (-1.0 is
    // the default value)
    if (d_regrid_structure_cfl_interval != -1.0)
        d_regrid_structure_cfl_estimate = d_ib_method_ops->getMaxPointDisplacement();

    if (d_enable_logging)
    {
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl << "\n";
        plog << d_object_name
             << "::postprocessIntegrateHierarchy(): Eulerian estimate of "
                "upper bound on IB point displacement since last regrid = "
             << d_regrid_fluid_cfl_estimate << "\n";

        if (d_regrid_structure_cfl_interval != -1.0)
        {
            plog << d_object_name
                 << "::postprocessIntegrateHierarchy(): Lagrangian estimate of "
                    "upper bound on IB point displacement since last regrid = "
                 << d_regrid_structure_cfl_estimate << "\n";
        }
    }

    // Keep track of the number of function evaluations.
    d_snes_f_evals_previous = d_snes_f_evals;

    // Deallocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_f_idx);
        if (d_f_current_idx != -1) level->deallocatePatchData(d_f_current_idx);
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
IBImplicitHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                             Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Register body force function with INSHierarchyIntegrator
    d_ins_hier_integrator->registerBodyForceFunction(new IBEulerianForceFunction(this));

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
    return;
} // initializeHierarchyIntegrator

int
IBImplicitHierarchyIntegrator::getNumberOfCycles() const
{
    return d_ins_hier_integrator->getNumberOfCycles();
} // getNumberOfCycles

/////////////////////////////// PROTECTED ////////////////////////////////////

double
IBImplicitHierarchyIntegrator::getMaximumTimeStepSizeSpecialized()
{
    auto dt_max = IBHierarchyIntegrator::getMaximumTimeStepSizeSpecialized();
    const bool initial_time = IBTK::rel_equal_eps(d_integrator_time, d_start_time);
    if (initial_time) return dt_max;
    if (d_snes_f_evals_previous < d_snes_f_evals_target_min)
    {
        pout << "BELOW TARGET!!!\n";
        dt_max = std::min(dt_max, std::max(1.0, d_dt_growth_factor) * d_dt_previous[0]);
    }
    if (d_snes_f_evals_previous > d_snes_f_evals_target_max)
    {
        pout << "ABOVE TARGET!!!\n";
        dt_max = std::min(dt_max, std::min(1.0, d_dt_shrink_factor) * d_dt_previous[0]);
    }
    return dt_max;
} // getMaximumTimeStepSizeSpecialized

void
IBImplicitHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    IBHierarchyIntegrator::putToDatabaseSpecialized(db);
    db->putInteger("IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION", IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBImplicitHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_IMPLICIT_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

void
IBImplicitHierarchyIntegrator::updateSolution(Vec Y, Vec R)
{
    Vec X = d_solve_for_position ? Y : nullptr;
    Vec F = d_solve_for_position ? nullptr : Y;

    const double current_time = d_current_time;
    const double new_time = d_new_time;
    const double half_time = current_time + 0.5 * (new_time - current_time);
    const int cycle_num = d_cycle_num;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getCurrentContext());
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    if (d_solve_for_position)
    {
        // Set the current position data.
        d_ib_implicit_ops->setUpdatedPosition(X);
    }
    else
    {
        // Set the current force data.
        d_ib_implicit_ops->setUpdatedForce(F);
    }

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
        if (d_solve_for_position)
        {
            if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
            d_ib_method_ops->computeLagrangianForce(d_new_time);
        }
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        break;
    case MIDPOINT_RULE:
        if (d_solve_for_position)
        {
            if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
            d_ib_method_ops->computeLagrangianForce(half_time);
        }
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        break;
    case TRAPEZOIDAL_RULE:
        if (d_solve_for_position)
        {
            if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
            d_ib_method_ops->computeLagrangianForce(new_time);
        }
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_hier_velocity_data_ops->linearSum(d_f_idx, 0.5, d_f_current_idx, 0.5, d_f_idx);
        break;
    default:
        TBOX_ERROR(
            d_object_name << "::integrateHierarchy():\n"
                          << "  unsupported time stepping type: "
                          << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                          << "  supported time stepping types are: BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Solve the incompressible Navier-Stokes equations.
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): solving the incompressible Navier-Stokes equations\n";
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    d_ins_hier_integrator->integrateHierarchy(current_time, new_time, d_ins_cycle_num);
    d_ins_cycle_num++;
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
    case TRAPEZOIDAL_RULE:
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        break;
    case MIDPOINT_RULE:
        d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             half_time);
        break;
    default:
        TBOX_ERROR(
            d_object_name << "::integrateHierarchy():\n"
                          << "  unsupported time stepping type: "
                          << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                          << "  supported time stepping types are: BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute the time stepping-based residual before updating the IB variables.
    //
    // (Once the IB variables are updated, we can no longer get the residual from the IB ops object.)
    if (R != nullptr && d_solve_for_position)
    {
        switch (d_time_stepping_type)
        {
        case BACKWARD_EULER:
        {
            d_ib_implicit_ops->computeResidualBackwardEuler(R);
            break;
        }
        case MIDPOINT_RULE:
        {
            d_ib_implicit_ops->computeResidualMidpointRule(R);
            break;
        }
        case TRAPEZOIDAL_RULE:
        {
            d_ib_implicit_ops->computeResidualTrapezoidalRule(R);
            break;
        }
        default:
            TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                     << "  unsupported time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                                     << "  supported time stepping types are: BACKWARD_EULER, "
                                        "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }

    // Update the IB variables.
    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
    {
        d_ib_implicit_ops->backwardEulerStep(current_time, new_time);
        break;
    }
    case MIDPOINT_RULE:
    {
        d_ib_implicit_ops->midpointStep(current_time, new_time);
        break;
    }
    case TRAPEZOIDAL_RULE:
    {
        d_ib_implicit_ops->trapezoidalStep(current_time, new_time);
        break;
    }
    default:
        TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                 << "  unsupported time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                                 << "  supported time stepping types are: BACKWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute the force-based residual after updating the IB state variables.
    if (R != nullptr && !d_solve_for_position)
    {
        switch (d_time_stepping_type)
        {
        case BACKWARD_EULER:
        case TRAPEZOIDAL_RULE:
        {
            d_ib_method_ops->computeLagrangianForce(new_time);
            d_ib_implicit_ops->getUpdatedForce(R);
            PetscErrorCode ierr = VecAXPBY(R, 1.0, -1.0, F);
            IBTK_CHKERRQ(ierr);
            break;
        }
        case MIDPOINT_RULE:
        {
            TBOX_ERROR("unimplemented!");
            break;
        }
        default:
            TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                     << "  unsupported time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                                     << "  supported time stepping types are: BACKWARD_EULER, "
                                        "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }
    return;
}

PetscErrorCode
IBImplicitHierarchyIntegrator::IBFunction(SNES /*snes*/, Vec Y, Vec R, void* ctx)
{
    auto integrator = static_cast<IBImplicitHierarchyIntegrator*>(ctx);
    integrator->updateSolution(Y, R);
    return 0;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
