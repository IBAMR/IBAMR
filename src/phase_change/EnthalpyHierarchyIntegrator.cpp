// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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
#include "ibamr/AdvDiffConservativeMassScalarTransportRKIntegrator.h"
#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/EnthalpyHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CCLaplaceOperator.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/PoissonSolver.h"

#include "BasePatchHierarchy.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellDataFactory.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchFaceDataOpsReal.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
} // namespace hier
} // namespace SAMRAI

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
// Version of INSHierarchyIntegrator restart file data.
static const int ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION = 4;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int NOGHOSTS = 0;

static const double H_LIM = 0.5;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

EnthalpyHierarchyIntegrator::EnthalpyHierarchyIntegrator(const std::string& object_name,
                                                         Pointer<Database> input_db,
                                                         bool register_for_restart)
    : PhaseChangeHierarchyIntegrator(object_name, input_db, register_for_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();

    getFromInput(input_db, from_restart);
    if (from_restart) getFromRestart();

    return;
} // EnthalpyHierarchyIntegrator

void
EnthalpyHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                           Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Perform hierarchy initialization operations common to all implementations
    // of PhaseChangeHierarchyIntegrator.
    PhaseChangeHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Register additional variables required for present time stepping algorithm.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> no_ghosts = NOGHOSTS;

    // Register specific enthalpy.
    d_h_var = new CellVariable<NDIM, double>(d_object_name + "::enthalpy");
    registerVariable(d_h_current_idx,
                     d_h_new_idx,
                     d_h_scratch_idx,
                     d_h_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");
    if (d_visit_writer) d_visit_writer->registerPlotQuantity("enthalpy", "SCALAR", d_h_current_idx);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_T_pre_var = new CellVariable<NDIM, double>(d_object_name + "::T_pre_var");
    d_T_pre_idx = var_db->registerVariableAndContext(d_T_pre_var, getCurrentContext());

    d_dh_dT_var = new CellVariable<NDIM, double>(d_object_name + "::dh_dT_var");
    d_dh_dT_scratch_idx = var_db->registerVariableAndContext(d_dh_dT_var, getCurrentContext(), no_ghosts);

    d_grad_T_var = new SideVariable<NDIM, double>(d_object_name + "::grad_T");
    d_grad_T_idx =
        var_db->registerVariableAndContext(d_grad_T_var, var_db->getContext(d_object_name + "grad_T::SCRATCH"));

    if (d_solve_mass_conservation)
    {
        // Set various objects with conservative time integrator.
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
        rho_p_cc_integrator->setCellCenteredTransportQuantityBoundaryConditions(d_h_bc_coef);
    }

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
EnthalpyHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                          const double new_time,
                                                          const int num_cycles)
{
    PhaseChangeHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_dh_dT_scratch_idx)) level->allocatePatchData(d_dh_dT_scratch_idx, current_time);
        if (!level->checkAllocated(d_T_pre_idx)) level->allocatePatchData(d_T_pre_idx, current_time);
        if (!level->checkAllocated(d_grad_T_idx)) level->allocatePatchData(d_grad_T_idx, current_time);
    }

    const int H_current_idx = var_db->mapVariableAndContextToIndex(d_H_var, getCurrentContext());

    // Initialize enthalpy h only at the start of the simulation.
    if (initial_time)
        computeEnthalpyBasedOnTemperature(
            d_h_current_idx, d_T_current_idx, d_rho_current_idx, d_lf_current_idx, H_current_idx);

    if (d_solve_mass_conservation)
    {
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;

        // Set the velocities used to update the density and the previous time step
        // size
        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ -1, /*current*/ d_u_adv_current_idx, /*new*/ -1);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*current*/ d_h_current_idx, /*new*/ -1);
        }
        else
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ d_U_old_current_idx, /*current*/ d_u_adv_current_idx, /*new*/ -1);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*current*/ d_h_current_idx, /*new*/ -1);
            d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        // Integrate density and convective term of energy equation.
        d_rho_p_integrator->integrate(dt);
    }
    // Setup the operators and solvers and compute the right-hand-side terms.
    // Setup the problem coefficients for the linear solver.
    double alpha = 0.0;
    switch (d_T_diffusion_time_stepping_type)
    {
    case BACKWARD_EULER:
        alpha = 1.0;
        break;
    case FORWARD_EULER:
        alpha = 0.0;
        break;
    case TRAPEZOIDAL_RULE:
        alpha = 0.5;
        break;
    default:
        TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                 << "  unsupported diffusion time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_T_diffusion_time_stepping_type) << " \n"
                                 << "  valid choices are: BACKWARD_EULER, "
                                    "FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }
    PoissonSpecifications T_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_T_var->getName());

    // There is no coefficient in the RHS for the enthalpy formulation.
    T_rhs_op_spec.setCZero();

    const double apply_time = current_time;
    for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
    {
        d_reset_kappa_fcns[k](d_T_diffusion_coef_cc_current_idx,
                              d_T_diffusion_coef_cc_var,
                              d_hier_math_ops,
                              -1 /*cycle_num*/,
                              apply_time,
                              current_time,
                              new_time,
                              d_reset_kappa_fcns_ctx[k]);
    }

    // Interpolate the cell-centered diffusion coef to side-centered.
    d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_scratch_idx, d_T_diffusion_coef_cc_current_idx);
    d_k_bdry_bc_fill_op->fillData(current_time);

    if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
    {
        interpolateCCToSCSimpleAveraging(d_T_diffusion_coef_current_idx, d_T_diffusion_coef_cc_scratch_idx);
    }
    else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
    {
        interpolateCCToSCHarmonicAveraging(d_T_diffusion_coef_current_idx, d_T_diffusion_coef_cc_scratch_idx);
    }
    else
    {
        TBOX_ERROR("this statement should not be reached");
    }

    // For plotting purpose.
    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(d_D_cc_new_idx,
                            d_D_cc_var,
                            d_T_diffusion_coef_current_idx,
                            d_T_diffusion_coef_var,
                            d_no_fill_op,
                            d_integrator_time,
                            synch_cf_interface);

    d_hier_sc_data_ops->scale(d_T_diffusion_coef_rhs_scratch_idx, (1.0 - alpha), d_T_diffusion_coef_current_idx);
    T_rhs_op_spec.setDPatchDataId(d_T_diffusion_coef_rhs_scratch_idx);

    // Initialize the RHS operator and compute the RHS vector for the temperature equation.
    Pointer<LaplaceOperator> T_rhs_op = d_T_rhs_op;
    T_rhs_op->setPoissonSpecifications(T_rhs_op_spec);
    T_rhs_op->setPhysicalBcCoef(d_T_bc_coef);
    T_rhs_op->setHomogeneousBc(false);
    T_rhs_op->setSolutionTime(current_time);
    T_rhs_op->setTimeInterval(current_time, new_time);
    if (d_T_rhs_op_needs_init)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the RHS operator for" << d_T_var->getName() << "\n";
        }
        T_rhs_op->initializeOperatorState(*d_T_sol, *d_T_rhs);
        d_T_rhs_op_needs_init = false;
    }
    d_hier_cc_data_ops->copyData(d_T_scratch_idx, d_T_current_idx, false);
    T_rhs_op->apply(*d_T_sol, *d_T_rhs);

    d_hier_cc_data_ops->copyData(d_h_new_idx, d_h_current_idx);
    d_hier_cc_data_ops->copyData(d_T_new_idx, d_T_current_idx);
    d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_new_idx, d_T_diffusion_coef_cc_current_idx);
    d_hier_cc_data_ops->copyData(d_lf_new_idx, d_lf_current_idx);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
EnthalpyHierarchyIntegrator::integrateHierarchySpecialized(const double current_time,
                                                           const double new_time,
                                                           const int cycle_num)
{
    AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchySpecialized(current_time, new_time, cycle_num);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE({
            pout << "EnthalpyHierarchyIntegrator::integrateHierarchy():\n"
                 << "  WARNING: num_cycles = " << d_current_num_cycles
                 << " but expected num_cycles = " << expected_num_cycles << ".\n";
        });
    }

    // Perform a single step of fixed point iteration.
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());

    // In the special case of conservative discretization, the updated
    // density is calculated by the mass and convective integrator.
    Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
    // Update N_idx if necessary
    if (cycle_num > 0 && d_solve_mass_conservation)
    {
        const double dt = new_time - current_time;
        const double half_time = current_time + 0.5 * dt;
        d_rho_p_integrator->setSolutionTime(half_time);

        // Set the cycle number
        d_rho_p_integrator->setCycleNumber(cycle_num);

        // Set the patch data index for convective derivative.
        d_rho_p_integrator->setConvectiveDerivativePatchDataIndex(d_T_N_scratch_idx);

        // Always set to current because we want to update rho^{n} to rho^{n+1}
        d_rho_p_integrator->setDensityPatchDataIndex(d_rho_current_idx);

        // Set the velocities used to update the density
        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ -1, /*current*/ d_u_adv_current_idx, /*new*/ d_u_adv_new_idx);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*current*/ d_h_current_idx, /*new*/ d_h_new_idx);
        }
        else
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ d_U_old_current_idx,
                /*current*/ d_u_adv_current_idx,
                /*new*/ d_u_adv_new_idx);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*current*/ d_h_current_idx,
                /*new*/ d_h_new_idx);

            d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        d_rho_p_integrator->integrate(dt);
    }

    d_updated_rho_idx = d_rho_p_integrator ? d_rho_p_integrator->getUpdatedDensityPatchDataIndex() : d_rho_new_idx;
    d_hier_cc_data_ops->copyData(d_rho_new_idx,
                                 d_updated_rho_idx,
                                 /*interior_only*/ true);

    // Account for the convective acceleration term N_full.
    if (d_u_adv_var) d_hier_cc_data_ops->axpy(d_T_rhs_scratch_idx, -1.0, d_T_N_scratch_idx, d_T_rhs_scratch_idx);

    PoissonSpecifications T_solver_spec(d_object_name + "::solver_spec::" + d_T_var->getName());

    double lf_relative_iteration_error = 1.0;
    double inner_iterations = 1.0;

    // Inner iterations for the Newton-Ralphson scheme.
    while (lf_relative_iteration_error >= d_lf_iteration_error_tolerance && inner_iterations <= d_max_inner_iterations)
    {
        // Setup the problem coefficients for the linear solve
        double alpha = 0.0;
        switch (d_T_diffusion_time_stepping_type)
        {
        case BACKWARD_EULER:
            alpha = 1.0;
            break;
        case FORWARD_EULER:
            alpha = 0.0;
            break;
        case TRAPEZOIDAL_RULE:
            alpha = 0.5;
            break;
        default:
            TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                     << "  unsupported diffusion time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_T_diffusion_time_stepping_type) << " \n"
                                     << "  valid choices are: BACKWARD_EULER, "
                                        "FORWARD_EULER, TRAPEZOIDAL_RULE\n");
        }

        // Interpolate the cell-centered diffusion coef to side-centered.
        d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_scratch_idx, d_T_diffusion_coef_cc_new_idx);
        d_k_bdry_bc_fill_op->fillData(new_time);

        if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
        {
            interpolateCCToSCSimpleAveraging(d_T_diffusion_coef_new_idx, d_T_diffusion_coef_cc_scratch_idx);
        }
        else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
        {
            interpolateCCToSCHarmonicAveraging(d_T_diffusion_coef_new_idx, d_T_diffusion_coef_cc_scratch_idx);
        }
        else
        {
            TBOX_ERROR("this statement should not be reached");
        }

        // For plotting purpose.
        static const bool synch_cf_interface = true;
        d_hier_math_ops->interp(d_D_cc_new_idx,
                                d_D_cc_var,
                                d_T_diffusion_coef_new_idx,
                                d_T_diffusion_coef_var,
                                d_no_fill_op,
                                d_integrator_time,
                                synch_cf_interface);

        d_hier_sc_data_ops->scale(d_T_diffusion_coef_scratch_idx, -alpha, d_T_diffusion_coef_new_idx);
        T_solver_spec.setDPatchDataId(d_T_diffusion_coef_scratch_idx);

        computeEnthalpyDerivative(d_dh_dT_scratch_idx, d_T_new_idx, H_new_idx);

        // Set rho*Cp/dt.
        d_hier_cc_data_ops->multiply(d_C_new_idx, d_rho_new_idx, d_dh_dT_scratch_idx);
        d_hier_cc_data_ops->copyData(d_T_C_idx, d_C_new_idx);
        d_hier_cc_data_ops->scale(d_T_C_idx, 1.0 / dt, d_T_C_idx);
        T_solver_spec.setCPatchDataId(d_T_C_idx);

        // Initialize the linear solver for temperature equation.
        Pointer<PoissonSolver> T_solver = d_T_solver;
        T_solver->setPoissonSpecifications(T_solver_spec);
        T_solver->setPhysicalBcCoef(d_T_bc_coef);
        T_solver->setHomogeneousBc(false);
        T_solver->setSolutionTime(new_time);
        T_solver->setTimeInterval(current_time, new_time);
        // Initialize solver each time because the coefficients are changing.
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the solvers for" << d_T_var->getName() << "\n";
        }
        T_solver->initializeSolverState(*d_T_sol, *d_T_rhs);

        // Account for forcing terms.
        if (d_T_F_fcn)
        {
            d_T_F_fcn->setDataOnPatchHierarchy(d_T_F_scratch_idx, d_T_F_var, d_hierarchy, half_time);
        }
        else
            d_hier_cc_data_ops->setToScalar(d_T_F_scratch_idx, 0.0);

        // Compute and add temporal and linearized terms to the RHS of the energy equation.
        addTemporalAndLinearTermstoRHSOfEnergyEquation(d_T_F_scratch_idx, dt);
        d_hier_cc_data_ops->axpy(d_T_rhs_scratch_idx, +1.0, d_T_F_scratch_idx, d_T_rhs_scratch_idx);

        // Storing T^n+1,m.
        d_hier_cc_data_ops->copyData(d_T_pre_idx, d_T_new_idx);

        // Solve for T(n+1, m+1).
        T_solver->solveSystem(*d_T_sol, *d_T_rhs);
        d_hier_cc_data_ops->copyData(d_T_new_idx, d_T_scratch_idx);

        if (d_enable_logging && d_enable_logging_solver_iterations)
            plog << d_object_name << ":" << d_T_var->getName()
                 << "::integrateHierarchy():diffusion solve number of iterations = " << T_solver->getNumIterations()
                 << "\n";
        if (d_enable_logging)
            plog << d_object_name << ":" << d_T_var->getName()
                 << "::integrateHierarchy():diffusion solve residual norm        = " << T_solver->getResidualNorm()
                 << "\n";
        if (T_solver->getNumIterations() == T_solver->getMaxIterations())
        {
            pout << d_object_name << ":" << d_T_var->getName()
                 << "::integrateHierarchy():WARNING: linear solver iterations == max iterations\n";
        }

        // Find h^n+1, m+1.
        updateEnthalpy(d_h_new_idx, d_T_new_idx, d_T_pre_idx);

        // Find T^n+1,m+1 based on h^n+1, m+1.
        computeTemperatureBasedOnEnthalpy(d_T_new_idx, d_h_new_idx, H_new_idx);

        // Find lf^n+1, m+1 based on h^n+1, m+1.
        d_hier_cc_data_ops->copyData(d_lf_pre_idx, d_lf_new_idx);
        computeLiquidFraction(d_lf_new_idx, d_h_new_idx, H_new_idx);

        // Update specific heat
        const double apply_time = new_time;
        for (unsigned k = 0; k < d_reset_specific_heat_fcns.size(); ++k)
        {
            d_reset_specific_heat_fcns[k](d_specific_heat_new_idx,
                                          d_specific_heat_var,
                                          d_hier_math_ops,
                                          -1 /*cycle_num*/,
                                          apply_time,
                                          current_time,
                                          new_time,
                                          d_reset_specific_heat_fcns_ctx[k]);
        }

        // Update conductivity
        for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
        {
            d_reset_kappa_fcns[k](d_T_diffusion_coef_cc_new_idx,
                                  d_T_diffusion_coef_cc_var,
                                  d_hier_math_ops,
                                  -1 /*cycle_num*/,
                                  apply_time,
                                  current_time,
                                  new_time,
                                  d_reset_kappa_fcns_ctx[k]);
        }

        // Finding L2 norm of lf^m-1 iteration.
        const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
        lf_relative_iteration_error = d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx);

        // Finding lf^m - lf^m-1.
        d_hier_cc_data_ops->subtract(d_lf_pre_idx, d_lf_new_idx, d_lf_pre_idx);

        plog << "liquid fraction relative error norms at Newton iteration: " << inner_iterations << "\n"
             << "L1 : " << d_hier_cc_data_ops->L1Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error)
             << " || "
             << "L2 : " << d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error)
             << " || "
             << "L_oo : " << d_hier_cc_data_ops->maxNorm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error)
             << "\n";

        lf_relative_iteration_error =
            d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error);
        inner_iterations++;

        d_hier_cc_data_ops->axpy(d_T_rhs_scratch_idx, -1.0, d_T_F_scratch_idx, d_T_rhs_scratch_idx);
        d_hier_cc_data_ops->copyData(d_T_F_new_idx, d_T_F_scratch_idx);
    }

    // Reset the right-hand side vector.
    if (d_u_adv_var) d_hier_cc_data_ops->axpy(d_T_rhs_scratch_idx, +1.0, d_T_N_scratch_idx, d_T_rhs_scratch_idx);

    // Compute the source term for Div U equation.
    computeDivergenceVelocitySourceTerm(d_Div_U_F_idx, new_time);

    return;
} // integrateHierarchy

void
EnthalpyHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                           const double new_time,
                                                           const bool skip_synchronize_new_state_data,
                                                           const int num_cycles)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Deallocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_T_pre_idx);
        level->deallocatePatchData(d_dh_dT_scratch_idx);
        level->deallocatePatchData(d_grad_T_idx);
    }

    PhaseChangeHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
EnthalpyHierarchyIntegrator::setEnthalpyBcCoef(RobinBcCoefStrategy<NDIM>* h_bc_coef)
{
    d_h_bc_coef = h_bc_coef;
    return;
} // setEnthalpyBcCoef

void
EnthalpyHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    db->putInteger("ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION", ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("specific_heat_liquid", d_specific_heat_liquid);
    db->putDouble("specific_heat_solid", d_specific_heat_solid);
    db->putDouble("specific_heat_gas", d_specific_heat_gas);
    db->putDouble("liquidus_temperature", d_liquidus_temperature);
    db->putDouble("solidus_temperature", d_solidus_temperature);
    db->putDouble("reference_temperature", d_reference_temperature);
    db->putDouble("specific_heat_mushy", d_specific_heat_mushy);
    db->putDouble("gas_liquid_fraction", d_gas_liquid_fraction);
    db->putInteger("max_inner_iterations", d_max_inner_iterations);
    db->putDouble("lf_iteration_error_tolerance", d_lf_iteration_error_tolerance);

    AdvDiffSemiImplicitHierarchyIntegrator::putToDatabaseSpecialized(db);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
EnthalpyHierarchyIntegrator::addTemporalAndLinearTermstoRHSOfEnergyEquation(int F_scratch_idx, const double dt)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_new_data = patch->getPatchData(d_T_new_idx);
            Pointer<CellData<NDIM, double> > h_new_data = patch->getPatchData(d_h_new_idx);
            Pointer<CellData<NDIM, double> > h_current_data = patch->getPatchData(d_h_current_idx);
            Pointer<CellData<NDIM, double> > rho_new_data = patch->getPatchData(d_rho_new_idx);
            Pointer<CellData<NDIM, double> > rho_current_data = patch->getPatchData(d_rho_current_idx);
            Pointer<CellData<NDIM, double> > dh_dT_data = patch->getPatchData(d_dh_dT_scratch_idx);
            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_scratch_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*F_data)(ci) += -1.0 / dt *
                                 ((*rho_new_data)(ci) * ((*h_new_data)(ci) - (*dh_dT_data)(ci) * (*T_new_data)(ci)) -
                                  (*rho_current_data)(ci) * (*h_current_data)(ci));
            }
        }
    }
    return;
} // addTemporalAndLinearTermstoRHSOfEnergyEquation

void
EnthalpyHierarchyIntegrator::computeDivergenceVelocitySourceTerm(int Div_U_F_idx, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());

    // Filling ghost cells for temperature.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> T_transaction_comps(1);
    T_transaction_comps[0] = InterpolationTransactionComponent(d_T_scratch_idx,
                                                               d_T_new_idx,
                                                               "CONSERVATIVE_LINEAR_REFINE",
                                                               false,
                                                               "CONSERVATIVE_COARSEN",
                                                               "LINEAR",
                                                               false,
                                                               d_T_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(T_transaction_comps, d_hierarchy);
    hier_bdry_fill->fillData(new_time);

    // Compute gradient of temperature.
    d_hier_math_ops->grad(d_grad_T_idx, d_grad_T_var, true, 1.0, d_T_scratch_idx, d_T_var, nullptr, new_time);

    // Compute k*grad_T.
    d_hier_sc_data_ops->multiply(d_grad_T_idx, d_grad_T_idx, d_T_diffusion_coef_new_idx);

    // Compute div(k*grad_T).
    d_hier_math_ops->div(Div_U_F_idx, d_Div_U_F_var, 1.0, d_grad_T_idx, d_grad_T_var, nullptr, new_time, false);

    const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(d_h_new_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(d_rho_new_idx);
            Pointer<CellData<NDIM, double> > Div_U_F_data = patch->getPatchData(Div_U_F_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                double material_derivative = 0.0;
                if ((*h_data)(ci) >= h_s && (*h_data)(ci) <= h_l && (*H_data)(ci) >= H_LIM)
                {
                    const double denominator = (*rho_data)(ci)*std::pow(
                        (*h_data)(ci) * (d_rho_liquid - d_rho_solid) - d_rho_liquid * h_l + d_rho_solid * h_s, 2.0);

                    material_derivative = (*Div_U_F_data)(ci)*d_rho_solid * d_rho_liquid * (h_l - h_s) /
                                          denominator; // div k grad T rho_s*rho_l (h_l - h_s) / denominator
                }

                (*Div_U_F_data)(ci) =
                    -(d_rho_liquid - d_rho_solid) * material_derivative * (*H_data)(ci) / (*rho_data)(ci);
            }
        }
    }

    return;
} // computeDivergenceVelocitySourceTerm

void
EnthalpyHierarchyIntegrator::computeEnthalpyBasedOnTemperature(int h_idx,
                                                               const int T_idx,
                                                               const int rho_idx,
                                                               const int lf_idx,
                                                               const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(h_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if ((*H_data)(ci) >= H_LIM)
                {
                    if ((*T_data)(ci) < d_solidus_temperature)
                    {
                        (*h_data)(ci) = d_specific_heat_solid * ((*T_data)(ci)-d_reference_temperature);
                    }
                    else if ((*T_data)(ci) >= d_solidus_temperature && (*T_data)(ci) <= d_liquidus_temperature)
                    {
                        (*h_data)(ci) = d_specific_heat_mushy * ((*T_data)(ci)-d_solidus_temperature) + h_s +
                                        (*lf_data)(ci)*d_rho_liquid * d_latent_heat / (*rho_data)(ci);
                    }
                    else
                    {
                        (*h_data)(ci) = d_specific_heat_liquid * ((*T_data)(ci)-d_liquidus_temperature) + h_l;
                    }
                }
                else
                {
                    (*h_data)(ci) = d_specific_heat_gas * ((*T_data)(ci)-d_reference_temperature);
                }
            }
        }
    }
    return;
} // computeEnthalpyBasedOnTemperature

void
EnthalpyHierarchyIntegrator::computeTemperatureBasedOnEnthalpy(int T_idx, const int h_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(h_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if ((*H_data)(ci) >= H_LIM)
                {
                    if ((*h_data)(ci) < h_s)
                    {
                        (*T_data)(ci) = (*h_data)(ci) / d_specific_heat_solid + d_reference_temperature;
                    }
                    else if ((*h_data)(ci) >= h_s && (*h_data)(ci) <= h_l)
                    {
                        (*T_data)(ci) = d_solidus_temperature + ((*h_data)(ci)-h_s) / (h_l - h_s) *
                                                                    (d_liquidus_temperature - d_solidus_temperature);
                    }
                    else
                    {
                        (*T_data)(ci) = d_liquidus_temperature + ((*h_data)(ci)-h_l) / d_specific_heat_liquid;
                    }
                }
                else
                {
                    (*T_data)(ci) = (*h_data)(ci) / d_specific_heat_gas + d_reference_temperature;
                }
            }
        }
    }
    return;
} // computeTemperatureBasedOnEnthalpy

void
EnthalpyHierarchyIntegrator::updateEnthalpy(int h_new_idx, const int T_new_idx, const int T_pre_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > h_new_data = patch->getPatchData(h_new_idx);
            Pointer<CellData<NDIM, double> > T_new_data = patch->getPatchData(T_new_idx);
            Pointer<CellData<NDIM, double> > T_pre_data = patch->getPatchData(T_pre_idx);
            Pointer<CellData<NDIM, double> > dh_dT_data = patch->getPatchData(d_dh_dT_scratch_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*h_new_data)(ci) += (*dh_dT_data)(ci) * ((*T_new_data)(ci) - (*T_pre_data)(ci));
            }
        }
    }
    return;
} // updateEnthalpy

void
EnthalpyHierarchyIntegrator::computeEnthalpyDerivative(int dh_dT_idx, const int T_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            Pointer<CellData<NDIM, double> > dh_dT_data = patch->getPatchData(dh_dT_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if ((*H_data)(ci) >= H_LIM)
                {
                    if ((*T_data)(ci) < d_solidus_temperature)
                    {
                        (*dh_dT_data)(ci) = d_specific_heat_solid;
                    }
                    else if ((*T_data)(ci) >= d_solidus_temperature && (*T_data)(ci) <= d_liquidus_temperature)
                    {
                        (*dh_dT_data)(ci) =
                            d_specific_heat_mushy + d_latent_heat / (d_liquidus_temperature - d_solidus_temperature);
                    }
                    else
                    {
                        (*dh_dT_data)(ci) = d_specific_heat_liquid;
                    }
                }
                else
                {
                    (*dh_dT_data)(ci) = d_specific_heat_gas;
                }
            }
        }
    }
    return;
} // computeEnthalpyDerivative

void
EnthalpyHierarchyIntegrator::computeLiquidFraction(int lf_idx, const int h_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(h_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if ((*H_data)(ci) >= H_LIM)
                {
                    if ((*h_data)(ci) < h_s)
                    {
                        (*lf_data)(ci) = 0.0;
                    }
                    else if ((*h_data)(ci) > h_l)
                    {
                        (*lf_data)(ci) = 1.0;
                    }
                    else
                    {
                        (*lf_data)(ci) =
                            d_rho_solid * (h_s - (*h_data)(ci)) /
                            ((d_rho_liquid - d_rho_solid) * (*h_data)(ci)-d_rho_liquid * h_l + d_rho_solid * h_s);
                    }
                }
                else
                {
                    (*lf_data)(ci) = d_gas_liquid_fraction;
                }
            }
        }
    }
    return;
} // computeLiquidFraction

void
EnthalpyHierarchyIntegrator::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db);
#endif
    // Read in data members from input database.
    if (!is_from_restart)
    {
        d_specific_heat_liquid = input_db->getDouble("specific_heat_liquid");
        d_specific_heat_solid = input_db->getDouble("specific_heat_solid");
        d_specific_heat_gas = input_db->getDouble("specific_heat_gas");
        d_liquidus_temperature = input_db->getDouble("liquidus_temperature");
        d_solidus_temperature = input_db->getDouble("solidus_temperature");
        d_reference_temperature = input_db->getDouble("solidus_temperature");

        if (input_db->keyExists("specific_heat_mushy"))
            d_specific_heat_mushy = input_db->getDouble("specific_heat_mushy");
        else
            d_specific_heat_mushy = 0.5 * (d_specific_heat_liquid + d_specific_heat_solid);
        if (input_db->keyExists("gas_liquid_fraction"))
            d_gas_liquid_fraction = input_db->getDouble("gas_liquid_fraction");
        if (input_db->keyExists("max_inner_iterations"))
            d_max_inner_iterations = input_db->getInteger("max_inner_iterations");
        if (input_db->keyExists("lf_iteration_error_tolerance"))
            d_lf_iteration_error_tolerance = input_db->getDouble("lf_iteration_error_tolerance");
    }
    return;
} // getFromInput

void
EnthalpyHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }

    d_specific_heat_liquid = db->getDouble("specific_heat_liquid");
    d_specific_heat_solid = db->getDouble("specific_heat_solid");
    d_specific_heat_gas = db->getDouble("specific_heat_gas");
    d_liquidus_temperature = db->getDouble("liquidus_temperature");
    d_solidus_temperature = db->getDouble("solidus_temperature");
    d_reference_temperature = db->getDouble("reference_temperature");
    d_specific_heat_mushy = db->getDouble("specific_heat_mushy");
    d_gas_liquid_fraction = db->getDouble("gas_liquid_fraction");
    d_max_inner_iterations = db->getInteger("max_inner_iterations");
    d_lf_iteration_error_tolerance = db->getDouble("lf_iteration_error_tolerance");

    return;
} // getFromRestart

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
