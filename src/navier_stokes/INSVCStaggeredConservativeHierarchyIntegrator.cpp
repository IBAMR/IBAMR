// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
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

#include "IBAMR_config.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/BrinkmanPenalizationStrategy.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredPressureBcCoef.h"
#include "ibamr/INSVCStaggeredVelocityBcCoef.h"
#include "ibamr/PETScKrylovStaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/VCStaggeredStokesOperator.h"
#include "ibamr/VCStaggeredStokesProjectionPreconditioner.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartSideDoubleDivPreservingRefine.h"
#include "ibtk/CartSideDoubleRT0Refine.h"
#include "ibtk/CartSideDoubleSpecializedLinearRefine.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/VCSCViscousOpPointRelaxationFACOperator.h"
#include "ibtk/VCSCViscousOperator.h"
#include "ibtk/VCSCViscousPETScLevelSolver.h"
#include "ibtk/ibtk_enums.h"

#include "ArrayData.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "EdgeData.h"
#include "EdgeVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "NodeData.h"
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "PoissonSpecifications.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghost cells used for each variable quantity.
static const int SIDEG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredConservativeHierarchyIntegrator::INSVCStaggeredConservativeHierarchyIntegrator(std::string object_name,
                                                                                             Pointer<Database> input_db,
                                                                                             bool register_for_restart)
    : INSVCStaggeredHierarchyIntegrator(std::move(object_name), input_db, register_for_restart),
      d_rho_sc_bc_coefs(NDIM, nullptr)
{
    if (!(d_convective_difference_form == CONSERVATIVE))
    {
        TBOX_ERROR(d_object_name << "::INSVCStaggeredConservativeHierarchyIntegrator():\n"
                                 << " variable coefficient discretization\n"
                                 << " requires CONSERVATIVE convective difference form\n");
    }

    if (d_rho_is_const)
    {
        TBOX_ERROR(d_object_name << "::INSVCStaggeredConservativeHierarchyIntegrator():\n"
                                 << " conservative variable coefficient discretization\n"
                                 << " requires non-constant density\n");
    }

    // Initialize conservative mass and momentum integrator.
    d_rho_p_integrator = new INSVCStaggeredConservativeMassMomentumIntegrator(
        "INSVCStaggeredConservativeHierarchyIntegrator::MassMomentumIntegrator",
        input_db->getDatabase("mass_momentum_integrator_db"));
    d_convective_op_type = "NONE";

    // Side centered state variable for density and interpolated density variable
    // for plotting.
    d_rho_sc_var = new SideVariable<NDIM, double>(d_object_name + "::rho_sc");
    d_rho_interp_cc_var = new CellVariable<NDIM, double>(d_object_name + "::rho_interp_cc", NDIM);

    return;
} // INSVCStaggeredConservativeHierarchyIntegrator

void
INSVCStaggeredConservativeHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    INSVCStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Get the density variable, which must be side-centered and maintained by the
    // INS integrator for
    // this conservative discretization form.
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> no_ghosts = 0;
    if (INSVCStaggeredHierarchyIntegrator::d_rho_var)
    {
        d_rho_sc_var = INSVCStaggeredHierarchyIntegrator::d_rho_var;
#if !defined(NDEBUG)
        TBOX_ASSERT(d_rho_sc_var);
        TBOX_ASSERT(d_rho_init_fcn || d_reset_rho_fcns.size() > 0);
#endif
        registerVariable(d_rho_sc_current_idx,
                         d_rho_sc_new_idx,
                         d_rho_sc_scratch_idx,
                         d_rho_sc_var,
                         side_ghosts,
                         d_rho_coarsen_type,
                         d_rho_refine_type,
                         d_rho_init_fcn);
    }
    else
    {
        TBOX_ERROR(
            "INSVCStaggeredConservativeHierarchyIntegrator::"
            "initializeHierarchyIntegrator():\n"
            << "  rho_is_const == false but no mass density variable has "
               "been registered.\n");
    }

    // Register variables for plotting.
    registerVariable(d_rho_interp_cc_idx, d_rho_interp_cc_var, no_ghosts, getCurrentContext());
    if (d_visit_writer)
    {
        if (d_output_rho)
        {
            d_visit_writer->registerPlotQuantity("rho_ins_cc", "VECTOR", d_rho_interp_cc_idx, 0, d_rho_scale);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0)
                    d_visit_writer->registerPlotQuantity("rho_x", "SCALAR", d_rho_interp_cc_idx, d, d_rho_scale);
                if (d == 1)
                    d_visit_writer->registerPlotQuantity("rho_y", "SCALAR", d_rho_interp_cc_idx, d, d_rho_scale);
                if (d == 2)
                    d_visit_writer->registerPlotQuantity("rho_z", "SCALAR", d_rho_interp_cc_idx, d, d_rho_scale);
            }
        }
    }

    // Set various objects with conservative time integrator.
    d_rho_p_integrator->setSideCenteredVelocityBoundaryConditions(d_U_bc_coefs);
    d_rho_p_integrator->setSideCenteredDensityBoundaryConditions(d_rho_sc_bc_coefs);
    if (d_S_fcn) d_rho_p_integrator->setMassDensitySourceTerm(d_S_fcn);

    return;
} // initializeHierarchyIntegrator

void
INSVCStaggeredConservativeHierarchyIntegrator::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                        Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    INSVCStaggeredHierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);
    return;
} // initializePatchHierarhcy

void
INSVCStaggeredConservativeHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                            const double new_time,
                                                                            const int num_cycles)
{
    INSVCStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);
    const double dt = new_time - current_time;

    // Get the current value of density and store it in scratch.  Note that
    // for conservative differencing, it is required that the side-centered
    // density is a state variable maintained by this integrator.
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_sc_var);
    TBOX_ASSERT(d_rho_sc_current_idx >= 0);
    TBOX_ASSERT(d_rho_sc_scratch_idx >= 0);
#endif

    // Note that we always reset current context of state variables here, if
    // necessary.
    const double apply_time = current_time;
    for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
    {
        d_reset_rho_fcns[k](d_rho_sc_current_idx,
                            d_rho_sc_var,
                            d_hier_math_ops,
                            -1 /*cycle_num*/,
                            apply_time,
                            current_time,
                            new_time,
                            d_reset_rho_fcns_ctx[k]);
    }
    d_hier_sc_data_ops->copyData(d_rho_sc_scratch_idx,
                                 d_rho_sc_current_idx,
                                 /*interior_only*/ true);

    if (!d_mu_is_const && d_mu_var)
    {
        for (unsigned k = 0; k < d_reset_mu_fcns.size(); ++k)
        {
            d_reset_mu_fcns[k](d_mu_current_idx,
                               d_mu_var,
                               d_hier_math_ops,
                               -1 /*cycle_num*/,
                               apply_time,
                               current_time,
                               new_time,
                               d_reset_mu_fcns_ctx[k]);
        }
    }

    // Get the current value of viscosity
    if (!d_mu_is_const)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int mu_current_idx;
        if (d_adv_diff_hier_integrator && d_mu_adv_diff_var)
        {
            mu_current_idx = var_db->mapVariableAndContextToIndex(d_mu_adv_diff_var,
                                                                  d_adv_diff_hier_integrator->getCurrentContext());
        }
        else
        {
            mu_current_idx = d_mu_current_idx;
        }

        d_hier_cc_data_ops->copyData(d_mu_scratch_idx,
                                     mu_current_idx,
                                     /*interior_only*/ true);
        d_mu_bdry_bc_fill_op->fillData(current_time);

        // Interpolate onto node or edge centers
        if (d_mu_vc_interp_type == VC_AVERAGE_INTERP)
        {
            d_hier_math_ops->interp(d_mu_interp_idx,
                                    d_mu_interp_var,
                                    /*dst_ghost_interp*/ true,
                                    d_mu_scratch_idx,
                                    d_mu_var,
                                    d_no_fill_op,
                                    current_time);
        }
        else if (d_mu_vc_interp_type == VC_HARMONIC_INTERP)
        {
            d_hier_math_ops->harmonic_interp(d_mu_interp_idx,
                                             d_mu_interp_var,
                                             /*dst_ghost_interp*/ true,
                                             d_mu_scratch_idx,
                                             d_mu_var,
                                             d_no_fill_op,
                                             current_time);
        }
        else
        {
            TBOX_ERROR("this statement should not be reached");
        }

        // Store the viscosities for later use
        d_hier_cc_data_ops->copyData(d_mu_linear_op_idx,
                                     d_mu_scratch_idx,
                                     /*interior_only*/ false);
#if (NDIM == 2)
        d_hier_nc_data_ops->copyData(d_mu_interp_linear_op_idx,
                                     d_mu_interp_idx,
                                     /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->copyData(d_mu_interp_linear_op_idx,
                                     d_mu_interp_idx,
                                     /*interior_only*/ false);
#endif
    }

    // Allocate solver vectors.
    d_U_rhs_vec->allocateVectorData(current_time);
    d_U_rhs_vec->setToScalar(0.0);
    d_P_rhs_vec->allocateVectorData(current_time);
    d_P_rhs_vec->setToScalar(0.0);
    if (!d_creeping_flow)
    {
        d_U_adv_vec->allocateVectorData(current_time);
        d_U_adv_vec->setToScalar(0.0);
        d_N_vec->allocateVectorData(current_time);
        d_N_vec->setToScalar(0.0);
    }

    // Cache BC data.
    d_bc_helper->cacheBcCoefData(d_bc_coefs, new_time, d_hierarchy);

    // Compute viscous right-hand side terms
    const double mu = d_mu_is_const ? d_problem_coefs.getMu() : -1.0;
    const double lambda = d_problem_coefs.getLambda();
    double K_rhs = 0.0;
    switch (d_viscous_time_stepping_type)
    {
    case BACKWARD_EULER:
        K_rhs = 0.0;
        break;
    case FORWARD_EULER:
        K_rhs = 1.0;
        break;
    case TRAPEZOIDAL_RULE:
        K_rhs = 0.5;
        break;
    default:
        TBOX_ERROR("this statement should not be reached");
    }

    PoissonSpecifications U_rhs_problem_coefs(d_object_name + "::U_rhs_problem_coefs");
    U_rhs_problem_coefs.setCConstant(-K_rhs * lambda);

    // rhs_D_{ec,nc} = K * mu
    if (d_mu_is_const)
    {
#if (NDIM == 2)
        d_hier_nc_data_ops->setToScalar(d_velocity_rhs_D_idx,
                                        +K_rhs * mu,
                                        /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->setToScalar(d_velocity_rhs_D_idx,
                                        +K_rhs * mu,
                                        /*interior_only*/ false);
#endif
    }
    else
    {
#if (NDIM == 2)
        d_hier_nc_data_ops->scale(d_velocity_rhs_D_idx,
                                  +K_rhs,
                                  d_mu_interp_idx,
                                  /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->scale(d_velocity_rhs_D_idx,
                                  +K_rhs,
                                  d_mu_interp_idx,
                                  /*interior_only*/ false);
#endif
    }
    U_rhs_problem_coefs.setDPatchDataId(d_velocity_rhs_D_idx);

    const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
    const Pointer<SideVariable<NDIM, double> > U_rhs_var = d_U_rhs_vec->getComponentVariable(0);
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_U_bc_coefs,
                                                              /*P_bc_coef*/ nullptr,
                                                              d_U_scratch_idx,
                                                              /*P_data_idx*/ -1,
                                                              /*homogeneous_bc*/ false);
    d_U_bdry_bc_fill_op->fillData(current_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs,
                                                              /*P_bc_coef*/ nullptr);
    d_bc_helper->enforceDivergenceFreeConditionAtBoundary(d_U_scratch_idx);
    // RHS^n = (C_rhs*I + L(D_rhs))*U^n
    d_hier_math_ops->vc_laplace(U_rhs_idx,
                                U_rhs_var,
                                1.0,
                                0.0,
                                U_rhs_problem_coefs.getDPatchDataId(),
#if (NDIM == 2)
                                Pointer<NodeVariable<NDIM, double> >(nullptr),
#elif (NDIM == 3)
                                Pointer<EdgeVariable<NDIM, double> >(nullptr),
#endif
                                d_U_scratch_idx,
                                d_U_var,
                                d_no_fill_op,
                                current_time,
                                d_mu_vc_interp_type);

    // Add the momentum portion of the RHS in the case of conservative
    // discretization form
    // RHS^n = RHS^n + 1/dt*(rho*U)^n
    d_hier_sc_data_ops->multiply(d_temp_sc_idx, d_rho_sc_scratch_idx, d_U_scratch_idx, /*interior_only*/ true);
    d_hier_sc_data_ops->axpy(U_rhs_idx,
                             1.0 / dt,
                             d_temp_sc_idx,
                             U_rhs_idx,
                             /*interior_only*/ true);
    d_hier_sc_data_ops->copyData(d_U_src_idx,
                                 d_U_scratch_idx,
                                 /*interior_only*/ false);

    // Set the initial guess.
    d_hier_sc_data_ops->copyData(d_U_new_idx, d_U_current_idx);
    d_hier_cc_data_ops->copyData(d_P_new_idx, d_P_current_idx);

    // Set up inhomogeneous BCs.
    d_stokes_solver->setHomogeneousBc(false);

    // Initialize any registered advection-diffusion solver.
    if (d_adv_diff_hier_integrator)
    {
        const int adv_diff_num_cycles = d_adv_diff_hier_integrator->getNumberOfCycles();
        if (adv_diff_num_cycles != d_current_num_cycles && d_current_num_cycles != 1)
        {
            TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                     << "  attempting to perform " << d_current_num_cycles
                                     << " cycles of fixed point iteration.\n"
                                     << "  number of cycles required by coupled advection-diffusion "
                                        "solver = "
                                     << adv_diff_num_cycles << ".\n"
                                     << "  current implementation requires either that both solvers use "
                                        "the same "
                                        "number of cycles,\n"
                                     << "  or that the Navier-Stokes solver use only a single cycle.\n");
        }
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_current_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        if (isAllocatedPatchData(U_adv_diff_current_idx))
        {
            INSVCStaggeredHierarchyIntegrator::copySideToFace(U_adv_diff_current_idx, d_U_current_idx, d_hierarchy);
        }
        d_adv_diff_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, adv_diff_num_cycles);
        const int U_adv_diff_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getScratchContext());
        if (isAllocatedPatchData(U_adv_diff_scratch_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_scratch_idx, U_adv_diff_current_idx);
        }
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_new_idx))
        {
            d_hier_fc_data_ops->copyData(U_adv_diff_new_idx, U_adv_diff_current_idx);
        }
    }

    // Account for the convective acceleration term.
    if (!d_creeping_flow)
    {
        d_rho_p_integrator->setSolutionTime(current_time);
        d_rho_p_integrator->setTimeInterval(current_time, new_time);

        // For conservative momentum discretization, an approximation to rho^{n+1}
        // will be computed from rho^{n}, which requires additional options to be
        // set.

        // Set the rho^{n} density
        d_rho_p_integrator->setSideCenteredDensityPatchDataIndex(d_rho_sc_current_idx);

        // Set the convective derivative patch data index.
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        d_rho_p_integrator->setSideCenteredConvectiveDerivativePatchDataIndex(N_idx);

        // Data for the conservative time integrator is for cycle 0
        const int ins_cycle_num = 0;
        d_rho_p_integrator->setCycleNumber(ins_cycle_num);

        // Set the velocities used to update the density and the previous time step
        // size
        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ -1, /*current*/ d_U_current_idx, /*new*/ -1);
        }
        else
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ d_U_old_current_idx, /*current*/ d_U_current_idx, /*new*/ -1);
            d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        // Integrate density and convective momentum.
        d_rho_p_integrator->integrate(dt);
    }

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
INSVCStaggeredConservativeHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                                  const double new_time,
                                                                  const int cycle_num)
{
    INSVCStaggeredHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    // Get the coarsest and finest level numbers.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Update the state variables of any linked advection-diffusion solver.
    // NOTE: This also updates rho and mu if they are maintained by adv-diff
    // integrator.
    if (d_adv_diff_hier_integrator)
    {
        d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
    }

    // Update viscosity if it is maintained by the fluid integrator.
    const double apply_time = new_time;
    if (!d_mu_is_const && d_mu_var)
    {
        for (unsigned k = 0; k < d_reset_mu_fcns.size(); ++k)
        {
            d_reset_mu_fcns[k](d_mu_new_idx,
                               d_mu_var,
                               d_hier_math_ops,
                               cycle_num,
                               apply_time,
                               current_time,
                               new_time,
                               d_reset_mu_fcns_ctx[k]);
        }
    }

    if (!d_mu_is_const)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int mu_new_idx;
        if (d_adv_diff_hier_integrator && d_mu_adv_diff_var)
        {
            mu_new_idx =
                var_db->mapVariableAndContextToIndex(d_mu_adv_diff_var, d_adv_diff_hier_integrator->getNewContext());
        }
        else
        {
            mu_new_idx = d_mu_new_idx;
        }
        d_hier_cc_data_ops->copyData(d_mu_scratch_idx,
                                     mu_new_idx,
                                     /*interior_only*/ true);
        d_mu_bdry_bc_fill_op->fillData(new_time);

        // Interpolate onto node or edge centers
        if (d_mu_vc_interp_type == VC_AVERAGE_INTERP)
        {
            d_hier_math_ops->interp(d_mu_interp_idx,
                                    d_mu_interp_var,
                                    /*dst_ghost_interp*/ true,
                                    d_mu_scratch_idx,
                                    d_mu_var,
                                    d_no_fill_op,
                                    new_time);
        }
        else if (d_mu_vc_interp_type == VC_HARMONIC_INTERP)
        {
            d_hier_math_ops->harmonic_interp(d_mu_interp_idx,
                                             d_mu_interp_var,
                                             /*dst_ghost_interp*/ true,
                                             d_mu_scratch_idx,
                                             d_mu_var,
                                             d_no_fill_op,
                                             new_time);
        }
        else
        {
            TBOX_ERROR("this statement should not be reached");
        }

        // Store the viscosities for later use
        d_hier_cc_data_ops->copyData(d_mu_linear_op_idx,
                                     d_mu_scratch_idx,
                                     /*interior_only*/ false);
#if (NDIM == 2)
        d_hier_nc_data_ops->copyData(d_mu_interp_linear_op_idx,
                                     d_mu_interp_idx,
                                     /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->copyData(d_mu_interp_linear_op_idx,
                                     d_mu_interp_idx,
                                     /*interior_only*/ false);
#endif
    }

    // In the special case of a conservative discretization form, the updated
    // density is calculated by application of the mass and convective
    // momentum integrator.
    if (!d_creeping_flow)
    {
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        // Update N_idx if necessary
        if (cycle_num > 0)
        {
            const double dt = new_time - current_time;
            const double half_time = current_time + 0.5 * dt;
            d_rho_p_integrator->setSolutionTime(half_time);

            // Set the cycle number
            d_rho_p_integrator->setCycleNumber(cycle_num);

            // Set the patch data index for convective derivative.
            d_rho_p_integrator->setSideCenteredConvectiveDerivativePatchDataIndex(N_idx);

            // Always set to current because we want to update rho^{n} to rho^{n+1}
            d_rho_p_integrator->setSideCenteredDensityPatchDataIndex(d_rho_sc_current_idx);

            // Set the velocities used to update the density
            if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ -1, /*current*/ d_U_current_idx, /*new*/ d_U_new_idx);
            }
            else
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ d_U_old_current_idx,
                    /*current*/ d_U_current_idx,
                    /*new*/ d_U_new_idx);
                d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
            }

            d_rho_p_integrator->integrate(dt);
        }

        // Set the full convective term depending on the time stepping type
        d_hier_sc_data_ops->copyData(d_N_full_idx, N_idx, /*interior_only*/ true);
    }
    const int rho_sc_new_idx = d_rho_p_integrator->getUpdatedSideCenteredDensityPatchDataIndex();
    d_hier_sc_data_ops->copyData(d_rho_sc_new_idx,
                                 rho_sc_new_idx,
                                 /*interior_only*/ true);

    // Copy new into scratch
    d_hier_sc_data_ops->copyData(d_rho_sc_scratch_idx, d_rho_sc_new_idx);

    // Synchronize newest density
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent rho_scratch_synch_transaction =
        SynchronizationTransactionComponent(d_rho_sc_scratch_idx, d_rho_coarsen_type);
    d_side_synch_op->resetTransactionComponent(rho_scratch_synch_transaction);
    d_side_synch_op->synchronizeData(d_integrator_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, d_U_coarsen_type);
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);

    // Store the density for later use
    d_hier_sc_data_ops->copyData(d_rho_linear_op_idx,
                                 d_rho_sc_scratch_idx,
                                 /*interior_only*/ true);

    // Compute the Brinkman contribution to the velocity operator.
    d_hier_sc_data_ops->setToScalar(d_velocity_L_idx, 0.0);
    for (auto& brinkman_force : d_brinkman_force)
    {
        brinkman_force->demarcateBrinkmanZone(d_velocity_L_idx, new_time, cycle_num);
    }

    // Synchronize Brinkman coefficients.
    SynchronizationTransactionComponent L_synch_transaction =
        SynchronizationTransactionComponent(d_velocity_L_idx, "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(L_synch_transaction);
    d_side_synch_op->synchronizeData(d_integrator_time);
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);

    // Update the solvers and operators to take into account new state variables
    updateOperatorsAndSolvers(current_time, new_time, cycle_num);

    // Compute the Brinkman velocity for the RHS vector.
    d_hier_sc_data_ops->setToScalar(d_velocity_L_idx, 0.0);
    for (auto& brinkman_force : d_brinkman_force)
    {
        brinkman_force->computeBrinkmanVelocity(d_velocity_L_idx, new_time, cycle_num);
    }
    d_side_synch_op->resetTransactionComponent(L_synch_transaction);
    d_side_synch_op->synchronizeData(d_integrator_time);
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);

    // Setup the solution and right-hand-side vectors.
    setupSolverVectors(d_sol_vec, d_rhs_vec, current_time, new_time, cycle_num);

    // Scale rhs if necessary
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hier_sc_data_ops->resetLevels(ln, ln);
        const double A_scale = d_A_scale[ln];
        if (!MathUtilities<double>::equalEps(A_scale, 1.0))
        {
            d_hier_sc_data_ops->scale(d_rhs_vec->getComponentDescriptorIndex(0),
                                      A_scale,
                                      d_rhs_vec->getComponentDescriptorIndex(0),
                                      /*interior_only*/ true);
        }
        d_hier_sc_data_ops->resetLevels(coarsest_ln, finest_ln);
    }

    // Solve for u(n+1), p(n+1/2).
    d_stokes_solver->solveSystem(*d_sol_vec, *d_rhs_vec);

    // Unscale rhs if necessary
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hier_sc_data_ops->resetLevels(ln, ln);
        const double A_scale = d_A_scale[ln];
        if (!MathUtilities<double>::equalEps(A_scale, 1.0))
        {
            d_hier_sc_data_ops->scale(d_rhs_vec->getComponentDescriptorIndex(0),
                                      1.0 / A_scale,
                                      d_rhs_vec->getComponentDescriptorIndex(0),
                                      /*interior_only*/ true);
        }
        d_hier_sc_data_ops->resetLevels(coarsest_ln, finest_ln);
    }

    if (d_enable_logging && d_enable_logging_solver_iterations)
        plog << d_object_name
             << "::integrateHierarchy(): stokes solve number of iterations = " << d_stokes_solver->getNumIterations()
             << "\n";
    if (d_enable_logging)
        plog << d_object_name
             << "::integrateHierarchy(): stokes solve residual norm        = " << d_stokes_solver->getResidualNorm()
             << "\n";
    if (d_explicitly_remove_nullspace) removeNullSpace(d_sol_vec);

    // Reset the solution and right-hand-side vectors.
    resetSolverVectors(d_sol_vec, d_rhs_vec, current_time, new_time, cycle_num);

    // Update the state variables of any linked advection-diffusion solver.
    if (d_adv_diff_hier_integrator)
    {
        // Update the advection velocities used by the advection-diffusion
        // solver.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int U_adv_diff_new_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getNewContext());
        if (isAllocatedPatchData(U_adv_diff_new_idx))
        {
            INSVCStaggeredHierarchyIntegrator::copySideToFace(U_adv_diff_new_idx, d_U_new_idx, d_hierarchy);
        }
        const int U_adv_diff_current_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getCurrentContext());
        const int U_adv_diff_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_U_adv_diff_var, d_adv_diff_hier_integrator->getScratchContext());
        if (isAllocatedPatchData(U_adv_diff_scratch_idx))
        {
            d_hier_fc_data_ops->linearSum(U_adv_diff_scratch_idx, 0.5, U_adv_diff_current_idx, 0.5, U_adv_diff_new_idx);
        }

        // Update the state variables maintained by the advection-diffusion
        // solver.
        //
        // NOTE: We already performed cycle 0 above.
        const int adv_diff_num_cycles = d_adv_diff_hier_integrator->getNumberOfCycles();
        if (d_current_num_cycles != adv_diff_num_cycles)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(d_current_num_cycles == 1);
#endif
            for (int adv_diff_cycle_num = 1; adv_diff_cycle_num < adv_diff_num_cycles; ++adv_diff_cycle_num)
            {
                d_adv_diff_hier_integrator->integrateHierarchy(current_time, new_time, adv_diff_cycle_num);
            }
        }
    }

    // Re-update viscosity if it is maintained by the integrator
    // using the newest available data from INS and advection-diffusion solvers
    if (d_current_num_cycles == cycle_num + 1)
    {
        if (!d_mu_is_const && d_mu_var)
        {
            for (unsigned k = 0; k < d_reset_mu_fcns.size(); ++k)
            {
                d_reset_mu_fcns[k](d_mu_new_idx,
                                   d_mu_var,
                                   d_hier_math_ops,
                                   cycle_num,
                                   apply_time,
                                   current_time,
                                   new_time,
                                   d_reset_mu_fcns_ctx[k]);
            }
        }
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
INSVCStaggeredConservativeHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                             const double new_time,
                                                                             const bool skip_synchronize_new_state_data,
                                                                             const int num_cycles)
{
    INSVCStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    return;
} // postprocessIntegrateHierarchy

void
INSVCStaggeredConservativeHierarchyIntegrator::removeNullSpace(const Pointer<SAMRAIVectorReal<NDIM, double> >& sol_vec)
{
    INSVCStaggeredHierarchyIntegrator::removeNullSpace(sol_vec);
    return;
} // removeNullSpace

void
INSVCStaggeredConservativeHierarchyIntegrator::registerMassDensityBoundaryConditions(
    RobinBcCoefStrategy<NDIM>* rho_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    std::vector<RobinBcCoefStrategy<NDIM>*> rho_sc_bc_coefs(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) rho_sc_bc_coefs[d] = rho_bc_coef;
    d_rho_sc_bc_coefs = rho_sc_bc_coefs;
    return;
} // registerMassDensityBoundaryConditions

void
INSVCStaggeredConservativeHierarchyIntegrator::registerMassDensityBoundaryConditions(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& rho_sc_bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(rho_sc_bc_coefs.size() == NDIM);
#endif
    d_rho_sc_bc_coefs = rho_sc_bc_coefs;
    return;
} // registerMassDensityBoundaryConditions

void
INSVCStaggeredConservativeHierarchyIntegrator::registerMassDensitySourceTerm(Pointer<CartGridFunction> S_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (!d_S_fcn)
    {
        d_S_fcn = S_fcn;
    }
    else
    {
        TBOX_ERROR(d_object_name << "::INSVCStaggeredConservativeHierarchyIntegrator():\n"
                                 << " present implementation allows for only one mass density source\n"
                                 << " term to be set. Consider combining source terms into single "
                                    "CartGridFunction.\n");
    }
    return;
} // registerMassDensitySourceTerm

int
INSVCStaggeredConservativeHierarchyIntegrator::getNumberOfCycles() const
{
    return d_num_cycles;
} // getNumberOfCycles

Pointer<ConvectiveOperator>
INSVCStaggeredConservativeHierarchyIntegrator::getConvectiveOperator()
{
    d_convective_op.setNull();
    d_convective_op_needs_init = false;
    return d_convective_op;
} // getConvectiveOperator

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSVCStaggeredConservativeHierarchyIntegrator::initializeLevelDataSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    INSVCStaggeredHierarchyIntegrator::initializeLevelDataSpecialized(
        base_hierarchy, level_number, init_data_time, can_be_refined, initial_time, base_old_level, allocate_data);
    return;
} // initializeLevelDataSpecialized

void
INSVCStaggeredConservativeHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    INSVCStaggeredHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
        base_hierarchy, coarsest_level, finest_level);
    d_rho_p_integrator->setHierarchyMathOps(d_hier_math_ops);
    d_rho_p_integrator->initializeTimeIntegrator(base_hierarchy);
    return;
} // resetHierarchyConfigurationSpecialized

void
INSVCStaggeredConservativeHierarchyIntegrator::applyGradientDetectorSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    INSVCStaggeredHierarchyIntegrator::applyGradientDetectorSpecialized(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    return;
} // applyGradientDetectorSpecialized

void
INSVCStaggeredConservativeHierarchyIntegrator::setupPlotDataSpecialized()
{
    INSVCStaggeredHierarchyIntegrator::setupPlotDataSpecialized();

    Pointer<VariableContext> ctx = getCurrentContext();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    static const bool synch_cf_interface = true;

    // Put interpolated rho onto cell centers for plotting purposes.
    if (d_output_rho)
    {
        const int rho_sc_idx = var_db->mapVariableAndContextToIndex(d_rho_sc_var, ctx);
        const int rho_cc_idx = var_db->mapVariableAndContextToIndex(d_rho_interp_cc_var, ctx);
        d_hier_math_ops->interp(rho_cc_idx,
                                d_rho_interp_cc_var,
                                rho_sc_idx,
                                d_rho_sc_var,
                                d_no_fill_op,
                                d_integrator_time,
                                synch_cf_interface);
    }

    return;
} // setupPlotDataSpecialized

void
INSVCStaggeredConservativeHierarchyIntegrator::regridProjection()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Setup the solver vectors.
    SAMRAIVectorReal<NDIM, double> sol_vec(d_object_name + "::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_P_var, d_P_scratch_idx, wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM, double> rhs_vec(d_object_name + "::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_Div_U_var, d_Div_U_idx, wgt_cc_idx, d_hier_cc_data_ops);

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_U_scratch_idx);
    scratch_idxs.setFlag(d_P_scratch_idx);
    scratch_idxs.setFlag(d_pressure_D_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, d_integrator_time);
    }

    // Setup the regrid Poisson solver.
    Pointer<PoissonSolver> regrid_projection_solver =
        CCPoissonSolverManager::getManager()->allocateSolver(d_regrid_projection_solver_type,
                                                             d_object_name + "::regrid_projection_solver",
                                                             d_regrid_projection_solver_db,
                                                             "regrid_projection_",
                                                             d_regrid_projection_precond_type,
                                                             d_object_name + "::regrid_projection_precond",
                                                             d_regrid_projection_precond_db,
                                                             "regrid_projection_pc_");
    PoissonSpecifications regrid_projection_spec(d_object_name + "::regrid_projection_spec");
    regrid_projection_spec.setCZero();

    // Get the current value of density
    for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
    {
        const double apply_time = d_integrator_time;
        d_reset_rho_fcns[k](d_rho_sc_current_idx,
                            d_rho_sc_var,
                            d_hier_math_ops,
                            -1 /*cycle_num*/,
                            apply_time,
                            d_integrator_time,
                            d_integrator_time,
                            d_reset_rho_fcns_ctx[k]);
    }
    d_hier_sc_data_ops->reciprocal(d_pressure_D_idx, d_rho_sc_current_idx);
    d_hier_sc_data_ops->scale(d_pressure_D_idx, -1.0, d_pressure_D_idx);

    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent p_coef_synch_transaction =
        SynchronizationTransactionComponent(d_pressure_D_idx, "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(p_coef_synch_transaction);
    d_side_synch_op->synchronizeData(d_integrator_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, d_U_coarsen_type);
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);

    regrid_projection_spec.setDPatchDataId(d_pressure_D_idx);

    LocationIndexRobinBcCoefs<NDIM> Phi_bc_coef;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        Phi_bc_coef.setBoundarySlope(2 * d, 0.0);
        Phi_bc_coef.setBoundarySlope(2 * d + 1, 0.0);
    }
    regrid_projection_solver->setPoissonSpecifications(regrid_projection_spec);
    regrid_projection_solver->setPhysicalBcCoef(&Phi_bc_coef);
    regrid_projection_solver->setHomogeneousBc(true);
    regrid_projection_solver->setSolutionTime(d_integrator_time);
    regrid_projection_solver->setTimeInterval(d_integrator_time, d_integrator_time);
    auto p_regrid_projection_solver = dynamic_cast<LinearSolver*>(regrid_projection_solver.getPointer());
    if (p_regrid_projection_solver)
    {
        p_regrid_projection_solver->setInitialGuessNonzero(false);
        p_regrid_projection_solver->setNullspace(true);
    }

    // Setup the right-hand-side vector for the projection-Poisson solve.
    d_hier_math_ops->div(d_Div_U_idx,
                         d_Div_U_var,
                         -1.0,
                         d_U_current_idx,
                         d_U_var,
                         d_no_fill_op,
                         d_integrator_time,
                         /*synch_cf_bdry*/ false,
                         +1.0,
                         d_Q_current_idx,
                         d_Q_var);
    const double Div_U_mean = (1.0 / volume) * d_hier_cc_data_ops->integral(d_Div_U_idx, wgt_cc_idx);
    d_hier_cc_data_ops->addScalar(d_Div_U_idx, d_Div_U_idx, -Div_U_mean);

    // Solve the projection pressure-Poisson problem.
    regrid_projection_solver->solveSystem(sol_vec, rhs_vec);
    if (d_enable_logging && d_enable_logging_solver_iterations)
        plog << d_object_name
             << "::regridProjection(): regrid projection solve "
                "number of iterations = "
             << regrid_projection_solver->getNumIterations() << "\n";
    if (d_enable_logging)
        plog << d_object_name
             << "::regridProjection(): regrid projection solve "
                "residual norm        = "
             << regrid_projection_solver->getResidualNorm() << "\n";

    // Fill ghost cells for Phi, compute Grad Phi, and set U := U - 1/rho * Grad
    // Phi
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_bc_component(d_P_scratch_idx,
                                                       DATA_REFINE_TYPE,
                                                       USE_CF_INTERPOLATION,
                                                       DATA_COARSEN_TYPE,
                                                       d_bdry_extrap_type, // TODO: update variable name
                                                       CONSISTENT_TYPE_2_BDRY,
                                                       &Phi_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> Phi_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    Phi_bdry_bc_fill_op->initializeOperatorState(Phi_bc_component, d_hierarchy);
    Phi_bdry_bc_fill_op->setHomogeneousBc(true);
    Phi_bdry_bc_fill_op->fillData(d_integrator_time);
    d_hier_math_ops->grad(d_U_current_idx,
                          d_U_var,
                          /*synch_cf_bdry*/ true,
                          d_pressure_D_idx,
                          d_pressure_D_var,
                          d_P_scratch_idx,
                          d_P_var,
                          d_no_fill_op,
                          d_integrator_time,
                          +1.0,
                          d_U_current_idx,
                          d_U_var);

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }

    // Synchronize data on the patch hierarchy.
    synchronizeHierarchyData(CURRENT_DATA);
    return;
} // regridProjection

/////////////////////////////// PRIVATE //////////////////////////////////////
void
INSVCStaggeredConservativeHierarchyIntegrator::updateOperatorsAndSolvers(const double current_time,
                                                                         const double new_time,
                                                                         const int cycle_num)
{
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    const double mu = d_mu_is_const ? d_problem_coefs.getMu() : -1.0;
    const double lambda = d_problem_coefs.getLambda();
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    double K = 0.0;
    switch (d_viscous_time_stepping_type)
    {
    case BACKWARD_EULER:
        K = 1.0;
        break;
    case FORWARD_EULER:
        K = 0.0;
        break;
    case TRAPEZOIDAL_RULE:
        K = 0.5;
        break;
    default:
        TBOX_ERROR("this statement should not be reached");
    }
    PoissonSpecifications U_problem_coefs(d_object_name + "::U_problem_coefs");
    PoissonSpecifications P_problem_coefs(d_object_name + "::P_problem_coefs");

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        // Operate on a single level at a time
        d_hier_cc_data_ops->resetLevels(ln, ln);
        d_hier_sc_data_ops->resetLevels(ln, ln);
#if (NDIM == 2)
        d_hier_nc_data_ops->resetLevels(ln, ln);
#elif (NDIM == 3)
        d_hier_ec_data_ops->resetLevels(ln, ln);
#endif
        // Get the condition number scaling on this level
        const double A_scale = d_A_scale[ln];

        // C_sc = (rho / dt) + K * lambda + L(x,t^n+1)
        d_hier_sc_data_ops->scale(d_velocity_C_idx, A_scale / dt, d_rho_sc_scratch_idx, /*interior_only*/ true);
        if (!MathUtilities<double>::equalEps(lambda, 0.0))
        {
            d_hier_sc_data_ops->addScalar(d_velocity_C_idx,
                                          d_velocity_C_idx,
                                          A_scale * K * lambda,
                                          /*interior_only*/ true);
        }
        d_hier_sc_data_ops->axpy(d_velocity_C_idx, A_scale, d_velocity_L_idx, d_velocity_C_idx);

        // D_{ec,nc} = -K * mu
        if (d_mu_is_const)
        {
#if (NDIM == 2)
            d_hier_nc_data_ops->setToScalar(d_velocity_D_idx,
                                            A_scale * (-K * mu),
                                            /*interior_only*/ false);
#elif (NDIM == 3)
            d_hier_ec_data_ops->setToScalar(d_velocity_D_idx,
                                            A_scale * (-K * mu),
                                            /*interior_only*/ false);
#endif
            d_hier_cc_data_ops->setToScalar(d_velocity_D_cc_idx,
                                            A_scale * (-K * mu),
                                            /*interior_only*/ false);
        }
        else
        {
#if (NDIM == 2)
            d_hier_nc_data_ops->scale(d_velocity_D_idx, A_scale * (-K), d_mu_interp_idx, /*interior_only*/ false);
#elif (NDIM == 3)
            d_hier_ec_data_ops->scale(d_velocity_D_idx, A_scale * (-K), d_mu_interp_idx, /*interior_only*/ false);
#endif
            d_hier_cc_data_ops->scale(d_velocity_D_cc_idx, A_scale * (-K), d_mu_scratch_idx, /*interior_only*/ false);
        }
    }
    U_problem_coefs.setCPatchDataId(d_velocity_C_idx);
    U_problem_coefs.setDPatchDataId(d_velocity_D_idx);

    // Ensure that hierarchy operators operate on all levels in the future
    d_hier_cc_data_ops->resetLevels(coarsest_ln, finest_ln);
    d_hier_sc_data_ops->resetLevels(coarsest_ln, finest_ln);
#if (NDIM == 2)
    d_hier_nc_data_ops->resetLevels(coarsest_ln, finest_ln);
#elif (NDIM == 3)
    d_hier_ec_data_ops->resetLevels(coarsest_ln, finest_ln);
#endif

    // D_sc = -1/rho
    P_problem_coefs.setCZero();
    d_hier_sc_data_ops->reciprocal(d_pressure_D_idx,
                                   d_rho_sc_scratch_idx,
                                   /*interior_only*/ false);
    d_hier_sc_data_ops->scale(d_pressure_D_idx,
                              -1.0,
                              d_pressure_D_idx,
                              /*interior_only*/ false);

    // Synchronize pressure patch data coefficient
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent p_coef_synch_transaction =
        SynchronizationTransactionComponent(d_pressure_D_idx, "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(p_coef_synch_transaction);
    d_side_synch_op->synchronizeData(d_integrator_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, d_U_coarsen_type);
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);
    P_problem_coefs.setDPatchDataId(d_pressure_D_idx);

    // Ensure that solver components are appropriately reinitialized at the
    // correct intervals or when the time step size changes. Subdomain solvers
    // are only reinitialized during the first cycle
    const bool dt_change = initial_time || !MathUtilities<double>::equalEps(dt, d_dt_previous[0]);
    const bool precond_reinit = d_integrator_step % d_precond_reinit_interval == 0;
    const bool first_cycle = cycle_num == 0;
    if (first_cycle && precond_reinit)
    {
        d_velocity_solver_needs_init = true;
        d_pressure_solver_needs_init = true;
        d_stokes_solver_needs_init = true;
    }
    else if (first_cycle && dt_change)
    {
        d_velocity_solver_needs_init = true;
        d_stokes_solver_needs_init = true;
    }

    // Setup subdomain solvers.
    const bool has_velocity_nullspace = d_normalize_velocity;
    const bool has_pressure_nullspace = d_normalize_pressure;

    if (d_velocity_solver)
    {
        d_velocity_solver->setPoissonSpecifications(U_problem_coefs);
        d_velocity_solver->setPhysicalBcCoefs(d_U_star_bc_coefs);
        d_velocity_solver->setSolutionTime(new_time);
        d_velocity_solver->setTimeInterval(current_time, new_time);
        if (d_velocity_solver_needs_init)
        {
            if (d_enable_logging)
                plog << d_object_name
                     << "::updateOperatorsAndSolvers`(): initializing "
                        "velocity subdomain solver"
                     << std::endl;
            auto p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
            if (p_velocity_solver)
            {
                p_velocity_solver->setInitialGuessNonzero(false);
                if (has_velocity_nullspace) p_velocity_solver->setNullspace(false, d_U_nul_vecs);
            }
            d_velocity_solver->initializeSolverState(*d_U_scratch_vec, *d_U_rhs_vec);
            d_velocity_solver_needs_init = false;
        }
    }

    if (d_pressure_solver)
    {
        d_pressure_solver->setPoissonSpecifications(P_problem_coefs);
        d_pressure_solver->setPhysicalBcCoef(d_Phi_bc_coef);
        d_pressure_solver->setSolutionTime(half_time);
        d_pressure_solver->setTimeInterval(current_time, new_time);
        if (d_pressure_solver_needs_init)
        {
            if (d_enable_logging)
                plog << d_object_name
                     << "::updateOperatorsAndSolvers(): initializing "
                        "pressure subdomain solver"
                     << std::endl;
            auto p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
            if (p_pressure_solver)
            {
                p_pressure_solver->setInitialGuessNonzero(false);
                if (has_pressure_nullspace) p_pressure_solver->setNullspace(true);
            }
            d_pressure_solver->initializeSolverState(*d_P_scratch_vec, *d_P_rhs_vec);
            d_pressure_solver_needs_init = false;
        }
    }

    // Setup Stokes solver.
    d_stokes_solver->setVelocityPoissonSpecifications(U_problem_coefs);
    d_stokes_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
    d_stokes_solver->setPhysicalBoundaryHelper(d_bc_helper);
    d_stokes_solver->setSolutionTime(new_time);
    d_stokes_solver->setTimeInterval(current_time, new_time);
    d_stokes_solver->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);

    auto p_stokes_linear_solver = dynamic_cast<LinearSolver*>(d_stokes_solver.getPointer());
    if (!p_stokes_linear_solver)
    {
        auto p_stokes_newton_solver = dynamic_cast<NewtonKrylovSolver*>(d_stokes_solver.getPointer());
        if (p_stokes_newton_solver) p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver().getPointer();
    }
    if (p_stokes_linear_solver)
    {
        auto p_stokes_block_pc = dynamic_cast<StaggeredStokesBlockPreconditioner*>(p_stokes_linear_solver);
        auto p_stokes_fac_pc = dynamic_cast<StaggeredStokesFACPreconditioner*>(p_stokes_linear_solver);
        auto p_vc_stokes_proj_pc = dynamic_cast<VCStaggeredStokesProjectionPreconditioner*>(p_stokes_linear_solver);
        if (!(p_stokes_block_pc || p_stokes_fac_pc))
        {
            auto p_stokes_krylov_solver = dynamic_cast<KrylovLinearSolver*>(p_stokes_linear_solver);
            if (p_stokes_krylov_solver)
            {
                p_stokes_block_pc = dynamic_cast<StaggeredStokesBlockPreconditioner*>(
                    p_stokes_krylov_solver->getPreconditioner().getPointer());

                p_stokes_fac_pc = dynamic_cast<StaggeredStokesFACPreconditioner*>(
                    p_stokes_krylov_solver->getPreconditioner().getPointer());

                p_vc_stokes_proj_pc = dynamic_cast<VCStaggeredStokesProjectionPreconditioner*>(
                    p_stokes_krylov_solver->getPreconditioner().getPointer());
            }
        }
        if (p_stokes_block_pc)
        {
            p_stokes_block_pc->setPressurePoissonSpecifications(P_problem_coefs);
            p_stokes_block_pc->setPhysicalBcCoefs(d_U_star_bc_coefs, d_Phi_bc_coef);
            p_stokes_block_pc->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
        }
        else if (p_stokes_fac_pc)
        {
            p_stokes_fac_pc->setPhysicalBcCoefs(d_U_star_bc_coefs, d_Phi_bc_coef);
            p_stokes_fac_pc->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
        }
        else
        {
            TBOX_WARNING("No special BCs set for the preconditioner \n");
        }

        if (p_vc_stokes_proj_pc)
        {
            p_vc_stokes_proj_pc->setVelocityCellCenteredDCoefficient(d_velocity_D_cc_idx);
        }
    }
    if (d_stokes_solver_needs_init)
    {
        if (d_enable_logging)
            plog << d_object_name
                 << "::updateOperatorsAndSolvers(): initializing "
                    "incompressible Stokes solver"
                 << std::endl;
        if (p_stokes_linear_solver)
        {
            p_stokes_linear_solver->setInitialGuessNonzero(true);
            if (has_velocity_nullspace || has_pressure_nullspace)
                p_stokes_linear_solver->setNullspace(false, d_nul_vecs);
        }
        d_stokes_solver->initializeSolverState(*d_sol_vec, *d_rhs_vec);
        d_stokes_solver_needs_init = false;
    }

} // updateOperatorsAndSolvers

void
INSVCStaggeredConservativeHierarchyIntegrator::setupSolverVectors(
    const Pointer<SAMRAIVectorReal<NDIM, double> >& sol_vec,
    const Pointer<SAMRAIVectorReal<NDIM, double> >& rhs_vec,
    const double current_time,
    const double new_time,
    const int /*cycle_num*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;

    if (rhs_vec->getComponentDescriptorIndex(0) != d_U_rhs_vec->getComponentDescriptorIndex(0))
    {
        d_hier_sc_data_ops->copyData(rhs_vec->getComponentDescriptorIndex(0),
                                     d_U_rhs_vec->getComponentDescriptorIndex(0));
    }
    if (rhs_vec->getComponentDescriptorIndex(1) != d_P_rhs_vec->getComponentDescriptorIndex(0))
    {
        d_hier_cc_data_ops->copyData(rhs_vec->getComponentDescriptorIndex(1),
                                     d_P_rhs_vec->getComponentDescriptorIndex(0));
    }

    // Account for the convective acceleration term N_full, which will contain the
    // rho scaling factor.
    if (!d_creeping_flow)
    {
        d_hier_sc_data_ops->axpy(
            rhs_vec->getComponentDescriptorIndex(0), -1.0, d_N_full_idx, rhs_vec->getComponentDescriptorIndex(0));
    }

    // Account for body forcing terms.
    if (d_F_fcn)
    {
        d_F_fcn->setDataOnPatchHierarchy(d_F_scratch_idx, d_F_var, d_hierarchy, half_time);
        d_hier_sc_data_ops->add(
            rhs_vec->getComponentDescriptorIndex(0), rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
    }

    // Add Brinkman penalized velocity term.
    d_hier_sc_data_ops->add(
        rhs_vec->getComponentDescriptorIndex(0), rhs_vec->getComponentDescriptorIndex(0), d_velocity_L_idx);

    // Account for internal source/sink distributions.
    if (d_Q_fcn)
    {
        TBOX_ERROR("Presently not supported for variable coefficient problems");
    }

    // Set solution components to equal most recent approximations to u(n+1) and
    // p(n+1/2).
    d_hier_sc_data_ops->copyData(sol_vec->getComponentDescriptorIndex(0), d_U_new_idx);
    d_hier_cc_data_ops->copyData(sol_vec->getComponentDescriptorIndex(1), d_P_new_idx);

    // Scale pressure approximation if necessary
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hier_cc_data_ops->resetLevels(ln, ln);
        const double A_scale = d_A_scale[ln];
        if (!MathUtilities<double>::equalEps(A_scale, 1.0))
        {
            d_hier_cc_data_ops->scale(sol_vec->getComponentDescriptorIndex(1),
                                      A_scale,
                                      d_sol_vec->getComponentDescriptorIndex(1),
                                      /*interior_only*/ true);
        }
        d_hier_cc_data_ops->resetLevels(coarsest_ln, finest_ln);
    }

    // Synchronize solution and right-hand-side data before solve.
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction =
        SynchronizationTransactionComponent(sol_vec->getComponentDescriptorIndex(0), d_U_coarsen_type);
    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent rhs_synch_transaction =
        SynchronizationTransactionComponent(rhs_vec->getComponentDescriptorIndex(0), d_F_coarsen_type);
    d_side_synch_op->resetTransactionComponent(rhs_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, d_U_coarsen_type);
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);

    return;
} // setupSolverVectors

void
INSVCStaggeredConservativeHierarchyIntegrator::resetSolverVectors(
    const Pointer<SAMRAIVectorReal<NDIM, double> >& sol_vec,
    const Pointer<SAMRAIVectorReal<NDIM, double> >& rhs_vec,
    const double current_time,
    const double /*new_time*/,
    const int /*cycle_num*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Synchronize solution data after solve.
    using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction =
        SynchronizationTransactionComponent(sol_vec->getComponentDescriptorIndex(0), d_U_coarsen_type);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, d_U_coarsen_type);
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);

    // Pull out solution components.
    d_hier_sc_data_ops->copyData(d_U_new_idx, sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P_new_idx, sol_vec->getComponentDescriptorIndex(1));

    // Scale pressure solution if necessary
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hier_cc_data_ops->resetLevels(ln, ln);
        const double A_scale = d_A_scale[ln];
        if (!MathUtilities<double>::equalEps(A_scale, 1.0))
        {
            d_hier_cc_data_ops->scale(d_P_new_idx,
                                      1.0 / A_scale,
                                      d_P_new_idx,
                                      /*interior_only*/ true);
        }
        d_hier_cc_data_ops->resetLevels(coarsest_ln, finest_ln);
    }

    // Reset the right-hand side vector.
    if (!d_creeping_flow)
    {
        d_hier_sc_data_ops->axpy(
            rhs_vec->getComponentDescriptorIndex(0), +1.0, d_N_full_idx, rhs_vec->getComponentDescriptorIndex(0));
    }

    if (d_F_fcn)
    {
        d_hier_sc_data_ops->subtract(
            rhs_vec->getComponentDescriptorIndex(0), rhs_vec->getComponentDescriptorIndex(0), d_F_scratch_idx);
        d_hier_sc_data_ops->copyData(d_F_new_idx, d_F_scratch_idx);
    }
    d_hier_sc_data_ops->subtract(
        rhs_vec->getComponentDescriptorIndex(0), rhs_vec->getComponentDescriptorIndex(0), d_velocity_L_idx);

    if (d_Q_fcn)
    {
        TBOX_ERROR("Presently not supported for variable coefficient problems");
    }
    return;
} // resetSolverVectors

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
