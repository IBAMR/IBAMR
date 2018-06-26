// Filename: INSVCStaggeredConservativeHierarchyIntegrator.cpp
// Created on 15 May 2018 by Nishant Nangia and Amneet Bhalla
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
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

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <ostream>
#include <stddef.h>
#include <string>
#include <vector>

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
#include "IBAMR_config.h"
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
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/INSStaggeredConvectiveOperatorManager.h"
#include "ibamr/INSVCStaggeredConservativeConvectiveOperator.h"
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
#include "ibtk/CartSideDoubleSpecializedConstantRefine.h"
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
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghost cells used for each variable quantity.
static const int SIDEG = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredConservativeHierarchyIntegrator::INSVCStaggeredConservativeHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : INSVCStaggeredHierarchyIntegrator(object_name, input_db, register_for_restart),
      d_rho_sc_bc_coefs(NDIM, NULL),
      d_S_fcn(NULL)
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

    switch (d_convective_time_stepping_type)
    {
    case FORWARD_EULER:
    case MIDPOINT_RULE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::INSVCStaggeredConservativeHierarchyIntegrator():\n"
                                 << "  unsupported convective time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type) << " \n"
                                 << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE\n");
    }

    // Check to see whether the convective operator type has been set.
    d_convective_op_type = INSStaggeredConvectiveOperatorManager::DEFAULT;
    if (input_db->keyExists("convective_op_type"))
        d_convective_op_type = input_db->getString("convective_op_type");
    else if (input_db->keyExists("convective_operator_type"))
        d_convective_op_type = input_db->getString("convective_operator_type");
    else if (input_db->keyExists("default_convective_op_type"))
        d_convective_op_type = input_db->getString("default_convective_op_type");
    else if (input_db->keyExists("default_convective_operator_type"))
        d_convective_op_type = input_db->getString("default_convective_operator_type");

    if (d_convective_op_type != "VC_CONSERVATIVE_OP")
    {
        TBOX_ERROR(d_object_name << "::INSVCStaggeredConservativeHierarchyIntegrator():\n"
                                 << " variable coefficient conservative discretization\n"
                                 << " requires VC_CONSERVATIVE_OP convective operator\n");
    }

    // Side centered state variable for density and interpolated density variable for plotting.
    d_rho_sc_var = new SideVariable<NDIM, double>(d_object_name + "::rho_sc");
    d_rho_interp_cc_var = new CellVariable<NDIM, double>(d_object_name + "::rho_interp_cc", NDIM);

    // Side centered state variable for old velocity, required for the convective operator.
    d_U_old_var = new SideVariable<NDIM, double>(d_object_name + "::U_old");

    return;
} // INSVCStaggeredConservativeHierarchyIntegrator

INSVCStaggeredConservativeHierarchyIntegrator::~INSVCStaggeredConservativeHierarchyIntegrator()
{
    // intentionally blank
    return;
} // ~INSVCStaggeredConservativeHierarchyIntegrator

void
INSVCStaggeredConservativeHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    INSVCStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Get the density variable, which must be side-centered and maintained by the INS integrator for
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
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_rho_init_fcn);
    }
    else
    {
        TBOX_ERROR("INSVCStaggeredConservativeHierarchyIntegrator::initializeHierarchyIntegrator():\n"
                   << "  rho_is_const == false but no mass density variable has been registered.\n");
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

    // The old velocity variable, required for conservative convective operator
    registerVariable(d_U_old_current_idx,
                     d_U_old_new_idx,
                     d_U_old_scratch_idx,
                     d_U_old_var,
                     side_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_U_init);

    // Set the optional density source function.
    INSVCStaggeredConservativeConvectiveOperator* p_vc_convective_op =
        dynamic_cast<INSVCStaggeredConservativeConvectiveOperator*>(d_convective_op.getPointer());
    if (p_vc_convective_op)
    {
        p_vc_convective_op->setSideCenteredDensityBoundaryConditions(d_rho_sc_bc_coefs);
        if (d_S_fcn) p_vc_convective_op->setMassDensitySourceTerm(d_S_fcn);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchyIntegrator():\n"
                                 << " variable coefficient conservative discretization form \n"
                                 << " presently requires VC_CONSERVATIVE_OP convective operator\n"
                                 << " this statement should not have been reached\n");
    }

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

    // Keep track of the number of cycles to be used for the present integration
    // step.
    if (!d_creeping_flow && (d_current_num_cycles == 1) && (d_convective_time_stepping_type == MIDPOINT_RULE))
    {
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                 << "  time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_convective_time_stepping_type)
                                 << " requires num_cycles > 1.\n"
                                 << "  at current time step, num_cycles = " << d_current_num_cycles << "\n");
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Get the current value of density and store it in scratch.
    // Note that for conservative differencing, it is required that
    // the side-centered density is a state variable maintained
    // by this integrator.
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_sc_var);
    TBOX_ASSERT(d_rho_sc_current_idx >= 0);
    TBOX_ASSERT(d_rho_sc_scratch_idx >= 0);
#endif
    
    // Note that we always reset current context of state variables here, if necessary.
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
    d_hier_sc_data_ops->copyData(d_rho_sc_scratch_idx, d_rho_sc_current_idx, /*interior_only*/ true);

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
        // Note that we always reset current context of state variables here, if necessary.
        for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
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
        
        d_hier_cc_data_ops->copyData(d_mu_scratch_idx, mu_current_idx, /*interior_only*/ true);
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
        d_hier_cc_data_ops->copyData(d_mu_linear_op_idx, d_mu_scratch_idx, /*interior_only*/ false);
#if (NDIM == 2)
        d_hier_nc_data_ops->copyData(d_mu_interp_linear_op_idx, d_mu_interp_idx, /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->copyData(d_mu_interp_linear_op_idx, d_mu_interp_idx, /*interior_only*/ false);
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
        d_hier_nc_data_ops->setToScalar(d_velocity_rhs_D_idx, +K_rhs * mu, /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->setToScalar(d_velocity_rhs_D_idx, +K_rhs * mu, /*interior_only*/ false);
#endif
    }
    else
    {
#if (NDIM == 2)
        d_hier_nc_data_ops->scale(d_velocity_rhs_D_idx, +K_rhs, d_mu_interp_idx, /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->scale(d_velocity_rhs_D_idx, +K_rhs, d_mu_interp_idx, /*interior_only*/ false);
#endif
    }
    U_rhs_problem_coefs.setDPatchDataId(d_velocity_rhs_D_idx);

    const int U_rhs_idx = d_U_rhs_vec->getComponentDescriptorIndex(0);
    const Pointer<SideVariable<NDIM, double> > U_rhs_var = d_U_rhs_vec->getComponentVariable(0);
    d_hier_sc_data_ops->copyData(d_U_scratch_idx, d_U_current_idx);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_U_bc_coefs,
                                                              /*P_bc_coef*/ NULL,
                                                              d_U_scratch_idx,
                                                              /*P_data_idx*/ -1,
                                                              /*homogeneous_bc*/ false);
    d_U_bdry_bc_fill_op->fillData(current_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs,
                                                              /*P_bc_coef*/ NULL);
    d_bc_helper->enforceDivergenceFreeConditionAtBoundary(d_U_scratch_idx);
    // RHS^n = (C_rhs*I + L(D_rhs))*U^n
    d_hier_math_ops->vc_laplace(U_rhs_idx,
                                U_rhs_var,
                                1.0,
                                0.0,
                                U_rhs_problem_coefs.getDPatchDataId(),
#if (NDIM == 2)
                                Pointer<NodeVariable<NDIM, double> >(NULL),
#elif (NDIM == 3)
                                Pointer<EdgeVariable<NDIM, double> >(NULL),
#endif
                                d_U_scratch_idx,
                                d_U_var,
                                d_no_fill_op,
                                current_time,
                                d_mu_vc_interp_type);

    // Add the momentum portion of the RHS in the case of conservative discretization form
    // RHS^n = RHS^n + 1/dt*(rho*U)^n
    d_hier_sc_data_ops->multiply(d_temp_sc_idx, d_rho_sc_scratch_idx, d_U_scratch_idx, /*interior_only*/ true);
    d_hier_sc_data_ops->scale(d_temp_sc_idx, 1.0 / dt, d_temp_sc_idx, /*interior_only*/ true);
    d_hier_sc_data_ops->add(U_rhs_idx, d_temp_sc_idx, U_rhs_idx, /*interior_only*/ true);
    d_hier_sc_data_ops->copyData(d_U_src_idx, d_U_scratch_idx, /*interior_only*/ false);

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
                                     << "  number of cycles required by coupled advection-diffusion solver = "
                                     << adv_diff_num_cycles << ".\n"
                                     << "  current implementation requires either that both solvers use the same "
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
        const int U_adv_idx = d_U_adv_vec->getComponentDescriptorIndex(0);
        d_hier_sc_data_ops->copyData(U_adv_idx, d_U_current_idx);
        for (int ln = finest_ln; ln > coarsest_ln; --ln)
        {
            Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            Pointer<CoarsenOperator<NDIM> > coarsen_op =
                grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
            coarsen_alg->registerCoarsen(U_adv_idx, U_adv_idx, coarsen_op);
            coarsen_alg->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
            getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]->coarsenData();
            getCoarsenAlgorithm(d_object_name + "::CONVECTIVE_OP")
                ->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
        }
        d_convective_op->setAdvectionVelocity(d_U_adv_vec->getComponentDescriptorIndex(0));
        d_convective_op->setSolutionTime(current_time);
        d_convective_op->setTimeInterval(current_time, new_time);

        // For conservative momentum discretization, an approximation to rho^{n+1}
        // will be computed from rho^{n}, which requires additional options to be set.
        INSVCStaggeredConservativeConvectiveOperator* p_vc_convective_op =
            dynamic_cast<INSVCStaggeredConservativeConvectiveOperator*>(d_convective_op.getPointer());
        if (p_vc_convective_op)
        {
            // Set the rho^{n} density
            p_vc_convective_op->setSideCenteredDensityPatchDataIndex(d_rho_sc_current_idx);

            // Data for the convective operator is for cycle 0
            const int ins_cycle_num = 0;
            p_vc_convective_op->setCycleNumber(ins_cycle_num);

            // Set the velocities used to update the density
            if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
            {
                p_vc_convective_op->setFluidVelocityPatchDataIndices(
                    /*old*/ -1, /*current*/ d_U_current_idx, /*new*/ -1);
            }
            else
            {
                p_vc_convective_op->setFluidVelocityPatchDataIndices(
                    /*old*/ d_U_old_current_idx, /*current*/ d_U_current_idx, /*new*/ -1);
            }
        }
        else
        {
            TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                     << " variable coefficient conservative discretization form \n"
                                     << " presently requires VC_CONSERVATIVE_OP convective operator\n"
                                     << " this statement should not have been reached\n");
        }

        // Compute the convective operator
        d_convective_op->apply(*d_U_adv_vec, *d_N_vec);

        // Keep track of the time-lagged velocity
        d_hier_sc_data_ops->copyData(d_U_old_new_idx, d_U_current_idx);
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

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE({
            pout << "INSVCStaggeredConservativeHierarchyIntegrator::integrateHierarchy():\n"
                 << "  WARNING: num_cycles = " << d_current_num_cycles
                 << " but expected num_cycles = " << expected_num_cycles << ".\n";
        });
    }

    // Update the state variables of any linked advection-diffusion solver.
    // NOTE: This also updates rho and mu if they are maintained by adv-diff integrator.
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

    // In the special case of a conservative discretization form, the updated density is previously calculated by
    // application of the convective operator
    INSVCStaggeredConservativeConvectiveOperator* p_vc_convective_op =
        dynamic_cast<INSVCStaggeredConservativeConvectiveOperator*>(d_convective_op.getPointer());
    if (p_vc_convective_op)
    {
        const int rho_sc_new_idx = p_vc_convective_op->getUpdatedSideCenteredDensityPatchDataIndex();
        d_hier_sc_data_ops->copyData(d_rho_sc_new_idx, rho_sc_new_idx, /*interior_only*/ true);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                 << " variable coefficient conservative discretization form \n"
                                 << " presently requires VC_CONSERVATIVE_OP convective operator\n"
                                 << " this statement should not have been reached\n");
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
        d_hier_cc_data_ops->copyData(d_mu_scratch_idx, mu_new_idx, /*interior_only*/ true);
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
        d_hier_cc_data_ops->copyData(d_mu_linear_op_idx, d_mu_scratch_idx, /*interior_only*/ false);
#if (NDIM == 2)
        d_hier_nc_data_ops->copyData(d_mu_interp_linear_op_idx, d_mu_interp_idx, /*interior_only*/ false);
#elif (NDIM == 3)
        d_hier_ec_data_ops->copyData(d_mu_interp_linear_op_idx, d_mu_interp_idx, /*interior_only*/ false);
#endif
    }

    // Copy new into scratch
    d_hier_sc_data_ops->copyData(d_rho_sc_scratch_idx, d_rho_sc_new_idx);

    // Store the density for later use
    d_hier_sc_data_ops->copyData(d_rho_linear_op_idx, d_rho_sc_scratch_idx, /*interior_only*/ true);

    // Update the solvers and operators to take into account new state variables
    updateOperatorsAndSolvers(current_time, new_time);

    // Setup the solution and right-hand-side vectors.
    setupSolverVectors(d_sol_vec, d_rhs_vec, current_time, new_time, cycle_num);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hier_sc_data_ops->resetLevels(ln, ln);
        const double A_scale = d_A_scale[ln];
        // Scale rhs if necessary
        if (!MathUtilities<double>::equalEps(A_scale, 1.0))
        {
            d_hier_sc_data_ops->scale(d_rhs_vec->getComponentDescriptorIndex(0),
                                      A_scale,
                                      d_rhs_vec->getComponentDescriptorIndex(0),
                                      /*interior_only*/ true);
        }
    }
        // Solve for u(n+1), p(n+1/2).
        d_stokes_solver->solveSystem(*d_sol_vec, *d_rhs_vec);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hier_sc_data_ops->resetLevels(ln, ln);
        const double A_scale = d_A_scale[ln];
        // Unscale rhs if necessary
        if (!MathUtilities<double>::equalEps(A_scale, 1.0))
        {
            d_hier_sc_data_ops->scale(d_rhs_vec->getComponentDescriptorIndex(0),
                                      1.0 / A_scale,
                                      d_rhs_vec->getComponentDescriptorIndex(0),
                                      /*interior_only*/ true);
        }
        d_hier_sc_data_ops->resetLevels(coarsest_ln, finest_ln);
    }

    if (d_enable_logging)
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

    // Re-update density and viscosity is they are maintained by the integrator
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
INSVCStaggeredConservativeHierarchyIntegrator::regridHierarchy()
{
    INSVCStaggeredHierarchyIntegrator::regridHierarchy();
    return;
} // regridHierarchy

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
                                 << " term to be set. Consider combining source terms into single CartGridFunction.\n");
    }
    return;
} // registerMassDensitySourceTerm

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
    INSVCStaggeredHierarchyIntegrator::regridProjection();
    return;
} // regridProjection

/////////////////////////////// PRIVATE //////////////////////////////////////
TimeSteppingType
INSVCStaggeredConservativeHierarchyIntegrator::getConvectiveTimeSteppingType(const int /*cycle_num*/)
{
    TimeSteppingType convective_time_stepping_type = d_convective_time_stepping_type;
    return convective_time_stepping_type;
} // getConvectiveTimeSteppingType

void
INSVCStaggeredConservativeHierarchyIntegrator::updateOperatorsAndSolvers(const double current_time,
                                                                         const double new_time)
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

        // C_sc = (rho / dt) + K * lambda
        d_hier_sc_data_ops->scale(d_velocity_C_idx, A_scale / dt, d_rho_sc_scratch_idx, /*interior_only*/ true);

        if (!MathUtilities<double>::equalEps(lambda, 0.0))
        {
            d_hier_sc_data_ops->addScalar(
                d_velocity_C_idx, d_velocity_C_idx, A_scale * K * lambda, /*interior_only*/ true);
        }
        U_problem_coefs.setCPatchDataId(d_velocity_C_idx);

        // D_{ec,nc} = -K * mu
        if (d_mu_is_const)
        {
#if (NDIM == 2)
            d_hier_nc_data_ops->setToScalar(d_velocity_D_idx, A_scale * (-K * mu), /*interior_only*/ false);
#elif (NDIM == 3)
            d_hier_ec_data_ops->setToScalar(d_velocity_D_idx, A_scale * (-K * mu), /*interior_only*/ false);
#endif
            d_hier_cc_data_ops->setToScalar(d_velocity_D_cc_idx, A_scale * (-K * mu), /*interior_only*/ false);
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
        U_problem_coefs.setDPatchDataId(d_velocity_D_idx);

        // D_sc = -1/rho
        P_problem_coefs.setCZero();
        d_hier_sc_data_ops->reciprocal(d_pressure_D_idx, d_rho_sc_scratch_idx, /*interior_only*/ false);
        d_hier_sc_data_ops->scale(d_pressure_D_idx, -1.0 / A_scale, d_pressure_D_idx, /*interior_only*/ false);
        P_problem_coefs.setDPatchDataId(d_pressure_D_idx);

        // Ensure that these objects will operate on all levels in the future
        d_hier_cc_data_ops->resetLevels(coarsest_ln, finest_ln);
        d_hier_sc_data_ops->resetLevels(coarsest_ln, finest_ln);
#if (NDIM == 2)
        d_hier_nc_data_ops->resetLevels(coarsest_ln, finest_ln);
#elif (NDIM == 3)
        d_hier_ec_data_ops->resetLevels(coarsest_ln, finest_ln);
#endif
    }

    // Ensure that solver components are appropriately reinitialized at the correct intervals or
    // when the time step size changes.
    const bool dt_change = initial_time || !MathUtilities<double>::equalEps(dt, d_dt_previous[0]);
    const bool precond_reinit = d_integrator_step % d_precond_reinit_interval == 0;
    if (precond_reinit)
    {
        d_velocity_solver_needs_init = true;
        d_pressure_solver_needs_init = true;
        d_stokes_solver_needs_init = true;
    }
    else if (dt_change)
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
            LinearSolver* p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
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
            LinearSolver* p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
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

    LinearSolver* p_stokes_linear_solver = dynamic_cast<LinearSolver*>(d_stokes_solver.getPointer());
    if (!p_stokes_linear_solver)
    {
        NewtonKrylovSolver* p_stokes_newton_solver = dynamic_cast<NewtonKrylovSolver*>(d_stokes_solver.getPointer());
        if (p_stokes_newton_solver) p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver().getPointer();
    }
    if (p_stokes_linear_solver)
    {
        StaggeredStokesBlockPreconditioner* p_stokes_block_pc =
            dynamic_cast<StaggeredStokesBlockPreconditioner*>(p_stokes_linear_solver);
        StaggeredStokesFACPreconditioner* p_stokes_fac_pc =
            dynamic_cast<StaggeredStokesFACPreconditioner*>(p_stokes_linear_solver);
        VCStaggeredStokesProjectionPreconditioner* p_vc_stokes_proj_pc =
            dynamic_cast<VCStaggeredStokesProjectionPreconditioner*>(p_stokes_linear_solver);
        if (!(p_stokes_block_pc || p_stokes_fac_pc))
        {
            KrylovLinearSolver* p_stokes_krylov_solver = dynamic_cast<KrylovLinearSolver*>(p_stokes_linear_solver);
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
    const int cycle_num)
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

    // Account for the convective acceleration term N_full, which will contain the rho scaling factor.
    if (!d_creeping_flow)
    {
        const int N_idx = d_N_vec->getComponentDescriptorIndex(0);
        TimeSteppingType convective_time_stepping_type = getConvectiveTimeSteppingType(cycle_num);

        // Update N_idx if necessary
        if (cycle_num > 0)
        {
            const int U_adv_idx = d_U_adv_vec->getComponentDescriptorIndex(0);
            double apply_time = std::numeric_limits<double>::quiet_NaN();
            if (convective_time_stepping_type == FORWARD_EULER)
            {
                d_hier_sc_data_ops->copyData(U_adv_idx, d_U_new_idx);
                apply_time = current_time;
            }
            else if (convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_sc_data_ops->linearSum(U_adv_idx, 0.5, d_U_current_idx, 0.5, d_U_new_idx);
                apply_time = half_time;
            }
            for (int ln = finest_ln; ln > coarsest_ln; --ln)
            {
                Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
                Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
                Pointer<CoarsenOperator<NDIM> > coarsen_op =
                    grid_geom->lookupCoarsenOperator(d_U_var, "CONSERVATIVE_COARSEN");
                coarsen_alg->registerCoarsen(U_adv_idx, U_adv_idx, coarsen_op);
                coarsen_alg->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
                getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]->coarsenData();
                getCoarsenAlgorithm(d_object_name + "::CONVECTIVE_OP")
                    ->resetSchedule(getCoarsenSchedules(d_object_name + "::CONVECTIVE_OP")[ln]);
            }
            d_convective_op->setAdvectionVelocity(d_U_adv_vec->getComponentDescriptorIndex(0));
            d_convective_op->setSolutionTime(apply_time);

            INSVCStaggeredConservativeConvectiveOperator* p_vc_convective_op =
                dynamic_cast<INSVCStaggeredConservativeConvectiveOperator*>(d_convective_op.getPointer());
            if (p_vc_convective_op)
            {
                // Always set to current because we want to update rho^{n} to rho^{n+1}
                p_vc_convective_op->setSideCenteredDensityPatchDataIndex(d_rho_sc_current_idx);

                // Set the cycle number
                p_vc_convective_op->setCycleNumber(cycle_num);

                // Set the velocities used to update the density
                if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
                {
                    p_vc_convective_op->setFluidVelocityPatchDataIndices(
                        /*old*/ -1, /*current*/ d_U_current_idx, /*new*/ d_U_new_idx);
                }
                else
                {
                    p_vc_convective_op->setFluidVelocityPatchDataIndices(
                        /*old*/ d_U_old_current_idx, /*current*/ d_U_current_idx, /*new*/ d_U_new_idx);
                }
            }
            else
            {
                TBOX_ERROR(d_object_name << "::setupSolverVectors():\n"
                                         << " variable coefficient conservative discretization form \n"
                                         << " presently requires VC_CONSERVATIVE_OP convective operator\n"
                                         << " this statement should not have been reached\n");
            }
            d_convective_op->apply(*d_U_adv_vec, *d_N_vec);
        }

        // Set the convective term depending on the time stepping type
        if (convective_time_stepping_type == FORWARD_EULER || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_sc_data_ops->copyData(d_N_full_idx, N_idx, /*interior_only*/ true);
        }
        else
        {
            TBOX_ERROR("this statement should not be reached");
        }

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
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction =
        SynchronizationTransactionComponent(sol_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(sol_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent rhs_synch_transaction =
        SynchronizationTransactionComponent(rhs_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(rhs_synch_transaction);
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, "CONSERVATIVE_COARSEN");
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
    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent;
    SynchronizationTransactionComponent sol_synch_transaction =
        SynchronizationTransactionComponent(sol_vec->getComponentDescriptorIndex(0), "CONSERVATIVE_COARSEN");
    d_side_synch_op->synchronizeData(current_time);
    SynchronizationTransactionComponent default_synch_transaction =
        SynchronizationTransactionComponent(d_U_scratch_idx, "CONSERVATIVE_COARSEN");
    d_side_synch_op->resetTransactionComponent(default_synch_transaction);

    // Pull out solution components.
    d_hier_sc_data_ops->copyData(d_U_new_idx, sol_vec->getComponentDescriptorIndex(0));
    d_hier_cc_data_ops->copyData(d_P_new_idx, sol_vec->getComponentDescriptorIndex(1));

    // Scale pressure solution if necessary
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hier_sc_data_ops->resetLevels(ln, ln);
        const double A_scale = d_A_scale[ln];
        if (!MathUtilities<double>::equalEps(A_scale, 1.0))
        {
            d_hier_cc_data_ops->scale(d_P_new_idx,
                                      1.0 / A_scale,
                                      d_P_new_idx,
                                      /*interior_only*/ true);
        }
        d_hier_sc_data_ops->resetLevels(coarsest_ln, finest_ln);
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
    if (d_Q_fcn)
    {
        TBOX_ERROR("Presently not supported for variable coefficient problems");
    }
    return;
} // resetSolverVectors

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
