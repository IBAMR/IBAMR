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

/////////////////////////////// INCLUDES /////////////////////////////////////
// #include "ibamr/AdvDiffCUIConservativeConvectiveOperator.h"
#include "ibamr/AdvDiffCUIConvectiveOperator.h"
#include "ibamr/AdvDiffConservativeMassScalarTransportRKIntegrator.h"
#include "ibamr/AllenCahnHierarchyIntegrator.h"
#include "ibamr/CellConvectiveOperator.h"
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
#include "MultiblockDataTranslator.h"
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

#include <algorithm>
#include <deque>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
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
static const int IEP_HIERARCHY_INTEGRATOR_VERSION = 4;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
static const int FACEG = 1;
static const int NOGHOSTS = 0;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

AllenCahnHierarchyIntegrator::AllenCahnHierarchyIntegrator(const std::string& object_name,
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
    if (from_restart) getFromRestart();
    getFromInput(input_db, from_restart);

    if (!(d_lf_convective_difference_form == CONSERVATIVE))
    {
        TBOX_ERROR(d_object_name << "::AllenCahnHierarchyIntegrator():\n"
                                 << " variable coefficient discretization\n"
                                 << " requires CONSERVATIVE convective difference form\n");
    }

    if (!(d_lf_convective_op_type == "CUI_CONSERVATIVE"))
    {
        TBOX_ERROR(d_object_name << "::AllenCahnHierarchyIntegrator():\n"
                                 << " current implementation supports only\n"
                                 << " CUI CONSERVATIVE limiter for Allen-Cahn equation\n");
    }

    if (input_db->keyExists("lf_solver_type"))
    {
        d_lf_solver_type = input_db->getString("lf_solver_type");
        if (input_db->keyExists("lf_solver_db")) d_lf_solver_db = input_db->getDatabase("lf_solver_db");
    }
    if (!d_lf_solver_db) d_lf_solver_db = new MemoryDatabase("lf_solver_db");

    if (input_db->keyExists("lf_precond_type"))
    {
        d_lf_precond_type = input_db->getString("lf_precond_type");
        if (input_db->keyExists("lf_precond_db")) d_lf_precond_db = input_db->getDatabase("lf_precond_db");
    }
    if (!d_lf_precond_db) d_lf_precond_db = new MemoryDatabase("lf_precond_db");

    return;
} // AllenCahnHierarchyIntegrator

void
AllenCahnHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                            Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Perform hierarchy initialization operations common to all implementations
    // of PhaseChangeHierarchyIntegrator.
    PhaseChangeHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Operators and solvers are maintained for each variable registered with the
    // integrator.
    if (d_lf_solver_type == CCPoissonSolverManager::UNDEFINED)
    {
        d_lf_solver_type = CCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER;
    }
    if (d_lf_precond_type == CCPoissonSolverManager::UNDEFINED)
    {
        const int max_levels = gridding_alg->getMaxLevels();
        if (max_levels == 1)
        {
            d_lf_precond_type = CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
        }
        else
        {
            d_lf_precond_type = CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
        }
        d_lf_precond_db->putInteger("max_iterations", 1);
    }

    d_lf_solver = getAllenCahnEquationHelmholtzSolver(d_lf_var);
    d_lf_rhs_op = getAllenCahnEquationHelmholtzRHSOperator(d_lf_var);

    // Register additional variables required for present time stepping algorithm.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> face_ghosts = FACEG;
    const IntVector<NDIM> no_ghosts = NOGHOSTS;

    int lf_F_current_idx, lf_F_scratch_idx, lf_F_new_idx;
    registerVariable(lf_F_current_idx,
                     lf_F_new_idx,
                     lf_F_scratch_idx,
                     d_lf_F_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    int lf_diff_coef_current_idx, lf_diff_coef_scratch_idx, lf_diff_coef_new_idx;
    if (d_lf_diffusion_coef_var)
        registerVariable(lf_diff_coef_current_idx,
                         lf_diff_coef_new_idx,
                         lf_diff_coef_scratch_idx,
                         d_lf_diffusion_coef_var,
                         side_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int lf_diffusion_coef_rhs_scratch_idx;
    d_lf_diffusion_coef_rhs_var = new SideVariable<NDIM, double>(d_lf_var->getName() + "::Diff");
    registerVariable(lf_diffusion_coef_rhs_scratch_idx, d_lf_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    int lf_rhs_scratch_idx;
    registerVariable(lf_rhs_scratch_idx, d_lf_rhs_var, cell_ghosts, getScratchContext());

    int lf_H_scratch_idx;
    d_lf_H_var = new CellVariable<NDIM, double>(d_object_name + "::lf_H");
    registerVariable(lf_H_scratch_idx, d_lf_H_var, cell_ghosts, getScratchContext());

    int Div_u_scratch_idx;
    d_Div_u_var = new CellVariable<NDIM, double>(d_object_name + "::Div_u");
    registerVariable(Div_u_scratch_idx, d_Div_u_var, no_ghosts, getScratchContext());

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_lf_C_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::C");
    d_lf_C_idx = var_db->registerVariableAndContext(d_lf_C_var, getCurrentContext(), no_ghosts);

    d_lf_temp_rhs_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::temp_rhs");
    d_lf_temp_rhs_idx = var_db->registerVariableAndContext(d_lf_temp_rhs_var, getCurrentContext(), no_ghosts);

    d_lf_N_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::N");
    int lf_N_scratch_idx;
    registerVariable(lf_N_scratch_idx, d_lf_N_var, cell_ghosts, getScratchContext());

    d_T_lf_N_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::lf_N");
    registerVariable(d_T_lf_N_scratch_idx, d_T_lf_N_var, cell_ghosts, getScratchContext());

    d_lf_N_old_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::lf_N_old");
    int lf_N_old_current_idx, lf_N_old_new_idx, lf_N_old_scratch_idx;
    registerVariable(lf_N_old_current_idx,
                     lf_N_old_new_idx,
                     lf_N_old_scratch_idx,
                     d_lf_N_old_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // Registering a temporary cell-centered vector variable to be used in the interpolation
    // function.
    d_g_firstder_var = new CellVariable<NDIM, double>("g_firstder_var");
    d_g_firstder_idx = var_db->registerVariableAndContext(d_g_firstder_var, getCurrentContext());

    d_g_secondder_var = new CellVariable<NDIM, double>("g_secondder_var");
    d_g_secondder_idx = var_db->registerVariableAndContext(d_g_secondder_var, getCurrentContext());

    d_q_firstder_var = new CellVariable<NDIM, double>("q_firstder_var");
    d_q_firstder_idx = var_db->registerVariableAndContext(d_q_firstder_var, getCurrentContext());

    d_chemical_potential_var = new CellVariable<NDIM, double>("chemical_potential_var");
    registerVariable(d_chemical_potential_idx, d_chemical_potential_var, no_ghosts, getCurrentContext());

    d_lf_interp_var = new FaceVariable<NDIM, double>(d_lf_var->getName() + "::interp");
    d_lf_interp_idx = var_db->registerVariableAndContext(d_lf_interp_var, getCurrentContext(), no_ghosts);

    d_H_interp_var = new FaceVariable<NDIM, double>(d_H_var->getName() + "::interp");
    d_H_interp_idx = var_db->registerVariableAndContext(d_H_interp_var, getCurrentContext(), no_ghosts);

    d_lf_flux_var = new FaceVariable<NDIM, double>(d_H_var->getName() + "::flux");
    d_lf_flux_idx = var_db->registerVariableAndContext(d_lf_flux_var, getCurrentContext(), no_ghosts);

    if (d_visit_writer)
    {
        d_visit_writer->registerPlotQuantity("zeta", "SCALAR", d_chemical_potential_idx);
    }

    d_grad_lf_var = new SideVariable<NDIM, double>(d_object_name + "::grad_lf");
    d_grad_lf_idx =
        var_db->registerVariableAndContext(d_grad_lf_var, var_db->getContext(d_object_name + "grad_lf::SCRATCH"));
    d_H_sc_idx = var_db->registerVariableAndContext(d_grad_lf_var, var_db->getContext(d_object_name + "H_sc::SCRATCH"));

    // Setup the convective operator. We use AdvDiffConservativeCUIOperator class and this class is not registered
    // with ConvectiveOperatorManager
    d_lf_convective_op = getAllenCahnEquationConvectiveOperator(d_lf_var, d_H_var);

    if (d_solve_mass_conservation)
    {
        // Set various objects with conservative time integrator.
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
        rho_p_cc_integrator->setMaterialPropertyVariable(d_Cp_var);
        rho_p_cc_integrator->setCellCenteredMaterialPropertyBoundaryConditions(d_Cp_bc_coef);
        rho_p_cc_integrator->setCellCenteredTransportQuantityBoundaryConditions(d_T_bc_coef);
    }

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
AllenCahnHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                           const double new_time,
                                                           const int num_cycles)
{
    PhaseChangeHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Indicate that all solvers need to be reinitialized if the current
    // timestep size is different from the previous one.
    const bool dt_change = initial_time || !MathUtilities<double>::equalEps(dt, d_dt_previous[0]);
    if (dt_change)
    {
        d_lf_solver_needs_init = true;
        d_lf_rhs_op_needs_init = true;
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_lf_C_idx)) level->allocatePatchData(d_lf_C_idx, current_time);
        if (!level->checkAllocated(d_lf_temp_rhs_idx)) level->allocatePatchData(d_lf_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_g_firstder_idx)) level->allocatePatchData(d_g_firstder_idx, current_time);
        if (!level->checkAllocated(d_g_secondder_idx)) level->allocatePatchData(d_g_secondder_idx, current_time);
        if (!level->checkAllocated(d_q_firstder_idx)) level->allocatePatchData(d_q_firstder_idx, current_time);
        if (!level->checkAllocated(d_grad_lf_idx)) level->allocatePatchData(d_grad_lf_idx, current_time);
        if (!level->checkAllocated(d_H_sc_idx)) level->allocatePatchData(d_H_sc_idx, current_time);
        if (!level->checkAllocated(d_lf_interp_idx)) level->allocatePatchData(d_lf_interp_idx, current_time);
        if (!level->checkAllocated(d_H_interp_idx)) level->allocatePatchData(d_H_interp_idx, current_time);
        if (!level->checkAllocated(d_lf_flux_idx)) level->allocatePatchData(d_lf_flux_idx, current_time);
    }

    if (d_u_adv_var)
    {
        // Update the advection velocity.
        const int u_current_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getCurrentContext());
        const int u_scratch_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getScratchContext());
        const int u_new_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getNewContext());

        if (d_u_fcn[d_u_adv_var])
        {
            d_u_fcn[d_u_adv_var]->setDataOnPatchHierarchy(u_current_idx, d_u_adv_var, d_hierarchy, current_time);
            d_u_fcn[d_u_adv_var]->setDataOnPatchHierarchy(u_new_idx, d_u_adv_var, d_hierarchy, new_time);
        }
        else
        {
            d_hier_fc_data_ops->copyData(u_new_idx, u_current_idx);
        }
        d_hier_fc_data_ops->linearSum(u_scratch_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    }
    // Setup the operators and solvers and compute the right-hand-side terms.
    const int lf_current_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getCurrentContext());
    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
    const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
    const int lf_diff_coef_current_idx =
        var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_var, getCurrentContext());
    const int lf_diff_coef_rhs_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_rhs_var, getScratchContext());
    const int lf_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_N_var, getScratchContext());
    const int H_current_idx = var_db->mapVariableAndContextToIndex(d_H_var, getCurrentContext());
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    const int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());

    // setting lf equation diffusion timestepping type
    double alpha = 0.0;
    switch (d_lf_diffusion_time_stepping_type)
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
        TBOX_ERROR(d_object_name << "::preprocessintegrateHierarchy():\n"
                                 << "  unsupported diffusion time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_lf_diffusion_time_stepping_type) << " \n"
                                 << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }
    PoissonSpecifications lf_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_lf_var->getName());

    // set C coefficients.
    computeDoubleWellPotential(d_g_firstder_idx, d_g_secondder_idx, lf_current_idx);

    d_hier_cc_data_ops->scale(d_lf_temp_rhs_idx, 1.0 / dt, H_current_idx);
    lf_rhs_op_spec.setCPatchDataId(d_lf_temp_rhs_idx);

    // Interpolate the cell-centered Heaviside to side-centered.
    d_hier_cc_data_ops->copyData(H_scratch_idx, H_current_idx);
    d_H_bdry_bc_fill_op->fillData(current_time);

    interpolateCCToSCSimpleAveraging(lf_diff_coef_current_idx, H_scratch_idx);

    // Add numerical diffusion to smoothly extend the liquid fraction.
    d_hier_sc_data_ops->addScalar(lf_diff_coef_current_idx, lf_diff_coef_current_idx, d_num_diffusion);
    d_hier_sc_data_ops->scale(lf_diff_coef_current_idx, d_M_lf * d_lambda_lf, lf_diff_coef_current_idx);

    d_hier_sc_data_ops->scale(lf_diff_coef_rhs_scratch_idx, (1.0 - alpha), lf_diff_coef_current_idx);
    lf_rhs_op_spec.setDPatchDataId(lf_diff_coef_rhs_scratch_idx);

    // Initialize the RHS operator and compute the RHS vector for lf equation.
    Pointer<LaplaceOperator> lf_rhs_op = d_lf_rhs_op;
    lf_rhs_op->setPoissonSpecifications(lf_rhs_op_spec);
    lf_rhs_op->setPhysicalBcCoef(d_lf_bc_coef);
    lf_rhs_op->setHomogeneousBc(false);
    lf_rhs_op->setSolutionTime(current_time);
    lf_rhs_op->setTimeInterval(current_time, new_time);
    if (d_lf_rhs_op_needs_init)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the RHS operator for" << d_lf_var->getName() << "\n";
        }
        lf_rhs_op->initializeOperatorState(*d_lf_sol, *d_lf_rhs);
        d_lf_rhs_op_needs_init = false;
    }
    d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_current_idx, false);
    lf_rhs_op->apply(*d_lf_sol, *d_lf_rhs);

    if (d_u_adv_var)
    {
        if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(d_lf_convective_time_stepping_type))
        {
            d_lf_convective_time_stepping_type = d_lf_init_convective_time_stepping_type;
        }
        if ((num_cycles == 1) && (d_lf_convective_time_stepping_type == MIDPOINT_RULE ||
                                  d_lf_convective_time_stepping_type == TRAPEZOIDAL_RULE))
        {
            TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                     << "  time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_lf_convective_time_stepping_type)
                                     << " requires num_cycles > 1.\n"
                                     << "  at current time step, num_cycles = " << num_cycles << "\n");
        }
        if (d_lf_convective_op_needs_init)
        {
            d_lf_convective_op->initializeOperatorState(*d_lf_sol, *d_lf_rhs);
            d_lf_convective_op_needs_init = false;
        }
        const int u_current_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getCurrentContext());
        d_lf_convective_op->setAdvectionVelocity(u_current_idx);
        d_lf_convective_op->setSolutionTime(current_time);
        d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_current_idx);
        d_hier_cc_data_ops->copyData(H_scratch_idx, H_current_idx);

        // Checking whether the CellConvectiveOperator can be used for div (H lf u)
        d_lf_convective_op->interpolateToFaceOnHierarchy(
            d_lf_interp_idx, lf_scratch_idx, u_current_idx, /*synch_cf_bdry*/ false);
        d_lf_convective_op->interpolateToFaceOnHierarchy(
            d_H_interp_idx, H_scratch_idx, u_current_idx, /*synch_cf_bdry*/ false);
        d_hier_fc_data_ops->multiply(d_lf_flux_idx, d_lf_interp_idx, d_H_interp_idx);
        d_hier_fc_data_ops->multiply(d_lf_flux_idx, d_lf_flux_idx, u_current_idx);
        d_lf_convective_op->computeConservativeDerivativeOnHierarchy(
            lf_N_scratch_idx, d_lf_flux_idx, /*synch_cf_bdry*/ true);

        // Pointer<AdvDiffCUIConservativeConvectiveOperator> lf_cui_conservative_convective_op = d_lf_convective_op;
        // lf_cui_conservative_convective_op->applyConvectiveOperator(lf_scratch_idx, H_scratch_idx, lf_N_scratch_idx);

        const int lf_N_old_new_idx = var_db->mapVariableAndContextToIndex(d_lf_N_old_var, getNewContext());
        d_hier_cc_data_ops->copyData(lf_N_old_new_idx, lf_N_scratch_idx);

        const int lf_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_rhs_var, getScratchContext());

        if (d_lf_convective_time_stepping_type == FORWARD_EULER)
        {
            d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, -1.0, lf_N_scratch_idx, lf_rhs_scratch_idx);
        }
        else if (d_lf_convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, -0.5, lf_N_scratch_idx, lf_rhs_scratch_idx);
        }
    }

    // Set the initial guess for liquid fraction and heaviside.
    d_hier_cc_data_ops->copyData(lf_new_idx, lf_current_idx);
    d_hier_cc_data_ops->copyData(H_new_idx, H_current_idx);
    d_hier_cc_data_ops->copyData(d_H_pre_idx, H_current_idx);

    if (d_solve_energy)
    {
        const int rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());
        const int Cp_current_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getCurrentContext());
        const int T_current_idx = var_db->mapVariableAndContextToIndex(d_T_var, getCurrentContext());
        const int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, getScratchContext());
        const int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());
        const int T_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_rhs_var, getScratchContext());
        const int T_diff_coef_current_idx =
            var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_var, getCurrentContext());
        const int T_diff_coef_rhs_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_rhs_var, getScratchContext());
        const int T_diff_coef_cc_current_idx =
            var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_cc_var, getCurrentContext());
        const int T_diff_coef_cc_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_cc_var, getScratchContext());

        if (d_solve_mass_conservation)
        {
            Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
            const int u_current_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getCurrentContext());
            // Set the velocities used to update the density and the previous time step
            // size
            if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ -1, /*current*/ u_current_idx, /*new*/ -1);
                rho_p_cc_integrator->setMaterialPropertyPatchDataIndices(
                    /*current*/ Cp_current_idx, /*new*/ -1);
                rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                    /*current*/ T_current_idx, /*new*/ -1);
            }
            else
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ d_U_old_current_idx, /*current*/ u_current_idx, /*new*/ -1);
                rho_p_cc_integrator->setMaterialPropertyPatchDataIndices(
                    /*current*/ Cp_current_idx, /*new*/ -1);
                rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                    /*current*/ T_current_idx, /*new*/ -1);
                d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
            }

            // Integrate density and convective term of energy equation.
            d_rho_p_integrator->integrate(dt);
        }
        // Setup the problem coefficients for the linear solve
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

        // set rho*Cp/dt
        d_hier_cc_data_ops->multiply(d_C_current_idx, rho_current_idx, Cp_current_idx);
        d_hier_cc_data_ops->scale(d_C_current_idx, 1.0 / dt, d_C_current_idx);

        d_hier_cc_data_ops->copyData(d_T_temp_rhs_idx, d_C_current_idx);
        T_rhs_op_spec.setCPatchDataId(d_T_temp_rhs_idx);

        const double apply_time = current_time;
        for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
        {
            d_reset_kappa_fcns[k](T_diff_coef_cc_current_idx,
                                  d_T_diffusion_coef_cc_var,
                                  d_hier_math_ops,
                                  -1 /*cycle_num*/,
                                  apply_time,
                                  current_time,
                                  new_time,
                                  d_reset_kappa_fcns_ctx[k]);
        }

        // Interpolate the cell-centered diffusion coef to side-centered.
        d_hier_cc_data_ops->copyData(T_diff_coef_cc_scratch_idx, T_diff_coef_cc_current_idx);
        d_k_bdry_bc_fill_op->fillData(current_time);

        if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
        {
            interpolateCCToSCSimpleAveraging(T_diff_coef_current_idx, T_diff_coef_cc_scratch_idx);
        }
        else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
        {
            interpolateCCToSCHarmonicAveraging(T_diff_coef_current_idx, T_diff_coef_cc_scratch_idx);
        }
        else
        {
            TBOX_ERROR("this statement should not be reached");
        }

        // for plotting purpose.
        static const bool synch_cf_interface = true;
        d_hier_math_ops->interp(d_D_cc_new_idx,
                                d_D_cc_var,
                                T_diff_coef_current_idx,
                                d_T_diffusion_coef_var,
                                d_no_fill_op,
                                d_integrator_time,
                                synch_cf_interface);

        d_hier_sc_data_ops->scale(T_diff_coef_rhs_scratch_idx, (1.0 - alpha), T_diff_coef_current_idx);
        T_rhs_op_spec.setDPatchDataId(T_diff_coef_rhs_scratch_idx);

        // Initialize the RHS operator and compute the RHS vector for temperature equation.
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
        d_hier_cc_data_ops->copyData(T_scratch_idx, T_current_idx, false);
        T_rhs_op->apply(*d_T_sol, *d_T_rhs);

        d_hier_cc_data_ops->copyData(T_new_idx, T_current_idx);

        if (d_u_adv_var)
        {
            // Add div (u H lf).
            d_hier_cc_data_ops->scale(d_T_lf_N_scratch_idx, d_rho_liquid * d_latent_heat, lf_N_scratch_idx);
            if (d_lf_convective_time_stepping_type == FORWARD_EULER)
            {
                d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -1.0, d_T_lf_N_scratch_idx, T_rhs_scratch_idx);
            }
            else if (d_lf_convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -0.5, d_T_lf_N_scratch_idx, T_rhs_scratch_idx);
            }
        }
    }
    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
AllenCahnHierarchyIntegrator::integrateHierarchy(const double current_time, const double new_time, const int cycle_num)
{
    PhaseChangeHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE({
            pout << "AllenCahnHierarchyIntegrator::integrateHierarchy():\n"
                 << "  WARNING: num_cycles = " << d_current_num_cycles
                 << " but expected num_cycles = " << expected_num_cycles << ".\n";
        });
    }

    // Perform a single step of fixed point iteration.
    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
    const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
    const int lf_current_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getCurrentContext());
    const int lf_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_rhs_var, getScratchContext());
    const int lf_diff_coef_new_idx = (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_var, getNewContext()));
    const int lf_diff_coef_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_var, getScratchContext()));
    const int lf_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_N_var, getScratchContext());
    const int H_current_idx = var_db->mapVariableAndContextToIndex(d_H_var, getCurrentContext());
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    const int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());

    // update C coefficients.
    if (cycle_num > 0) computeDoubleWellPotential(d_g_firstder_idx, d_g_secondder_idx, lf_new_idx);

    // Add H/dt + M lambda / eps^2 H g''. Double well potential term is linearized and added in both LHS and RHS.
    d_hier_cc_data_ops->scale(d_lf_C_idx, 1.0 / dt, H_new_idx);
    d_hier_cc_data_ops->multiply(d_C_scratch_idx, d_g_secondder_idx, H_new_idx);
    d_hier_cc_data_ops->scale(d_C_scratch_idx, d_M_lf * d_lambda_lf / std::pow(d_eps_lf, 2.0), d_C_scratch_idx);
    d_hier_cc_data_ops->add(d_lf_C_idx, d_C_scratch_idx, d_lf_C_idx);

    PoissonSpecifications lf_solver_spec(d_object_name + "::solver_spec::" + d_lf_var->getName());
    lf_solver_spec.setCPatchDataId(d_lf_C_idx);

    // update D coefficients.
    // setting Allen-Cahn equation diffusion timestepping type
    double alpha = 0.0;
    switch (d_lf_diffusion_time_stepping_type)
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
        TBOX_ERROR(d_object_name << "::preprocessintegrateHierarchy():\n"
                                 << "  unsupported diffusion time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_lf_diffusion_time_stepping_type) << " \n"
                                 << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }

    // Interpolate the cell-centered Heaviside to side-centered.
    d_hier_cc_data_ops->copyData(H_scratch_idx, H_new_idx);
    d_H_bdry_bc_fill_op->fillData(new_time);

    interpolateCCToSCSimpleAveraging(lf_diff_coef_new_idx, H_scratch_idx);

    // Add numerical diffusion to smoothly extend the liquid fraction.
    d_hier_sc_data_ops->addScalar(lf_diff_coef_new_idx, lf_diff_coef_new_idx, d_num_diffusion);
    d_hier_sc_data_ops->copyData(d_H_sc_idx, lf_diff_coef_new_idx);
    d_hier_sc_data_ops->scale(lf_diff_coef_new_idx, d_M_lf * d_lambda_lf, lf_diff_coef_new_idx);

    d_hier_sc_data_ops->scale(lf_diff_coef_scratch_idx, -alpha, lf_diff_coef_new_idx);
    lf_solver_spec.setDPatchDataId(lf_diff_coef_scratch_idx);

    // Initialize the linear solver for lf equation.
    Pointer<PoissonSolver> lf_solver = d_lf_solver;
    lf_solver->setPoissonSpecifications(lf_solver_spec);
    lf_solver->setPhysicalBcCoef(d_lf_bc_coef);
    lf_solver->setHomogeneousBc(false);
    lf_solver->setSolutionTime(new_time);
    lf_solver->setTimeInterval(current_time, new_time);
    // Initializing solver every time step.
    if (d_enable_logging)
    {
        plog << d_object_name << ": "
             << "Initializing the solvers for" << d_lf_var->getName() << "\n";
    }
    lf_solver->initializeSolverState(*d_lf_sol, *d_lf_rhs);
    d_lf_solver_needs_init = false;

    if (d_u_adv_var)
    {
        if (cycle_num > 0)
        {
            // Update the advection velocity for lf.
            const int u_current_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getCurrentContext());
            const int u_scratch_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getScratchContext());
            const int u_new_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getNewContext());
            Pointer<IBTK::CartGridFunction> u_fcn = d_u_fcn[d_u_adv_var];
            if (u_fcn)
            {
                u_fcn->setDataOnPatchHierarchy(u_new_idx, d_u_adv_var, d_hierarchy, new_time);
            }
            d_hier_fc_data_ops->linearSum(u_scratch_idx, 0.5, u_current_idx, 0.5, u_new_idx);
        }
    }

    // Account for the convective difference term.
    const int lf_H_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_H_var, getScratchContext());
    TimeSteppingType convective_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
    if (d_u_adv_var)
    {
        convective_time_stepping_type = d_lf_convective_time_stepping_type;
        if (is_multistep_time_stepping_type(convective_time_stepping_type))
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(convective_time_stepping_type == ADAMS_BASHFORTH);
#endif
            if (getIntegratorStep() == 0)
            {
                convective_time_stepping_type = d_lf_init_convective_time_stepping_type;
            }
            else if (cycle_num > 0)
            {
                convective_time_stepping_type = MIDPOINT_RULE;
                IBAMR_DO_ONCE({
                    pout << "AllenCahnHierarchyIntegrator::"
                            "integrateHierarchy():"
                            "\n"
                         << "  WARNING: convective_time_stepping_type = "
                         << enum_to_string<TimeSteppingType>(d_lf_convective_time_stepping_type)
                         << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                         << "           using " << enum_to_string<TimeSteppingType>(d_lf_convective_time_stepping_type)
                         << " only for the first cycle in each time step;\n"
                         << "           using " << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                         << " for subsequent cycles.\n";
                });
            }
        }

        if (cycle_num > 0)
        {
            if (convective_time_stepping_type == MIDPOINT_RULE)
            {
                const int u_scratch_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getScratchContext());
                d_lf_convective_op->setAdvectionVelocity(u_scratch_idx);
                d_hier_cc_data_ops->linearSum(lf_scratch_idx, 0.5, lf_current_idx, 0.5, lf_new_idx);
                // H_pre_idx is used to use the same flux used in the advection of H.
                d_hier_cc_data_ops->linearSum(H_scratch_idx, 0.5, H_current_idx, 0.5, d_H_pre_idx);
                d_lf_convective_op->setSolutionTime(half_time);

                // Checking whether the CellConvectiveOperator can be used for div (H lf u)
                d_lf_convective_op->interpolateToFaceOnHierarchy(
                    d_lf_interp_idx, lf_scratch_idx, u_scratch_idx, /*synch_cf_bdry*/ false);
                d_lf_convective_op->interpolateToFaceOnHierarchy(
                    d_H_interp_idx, H_scratch_idx, u_scratch_idx, /*synch_cf_bdry*/ false);
                d_hier_fc_data_ops->multiply(d_lf_flux_idx, d_lf_interp_idx, d_H_interp_idx);
                d_hier_fc_data_ops->multiply(d_lf_flux_idx, d_lf_flux_idx, u_scratch_idx);
                d_lf_convective_op->computeConservativeDerivativeOnHierarchy(
                    lf_N_scratch_idx, d_lf_flux_idx, /*synch_cf_bdry*/ true);

                // Pointer<AdvDiffCUIConservativeConvectiveOperator> lf_cui_conservative_convective_op =
                //     d_lf_convective_op;
                // lf_cui_conservative_convective_op->applyConvectiveOperator(
                //     lf_scratch_idx, H_scratch_idx, lf_N_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                const int u_new_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getNewContext());
                d_lf_convective_op->setAdvectionVelocity(u_new_idx);
                d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_new_idx);
                // H_pre_idx is used to use the same flux used in the advection of H.
                d_hier_cc_data_ops->copyData(H_scratch_idx, d_H_pre_idx);
                d_lf_convective_op->setSolutionTime(new_time);

                // Checking whether the CellConvectiveOperator can be used for div (H lf u)
                d_lf_convective_op->interpolateToFaceOnHierarchy(
                    d_lf_interp_idx, lf_scratch_idx, u_new_idx, /*synch_cf_bdry*/ false);
                d_lf_convective_op->interpolateToFaceOnHierarchy(
                    d_H_interp_idx, H_scratch_idx, u_new_idx, /*synch_cf_bdry*/ false);
                d_hier_fc_data_ops->multiply(d_lf_flux_idx, d_lf_interp_idx, d_H_interp_idx);
                d_hier_fc_data_ops->multiply(d_lf_flux_idx, d_lf_flux_idx, u_new_idx);
                d_lf_convective_op->computeConservativeDerivativeOnHierarchy(
                    lf_N_scratch_idx, d_lf_flux_idx, /*synch_cf_bdry*/ true);

                // Pointer<AdvDiffCUIConservativeConvectiveOperator> lf_cui_conservative_convective_op =
                //     d_lf_convective_op;
                // lf_cui_conservative_convective_op->applyConvectiveOperator(
                //     lf_scratch_idx, H_scratch_idx, lf_N_scratch_idx);
            }
        }
        if (convective_time_stepping_type == ADAMS_BASHFORTH)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(cycle_num == 0);
#endif
            const int lf_N_old_current_idx = var_db->mapVariableAndContextToIndex(d_lf_N_old_var, getCurrentContext());
            const double omega = dt / d_dt_previous[0];
            d_hier_cc_data_ops->linearSum(
                lf_N_scratch_idx, 1.0 + 0.5 * omega, lf_N_scratch_idx, -0.5 * omega, lf_N_old_current_idx);
        }

        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, -1.0, lf_N_scratch_idx, lf_rhs_scratch_idx);
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, -0.5, lf_N_scratch_idx, lf_rhs_scratch_idx);
        }
    }

    // Account for forcing terms.
    const int lf_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_F_var, getScratchContext());
    const int lf_F_new_idx = var_db->mapVariableAndContextToIndex(d_lf_F_var, getNewContext());

    computeInterpolationFunction(d_q_firstder_idx, lf_new_idx, d_T_new_idx);
    computeLiquidFractionSourceTerm(lf_F_scratch_idx);

    if (d_u_adv_var)
    {
        // compute H \varphi Div u.
        int Div_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_Div_u_var, getScratchContext());
        const int u_new_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getNewContext());
        d_hier_math_ops->div(Div_u_scratch_idx,
                             d_Div_u_var,
                             1.0,
                             u_new_idx,
                             d_u_adv_var,
                             d_no_fill_op,
                             d_integrator_time,
                             /*synch_cf_bdry*/ false);

        d_hier_cc_data_ops->multiply(lf_H_scratch_idx, lf_new_idx, H_new_idx);
        d_hier_cc_data_ops->multiply(Div_u_scratch_idx, lf_H_scratch_idx, Div_u_scratch_idx);
        d_hier_cc_data_ops->axpy(lf_F_scratch_idx, +1.0, Div_u_scratch_idx, lf_F_scratch_idx);
    }
    d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, +1.0, lf_F_scratch_idx, lf_rhs_scratch_idx);

    // This will be used while computing chemical potential.
    d_hier_cc_data_ops->copyData(d_lf_pre_idx, lf_new_idx);

    // Solve for lf(n+1).
    lf_solver->solveSystem(*d_lf_sol, *d_lf_rhs);
    d_hier_cc_data_ops->copyData(lf_new_idx, lf_scratch_idx);
    if (d_enable_logging && d_enable_logging_solver_iterations)
        plog << d_object_name << ":" << d_lf_var->getName()
             << "::integrateHierarchy():diffusion solve number of iterations = " << lf_solver->getNumIterations()
             << "\n";
    if (d_enable_logging)
        plog << d_object_name << ":" << d_lf_var->getName()
             << "::integrateHierarchy():diffusion solve residual norm        = " << lf_solver->getResidualNorm()
             << "\n";
    if (lf_solver->getNumIterations() == lf_solver->getMaxIterations())
    {
        pout << d_object_name << ":" << d_lf_var->getName()
             << "::integrateHierarchy():WARNING: linear solver iterations == max iterations\n";
    }

    // Reset the right-hand side vector.
    if (d_u_adv_var)
    {
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, +1.0, lf_N_scratch_idx, lf_rhs_scratch_idx);
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, +0.5, lf_N_scratch_idx, lf_rhs_scratch_idx);
        }
    }

    if (d_lf_F_var)
    {
        d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, -1.0, lf_F_scratch_idx, lf_rhs_scratch_idx);
        d_hier_cc_data_ops->copyData(lf_F_new_idx, lf_F_scratch_idx);
    }

    // To use same H used in advection of H in the convective term.
    d_hier_cc_data_ops->copyData(d_H_pre_idx, H_new_idx);

    // Bound the updated liquid fraction.
    boundLiquidFraction(lf_new_idx);

    if (d_solve_energy)
    {
        const int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());
        const int rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());

        const int Cp_new_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getNewContext());
        const int Cp_current_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getCurrentContext());

        const int T_current_idx = var_db->mapVariableAndContextToIndex(d_T_var, getCurrentContext());
        const int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, getScratchContext());
        const int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());
        const int T_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_rhs_var, getScratchContext());
        const int T_diff_coef_new_idx = (var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_var, getNewContext()));
        const int T_diff_coef_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_var, getScratchContext());
        const int T_diff_coef_cc_new_idx =
            var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_cc_var, getNewContext());
        const int T_diff_coef_cc_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_cc_var, getScratchContext());

        // update specific heat.
        double apply_time = new_time;
        for (unsigned k = 0; k < d_reset_Cp_fcns.size(); ++k)
        {
            d_reset_Cp_fcns[k](Cp_new_idx,
                               d_Cp_var,
                               d_hier_math_ops,
                               -1 /*cycle_num*/,
                               apply_time,
                               current_time,
                               new_time,
                               d_reset_Cp_fcns_ctx[k]);
        }

        // In the special case of a conservative discretization form, the updated
        // density is calculated by application of the mass and convective
        // momentum integrator.
        const int T_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_N_var, getScratchContext());

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
            d_rho_p_integrator->setConvectiveDerivativePatchDataIndex(T_N_scratch_idx);

            // Always set to current because we want to update rho^{n} to rho^{n+1}
            d_rho_p_integrator->setDensityPatchDataIndex(rho_current_idx);

            const int u_current_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getCurrentContext());
            const int u_new_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getNewContext());

            // Set the velocities used to update the density
            if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ -1, /*current*/ u_current_idx, /*new*/ u_new_idx);
                rho_p_cc_integrator->setMaterialPropertyPatchDataIndices(
                    /*current*/ Cp_current_idx, /*new*/ Cp_new_idx);
                rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                    /*current*/ T_current_idx, /*new*/ T_new_idx);
            }
            else
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ d_U_old_current_idx,
                    /*current*/ u_current_idx,
                    /*new*/ u_new_idx);
                rho_p_cc_integrator->setMaterialPropertyPatchDataIndices(
                    /*current*/ Cp_current_idx,
                    /*new*/ Cp_new_idx);
                rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                    /*current*/ T_current_idx,
                    /*new*/ T_new_idx);

                d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
            }

            d_rho_p_integrator->integrate(dt);
        }

        d_updated_rho_idx = d_rho_p_integrator ? d_rho_p_integrator->getUpdatedDensityPatchDataIndex() : rho_new_idx;
        d_hier_cc_data_ops->copyData(rho_new_idx,
                                     d_updated_rho_idx,
                                     /*interior_only*/ true);

        PoissonSpecifications T_solver_spec(d_object_name + "::solver_spec::" + d_T_var->getName());
        // set rho*Cp/dt
        d_hier_cc_data_ops->multiply(d_C_new_idx, rho_new_idx, Cp_new_idx);
        d_hier_cc_data_ops->scale(d_C_new_idx, 1.0 / dt, d_C_new_idx);
        d_hier_cc_data_ops->copyData(d_T_C_idx, d_C_new_idx);
        T_solver_spec.setCPatchDataId(d_T_C_idx);

        // Setup the problem coefficients for the linear solve
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

        apply_time = new_time;
        for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
        {
            d_reset_kappa_fcns[k](T_diff_coef_cc_new_idx,
                                  d_T_diffusion_coef_cc_var,
                                  d_hier_math_ops,
                                  -1 /*cycle_num*/,
                                  apply_time,
                                  current_time,
                                  new_time,
                                  d_reset_kappa_fcns_ctx[k]);
        }

        // Interpolate the cell-centered diffusion coef to side-centered.
        d_hier_cc_data_ops->copyData(T_diff_coef_cc_scratch_idx, T_diff_coef_cc_new_idx);
        d_k_bdry_bc_fill_op->fillData(new_time);

        if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
        {
            interpolateCCToSCSimpleAveraging(T_diff_coef_new_idx, T_diff_coef_cc_scratch_idx);
        }
        else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
        {
            interpolateCCToSCHarmonicAveraging(T_diff_coef_new_idx, T_diff_coef_cc_scratch_idx);
        }
        else
        {
            TBOX_ERROR("this statement should not be reached");
        }

        // for plotting purpose.
        static const bool synch_cf_interface = true;
        d_hier_math_ops->interp(d_D_cc_new_idx,
                                d_D_cc_var,
                                T_diff_coef_new_idx,
                                d_T_diffusion_coef_var,
                                d_no_fill_op,
                                d_integrator_time,
                                synch_cf_interface);

        d_hier_sc_data_ops->scale(T_diff_coef_scratch_idx, -alpha, T_diff_coef_new_idx);
        T_solver_spec.setDPatchDataId(T_diff_coef_scratch_idx);

        // Initialize the linear solver for temperature equation.
        Pointer<PoissonSolver> T_solver = d_T_solver;
        T_solver->setPoissonSpecifications(T_solver_spec);
        T_solver->setPhysicalBcCoef(d_T_bc_coef);
        T_solver->setHomogeneousBc(false);
        T_solver->setSolutionTime(new_time);
        T_solver->setTimeInterval(current_time, new_time);
        // Initializing solver every time step.
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the solvers for" << d_T_var->getName() << "\n";
        }
        T_solver->initializeSolverState(*d_T_sol, *d_T_rhs);
        d_T_solver_needs_init = true;

        if (d_u_adv_var)
        {
            // Account for the convective term computed from STSMassFluxIntegrator class.
            d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -1.0, T_N_scratch_idx, T_rhs_scratch_idx);

            // Add div (u H lf).
            d_hier_cc_data_ops->scale(d_T_lf_N_scratch_idx, d_rho_liquid * d_latent_heat, lf_N_scratch_idx);
            if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -1.0, d_T_lf_N_scratch_idx, T_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -0.5, d_T_lf_N_scratch_idx, T_rhs_scratch_idx);
            }
        }

        // Account for forcing terms.
        const int T_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_F_var, getScratchContext());
        const int T_F_new_idx = var_db->mapVariableAndContextToIndex(d_T_F_var, getNewContext());
        if (d_T_F_fcn)
        {
            d_T_F_fcn->setDataOnPatchHierarchy(T_F_scratch_idx, d_T_F_var, d_hierarchy, half_time);
        }
        else
            d_hier_cc_data_ops->setToScalar(T_F_scratch_idx, 0.0);

        // Add Allen-Cahn temporal term.
        addTemporalAndLinearTermstoRHSOfEnergyEquation(T_F_scratch_idx, dt);
        d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, +1.0, T_F_scratch_idx, T_rhs_scratch_idx);

        // Solve for T(n+1).
        T_solver->solveSystem(*d_T_sol, *d_T_rhs);
        d_hier_cc_data_ops->copyData(T_new_idx, T_scratch_idx);
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

        // Compute the source term for a Div U equation
        computeDivergenceVelocitySourceTerm(d_Div_U_F_idx, new_time);

        // Reset the right-hand side vector.
        if (d_u_adv_var)
        {
            d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, +1.0, T_N_scratch_idx, T_rhs_scratch_idx);

            if (d_lf_convective_time_stepping_type == ADAMS_BASHFORTH ||
                d_lf_convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, +1.0, d_T_lf_N_scratch_idx, T_rhs_scratch_idx);
            }
            else if (d_lf_convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, +0.5, d_T_lf_N_scratch_idx, T_rhs_scratch_idx);
            }
        }

        if (d_T_F_var)
        {
            d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -1.0, T_F_scratch_idx, T_rhs_scratch_idx);
            d_hier_cc_data_ops->copyData(T_F_new_idx, T_F_scratch_idx);
        }
    }
    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
AllenCahnHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
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
        level->deallocatePatchData(d_lf_C_idx);
        level->deallocatePatchData(d_lf_temp_rhs_idx);
        level->deallocatePatchData(d_g_firstder_idx);
        level->deallocatePatchData(d_g_secondder_idx);
        level->deallocatePatchData(d_q_firstder_idx);
        level->deallocatePatchData(d_grad_lf_idx);
        level->deallocatePatchData(d_H_sc_idx);
        level->deallocatePatchData(d_lf_interp_idx);
        level->deallocatePatchData(d_H_interp_idx);
        level->deallocatePatchData(d_lf_flux_idx);
    }

    PhaseChangeHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
AllenCahnHierarchyIntegrator::registerLiquidFractionVariable(Pointer<CellVariable<NDIM, double> > lf_var,
                                                             const bool output_lf_var)
{
    d_lf_var = lf_var;
    d_output_lf = output_lf_var;

    Pointer<CellDataFactory<NDIM, double> > lf_factory = lf_var->getPatchDataFactory();
    const int lf_depth = lf_factory->getDefaultDepth();
    Pointer<CellVariable<NDIM, double> > lf_rhs_var =
        new CellVariable<NDIM, double>(lf_var->getName() + "::lf_rhs", lf_depth);
    Pointer<CellVariable<NDIM, double> > lf_F_var = new CellVariable<NDIM, double>(lf_var->getName() + "::F", lf_depth);
    Pointer<SideVariable<NDIM, double> > lf_diff_coef_var =
        new SideVariable<NDIM, double>(lf_var->getName() + "::diff_coef", lf_depth);

    // Set default values.
    d_u_adv_var = nullptr;
    d_lf_F_var = lf_F_var;
    d_lf_rhs_var = lf_rhs_var;
    d_lf_diffusion_coef_var = lf_diff_coef_var;
    d_lf_init = nullptr;
    d_lf_F_fcn = nullptr;
    d_lf_bc_coef = nullptr;
    return;
} // registerLiquidFractionVariable

Pointer<PoissonSolver>
AllenCahnHierarchyIntegrator::getAllenCahnEquationHelmholtzSolver(Pointer<CellVariable<NDIM, double> > lf_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_lf_var);
#endif
    if (!d_lf_solver)
    {
        const std::string& name = lf_var->getName();
        d_lf_solver =
            CCPoissonSolverManager::getManager()->allocateSolver(d_lf_solver_type,
                                                                 d_object_name + "::helmholtz_solver::" + name,
                                                                 d_lf_solver_db,
                                                                 "liquid_fraction_",
                                                                 d_lf_precond_type,
                                                                 d_object_name + "::helmholtz_precond::" + name,
                                                                 d_lf_precond_db,
                                                                 "liquid_fraction_pc_");
        d_lf_solver_needs_init = true;
    }
    return d_lf_solver;
} // getAllenCahnEquationHelmholtzSolver

Pointer<LaplaceOperator>
AllenCahnHierarchyIntegrator::getAllenCahnEquationHelmholtzRHSOperator(Pointer<CellVariable<NDIM, double> > lf_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_lf_var);
#endif
    const std::string& name = lf_var->getName();
    if (!d_lf_rhs_op)
    {
        d_lf_rhs_op = new CCLaplaceOperator(d_object_name + "::helmholtz_rhs_op::" + name, /*homogeneous_bc*/ false);
        d_lf_rhs_op_needs_init = true;
    }
    return d_lf_rhs_op;
} // getAllenCahnEquationHelmholtzRHSOperator

void
AllenCahnHierarchyIntegrator::setLiquidFractionPhysicalBcCoef(Pointer<CellVariable<NDIM, double> > lf_var,
                                                              RobinBcCoefStrategy<NDIM>* lf_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(lf_var);
#endif
    d_lf_bc_coef = lf_bc_coef;
    return;
} // setLiquidFractionPhysicalBcCoef

RobinBcCoefStrategy<NDIM>*
AllenCahnHierarchyIntegrator::getLiquidFractionPhysicalBcCoef()
{
    return d_lf_bc_coef;
} // getLiquidFractionPhysicalBcCoef

Pointer<CellConvectiveOperator>
AllenCahnHierarchyIntegrator::getAllenCahnEquationConvectiveOperator(Pointer<CellVariable<NDIM, double> > lf_var,
                                                                     Pointer<CellVariable<NDIM, double> > H_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(lf_var);
#endif
    if (!d_lf_convective_op)
    {
        std::vector<RobinBcCoefStrategy<NDIM>*> lf_bc_coef(1, d_lf_bc_coef);
        std::vector<RobinBcCoefStrategy<NDIM>*> H_bc_coef = getPhysicalBcCoefs(H_var);

        // Since the Allen-Cahn equation requires the convective derivative div(lf*H*u),
        // we use AdvDiffCUIConservativeConvectiveOperator class.
        d_lf_convective_op = new AdvDiffCUIConvectiveOperator(d_object_name + "::lfConvectiveOperator",
                                                              lf_var,
                                                              d_lf_convective_op_input_db,
                                                              d_lf_convective_difference_form,
                                                              lf_bc_coef);
        d_lf_convective_op_needs_init = true;
    }
    return d_lf_convective_op;
} // getAllenCahnEquationConvectiveOperator

void
AllenCahnHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    db->putInteger("IEP_HIERARCHY_INTEGRATOR_VERSION", IEP_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_lf_diffusion_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_lf_diffusion_time_stepping_type));
    db->putString("d_lf_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_lf_convective_time_stepping_type));
    db->putString("d_lf_init_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_lf_init_convective_time_stepping_type));
    db->putString("d_lf_convective_difference_form",
                  enum_to_string<ConvectiveDifferencingType>(d_lf_convective_difference_form));
    db->putString("d_lf_convective_op_type", d_lf_convective_op_type);
    d_lf_convective_op_input_db = db->putDatabase("d_lf_convective_op_db");

    db->putDouble("d_M_lf", d_M_lf);
    db->putDouble("d_lambda_lf", d_lambda_lf);
    db->putDouble("d_eps_lf", d_eps_lf);
    db->putBool("d_solve_energy", d_solve_energy);
    db->putDouble("d_num_diffusion", d_num_diffusion);
    db->putString("d_interpolation_function_profile", d_interpolation_function_profile);

    AdvDiffSemiImplicitHierarchyIntegrator::putToDatabaseSpecialized(db);
    return;
} // putToDatabaseSpecialized

void
AllenCahnHierarchyIntegrator::addTemporalAndLinearTermstoRHSOfEnergyEquation(int F_scratch_idx, const double dt)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    int H_current_idx = var_db->mapVariableAndContextToIndex(d_H_var, getCurrentContext());
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > lf_new_data = patch->getPatchData(d_lf_new_idx);
            Pointer<CellData<NDIM, double> > H_new_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double> > lf_current_data = patch->getPatchData(d_lf_current_idx);
            Pointer<CellData<NDIM, double> > H_current_data = patch->getPatchData(H_current_idx);
            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_scratch_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*F_data)(ci) +=
                    -d_rho_liquid * d_latent_heat *
                    ((((*H_new_data)(ci) * (*lf_new_data)(ci)) - ((*H_current_data)(ci) * (*lf_current_data)(ci))) /
                     dt);
            }
        }
    }
    return;
} // addTemporalAndLinearTermstoRHSOfEnergyEquation

void
AllenCahnHierarchyIntegrator::computeDivergenceVelocitySourceTerm(int Div_U_F_idx, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
    const int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    const int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());

    // Filling ghost cells for liquid fraction.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> lf_transaction_comps(1);
    lf_transaction_comps[0] = InterpolationTransactionComponent(lf_scratch_idx,
                                                                lf_new_idx,
                                                                "CONSERVATIVE_LINEAR_REFINE",
                                                                false,
                                                                "CONSERVATIVE_COARSEN",
                                                                "LINEAR",
                                                                false,
                                                                d_lf_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(lf_transaction_comps, d_hierarchy);
    hier_bdry_fill->fillData(new_time);

    // perform gradient of liquid fraction.
    d_hier_math_ops->grad(d_grad_lf_idx, d_grad_lf_var, true, 1.0, lf_scratch_idx, d_lf_var, nullptr, new_time);

    // compute H*grad_lf.
    d_hier_sc_data_ops->multiply(d_grad_lf_idx, d_grad_lf_idx, d_H_sc_idx);

    // compute div(H*grad_lf).
    d_hier_math_ops->div(d_chemical_potential_idx,
                         d_chemical_potential_var,
                         1.0,
                         d_grad_lf_idx,
                         d_grad_lf_var,
                         nullptr,
                         new_time,
                         false);

    // update q' and g'.
    computeInterpolationFunction(d_q_firstder_idx, lf_new_idx, T_new_idx);
    computeDoubleWellPotential(d_g_firstder_idx, d_g_secondder_idx, lf_new_idx);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_new_idx);
            Pointer<CellData<NDIM, double> > q_firstder_data = patch->getPatchData(d_q_firstder_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_new_idx);
            Pointer<CellData<NDIM, double> > g_firstder_data = patch->getPatchData(d_g_firstder_idx);
            Pointer<CellData<NDIM, double> > g_secondder_data = patch->getPatchData(d_g_secondder_idx);
            Pointer<CellData<NDIM, double> > chemical_potential_data = patch->getPatchData(d_chemical_potential_idx);
            Pointer<CellData<NDIM, double> > Div_U_F_data = patch->getPatchData(Div_U_F_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                double F = d_rho_liquid * d_latent_heat * (*H_data)(ci) * (*q_firstder_data)(ci) *
                           (d_T_melt - (*T_data)(ci)) / d_T_melt;

                // g' is linearized while solving AC eqn. But for Div U eqn, it is computed based on
                // varphi^n+1,m+1.
                (*chemical_potential_data)(ci) =
                    -d_lambda_lf * (*chemical_potential_data)(ci) +
                    (d_lambda_lf * (*H_data)(ci) / std::pow(d_eps_lf, 2.0) * (*g_firstder_data)(ci)) + F;

                // Div_U_F = M_lf * zeta_lf * (rho_liquid - rho_solid) / rho
                (*Div_U_F_data)(ci) =
                    d_M_lf * (*chemical_potential_data)(ci) * (d_rho_liquid - d_rho_solid) / (*rho_data)(ci);
            }
        }
    }

    return;
} // computeDivergenceVelocitySourceTerm

/////////////////////////////// PROTECTED ////////////////////////////////////

void
AllenCahnHierarchyIntegrator::regridHierarchyBeginSpecialized()
{
    PhaseChangeHierarchyIntegrator::regridHierarchyBeginSpecialized();

    d_lf_rhs_op->deallocateOperatorState();
    d_lf_solver->deallocateSolverState();

    d_lf_solver_needs_init = true;
    d_lf_rhs_op_needs_init = true;
    d_lf_convective_op_needs_init = true;

    return;
} // regridHierarchyBeginSpecialized

void
AllenCahnHierarchyIntegrator::regridHierarchyEndSpecialized()
{
    PhaseChangeHierarchyIntegrator::regridHierarchyEndSpecialized();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    std::vector<RobinBcCoefStrategy<NDIM>*> H_bc_coef = getPhysicalBcCoefs(d_H_var);

    // Setup the patch boundary filling objects.
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();

    // Reset the solution and rhs vectors.
    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
    d_lf_sol = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::sol_vec::" + d_lf_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_lf_sol->addComponent(d_lf_var, lf_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    const int lf_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_rhs_var, getScratchContext());
    d_lf_rhs = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::rhs_vec::" + d_lf_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_lf_rhs->addComponent(d_lf_rhs_var, lf_rhs_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    d_lf_solver_needs_init = true;
    d_lf_rhs_op_needs_init = true;
    d_lf_convective_op_needs_init = true;
    return;
} // regridHierarchyEndSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AllenCahnHierarchyIntegrator::computeDoubleWellPotential(int g_firstder_idx,
                                                         int g_secondder_idx,
                                                         const int liquid_fraction_idx)
{
    const bool use_convex_splitting = true;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(liquid_fraction_idx);
            Pointer<CellData<NDIM, double> > g_firstder_data = patch->getPatchData(g_firstder_idx);
            Pointer<CellData<NDIM, double> > g_secondder_data = patch->getPatchData(g_secondder_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double lf = (*lf_data)(ci);
                if (use_convex_splitting)
                {
                    if (lf < 0.0)
                        (*g_firstder_data)(ci) = 2.0 * lf;
                    else if (lf > 1.0)
                        (*g_firstder_data)(ci) = 2.0 * (lf - 1.0);
                    else
                        (*g_firstder_data)(ci) = 2.0 * (lf - 1.0) * lf * (2.0 * lf - 1.0);

                    (*g_secondder_data)(ci) = 2.0;
                }
                else
                {
                    (*g_firstder_data)(ci) = 2.0 * (lf - 1.0) * lf * (2.0 * lf - 1.0);
                    (*g_secondder_data)(ci) = 12.0 * lf * lf - 12.0 * lf + 2.0;
                }
            }
        }
    }
    return;
} //  computeDoubleWellPotential

void
AllenCahnHierarchyIntegrator::computeInterpolationFunction(int q_firstder_idx,
                                                           const int liquid_fraction_idx,
                                                           const int T_idx)
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

            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(liquid_fraction_idx);
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double> > q_firstder_data = patch->getPatchData(q_firstder_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double lf = (*lf_data)(ci);
                const double T = (*T_data)(ci);

                if (d_interpolation_function_profile == "QUARTIC")
                {
                    // Ziyang's profile
                    (*q_firstder_data)(ci) = 30.0 * std::pow(lf, 4.0) - 60.0 * std::pow(lf, 3.0) + 30.0 * lf * lf;
                }
                else if (d_interpolation_function_profile == "QUADRATIC")
                {
                    // Li's profile
                    (*q_firstder_data)(ci) = 6.0 * lf - 6.0 * std::pow(lf, 2.0);
                }
                else if (d_interpolation_function_profile == "LINEAR_3")
                {
                    // Hybrid profile
                    // 3lf and 3-3lf
                    if (lf <= 0.133048682404023)
                        (*q_firstder_data)(ci) = 3.0 * lf;
                    else if (lf > 0.133048682404023 && lf <= 0.866951317595975)
                        (*q_firstder_data)(ci) = 30.0 * std::pow(lf, 4.0) - 60.0 * std::pow(lf, 3.0) + 30.0 * lf * lf;
                    else
                        (*q_firstder_data)(ci) = 3.0 - 3.0 * lf;
                }
                else if (d_interpolation_function_profile == "LINEAR_4")
                {
                    // 4lf and 4-4lf
                    if (lf <= 0.218078018145755)
                        (*q_firstder_data)(ci) = 4.0 * lf;
                    else if (lf > 0.218078018145755 && lf <= 0.781921981854249)
                        (*q_firstder_data)(ci) = 30.0 * std::pow(lf, 4.0) - 60.0 * std::pow(lf, 3.0) + 30.0 * lf * lf;
                    else
                        (*q_firstder_data)(ci) = 4.0 - 4.0 * lf;
                }
                else if (d_interpolation_function_profile == "LINEAR_2")
                {
                    // 2lf and 2-2lf
                    if (lf <= 0.0785105470229893)
                        (*q_firstder_data)(ci) = 2.0 * lf;
                    else if (lf > 0.0785105470229893 && lf <= 0.921489452977004)
                        (*q_firstder_data)(ci) = 30.0 * std::pow(lf, 4.0) - 60.0 * std::pow(lf, 3.0) + 30.0 * lf * lf;
                    else
                        (*q_firstder_data)(ci) = 2.0 - 2.0 * lf;
                }
                else if (d_interpolation_function_profile == "LINEAR_1")
                {
                    // lf and 1-lf
                    if (lf <= 0.0358589525337265)
                        (*q_firstder_data)(ci) = lf;
                    else if (lf > 0.0358589525337265 && lf <= 0.964141047466262)
                        (*q_firstder_data)(ci) = 30.0 * std::pow(lf, 4.0) - 60.0 * std::pow(lf, 3.0) + 30.0 * lf * lf;
                    else
                        (*q_firstder_data)(ci) = 1.0 - lf;
                }
                else if (d_interpolation_function_profile == "LINEAR_0.5")
                {
                    // 0.5*lf and 0.5-0.5*lf
                    if (lf <= 0.017257145471902)
                        (*q_firstder_data)(ci) = 0.5 * lf;
                    else if (lf > 0.017257145471902 && lf <= 0.982742854528096)
                        (*q_firstder_data)(ci) = 30.0 * std::pow(lf, 4.0) - 60.0 * std::pow(lf, 3.0) + 30.0 * lf * lf;
                    else
                        (*q_firstder_data)(ci) = 0.5 - 0.5 * lf;
                }
                else
                    TBOX_ERROR(
                        "Interpolation function valid options are: LINEAR_0.5, LINEAR_1, LINEAR_2, LINEAR_3, LINEAR_4, "
                        "QUADRATIC, QUARTIC");

                if (lf >= 1.0 - 1e-10 && T <= d_T_melt)
                {
                    (*q_firstder_data)(ci) = 1.0;
                }

                if (lf <= 1e-10 && T >= d_T_melt)
                {
                    (*q_firstder_data)(ci) = 1.0;
                }
            }
        }
    }
    return;
} // computeInterpolationFunction

void
AllenCahnHierarchyIntegrator::computeLiquidFractionSourceTerm(int F_scratch_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
    int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > lf_new_data = patch->getPatchData(lf_new_idx);
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_new_idx);
            Pointer<CellData<NDIM, double> > q_firstder_data = patch->getPatchData(d_q_firstder_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_scratch_idx);
            Pointer<CellData<NDIM, double> > g_firstder_data = patch->getPatchData(d_g_firstder_idx);
            Pointer<CellData<NDIM, double> > g_secondder_data = patch->getPatchData(d_g_secondder_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                double F = -d_M_lf * d_rho_liquid * d_latent_heat * (*H_data)(ci) * (*q_firstder_data)(ci) *
                           (d_T_melt - (*T_data)(ci)) / d_T_melt;
                (*F_data)(ci) = F - (d_M_lf * d_lambda_lf * (*H_data)(ci) / std::pow(d_eps_lf, 2.0) *
                                     ((*g_firstder_data)(ci) - ((*g_secondder_data)(ci) * (*lf_new_data)(ci))));
            }
        }
    }
    return;
} // computeLiquidFractionSourceTerm

void
AllenCahnHierarchyIntegrator::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db);
#endif
    // Read in data members from input database.
    if (!is_from_restart)
    {
        d_M_lf = input_db->getDouble("M_lf");
        d_lambda_lf = input_db->getDouble("lambda_lf");
        d_eps_lf = input_db->getDouble("eps_lf");

        if (input_db->keyExists("solve_energy")) d_solve_energy = input_db->getBool("solve_energy");
        if (input_db->keyExists("lf_diffusion_time_stepping_type"))
            d_lf_diffusion_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("lf_diffusion_time_stepping_type"));
        if (input_db->keyExists("lf_convective_difference_form"))
            d_lf_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("lf_convective_difference_form"));
        else if (input_db->keyExists("lf_convective_difference_type"))
            d_lf_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("lf_convective_difference_type"));
        else if (input_db->keyExists("default_lf_convective_difference_form"))
            d_lf_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(
                input_db->getString("default_lf_convective_difference_form"));
        else if (input_db->keyExists("default_lf_convective_difference_type"))
            d_lf_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(
                input_db->getString("default_lf_convective_difference_type"));

        if (input_db->keyExists("lf_convective_op_type"))
            d_lf_convective_op_type = input_db->getString("lf_convective_op_type");
        else if (input_db->keyExists("lf_convective_operator_type"))
            d_lf_convective_op_type = input_db->getString("lf_convective_operator_type");
        else if (input_db->keyExists("default_lf_convective_op_type"))
            d_lf_convective_op_type = input_db->getString("default_lf_convective_op_type");
        else if (input_db->keyExists("default_lf_convective_operator_type"))
            d_lf_convective_op_type = input_db->getString("default_lf_convective_operator_type");

        if (input_db->keyExists("lf_convective_op_db"))
            d_lf_convective_op_input_db = input_db->getDatabase("lf_convective_op_db");
        else if (input_db->keyExists("default_lf_convective_op_db"))
            d_lf_convective_op_input_db = input_db->getDatabase("default_lf_convective_op_db");

        if (input_db->keyExists("lf_convective_time_stepping_type"))
            d_lf_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("lf_convective_time_stepping_type"));
        if (input_db->keyExists("T_convective_time_stepping_type"))
            d_T_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("T_convective_time_stepping_type"));

        if (input_db->keyExists("num_diffusion")) d_num_diffusion = input_db->getDouble("num_diffusion");
        if (input_db->keyExists("interpolation_function_profile"))
            d_interpolation_function_profile = input_db->getString("interpolation_function_profile");
    }
}

void
AllenCahnHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("IEP_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IEP_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }

    d_lf_diffusion_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_lf_diffusion_time_stepping_type"));
    d_lf_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_lf_convective_time_stepping_type"));
    d_lf_init_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_lf_init_convective_time_stepping_type"));
    d_lf_convective_difference_form =
        string_to_enum<ConvectiveDifferencingType>(db->getString("d_lf_convective_difference_form"));
    d_lf_convective_op_type = db->getString("d_lf_convective_op_type");
    d_lf_convective_op_input_db = db->getDatabase("d_lf_convective_op_db");

    d_M_lf = db->getDouble("d_M_lf");
    d_lambda_lf = db->getDouble("d_lambda_lf");
    d_eps_lf = db->getDouble("d_eps_lf");
    d_solve_energy = db->getBool("d_solve_energy");

    d_num_diffusion = db->getDouble("d_num_diffusion");
    d_interpolation_function_profile = db->getString("d_interpolation_function_profile");
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
