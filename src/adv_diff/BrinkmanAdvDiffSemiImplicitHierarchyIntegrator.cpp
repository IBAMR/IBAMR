// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
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

#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/BrinkmanAdvDiffBcHelper.h"
#include "ibamr/BrinkmanAdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

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

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

BrinkmanAdvDiffSemiImplicitHierarchyIntegrator::BrinkmanAdvDiffSemiImplicitHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : AdvDiffSemiImplicitHierarchyIntegrator(object_name, input_db, register_for_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    if (input_db->keyExists("brinkman_time_independent"))
        d_brinkman_time_independent = input_db->getBool("brinkman_time_independent");

    return;
} // BrinkmanAdvDiffSemiImplicitHierarchyIntegrator

void
BrinkmanAdvDiffSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_fc_data_ops =
        hier_ops_manager->getOperationsDouble(new FaceVariable<NDIM, double>("fc_var"), hierarchy, true);

    // Register additional variables required for present time stepping algorithm.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> no_ghosts = 0;
    for (const auto& Q_var : d_Q_var)
    {
        Pointer<CellDataFactory<NDIM, double> > Q_factory = Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        // Variables required for Brinkman penalization
        const bool apply_brinkman =
            d_brinkman_penalization && d_brinkman_penalization->hasBrinkmanBoundaryCondition(Q_var);
        if (apply_brinkman)
        {
            Pointer<CellVariable<NDIM, double> > Cb_var =
                new CellVariable<NDIM, double>(Q_var->getName() + "::Cb", Q_depth);
            Pointer<CellVariable<NDIM, double> > Cb_rhs_var =
                new CellVariable<NDIM, double>(Q_var->getName() + "::Cb_RHS", Q_depth);
            Pointer<SideVariable<NDIM, double> > Db_var =
                new SideVariable<NDIM, double>(Q_var->getName() + "::Db", Q_depth);
            Pointer<SideVariable<NDIM, double> > Db_rhs_var =
                new SideVariable<NDIM, double>(Q_var->getName() + "::Db_RHS", Q_depth);
            Pointer<CellVariable<NDIM, double> > Fb_var =
                new CellVariable<NDIM, double>(Q_var->getName() + "::Fb", Q_depth);
            d_Q_Cb_map[Q_var] = Cb_var;
            d_Q_Cb_rhs_map[Q_var] = Cb_rhs_var;
            d_Q_Db_map[Q_var] = Db_var;
            d_Q_Db_rhs_map[Q_var] = Db_rhs_var;
            d_Q_Fb_map[Q_var] = Fb_var;
            int Cb_current_idx, Cb_scratch_idx, Cb_rhs_scratch_idx, Db_current_idx, Db_scratch_idx, Db_rhs_scratch_idx,
                Fb_scratch_idx;
            registerVariable(Cb_current_idx, Cb_var, no_ghosts, getCurrentContext());
            registerVariable(Cb_scratch_idx, Cb_var, no_ghosts, getScratchContext());
            registerVariable(Cb_rhs_scratch_idx, Cb_rhs_var, no_ghosts, getScratchContext());
            registerVariable(Db_current_idx, Db_var, no_ghosts, getCurrentContext());
            registerVariable(Db_scratch_idx, Db_var, side_ghosts, getScratchContext());
            registerVariable(Db_rhs_scratch_idx, Db_rhs_var, side_ghosts, getScratchContext());
            registerVariable(Fb_scratch_idx, Fb_var, no_ghosts, getScratchContext());
        }
    }

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffSemiImplicitHierarchyIntegrator.
    AdvDiffSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
BrinkmanAdvDiffSemiImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                             const double new_time,
                                                                             const int num_cycles)
{
    AdvDiffHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    d_brinkman_penalization->setTimeInterval(current_time, new_time);
    d_brinkman_penalization->preprocessBrinkmanAdvDiffBcHelper(current_time, new_time, num_cycles);

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
        std::fill(d_helmholtz_solvers_need_init.begin(), d_helmholtz_solvers_need_init.end(), true);
        std::fill(d_helmholtz_rhs_ops_need_init.begin(), d_helmholtz_rhs_ops_need_init.end(), true);
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Update the advection velocity.
    for (const auto& u_var : d_u_var)
    {
        const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
        const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, getScratchContext());
        const int u_new_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
        if (d_u_fcn[u_var])
        {
            d_u_fcn[u_var]->setDataOnPatchHierarchy(u_current_idx, u_var, d_hierarchy, current_time);
            d_u_fcn[u_var]->setDataOnPatchHierarchy(u_new_idx, u_var, d_hierarchy, new_time);
        }
        else
        {
            d_hier_fc_data_ops->copyData(u_new_idx, u_current_idx);
        }
        d_hier_fc_data_ops->linearSum(u_scratch_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    }

    // Update the diffusion coefficient
    for (const auto& D_var : d_diffusion_coef_var)
    {
        Pointer<CartGridFunction> D_fcn = d_diffusion_coef_fcn[D_var];
        if (D_fcn)
        {
            const int D_current_idx = var_db->mapVariableAndContextToIndex(D_var, getCurrentContext());
            D_fcn->setDataOnPatchHierarchy(D_current_idx, D_var, d_hierarchy, current_time);
        }
    }

    // Setup the operators and solvers and compute the right-hand-side terms.
    unsigned int l = 0;
    for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM, double> > Q_var = *cit;
        Pointer<CellVariable<NDIM, double> > Q_rhs_var = d_Q_Q_rhs_map[Q_var];
        Pointer<SideVariable<NDIM, double> > D_var = d_Q_diffusion_coef_variable[Q_var];
        Pointer<SideVariable<NDIM, double> > D_rhs_var = d_diffusion_coef_rhs_map[D_var];
        TimeSteppingType diffusion_time_stepping_type = d_Q_diffusion_time_stepping_type[Q_var];
        const double lambda = d_Q_damping_coef[Q_var];
        const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coef = d_Q_bc_coef[Q_var];

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int Q_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getScratchContext());
        const int D_current_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getCurrentContext()) : -1);
        const int D_scratch_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getScratchContext()) : -1);
        const int D_rhs_scratch_idx =
            (D_rhs_var ? var_db->mapVariableAndContextToIndex(D_rhs_var, getScratchContext()) : -1);

        // Note that when Brinkman penalization is applied, linear operators are reset on a cycle-by-cycle basis
        // in integrateHierarchy().
        const bool apply_brinkman =
            d_brinkman_penalization && d_brinkman_penalization->hasBrinkmanBoundaryCondition(Q_var);
        if (apply_brinkman)
        {
            if (lambda != 0.0)
            {
                TBOX_ERROR(d_object_name
                           << "::preprocessIntegrateHierarchy():\n"
                           << "  physical damping coefficient lambda must be 0.0 to apply Brinkman penalization BCs\n");
            }
            if (d_enable_logging)
            {
                plog << d_object_name << ": "
                     << "Applying Brinkman penalization for variable number " << l << "\n";
            }
        }
        else
        {
            // Setup the problem coefficients for the linear solve for Q(n+1).
            double K = 0.0;
            switch (diffusion_time_stepping_type)
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
                TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                         << "  unsupported diffusion time stepping type: "
                                         << enum_to_string<TimeSteppingType>(diffusion_time_stepping_type) << " \n"
                                         << "  valid choices are: BACKWARD_EULER, "
                                            "FORWARD_EULER, TRAPEZOIDAL_RULE\n");
            }
            PoissonSpecifications solver_spec(d_object_name + "::solver_spec::" + Q_var->getName());
            PoissonSpecifications rhs_op_spec(d_object_name + "::rhs_op_spec::" + Q_var->getName());
            solver_spec.setCConstant(1.0 / dt + K * lambda);
            rhs_op_spec.setCConstant(1.0 / dt - (1.0 - K) * lambda);

            if (isDiffusionCoefficientVariable(Q_var))
            {
                // set -K*kappa in solver_spec
                d_hier_sc_data_ops->scale(D_scratch_idx, -K, D_current_idx);
                solver_spec.setDPatchDataId(D_scratch_idx);
                // set (1.0-K)*kappa in rhs_op_spec
                d_hier_sc_data_ops->scale(D_rhs_scratch_idx, (1.0 - K), D_current_idx);
                rhs_op_spec.setDPatchDataId(D_rhs_scratch_idx);
            }
            else
            {
                const double kappa = d_Q_diffusion_coef[Q_var];
                solver_spec.setDConstant(-K * kappa);
                rhs_op_spec.setDConstant(+(1.0 - K) * kappa);
            }

            // Initialize the RHS operator and compute the RHS vector.
            Pointer<LaplaceOperator> helmholtz_rhs_op = d_helmholtz_rhs_ops[l];
            helmholtz_rhs_op->setPoissonSpecifications(rhs_op_spec);
            helmholtz_rhs_op->setPhysicalBcCoefs(Q_bc_coef);
            helmholtz_rhs_op->setHomogeneousBc(false);
            helmholtz_rhs_op->setSolutionTime(current_time);
            helmholtz_rhs_op->setTimeInterval(current_time, new_time);
            if (d_helmholtz_rhs_ops_need_init[l])
            {
                if (d_enable_logging)
                {
                    plog << d_object_name << ": "
                         << "Initializing Helmholtz RHS operator for variable number " << l << "\n";
                }
                helmholtz_rhs_op->initializeOperatorState(*d_sol_vecs[l], *d_rhs_vecs[l]);
                d_helmholtz_rhs_ops_need_init[l] = false;
            }
            d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx, false);
            helmholtz_rhs_op->apply(*d_sol_vecs[l], *d_rhs_vecs[l]);

            // Initialize the linear solver.
            Pointer<PoissonSolver> helmholtz_solver = d_helmholtz_solvers[l];
            helmholtz_solver->setPoissonSpecifications(solver_spec);
            helmholtz_solver->setPhysicalBcCoefs(Q_bc_coef);
            helmholtz_solver->setHomogeneousBc(false);
            helmholtz_solver->setSolutionTime(new_time);
            helmholtz_solver->setTimeInterval(current_time, new_time);
            if (d_helmholtz_solvers_need_init[l])
            {
                if (d_enable_logging)
                {
                    plog << d_object_name << ": "
                         << "Initializing Helmholtz solvers for variable number " << l << "\n";
                }
                helmholtz_solver->initializeSolverState(*d_sol_vecs[l], *d_rhs_vecs[l]);
                d_helmholtz_solvers_need_init[l] = false;
            }
        }

        // Account for the convective difference term.
        Pointer<FaceVariable<NDIM, double> > u_var = d_Q_u_map[Q_var];
        if (u_var)
        {
            Pointer<CellVariable<NDIM, double> > N_var = d_Q_N_map[Q_var];
            Pointer<CellVariable<NDIM, double> > N_old_var = d_Q_N_old_map[Q_var];
            TimeSteppingType convective_time_stepping_type = d_Q_convective_time_stepping_type[Q_var];
            if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(convective_time_stepping_type))
            {
                convective_time_stepping_type = d_Q_init_convective_time_stepping_type[Q_var];
            }
            if ((num_cycles == 1) &&
                (convective_time_stepping_type == MIDPOINT_RULE || convective_time_stepping_type == TRAPEZOIDAL_RULE))
            {
                TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                         << "  time stepping type: "
                                         << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                                         << " requires num_cycles > 1.\n"
                                         << "  at current time step, num_cycles = " << num_cycles << "\n");
            }
            if (d_Q_convective_op_needs_init[Q_var])
            {
                d_Q_convective_op[Q_var]->initializeOperatorState(*d_sol_vecs[l], *d_rhs_vecs[l]);
                d_Q_convective_op_needs_init[Q_var] = false;
            }
            const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
            d_Q_convective_op[Q_var]->setAdvectionVelocity(u_current_idx);
            const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
            const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx);
            d_Q_convective_op[Q_var]->setSolutionTime(current_time);
            d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
            const int N_old_new_idx = var_db->mapVariableAndContextToIndex(N_old_var, getNewContext());
            d_hier_cc_data_ops->copyData(N_old_new_idx, N_scratch_idx);

            if (!apply_brinkman)
            {
                if (convective_time_stepping_type == FORWARD_EULER)
                {
                    d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, N_scratch_idx, Q_rhs_scratch_idx);
                }
                else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
                {
                    d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -0.5, N_scratch_idx, Q_rhs_scratch_idx);
                }
            }
        }

        // Set the initial guess.
        d_hier_cc_data_ops->copyData(Q_new_idx, Q_current_idx);
    }

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
BrinkmanAdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                                   const double new_time,
                                                                   const int cycle_num)
{
    AdvDiffHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE(
            {
                pout << "BrinkmanAdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy():\n"
                     << "  WARNING: num_cycles = " << d_current_num_cycles
                     << " but expected num_cycles = " << expected_num_cycles << ".\n";
            });
    }

    // Perform a single step of fixed point iteration.
    unsigned int l = 0;
    for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM, double> > Q_var = *cit;
        Pointer<CellVariable<NDIM, double> > F_var = d_Q_F_map[Q_var];
        Pointer<CellVariable<NDIM, double> > Q_rhs_var = d_Q_Q_rhs_map[Q_var];

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int F_scratch_idx =
            d_F_fcn[F_var] ? var_db->mapVariableAndContextToIndex(F_var, getScratchContext()) : -1;
        const int F_new_idx = d_F_fcn[F_var] ? var_db->mapVariableAndContextToIndex(F_var, getNewContext()) : -1;
        const int Q_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getScratchContext());

        // Get Brinkman penalization variables and data indices if necessary
        int Fb_scratch_idx = -1;
        const bool apply_brinkman =
            d_brinkman_penalization && d_brinkman_penalization->hasBrinkmanBoundaryCondition(Q_var);
        if (apply_brinkman)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(d_Q_damping_coef[Q_var] == 0.0);
#endif
            TimeSteppingType diffusion_time_stepping_type = d_Q_diffusion_time_stepping_type[Q_var];
            const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coef = d_Q_bc_coef[Q_var];
            Pointer<SideVariable<NDIM, double> > D_var = d_Q_diffusion_coef_variable[Q_var];
            Pointer<CellVariable<NDIM, double> > Cb_var = d_Q_Cb_map[Q_var];
            Pointer<CellVariable<NDIM, double> > Cb_rhs_var = d_Q_Cb_rhs_map[Q_var];
            Pointer<SideVariable<NDIM, double> > Db_var = d_Q_Db_map[Q_var];
            Pointer<SideVariable<NDIM, double> > Db_rhs_var = d_Q_Db_rhs_map[Q_var];
            Pointer<CellVariable<NDIM, double> > Fb_var = d_Q_Fb_map[Q_var];
            int Cb_current_idx = -1, Cb_scratch_idx = -1, Cb_rhs_scratch_idx = -1;
            int Db_current_idx = -1, Db_scratch_idx = -1, Db_rhs_scratch_idx = -1;

            const int D_current_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getCurrentContext()) : -1);
            Cb_current_idx = var_db->mapVariableAndContextToIndex(Cb_var, getCurrentContext());
            Cb_scratch_idx = var_db->mapVariableAndContextToIndex(Cb_var, getScratchContext());
            Cb_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Cb_rhs_var, getScratchContext());
            Db_current_idx = var_db->mapVariableAndContextToIndex(Db_var, getCurrentContext());
            Db_scratch_idx = var_db->mapVariableAndContextToIndex(Db_var, getScratchContext());
            Db_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Db_rhs_var, getScratchContext());
            Fb_scratch_idx = var_db->mapVariableAndContextToIndex(Fb_var, getScratchContext());

            // Compute the Brinkman penalization contributions to the linear operators and RHS for Q_var
            d_brinkman_penalization->computeDampingCoefficient(Cb_current_idx, Q_var);
            d_brinkman_penalization->computeDiffusionCoefficient(
                Db_current_idx, Q_var, D_current_idx, d_Q_diffusion_coef[Q_var]);

            // Setup the problem coefficients for the linear solve
            double K = 0.0;
            switch (diffusion_time_stepping_type)
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
                TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                         << "  unsupported diffusion time stepping type: "
                                         << enum_to_string<TimeSteppingType>(diffusion_time_stepping_type) << " \n"
                                         << "  valid choices are: BACKWARD_EULER, "
                                            "FORWARD_EULER, TRAPEZOIDAL_RULE\n");
            }
            PoissonSpecifications solver_spec(d_object_name + "::solver_spec::" + Q_var->getName());
            PoissonSpecifications rhs_op_spec(d_object_name + "::rhs_op_spec::" + Q_var->getName());

            // set 1.0/dt + K*Cb in solver_spec (lambda is included within Cb).
            d_hier_cc_data_ops->scale(Cb_scratch_idx, K, Cb_current_idx);

            // In the case of Poisson applications zero out the temporal contribution.
            if (d_brinkman_time_independent)
            {
                solver_spec.setCPatchDataId(Cb_scratch_idx);
                // solver_spec.setCConstant(0.0);
            }
            else
            {
                d_hier_cc_data_ops->addScalar(Cb_scratch_idx, Cb_scratch_idx, 1.0 / dt);
                solver_spec.setCPatchDataId(Cb_scratch_idx);
            }

            // set 1.0/dt - (1.0-K)*Cb in rhs_op_spec, which is applied to Q later
            d_hier_cc_data_ops->scale(Cb_rhs_scratch_idx, -(1.0 - K), Cb_current_idx);
            // In the case of Poisson applications zero out the temporal contribution.
            if (d_brinkman_time_independent)
            {
                // rhs_op_spec.setCConstant(0.0);
                rhs_op_spec.setCPatchDataId(Cb_rhs_scratch_idx);
            }
            else
            {
                d_hier_cc_data_ops->addScalar(Cb_rhs_scratch_idx, Cb_rhs_scratch_idx, 1.0 / dt);
                rhs_op_spec.setCPatchDataId(Cb_rhs_scratch_idx);
            }

            // set -K*Db in solver_spec (kappa is included within Db)
            d_hier_sc_data_ops->scale(Db_scratch_idx, -K, Db_current_idx);
            solver_spec.setDPatchDataId(Db_scratch_idx);
            // set (1.0-K)*Db in rhs_op_spec
            d_hier_sc_data_ops->scale(Db_rhs_scratch_idx, (1.0 - K), Db_current_idx);
            rhs_op_spec.setDPatchDataId(Db_rhs_scratch_idx);

            // Initialize the RHS operator and compute the RHS vector.
            Pointer<LaplaceOperator> helmholtz_rhs_op = d_helmholtz_rhs_ops[l];
            helmholtz_rhs_op->setPoissonSpecifications(rhs_op_spec);
            helmholtz_rhs_op->setPhysicalBcCoefs(Q_bc_coef);
            helmholtz_rhs_op->setHomogeneousBc(false);
            helmholtz_rhs_op->setSolutionTime(current_time);
            helmholtz_rhs_op->setTimeInterval(current_time, new_time);

            if (d_helmholtz_rhs_ops_need_init[l]) // TODO: Does this need to happen every cycle?
            {
                if (d_enable_logging)
                {
                    plog << d_object_name << ": "
                         << "Initializing Helmholtz RHS operator for variable number " << l << "\n";
                }
                helmholtz_rhs_op->initializeOperatorState(*d_sol_vecs[l], *d_rhs_vecs[l]);
                d_helmholtz_rhs_ops_need_init[l] = false;
            }
            d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx, false);
            helmholtz_rhs_op->apply(*d_sol_vecs[l], *d_rhs_vecs[l]);

            // Initialize the linear solver.
            Pointer<PoissonSolver> helmholtz_solver = d_helmholtz_solvers[l];
            helmholtz_solver->setPoissonSpecifications(solver_spec);
            helmholtz_solver->setPhysicalBcCoefs(Q_bc_coef);
            helmholtz_solver->setHomogeneousBc(false);
            helmholtz_solver->setSolutionTime(new_time);
            helmholtz_solver->setTimeInterval(current_time, new_time);

            if (d_helmholtz_solvers_need_init[l]) // TODO: Does this need to happen every cycle?
            {
                if (d_enable_logging)
                {
                    plog << d_object_name << ": "
                         << "Initializing Helmholtz solvers for variable number " << l << "\n";
                }
                helmholtz_solver->initializeSolverState(*d_sol_vecs[l], *d_rhs_vecs[l]);
                d_helmholtz_solvers_need_init[l] = false;
            }
        }

        // Update the advection velocity.
        if (cycle_num > 0)
        {
            for (const auto& u_var : d_u_var)
            {
                const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
                const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, getScratchContext());
                const int u_new_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
                if (d_u_fcn[u_var])
                {
                    d_u_fcn[u_var]->setDataOnPatchHierarchy(u_new_idx, u_var, d_hierarchy, new_time);
                }
                d_hier_fc_data_ops->linearSum(u_scratch_idx, 0.5, u_current_idx, 0.5, u_new_idx);
            }
        }

        // Account for the convective difference term.
        Pointer<FaceVariable<NDIM, double> > u_var = d_Q_u_map[Q_var];
        Pointer<CellVariable<NDIM, double> > N_var = d_Q_N_map[Q_var];
        Pointer<CellVariable<NDIM, double> > N_old_var = d_Q_N_old_map[Q_var];
        TimeSteppingType convective_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
        if (u_var)
        {
            convective_time_stepping_type = d_Q_convective_time_stepping_type[Q_var];
            if (is_multistep_time_stepping_type(convective_time_stepping_type))
            {
#if !defined(NDEBUG)
                TBOX_ASSERT(convective_time_stepping_type == ADAMS_BASHFORTH);
#endif
                if (getIntegratorStep() == 0)
                {
                    convective_time_stepping_type = d_Q_init_convective_time_stepping_type[Q_var];
                }
                else if (cycle_num > 0)
                {
                    convective_time_stepping_type = MIDPOINT_RULE;
                    IBAMR_DO_ONCE(
                        {
                            pout << "BrinkmanAdvDiffSemiImplicitHierarchyIntegrator::"
                                    "integrateHierarchy():"
                                    "\n"
                                 << "  WARNING: convective_time_stepping_type = "
                                 << enum_to_string<TimeSteppingType>(d_Q_convective_time_stepping_type[Q_var])
                                 << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                                 << "           using "
                                 << enum_to_string<TimeSteppingType>(d_Q_convective_time_stepping_type[Q_var])
                                 << " only for the first cycle in each time step;\n"
                                 << "           using "
                                 << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                                 << " for subsequent cycles.\n";
                        });
                }
            }
            // Account for masked convective operator for Brinkman penalization
            if (apply_brinkman)
            {
                const int N_old_new_idx = var_db->mapVariableAndContextToIndex(N_old_var, getNewContext());
                const int N_old_scratch_idx = var_db->mapVariableAndContextToIndex(N_old_var, getScratchContext());
                d_hier_cc_data_ops->copyData(N_old_scratch_idx, N_old_new_idx);
                d_brinkman_penalization->maskForcingTerm(N_old_scratch_idx, Q_var, true /*mask_smeared_region*/);

                if (convective_time_stepping_type == FORWARD_EULER)
                {
                    d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, N_old_scratch_idx, Q_rhs_scratch_idx);
                }
                else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
                {
                    d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -0.5, N_old_scratch_idx, Q_rhs_scratch_idx);
                }
            }

            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            if (cycle_num > 0)
            {
                if (convective_time_stepping_type == MIDPOINT_RULE)
                {
                    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, getScratchContext());
                    d_Q_convective_op[Q_var]->setAdvectionVelocity(u_scratch_idx);
                    const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
                    const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
                    const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
                    d_hier_cc_data_ops->linearSum(Q_scratch_idx, 0.5, Q_current_idx, 0.5, Q_new_idx);
                    d_Q_convective_op[Q_var]->setSolutionTime(half_time);
                    d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
                }
                else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
                {
                    const int u_new_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
                    d_Q_convective_op[Q_var]->setAdvectionVelocity(u_new_idx);
                    const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
                    const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
                    d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_new_idx);
                    d_Q_convective_op[Q_var]->setSolutionTime(new_time);
                    d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
                }
            }
            if (convective_time_stepping_type == ADAMS_BASHFORTH)
            {
#if !defined(NDEBUG)
                TBOX_ASSERT(cycle_num == 0);
#endif
                const int N_old_current_idx = var_db->mapVariableAndContextToIndex(N_old_var, getCurrentContext());
                const double omega = dt / d_dt_previous[0];
                d_hier_cc_data_ops->linearSum(
                    N_scratch_idx, 1.0 + 0.5 * omega, N_scratch_idx, -0.5 * omega, N_old_current_idx);
            }

            // Mask the convective operator due to Brinkman penalization
            if (apply_brinkman)
                d_brinkman_penalization->maskForcingTerm(N_scratch_idx, Q_var, true /*mask_smeared_region*/);

            if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -0.5, N_scratch_idx, Q_rhs_scratch_idx);
            }
        }

        // Account for forcing terms.
        if (d_F_fcn[F_var])
        {
            d_F_fcn[F_var]->setDataOnPatchHierarchy(F_scratch_idx, F_var, d_hierarchy, half_time);

            // Mask the body forcing due to Brinkman penalization
            if (apply_brinkman) d_brinkman_penalization->maskForcingTerm(F_scratch_idx, Q_var);

            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, 1.0, F_scratch_idx, Q_rhs_scratch_idx);
        }

        // Account for optional Brinkman RHS forcing terms
        if (apply_brinkman)
        {
            d_brinkman_penalization->computeForcing(Fb_scratch_idx, Q_var);
            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, 1.0, Fb_scratch_idx, Q_rhs_scratch_idx);
        }

        if (isDiffusionCoefficientVariable(Q_var) || (d_Q_diffusion_coef[Q_var] != 0.0) || apply_brinkman)
        {
            // Solve for Q(n+1).
            Pointer<PoissonSolver> helmholtz_solver = d_helmholtz_solvers[l];
            helmholtz_solver->solveSystem(*d_sol_vecs[l], *d_rhs_vecs[l]);
            d_hier_cc_data_ops->copyData(Q_new_idx, Q_scratch_idx);
            if (d_enable_logging && d_enable_logging_solver_iterations)
                plog << d_object_name << "::integrateHierarchy(): diffusion solve number of iterations = "
                     << helmholtz_solver->getNumIterations() << "\n";
            if (d_enable_logging)
                plog << d_object_name << "::integrateHierarchy(): diffusion solve residual norm        = "
                     << helmholtz_solver->getResidualNorm() << "\n";
            if (helmholtz_solver->getNumIterations() == helmholtz_solver->getMaxIterations())
            {
                pout << d_object_name << "::integrateHierarchy():"
                     << "  WARNING: linear solver iterations == max iterations\n";
            }
        }
        else
        {
            // No solve needed for Q(n+1)
            d_hier_cc_data_ops->scale(Q_new_idx, dt, Q_rhs_scratch_idx);
            if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): completed solution update.\n";
        }

        // Reset the right-hand side vector.
        if (u_var)
        {
            if (apply_brinkman)
            {
                const int N_old_scratch_idx = var_db->mapVariableAndContextToIndex(N_old_var, getScratchContext());
                if (convective_time_stepping_type == FORWARD_EULER)
                {
                    d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, +1.0, N_old_scratch_idx, Q_rhs_scratch_idx);
                }
                else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
                {
                    d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, +0.5, N_old_scratch_idx, Q_rhs_scratch_idx);
                }
            }
            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, +1.0, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, +0.5, N_scratch_idx, Q_rhs_scratch_idx);
            }
        }
        if (d_F_fcn[F_var])
        {
            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, F_scratch_idx, Q_rhs_scratch_idx);
            d_hier_cc_data_ops->copyData(F_new_idx, F_scratch_idx);
        }

        if (apply_brinkman)
        {
            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, Fb_scratch_idx, Q_rhs_scratch_idx);
        }
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
BrinkmanAdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const bool skip_synchronize_new_state_data,
    const int num_cycles)
{
    AdvDiffHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    d_brinkman_penalization->postprocessBrinkmanAdvDiffBcHelper(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Determine the CFL number.
    double cfl_max = 0.0;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    for (const auto& u_var : d_u_var)
    {
        const int u_new_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
        PatchFaceDataOpsReal<NDIM, double> patch_fc_ops;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double dx_min = *(std::min_element(dx, dx + NDIM));
                Pointer<FaceData<NDIM, double> > u_fc_new_data = patch->getPatchData(u_new_idx);
                double u_max = 0.0;
                u_max = patch_fc_ops.maxNorm(u_fc_new_data, patch_box);
                cfl_max = std::max(cfl_max, u_max * dt / dx_min);
            }
        }
    }
    cfl_max = IBTK_MPI::maxReduction(cfl_max);
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
BrinkmanAdvDiffSemiImplicitHierarchyIntegrator::registerBrinkmanAdvDiffBcHelper(
    Pointer<BrinkmanAdvDiffBcHelper> brinkman_penalization)
{
    d_brinkman_penalization = brinkman_penalization;
    return;
} // registerBrinkmanAdvDiffBcHelper

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
