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

#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/TemperatureSemiImplicitHierarchyIntegrator.h"
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
namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
static const int NOGHOSTS = 0;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

TemperatureSemiImplicitHierarchyIntegrator::TemperatureSemiImplicitHierarchyIntegrator(const std::string& object_name,
                                                                                       Pointer<Database> input_db,
                                                                                       bool register_for_restart)
    : AdvDiffSemiImplicitHierarchyIntegrator(object_name, input_db, register_for_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    return;
} // TemperatureSemiImplicitHierarchyIntegrator

void
TemperatureSemiImplicitHierarchyIntegrator::registerSpecificHeatVariable(Pointer<CellVariable<NDIM, double> > Cp_var,
                                                                         const bool output_Cp)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Cp_var);
    TBOX_ASSERT(std::find(d_Cp_var.begin(), d_Cp_var.end(), Cp_var) == d_Cp_var.end());
#endif
    d_Cp_var.push_back(Cp_var);
    d_Cp_output[Cp_var] = output_Cp;

    return;
} // registerSpecificHeatVariable

void
TemperatureSemiImplicitHierarchyIntegrator::registerDensityVariable(Pointer<CellVariable<NDIM, double> > rho_var,
                                                                    const bool output_rho)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(rho_var);
    TBOX_ASSERT(std::find(d_rho_var.begin(), d_rho_var.end(), rho_var) == d_rho_var.end());
#endif
    d_rho_var.push_back(rho_var);
    d_rho_output[rho_var] = output_rho;

    return;
} // registerDensityVariable

void
TemperatureSemiImplicitHierarchyIntegrator::registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_rho_fcns.push_back(callback);
    d_reset_rho_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidDensityFcn

void
TemperatureSemiImplicitHierarchyIntegrator::registerResetSpecificHeatFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_Cp_fcns.push_back(callback);
    d_reset_Cp_fcns_ctx.push_back(ctx);
    return;
} // registerResetSpecificHeatFcn

void
TemperatureSemiImplicitHierarchyIntegrator::registerResetDiffusionCoefficientFcn(ResetFluidPropertiesFcnPtr callback,
                                                                                 void* ctx)
{
    d_reset_kappa_fcns.push_back(callback);
    d_reset_kappa_fcns_ctx.push_back(ctx);
    return;
} // registerResetDiffusionCoefficientFcn

void
TemperatureSemiImplicitHierarchyIntegrator::setDensityVariable(Pointer<CellVariable<NDIM, double> > Q_var,
                                                               Pointer<CellVariable<NDIM, double> > rho_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    TBOX_ASSERT(std::find(d_rho_var.begin(), d_rho_var.end(), rho_var) != d_rho_var.end());
#endif
    d_Q_rho_map[Q_var] = rho_var;
    return;
} // setDensityVariable

void
TemperatureSemiImplicitHierarchyIntegrator::registerINSVCStaggeredHierarchyIntegrator(
    Pointer<INSVCStaggeredHierarchyIntegrator> ins_cons_hier_integrator)
{
    d_ins_hierarchy_integrator = ins_cons_hier_integrator;
    return;
} // registerINSVCStaggeredHierarchyIntegrator

void
TemperatureSemiImplicitHierarchyIntegrator::interpolateSCMassDensityToCC(Pointer<CellVariable<NDIM, double> > Q_var,
                                                                         Pointer<VariableContext> ctx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > rho_var = d_ins_hierarchy_integrator->getMassDensityVariable();
    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    Pointer<CellVariable<NDIM, double> > rho_cc_var = rho_var;

    if (rho_sc_var)
    {
        int rho_sc_idx = 0;
        if (ctx == getNewContext())
        {
            rho_sc_idx = var_db->mapVariableAndContextToIndex(rho_sc_var, d_ins_hierarchy_integrator->getNewContext());
        }
        else if (ctx == getCurrentContext())
        {
            rho_sc_idx =
                var_db->mapVariableAndContextToIndex(rho_sc_var, d_ins_hierarchy_integrator->getCurrentContext());
        }

        Pointer<CellVariable<NDIM, double> > rho_cc_var = d_Q_rho_map[Q_var];
        const int rho_cc_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, ctx);
        const int rho_vec_cc_idx = var_db->mapVariableAndContextToIndex(d_rho_vec_cc_var, ctx);
        static const bool synch_cf_interface = true;

        d_hier_math_ops->interp(rho_vec_cc_idx,
                                d_rho_vec_cc_var,
                                rho_sc_idx,
                                rho_sc_var,
                                d_no_fill_op,
                                d_integrator_time,
                                synch_cf_interface);

        for (int ln = coarsest_ln; ln <= finest_ln; ln++)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            for (PatchLevelIterator<NDIM> p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > rho_cc_data = patch->getPatchData(rho_cc_idx);
                Pointer<CellData<NDIM, double> > rho_vec_data = patch->getPatchData(rho_vec_cc_idx);
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    double sum = 0.0;
                    for (int d = 0; d < NDIM; d++)
                    {
                        sum += (*rho_vec_data)(ci, d);
                    }
#if (NDIM == 2)
                    (*rho_cc_data)(ci) = 0.5 * sum;
#elif (NDIM == 3)
                    (*rho_cc_data)(ci) = 1.0 / 3.0 * sum;
#endif
                }
            }
        }
    }
    return;
} // interpolateSCMassDensityToCC

void
TemperatureSemiImplicitHierarchyIntegrator::setSpecificHeatVariable(Pointer<CellVariable<NDIM, double> > Q_var,
                                                                    Pointer<CellVariable<NDIM, double> > Cp_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    TBOX_ASSERT(std::find(d_Cp_var.begin(), d_Cp_var.end(), Cp_var) != d_Cp_var.end());
#endif
    d_Q_Cp_map[Q_var] = Cp_var;
    return;
} // setSpecificHeatVariable

void
TemperatureSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(
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
    const IntVector<NDIM> no_ghosts = NOGHOSTS;

    for (const auto& rho_var : d_rho_var)
    {
        int rho_current_idx, rho_new_idx, rho_scratch_idx;
        registerVariable(rho_current_idx,
                         rho_new_idx,
                         rho_scratch_idx,
                         rho_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
        const bool rho_data_output = d_rho_output[rho_var];
        if (d_visit_writer && rho_data_output)
            d_visit_writer->registerPlotQuantity(rho_var->getName(), "SCALAR", rho_current_idx);
    }

    for (const auto& Cp_var : d_Cp_var)
    {
        int Cp_current_idx, Cp_new_idx, Cp_scratch_idx;
        registerVariable(Cp_current_idx,
                         Cp_new_idx,
                         Cp_scratch_idx,
                         Cp_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
        const bool Cp_data_output = d_Cp_output[Cp_var];
        if (d_visit_writer && Cp_data_output)
            d_visit_writer->registerPlotQuantity(Cp_var->getName(), "SCALAR", Cp_current_idx);
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_C_var = new CellVariable<NDIM, double>("C_var");
    registerVariable(d_C_current_idx,
                     d_C_new_idx,
                     d_C_scratch_idx,
                     d_C_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_C_rhs_scratch_idx = var_db->registerVariableAndContext(d_C_var, var_db->getContext("C_rhs"));

    // Registering a temporary cell-centered vector variable to be used in the interpolation
    // function.
    d_rho_vec_cc_var = new CellVariable<NDIM, double>("rho_vec_cc", NDIM);
    registerVariable(d_rho_vec_cc_current_idx,
                     d_rho_vec_cc_new_idx,
                     d_rho_vec_cc_scratch_idx,
                     d_rho_vec_cc_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffSemiImplicitHierarchyIntegrator.
    AdvDiffSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
TemperatureSemiImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                         const double new_time,
                                                                         const int num_cycles)
{
    AdvDiffHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

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
        level->allocatePatchData(d_C_rhs_scratch_idx, current_time);
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

    // Update the diffusion coefficient.
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
        TimeSteppingType diffusion_time_stepping_type = d_Q_diffusion_time_stepping_type[Q_var];
        const double lambda = d_Q_damping_coef[Q_var];
        const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coef = d_Q_bc_coef[Q_var];

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int Q_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getScratchContext());

        const double apply_time = current_time;
        Pointer<CellVariable<NDIM, double> > rho_cc_var = d_Q_rho_map[Q_var];
        int rho_current_idx;
        if (rho_cc_var)
        {
            rho_current_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, getCurrentContext());
            for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
            {
                d_reset_rho_fcns[k](rho_current_idx,
                                    rho_cc_var,
                                    d_hier_math_ops,
                                    -1 /*cycle_num*/,
                                    apply_time,
                                    current_time,
                                    new_time,
                                    d_reset_rho_fcns_ctx[k]);
            }
        }
        //        if(rho_cc_var)
        //        {
        //            rho_current_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, getCurrentContext());
        //            interpolateSCMassDensityToCC(Q_var, getCurrentContext());
        //        }
        //
        Pointer<CellVariable<NDIM, double> > Cp_var = d_Q_Cp_map[Q_var];
        int Cp_current_idx;
        if (Cp_var)
        {
            Cp_current_idx = var_db->mapVariableAndContextToIndex(Cp_var, getCurrentContext());
            for (unsigned k = 0; k < d_reset_Cp_fcns.size(); ++k)
            {
                d_reset_Cp_fcns[k](Cp_current_idx,
                                   Cp_var,
                                   d_hier_math_ops,
                                   -1 /*cycle_num*/,
                                   apply_time,
                                   current_time,
                                   new_time,
                                   d_reset_Cp_fcns_ctx[k]);
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

            // Multiply rho*cp to div(uQ).
            if (rho_cc_var)
            {
                d_hier_cc_data_ops->multiply(d_C_current_idx, rho_current_idx, Cp_current_idx);
                d_hier_cc_data_ops->multiply(N_scratch_idx, N_scratch_idx, d_C_current_idx);
            }

            if (convective_time_stepping_type == FORWARD_EULER)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -0.5, N_scratch_idx, Q_rhs_scratch_idx);
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
TemperatureSemiImplicitHierarchyIntegrator::integrateHierarchy(const double current_time,
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
        IBAMR_DO_ONCE({
            pout << "TemperatureSemiImplicitHierarchyIntegrator::integrateHierarchy():\n"
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

        TimeSteppingType diffusion_time_stepping_type = d_Q_diffusion_time_stepping_type[Q_var];
        const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coef = d_Q_bc_coef[Q_var];
        Pointer<SideVariable<NDIM, double> > D_var = d_Q_diffusion_coef_variable[Q_var];
        Pointer<SideVariable<NDIM, double> > D_rhs_var = d_diffusion_coef_rhs_map[D_var];
        const int D_current_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getCurrentContext()) : -1);
        const int D_scratch_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getScratchContext()) : -1);
        const int D_new_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getNewContext()) : -1);
        const int D_rhs_scratch_idx =
            (D_rhs_var ? var_db->mapVariableAndContextToIndex(D_rhs_var, getScratchContext()) : -1);

        double apply_time = new_time;
        Pointer<CellVariable<NDIM, double> > rho_cc_var = d_Q_rho_map[Q_var];
        int rho_new_idx, rho_scratch_idx, rho_current_idx;
        if (rho_cc_var)
        {
            rho_new_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, getNewContext());
            rho_scratch_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, getScratchContext());
            for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
            {
                d_reset_rho_fcns[k](rho_new_idx,
                                    rho_cc_var,
                                    d_hier_math_ops,
                                    -1 /*cycle_num*/,
                                    apply_time,
                                    current_time,
                                    new_time,
                                    d_reset_rho_fcns_ctx[k]);
            }
            // interpolateSCMassDensityToCC(Q_var, getNewContext());
        }
        Pointer<CellVariable<NDIM, double> > Cp_var = d_Q_Cp_map[Q_var];
        int Cp_new_idx, Cp_scratch_idx, Cp_current_idx;
        if (Cp_var)
        {
            Cp_new_idx = var_db->mapVariableAndContextToIndex(Cp_var, getNewContext());
            Cp_scratch_idx = var_db->mapVariableAndContextToIndex(Cp_var, getScratchContext());
            for (unsigned k = 0; k < d_reset_Cp_fcns.size(); ++k)
            {
                d_reset_Cp_fcns[k](Cp_new_idx,
                                   Cp_var,
                                   d_hier_math_ops,
                                   -1 /*cycle_num*/,
                                   apply_time,
                                   current_time,
                                   new_time,
                                   d_reset_Cp_fcns_ctx[k]);
            }
        }

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

        // set rho*Cp/dt + K*lambda.
        const double lambda = d_Q_damping_coef[Q_var];
        if (rho_cc_var || Cp_var)
        {
            d_hier_cc_data_ops->multiply(d_C_new_idx, rho_new_idx, Cp_new_idx);
            d_hier_cc_data_ops->scale(d_C_new_idx, 1.0 / dt, d_C_new_idx);

            d_hier_cc_data_ops->addScalar(d_C_scratch_idx, d_C_new_idx, K * lambda);
            solver_spec.setCPatchDataId(d_C_scratch_idx);
        }
        else
        {
            solver_spec.setCConstant(1.0 / dt + K * lambda);
        }

        // set rho*cp/dt - (1.0-K)*lambda in rhs_op_spec, which is applied to Q later.
        if (rho_cc_var || Cp_var)
        {
            rho_current_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, getCurrentContext());
            Cp_current_idx = var_db->mapVariableAndContextToIndex(Cp_var, getCurrentContext());
            d_hier_cc_data_ops->multiply(d_C_current_idx, rho_current_idx, Cp_current_idx);
            d_hier_cc_data_ops->scale(d_C_current_idx, 1.0 / dt, d_C_current_idx);

            d_hier_cc_data_ops->addScalar(d_C_rhs_scratch_idx, d_C_current_idx, -(1.0 - K) * lambda);
            rhs_op_spec.setCPatchDataId(d_C_rhs_scratch_idx);
        }
        else
        {
            rhs_op_spec.setCConstant(1.0 / dt - (1.0 - K) * lambda);
        }

        if (rho_cc_var && isDiffusionCoefficientVariable(Q_var))
        {
            for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
            {
                d_reset_kappa_fcns[k](D_new_idx,
                                      D_var,
                                      d_hier_math_ops,
                                      -1 /*cycle_num*/,
                                      apply_time,
                                      current_time,
                                      new_time,
                                      d_reset_kappa_fcns_ctx[k]);
            }

            apply_time = current_time;
            for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
            {
                d_reset_kappa_fcns[k](D_current_idx,
                                      D_var,
                                      d_hier_math_ops,
                                      -1 /*cycle_num*/,
                                      apply_time,
                                      current_time,
                                      new_time,
                                      d_reset_kappa_fcns_ctx[k]);
            }
        }

        if (isDiffusionCoefficientVariable(Q_var))
        {
            // set -K*D in solver_spec (kappa is included within Db)
            d_hier_sc_data_ops->scale(D_scratch_idx, -K, D_new_idx);
            solver_spec.setDPatchDataId(D_scratch_idx);

            // set (1.0-K)*Db in rhs_op_spec
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
                    IBAMR_DO_ONCE({
                        pout << "TemperatureSemiImplicitHierarchyIntegrator::"
                                "integrateHierarchy():"
                                "\n"
                             << "  WARNING: convective_time_stepping_type = "
                             << enum_to_string<TimeSteppingType>(d_Q_convective_time_stepping_type[Q_var])
                             << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                             << "           using "
                             << enum_to_string<TimeSteppingType>(d_Q_convective_time_stepping_type[Q_var])
                             << " only for the first cycle in each time step;\n"
                             << "           using " << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                             << " for subsequent cycles.\n";
                    });
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

                    // compute rho^n+1/2 and Cp^n+1/2
                    if (rho_cc_var || Cp_var)
                    {
                        d_hier_cc_data_ops->linearSum(rho_scratch_idx, 0.5, rho_current_idx, 0.5, rho_new_idx);
                        d_hier_cc_data_ops->linearSum(Cp_scratch_idx, 0.5, Cp_current_idx, 0.5, Cp_new_idx);
                        d_hier_cc_data_ops->multiply(N_scratch_idx, N_scratch_idx, rho_scratch_idx);
                        d_hier_cc_data_ops->multiply(N_scratch_idx, N_scratch_idx, Cp_scratch_idx);
                    }
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

                    if (rho_cc_var || Cp_var)
                    {
                        d_hier_cc_data_ops->multiply(N_scratch_idx, N_scratch_idx, rho_new_idx);
                        d_hier_cc_data_ops->multiply(N_scratch_idx, N_scratch_idx, Cp_new_idx);
                    }
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

            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, 1.0, F_scratch_idx, Q_rhs_scratch_idx);
        }

        if (isDiffusionCoefficientVariable(Q_var) || (d_Q_diffusion_coef[Q_var] != 0.0))
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
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
TemperatureSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                          const double new_time,
                                                                          const bool skip_synchronize_new_state_data,
                                                                          const int num_cycles)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Deallocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_C_rhs_scratch_idx);
    }

    AdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
