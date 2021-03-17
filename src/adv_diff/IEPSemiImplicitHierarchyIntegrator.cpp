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
#include "ibamr/IEPSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CCLaplaceOperator.h"
#include "ibtk/CCPoissonSolverManager.h"
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

// FORTRAN ROUTINES
#if (NDIM == 2)
#define C_TO_S_CWISE_INTERP_FC IBTK_FC_FUNC(ctoscwiseinterp2nd2d, CTOSCWISEINTERP2ND2D)
#endif // if (NDIM == 2)

#if (NDIM == 3)
#define C_TO_S_CWISE_INTERP_FC IBTK_FC_FUNC(ctoscwiseinterp2nd3d, CTOSCWISEINTERP2ND3D)
#endif // if (NDIM == 3)

extern "C"
{
    void C_TO_S_CWISE_INTERP_FC(double* u0,
                                double* u1,
#if (NDIM == 3)
                                double* u2,
#endif
                                const int& u_gcw,
                                const double* V,
                                const int& V_gcw,
                                const int& ilower0,
                                const int& iupper0,
                                const int& ilower1,
                                const int& iupper1
#if (NDIM == 3)
                                ,
                                const int& ilower2,
                                const int& iupper2
#endif
    );
}

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
static const int FACEG = 1;
static const int NOGHOSTS = 0;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IEPSemiImplicitHierarchyIntegrator::IEPSemiImplicitHierarchyIntegrator(const std::string& object_name,
                                                                       Pointer<Database> input_db,
                                                                       bool register_for_restart)
    : AdvDiffSemiImplicitHierarchyIntegrator(object_name, input_db, register_for_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();

    getFromInput(input_db, from_restart);

    return;
} // IEPSemiImplicitHierarchyIntegrator

void
IEPSemiImplicitHierarchyIntegrator::registerSpecificHeatVariable(Pointer<CellVariable<NDIM, double> > Cp_var,
                                                                 const bool output_Cp)
{
    d_Cp_var = Cp_var;
    d_Cp_output = output_Cp;

    return;
} // registerSpecificHeatVariable

void
IEPSemiImplicitHierarchyIntegrator::registerDensityVariable(Pointer<CellVariable<NDIM, double> > rho_var,
                                                            const bool output_rho)
{
    d_rho_var = rho_var;
    d_rho_output = output_rho;

    return;
} // registerDensityVariable

void
IEPSemiImplicitHierarchyIntegrator::registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_rho_fcns.push_back(callback);
    d_reset_rho_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidDensityFcn

void
IEPSemiImplicitHierarchyIntegrator::registerResetSpecificHeatFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_Cp_fcns.push_back(callback);
    d_reset_Cp_fcns_ctx.push_back(ctx);
    return;
} // registerResetSpecificHeatFcn

void
IEPSemiImplicitHierarchyIntegrator::registerResetDiffusionCoefficientFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_kappa_fcns.push_back(callback);
    d_reset_kappa_fcns_ctx.push_back(ctx);
    return;
} // registerResetDiffusionCoefficientFcn

// void
// IEPSemiImplicitHierarchyIntegrator::registerINSVCStaggeredHierarchyIntegrator(
//    Pointer<INSVCStaggeredHierarchyIntegrator> ins_cons_hier_integrator)
//{
//    d_ins_hierarchy_integrator = ins_cons_hier_integrator;
//    return;
//} // registerINSVCStaggeredHierarchyIntegrator

void
IEPSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
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

    if (d_T_solver_type == CCPoissonSolverManager::UNDEFINED)
    {
        d_T_solver_type = CCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER;
    }
    if (d_T_precond_type == CCPoissonSolverManager::UNDEFINED)
    {
        const int max_levels = gridding_alg->getMaxLevels();
        if (max_levels == 1)
        {
            d_T_precond_type = CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
        }
        else
        {
            d_T_precond_type = CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
        }
        d_T_precond_db->putInteger("max_iterations", 1);
    }

    d_lf_solver = getHelmholtzSolverLiquidFractionEquation(d_lf_var);
    d_T_solver = getHelmholtzSolverTemperatureEquation(d_T_var);
    d_lf_rhs_op = getHelmholtzRHSOperatorLiquidFractionEquation(d_lf_var);
    d_T_rhs_op = getHelmholtzRHSOperatorTemperatureEquation(d_T_var);

    // Register additional variables required for present time stepping algorithm.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> face_ghosts = FACEG;
    const IntVector<NDIM> no_ghosts = NOGHOSTS;

    if (d_ls_var)
        registerVariable(d_ls_current_idx,
                         d_ls_new_idx,
                         d_ls_scratch_idx,
                         d_ls_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_ls_init);

    if (d_lf_var)
        registerVariable(d_lf_current_idx,
                         d_lf_new_idx,
                         d_lf_scratch_idx,
                         d_lf_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_lf_init);

    if (d_T_var)
        registerVariable(d_T_current_idx,
                         d_T_new_idx,
                         d_T_scratch_idx,
                         d_T_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_T_init);

    int lf_F_current_idx, lf_F_scratch_idx, lf_F_new_idx;
    if (d_lf_F_var)
        registerVariable(lf_F_current_idx,
                         lf_F_new_idx,
                         lf_F_scratch_idx,
                         d_lf_F_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int T_F_current_idx, T_F_scratch_idx, T_F_new_idx;
    if (d_T_F_var)
        registerVariable(T_F_current_idx,
                         T_F_new_idx,
                         T_F_scratch_idx,
                         d_T_F_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int lf_diff_coef_current_idx, lf_diff_coef_scratch_idx, lf_diff_coef_new_idx;
    if (d_lf_diffusion_coef_var)
        registerVariable(lf_diff_coef_current_idx,
                         lf_diff_coef_new_idx,
                         lf_diff_coef_scratch_idx,
                         d_lf_diffusion_coef_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int T_diff_coef_current_idx, T_diff_coef_scratch_idx, T_diff_coef_new_idx;
    if (d_T_diffusion_coef_var)
        registerVariable(T_diff_coef_current_idx,
                         T_diff_coef_new_idx,
                         T_diff_coef_scratch_idx,
                         d_T_diffusion_coef_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int lf_diffusion_coef_rhs_scratch_idx, T_diffusion_coef_rhs_scratch_idx;
    d_lf_diffusion_coef_rhs_var = new SideVariable<NDIM, double>("lf_diff_var");
    registerVariable(lf_diffusion_coef_rhs_scratch_idx, d_lf_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    d_T_diffusion_coef_rhs_var = new SideVariable<NDIM, double>("T_diff_var");
    registerVariable(T_diffusion_coef_rhs_scratch_idx, d_T_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    int lf_rhs_scratch_idx, T_rhs_scratch_idx;
    registerVariable(lf_rhs_scratch_idx, d_lf_rhs_var, cell_ghosts, getScratchContext());
    registerVariable(T_rhs_scratch_idx, d_T_rhs_var, cell_ghosts, getScratchContext());

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_lf_C_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::C");
    d_lf_C_idx = var_db->registerVariableAndContext(d_lf_C_var, getCurrentContext(), no_ghosts);

    d_T_C_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::C");
    d_T_C_idx = var_db->registerVariableAndContext(d_T_C_var, getCurrentContext(), no_ghosts);

    d_lf_temp_rhs_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::temp_rhs");
    d_lf_temp_rhs_idx = var_db->registerVariableAndContext(d_lf_temp_rhs_var, getCurrentContext(), no_ghosts);

    d_T_temp_rhs_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::temp_rhs");
    d_T_temp_rhs_idx = var_db->registerVariableAndContext(d_T_temp_rhs_var, getCurrentContext(), no_ghosts);

    int rho_current_idx, rho_new_idx, rho_scratch_idx;
    registerVariable(rho_current_idx,
                     rho_new_idx,
                     rho_scratch_idx,
                     d_rho_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_rho_init);

    if (d_visit_writer) d_visit_writer->registerPlotQuantity(d_rho_var->getName(), "SCALAR", rho_current_idx);

    int Cp_current_idx, Cp_new_idx, Cp_scratch_idx;
    registerVariable(Cp_current_idx,
                     Cp_new_idx,
                     Cp_scratch_idx,
                     d_Cp_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");
    if (d_visit_writer) d_visit_writer->registerPlotQuantity(d_Cp_var->getName(), "SCALAR", Cp_current_idx);

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

    d_H_var = new CellVariable<NDIM, double>("H_var");
    registerVariable(d_H_current_idx,
                     d_H_new_idx,
                     d_H_scratch_idx,
                     d_H_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_g_firstder_var = new CellVariable<NDIM, double>("g_firstder_var");
    d_g_firstder_idx = var_db->registerVariableAndContext(d_g_firstder_var, getCurrentContext());

    d_g_secondder_var = new CellVariable<NDIM, double>("g_secondder_var");
    d_g_secondder_idx = var_db->registerVariableAndContext(d_g_secondder_var, getCurrentContext());

    d_p_firstder_var = new CellVariable<NDIM, double>("p_firstder_var");
    d_p_firstder_idx = var_db->registerVariableAndContext(d_p_firstder_var, getCurrentContext());

    d_D_cc_var = new CellVariable<NDIM, double>("D_cc", NDIM);
    registerVariable(d_D_cc_current_idx,
                     d_D_cc_new_idx,
                     d_D_cc_scratch_idx,
                     d_D_cc_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");
    if (d_visit_writer) d_visit_writer->registerPlotQuantity(d_D_cc_var->getName(), "SCALAR", d_D_cc_current_idx);

    // Register variables for plotting.
    if (d_visit_writer)
    {
        if (d_output_ls)
        {
            d_visit_writer->registerPlotQuantity("level_set", "SCALAR", d_ls_current_idx, 0);
        }

        if (d_output_lf)
        {
            d_visit_writer->registerPlotQuantity("liquid_fraction", "SCALAR", d_lf_current_idx, 0);
        }

        if (d_output_T)
        {
            d_visit_writer->registerPlotQuantity("Temperature", "SCALAR", d_T_current_idx, 0);
        }
        d_visit_writer->registerPlotQuantity("Heaviside", "SCALAR", d_H_current_idx, 0);
    }

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffSemiImplicitHierarchyIntegrator.
    AdvDiffSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
IEPSemiImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                 const double new_time,
                                                                 const int num_cycles)
{
    AdvDiffSemiImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

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
        d_T_solver_needs_init = true;
        d_T_rhs_op_needs_init = true;
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_C_rhs_scratch_idx, current_time);
        if (!level->checkAllocated(d_lf_C_idx)) level->allocatePatchData(d_lf_C_idx, current_time);
        if (!level->checkAllocated(d_T_C_idx)) level->allocatePatchData(d_T_C_idx, current_time);
        if (!level->checkAllocated(d_lf_temp_rhs_idx)) level->allocatePatchData(d_lf_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_T_temp_rhs_idx)) level->allocatePatchData(d_T_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_g_firstder_idx))
            level->allocatePatchData(d_g_firstder_idx, current_time); // should this be new time?
        if (!level->checkAllocated(d_g_secondder_idx)) level->allocatePatchData(d_g_secondder_idx, current_time);
        if (!level->checkAllocated(d_p_firstder_idx)) level->allocatePatchData(d_p_firstder_idx, current_time);
    }

    int ls_current_idx;
    for (auto Q_var : d_Q_var)
    {
        ls_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
    }
    computeHeavisideFunction(d_H_current_idx, ls_current_idx);
    d_hier_cc_data_ops->copyData(d_H_new_idx, d_H_current_idx);

    // Set the initial guess for liquid fraction, temperature and density.
    d_hier_cc_data_ops->copyData(d_lf_new_idx, d_lf_current_idx);
    d_hier_cc_data_ops->copyData(d_T_new_idx, d_T_current_idx);

    int rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());
    int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());
    d_hier_cc_data_ops->copyData(rho_new_idx, rho_current_idx);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
IEPSemiImplicitHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                       const double new_time,
                                                       const int cycle_num)
{
    AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE({
            pout << "IEPSemiImplicitHierarchyIntegrator::integrateHierarchy():\n"
                 << "  WARNING: num_cycles = " << d_current_num_cycles
                 << " but expected num_cycles = " << expected_num_cycles << ".\n";
        });
    }

    // Perform a single step of fixed point iteration.
    const int lf_current_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getCurrentContext());
    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
    const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
    const int lf_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_rhs_var, getScratchContext());
    const int lf_diff_coef_new_idx = (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_var, getNewContext()));
    const int lf_diff_coef_current_idx =
        (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_var, getCurrentContext()));
    const int lf_diff_coef_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_var, getScratchContext()));
    const int lf_diff_coef_rhs_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_rhs_var, getScratchContext()));

    // setting k equation diffusion timestepping type
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
    PoissonSpecifications lf_solver_spec(d_object_name + "::solver_spec::" + d_lf_var->getName());
    PoissonSpecifications lf_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_lf_var->getName());

    // set C coefficients.
    // compute heaviside based on phi^n+1.
    int ls_new_idx, ls_scratch_idx;
    for (auto Q_var : d_Q_var)
    {
        ls_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        ls_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
    }
    computeHeavisideFunction(d_H_new_idx, ls_new_idx);
    computeDoubleWellPotential(d_g_firstder_idx, d_g_secondder_idx, d_lf_new_idx);

    d_hier_cc_data_ops->scale(d_lf_C_idx, 1.0 / dt, d_H_new_idx);
    d_hier_cc_data_ops->multiply(d_C_rhs_scratch_idx, d_g_secondder_idx, d_H_new_idx);
    d_hier_cc_data_ops->scale(d_C_rhs_scratch_idx, d_M_lf * d_lambda_lf / std::pow(d_eta_lf, 2.0), d_C_rhs_scratch_idx);
    d_hier_cc_data_ops->add(d_lf_C_idx, d_C_rhs_scratch_idx, d_lf_C_idx);
    lf_solver_spec.setCPatchDataId(d_lf_C_idx);

    d_hier_cc_data_ops->scale(d_lf_temp_rhs_idx, 1.0 / dt, d_H_current_idx);
    lf_rhs_op_spec.setCPatchDataId(d_lf_temp_rhs_idx);

    // set D coefficients.
    // Interpolate the cell-centered Heaviside to side-centered.
    d_hier_cc_data_ops->copyData(d_H_scratch_idx, d_H_new_idx);
    d_H_bdry_bc_fill_op->fillData(new_time);

    interpolateCCHeaviside(lf_diff_coef_new_idx, d_H_scratch_idx);
    //    d_hier_math_ops->interp(lf_diff_coef_new_idx,
    //                            d_lf_diffusion_coef_var,
    //                            true,
    //                            d_H_scratch_idx,
    //                            d_H_var,
    //                            d_no_fill_op,
    //                            d_integrator_time);

    d_hier_sc_data_ops->scale(lf_diff_coef_new_idx, d_M_lf * d_lambda_lf, lf_diff_coef_new_idx);
    d_hier_sc_data_ops->scale(lf_diff_coef_scratch_idx, -alpha, lf_diff_coef_new_idx);
    lf_solver_spec.setDPatchDataId(lf_diff_coef_scratch_idx);

    d_hier_cc_data_ops->copyData(d_H_scratch_idx, d_H_current_idx);
    d_H_bdry_bc_fill_op->fillData(current_time);

    interpolateCCHeaviside(lf_diff_coef_current_idx, d_H_scratch_idx);
    //    d_hier_math_ops->interp(lf_diff_coef_current_idx,
    //                            d_lf_diffusion_coef_var,
    //                            true,
    //                            d_H_scratch_idx,
    //                            d_H_var,
    //                            d_no_fill_op,
    //                            d_integrator_time);

    d_hier_sc_data_ops->scale(lf_diff_coef_current_idx, d_M_lf * d_lambda_lf, lf_diff_coef_current_idx);
    d_hier_sc_data_ops->scale(lf_diff_coef_rhs_scratch_idx, (1.0 - alpha), lf_diff_coef_current_idx);
    lf_rhs_op_spec.setDPatchDataId(lf_diff_coef_rhs_scratch_idx);

    // Initialize the RHS operator and compute the RHS vector for k equation.
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

    // Initialize the linear solver for k equation.
    Pointer<PoissonSolver> lf_solver = d_lf_solver;
    lf_solver->setPoissonSpecifications(lf_solver_spec);
    lf_solver->setPhysicalBcCoef(d_lf_bc_coef);
    lf_solver->setHomogeneousBc(false);
    lf_solver->setSolutionTime(new_time);
    lf_solver->setTimeInterval(current_time, new_time);
    if (d_lf_solver_needs_init)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the solvers for" << d_lf_var->getName() << "\n";
        }
        lf_solver->initializeSolverState(*d_lf_sol, *d_lf_rhs);
        d_lf_solver_needs_init = false;
    }

    // Account for forcing terms.
    const int lf_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_F_var, getScratchContext());
    const int lf_F_new_idx = var_db->mapVariableAndContextToIndex(d_lf_F_var, getNewContext());

    if (d_lf_F_var)
    {
        computeInterpolationFunction(d_p_firstder_idx, d_lf_new_idx, d_T_new_idx);
        computeLiquidFractionSourceTerm(lf_F_scratch_idx);
        d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, +1.0, lf_F_scratch_idx, lf_rhs_scratch_idx);
        std::cout << "L2 norm of lf_rhs_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(lf_rhs_scratch_idx) << std::endl;
    }

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
    //    if (d_k_u_var)
    //    {
    //        const int k_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_N_var, getScratchContext());
    //        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
    //        {
    //            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, +1.0, k_N_scratch_idx, k_rhs_scratch_idx);
    //        }
    //        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
    //        {
    //            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, +0.5, k_N_scratch_idx, k_rhs_scratch_idx);
    //        }
    //    }
    if (d_lf_F_var)
    {
        d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, -1.0, lf_F_scratch_idx, lf_rhs_scratch_idx);
        d_hier_cc_data_ops->copyData(lf_F_new_idx, lf_F_scratch_idx);
    }

    if (d_solve_energy)
    {
        // Stuff related to temperature equation.
        int rho_new_idx, rho_scratch_idx, rho_current_idx;
        rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());
        rho_scratch_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getScratchContext());
        rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());

        int Cp_new_idx, Cp_scratch_idx, Cp_current_idx;
        Cp_new_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getNewContext());
        Cp_scratch_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getScratchContext());
        Cp_current_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getCurrentContext());

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

        const int T_current_idx = var_db->mapVariableAndContextToIndex(d_T_var, getCurrentContext());
        const int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, getScratchContext());
        const int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());
        const int T_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_rhs_var, getScratchContext());
        const int T_diff_coef_new_idx = (var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_var, getNewContext()));
        const int T_diff_coef_current_idx =
            (var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_var, getCurrentContext()));
        const int T_diff_coef_scratch_idx =
            (var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_var, getScratchContext()));
        const int T_diff_coef_rhs_scratch_idx =
            (var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_rhs_var, getScratchContext()));

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
        PoissonSpecifications T_solver_spec(d_object_name + "::solver_spec::" + d_lf_var->getName());
        PoissonSpecifications T_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_lf_var->getName());

        // set rho*Cp/dt + K*lambda.
        const double lambda = 0.0;
        d_hier_cc_data_ops->multiply(d_C_new_idx, rho_new_idx, Cp_new_idx);
        d_hier_cc_data_ops->scale(d_C_new_idx, 1.0 / dt, d_C_new_idx);
        d_hier_cc_data_ops->copyData(d_T_C_idx, d_C_new_idx);
        T_solver_spec.setCPatchDataId(d_T_C_idx);

        apply_time = current_time;
        for (unsigned k = 0; k < d_reset_Cp_fcns.size(); ++k)
        {
            d_reset_Cp_fcns[k](Cp_current_idx,
                               d_Cp_var,
                               d_hier_math_ops,
                               -1 /*cycle_num*/,
                               apply_time,
                               current_time,
                               new_time,
                               d_reset_Cp_fcns_ctx[k]);
        }

        d_hier_cc_data_ops->multiply(d_C_current_idx, rho_current_idx, Cp_current_idx);
        d_hier_cc_data_ops->scale(d_C_current_idx, 1.0 / dt, d_C_current_idx);
        d_hier_cc_data_ops->copyData(d_T_temp_rhs_idx, d_C_current_idx);
        T_rhs_op_spec.setCPatchDataId(d_T_temp_rhs_idx);

        apply_time = new_time;
        for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
        {
            d_reset_kappa_fcns[k](T_diff_coef_new_idx,
                                  d_T_diffusion_coef_var,
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
            d_reset_kappa_fcns[k](T_diff_coef_current_idx,
                                  d_T_diffusion_coef_var,
                                  d_hier_math_ops,
                                  -1 /*cycle_num*/,
                                  apply_time,
                                  current_time,
                                  new_time,
                                  d_reset_kappa_fcns_ctx[k]);

            static const bool synch_cf_interface = true;
            d_hier_math_ops->interp(d_D_cc_new_idx,
                                    d_D_cc_var,
                                    T_diff_coef_new_idx,
                                    d_T_diffusion_coef_var,
                                    d_no_fill_op,
                                    d_integrator_time,
                                    synch_cf_interface);
        }

        d_hier_sc_data_ops->scale(T_diff_coef_scratch_idx, -alpha, T_diff_coef_new_idx);
        T_solver_spec.setDPatchDataId(T_diff_coef_scratch_idx);

        d_hier_sc_data_ops->scale(T_diff_coef_rhs_scratch_idx, (1.0 - alpha), T_diff_coef_current_idx);
        T_rhs_op_spec.setDPatchDataId(T_diff_coef_rhs_scratch_idx);

        // Initialize the RHS operator and compute the RHS vector for k equation.
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

        // Initialize the linear solver for temperature equation.
        Pointer<PoissonSolver> T_solver = d_T_solver;
        T_solver->setPoissonSpecifications(T_solver_spec);
        T_solver->setPhysicalBcCoef(d_T_bc_coef);
        T_solver->setHomogeneousBc(false);
        T_solver->setSolutionTime(new_time);
        T_solver->setTimeInterval(current_time, new_time);
        if (d_T_solver_needs_init)
        {
            if (d_enable_logging)
            {
                plog << d_object_name << ": "
                     << "Initializing the solvers for" << d_T_var->getName() << "\n";
            }
            T_solver->initializeSolverState(*d_T_sol, *d_T_rhs);
            d_T_solver_needs_init = false;
        }

        // Account for forcing terms.
        const int T_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_F_var, getScratchContext());
        const int T_F_new_idx = var_db->mapVariableAndContextToIndex(d_T_F_var, getNewContext());

        if (d_T_F_var)
        {
            std::cout << "L2 norm of T_rhs_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(T_rhs_scratch_idx) << std::endl;
            computeTemperatureSourceTerm(T_F_scratch_idx, dt);
            d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, +1.0, T_F_scratch_idx, T_rhs_scratch_idx);
            std::cout << "L2 norm of T_rhs_scratch_idx after F\t" << d_hier_cc_data_ops->L2Norm(T_rhs_scratch_idx)
                      << std::endl;
        }

        // Solve for lf(n+1).
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
        // Reset the right-hand side vector.
        //    if (d_k_u_var)
        //    {
        //        const int k_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_N_var, getScratchContext());
        //        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type ==
        //        MIDPOINT_RULE)
        //        {
        //            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, +1.0, k_N_scratch_idx, k_rhs_scratch_idx);
        //        }
        //        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        //        {
        //            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, +0.5, k_N_scratch_idx, k_rhs_scratch_idx);
        //        }
        //    }
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
IEPSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
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
        level->deallocatePatchData(d_lf_C_idx);
        level->deallocatePatchData(d_T_C_idx);
        level->deallocatePatchData(d_lf_temp_rhs_idx);
        level->deallocatePatchData(d_T_temp_rhs_idx);
    }

    AdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

Pointer<CellVariable<NDIM, double> >
IEPSemiImplicitHierarchyIntegrator::getHeavisideVariable() const
{
    return d_H_var;
} // getHeavisideVariable

void
IEPSemiImplicitHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    AdvDiffSemiImplicitHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
        base_hierarchy, coarsest_level, finest_level);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());

    // Setup the patch boundary filling objects.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent H_bc_component(d_H_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     BDRY_EXTRAP_TYPE,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_lf_bc_coef); // liquid fraction boundary condition is used.
    d_H_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_H_bdry_bc_fill_op->initializeOperatorState(H_bc_component, d_hierarchy);

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

    const int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, getScratchContext());
    d_T_sol = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::sol_vec::" + d_T_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_T_sol->addComponent(d_T_var, T_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    const int T_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_rhs_var, getScratchContext());
    d_T_rhs = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::rhs_vec::" + d_T_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_T_rhs->addComponent(d_T_rhs_var, T_rhs_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    d_lf_solver_needs_init = true;
    d_T_solver_needs_init = true;

    //   d_lf_convective_op_needs_init = true;
    //    d_T_convective_op_needs_init = true;
    return;
}

void
IEPSemiImplicitHierarchyIntegrator::registerLevelSetVariable(Pointer<CellVariable<NDIM, double> > ls_var,
                                                             const bool output_ls_var)
{
    d_ls_var = ls_var;
    d_output_ls = output_ls_var;
    return;
} // registerLevelSetVariable

void
IEPSemiImplicitHierarchyIntegrator::registerLiquidFractionVariable(Pointer<CellVariable<NDIM, double> > lf_var,
                                                                   const bool output_lf_var)
{
    d_lf_var = lf_var;
    d_output_lf = output_lf_var;

    Pointer<CellDataFactory<NDIM, double> > lf_factory = lf_var->getPatchDataFactory();
    const int lf_depth = lf_factory->getDefaultDepth();
    // Pointer<CellVariable<NDIM, double> > k_u_var = new CellVariable<NDIM, double>(k_var->getName() + "::u", k_depth);
    Pointer<CellVariable<NDIM, double> > lf_rhs_var =
        new CellVariable<NDIM, double>(lf_var->getName() + "::lf_rhs", lf_depth);
    Pointer<CellVariable<NDIM, double> > lf_F_var = new CellVariable<NDIM, double>(lf_var->getName() + "::F", lf_depth);
    Pointer<SideVariable<NDIM, double> > lf_diff_coef_var =
        new SideVariable<NDIM, double>(lf_var->getName() + "::diff_coef", lf_depth);

    // Set default values.
    // d_k_u_var = nullptr;
    d_lf_F_var = lf_F_var;
    d_lf_rhs_var = lf_rhs_var;
    d_lf_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    // d_lf_convective_difference_form = d_default_convective_difference_form;
    // if(!d_k_convective_time_stepping_type) d_k_convective_time_stepping_type =
    // d_default_convective_time_stepping_type;
    d_lf_diffusion_coef_var = lf_diff_coef_var;
    d_lf_init = nullptr;
    d_lf_F_fcn = nullptr;
    d_lf_bc_coef = nullptr;
    return;
} // registerLiquidFractionVariable

void
IEPSemiImplicitHierarchyIntegrator::registerTemperatureVariable(Pointer<CellVariable<NDIM, double> > T_var,
                                                                const bool output_T_var)
{
    d_T_var = T_var;
    d_output_T = output_T_var;

    Pointer<CellDataFactory<NDIM, double> > T_factory = T_var->getPatchDataFactory();
    const int T_depth = T_factory->getDefaultDepth();
    // Pointer<CellVariable<NDIM, double> > k_u_var = new CellVariable<NDIM, double>(k_var->getName() + "::u", k_depth);
    Pointer<CellVariable<NDIM, double> > T_rhs_var =
        new CellVariable<NDIM, double>(T_var->getName() + "::T_rhs", T_depth);
    Pointer<CellVariable<NDIM, double> > T_F_var = new CellVariable<NDIM, double>(T_var->getName() + "::F", T_depth);
    Pointer<SideVariable<NDIM, double> > T_diff_coef_var =
        new SideVariable<NDIM, double>(T_var->getName() + "::diff_coef", T_depth);

    // Set default values.
    // d_k_u_var = nullptr;
    d_T_F_var = T_F_var;
    d_T_rhs_var = T_rhs_var;
    d_T_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    // d_T_convective_difference_form = d_default_convective_difference_form;
    // if(!d_k_convective_time_stepping_type) d_k_convective_time_stepping_type =
    // d_default_convective_time_stepping_type;
    d_T_diffusion_coef_var = T_diff_coef_var;
    d_T_init = nullptr;
    d_T_F_fcn = nullptr;
    d_T_bc_coef = nullptr;
    return;
} // registerTemperatureVariable

Pointer<PoissonSolver>
IEPSemiImplicitHierarchyIntegrator::getHelmholtzSolverLiquidFractionEquation(
    Pointer<CellVariable<NDIM, double> > lf_var)
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
} // getHelmholtzSolverLiquidFractionEquation

Pointer<PoissonSolver>
IEPSemiImplicitHierarchyIntegrator::getHelmholtzSolverTemperatureEquation(Pointer<CellVariable<NDIM, double> > T_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_T_var);
#endif
    if (!d_T_solver)
    {
        const std::string& name = T_var->getName();
        d_T_solver =
            CCPoissonSolverManager::getManager()->allocateSolver(d_T_solver_type,
                                                                 d_object_name + "::helmholtz_solver::" + name,
                                                                 d_T_solver_db,
                                                                 "temperature_",
                                                                 d_T_precond_type,
                                                                 d_object_name + "::helmholtz_precond::" + name,
                                                                 d_T_precond_db,
                                                                 "temperature_pc_");
        d_T_solver_needs_init = true;
    }
    return d_T_solver;
} // getHelmholtzSolverTemperatureEquation

Pointer<LaplaceOperator>
IEPSemiImplicitHierarchyIntegrator::getHelmholtzRHSOperatorLiquidFractionEquation(
    Pointer<CellVariable<NDIM, double> > lf_var)
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
} // getHelmholtzRHSOperatorLiquidFractionEquation

Pointer<LaplaceOperator>
IEPSemiImplicitHierarchyIntegrator::getHelmholtzRHSOperatorTemperatureEquation(
    Pointer<CellVariable<NDIM, double> > T_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_T_var);
#endif
    const std::string& name = T_var->getName();
    if (!d_T_rhs_op)
    {
        d_T_rhs_op = new CCLaplaceOperator(d_object_name + "::helmholtz_rhs_op::" + name, /*homogeneous_bc*/ false);
        d_T_rhs_op_needs_init = true;
    }
    return d_T_rhs_op;
} // getHelmholtzRHSOperatorTemperatureEquation

void
IEPSemiImplicitHierarchyIntegrator::setInitialConditionsLiquidFractionEquation(
    Pointer<CellVariable<NDIM, double> > lf_var,
    Pointer<IBTK::CartGridFunction> lf_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_lf_var);
#endif
    d_lf_init = lf_init;
    return;
} // setInitialConditionsLiquidFractionEquation

void
IEPSemiImplicitHierarchyIntegrator::setInitialConditionsTemperatureEquation(Pointer<CellVariable<NDIM, double> > T_var,
                                                                            Pointer<IBTK::CartGridFunction> T_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_T_var);
#endif
    d_T_init = T_init;
    return;
} // setInitialConditionsTemperatureEquation

void
IEPSemiImplicitHierarchyIntegrator::setDensityInitialCondition(Pointer<CellVariable<NDIM, double> > rho_var,
                                                               Pointer<IBTK::CartGridFunction> rho_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_var);
#endif
    d_rho_init = rho_init;
    return;
} // setDensityInitialCondition

void
IEPSemiImplicitHierarchyIntegrator::setPhysicalBcCoefLiquidFractionEquation(Pointer<CellVariable<NDIM, double> > lf_var,
                                                                            RobinBcCoefStrategy<NDIM>* lf_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_lf_var);
#endif
    d_lf_bc_coef = lf_bc_coef;
    return;
} // setPhysicalBcCoefLiquidFractionEquation

RobinBcCoefStrategy<NDIM>*
IEPSemiImplicitHierarchyIntegrator::getPhysicalBcCoefLiquidFractionEquation()
{
    return d_lf_bc_coef;
} // getPhysicalBcCoefLiquidFractionEquation

void
IEPSemiImplicitHierarchyIntegrator::setPhysicalBcCoefTemperatureEquation(Pointer<CellVariable<NDIM, double> > T_var,
                                                                         RobinBcCoefStrategy<NDIM>* T_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_T_var);
#endif
    d_T_bc_coef = T_bc_coef;
    return;
} // setPhysicalBcCoefTemperatureEquation

RobinBcCoefStrategy<NDIM>*
IEPSemiImplicitHierarchyIntegrator::getPhysicalBcCoefTemperatureEquation()
{
    return d_T_bc_coef;
} // getPhysicalBcCoefTemperatureEquation

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IEPSemiImplicitHierarchyIntegrator::computeHeavisideFunction(int H_idx, const int phi_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const double* patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*phi_data)(ci);
                (*H_data)(ci) = IBTK::smooth_heaviside(phi, alpha);
            }
        }
    }

    return;
} // computeHeavisideFunction

void
IEPSemiImplicitHierarchyIntegrator::computeDoubleWellPotential(int g_firstder_idx,
                                                               int g_secondder_idx,
                                                               const int liquid_fraction_idx)
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
            Pointer<CellData<NDIM, double> > g_firstder_data = patch->getPatchData(g_firstder_idx);
            Pointer<CellData<NDIM, double> > g_secondder_data = patch->getPatchData(g_secondder_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double lf = (*lf_data)(ci);
                (*g_firstder_data)(ci) = 2.0 * (lf - 1.0) * lf * (2.0 * lf - 1.0);
                (*g_secondder_data)(ci) = 12.0 * lf * lf - 12.0 * lf + 2.0;
            }
        }
    }
    pout << "L2 norm of g'\t" << d_hier_cc_data_ops->L2Norm(g_firstder_idx) << "\n";
    pout << "L2 norm of g''\t" << d_hier_cc_data_ops->L2Norm(g_secondder_idx) << "\n";
    return;
} //  computeDoubleWellPotential

void
IEPSemiImplicitHierarchyIntegrator::computeInterpolationFunction(int p_firstder_idx,
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
            Pointer<CellData<NDIM, double> > p_firstder_data = patch->getPatchData(p_firstder_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double lf = (*lf_data)(ci);
                const double T = (*T_data)(ci);
                // std::cout << "T value is\t" << T << "\tand lf value is\t" << lf  << "\tT_ref is\t" << d_T_ref <<
                // std::endl;
                if (lf == 1.0 && T <= d_T_ref)
                {
                    (*p_firstder_data)(ci) = 1.0;
                }
                else if (lf == 0.0 && T >= d_T_ref)
                {
                    (*p_firstder_data)(ci) = 1.0;
                }
                else
                {
                    (*p_firstder_data)(ci) = 30.0 * std::pow(lf, 4.0) - 60.0 * std::pow(lf, 3.0) + 30.0 * lf * lf;
                    // (*p_firstder_data)(ci) = 6.0*lf - 6.0*lf*lf;
                }
                // std::cout << "p' value is\t" << (*p_firstder_data)(ci) <<  std::endl;
            }
        }
    }
    pout << "L2 norm of p'\t" << d_hier_cc_data_ops->L2Norm(p_firstder_idx) << "\n";
    return;
} // computeInterpolationFunction

void
IEPSemiImplicitHierarchyIntegrator::interpolateCCHeaviside(int lf_diff_coef_idx, const int H_idx)
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

            Pointer<SideData<NDIM, double> > diff_coef_data = patch->getPatchData(lf_diff_coef_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);

            C_TO_S_CWISE_INTERP_FC(diff_coef_data->getPointer(0),
                                   diff_coef_data->getPointer(1),
#if (NDIM == 3)
                                   diff_coef_data->getPointer(2),
#endif
                                   diff_coef_data->getGhostCellWidth().max(),
                                   H_data->getPointer(),
                                   H_data->getGhostCellWidth().max(),
                                   patch_box.lower(0),
                                   patch_box.upper(0),
                                   patch_box.lower(1),
                                   patch_box.upper(1)
#if (NDIM == 3)
                                       ,
                                   patch_box.lower(2),
                                   patch_box.upper(2)
#endif
            );
        }
    }

    return;
} // interpolateCCHeaviside

void
IEPSemiImplicitHierarchyIntegrator::computeLiquidFractionSourceTerm(int F_scratch_idx)
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

            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(d_lf_new_idx);
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(d_T_new_idx);
            Pointer<CellData<NDIM, double> > p_firstder_data = patch->getPatchData(d_p_firstder_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(d_H_new_idx);
            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_scratch_idx);
            Pointer<CellData<NDIM, double> > g_firstder_data = patch->getPatchData(d_g_firstder_idx);
            Pointer<CellData<NDIM, double> > g_secondder_data = patch->getPatchData(d_g_secondder_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                double F = -d_M_lf * d_rho_liquid * d_latent_heat * (*H_data)(ci) * (*p_firstder_data)(ci) *
                           (d_T_ref - (*T_data)(ci)) / d_T_ref;
                (*F_data)(ci) = F - (d_M_lf * d_lambda_lf * (*H_data)(ci) / std::pow(d_eta_lf, 2.0) *
                                     ((*g_firstder_data)(ci) - ((*g_secondder_data)(ci) * (*lf_data)(ci))));

                // g' is evaluated at previous cycle.
                // (*F_data)(ci) = F - (d_M_lf * d_lambda_lf * (*H_data)(ci) / std::pow(d_eta_lf, 2.0) *
                // (*g_firstder_data)(ci));
            }
        }
    }
    return;
} // computeLiquidFractionSourceTerm

void
IEPSemiImplicitHierarchyIntegrator::computeTemperatureSourceTerm(int F_scratch_idx, const double dt)
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
            Pointer<CellData<NDIM, double> > lf_new_data = patch->getPatchData(d_lf_new_idx);
            Pointer<CellData<NDIM, double> > H_new_data = patch->getPatchData(d_H_new_idx);
            Pointer<CellData<NDIM, double> > lf_current_data = patch->getPatchData(d_lf_current_idx);
            Pointer<CellData<NDIM, double> > H_current_data = patch->getPatchData(d_H_current_idx);
            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_scratch_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*F_data)(ci) =
                    -d_rho_liquid * d_latent_heat *
                    (((*H_new_data)(ci) * (*lf_new_data)(ci)) - ((*H_current_data)(ci) * (*lf_current_data)(ci))) / dt;
            }
        }
    }
    return;
} // computeTemperatureSourceTerm

void
IEPSemiImplicitHierarchyIntegrator::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db);
#endif
    // Read in data members from input database.
    if (!is_from_restart)
    {
        d_lf_solver_type = CCPoissonSolverManager::UNDEFINED;
        d_lf_precond_type = CCPoissonSolverManager::UNDEFINED;
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

        d_T_solver_type = CCPoissonSolverManager::UNDEFINED;
        d_T_precond_type = CCPoissonSolverManager::UNDEFINED;
        if (input_db->keyExists("T_solver_type"))
        {
            d_T_solver_type = input_db->getString("T_solver_type");
            if (input_db->keyExists("T_solver_db")) d_T_solver_db = input_db->getDatabase("T_solver_db");
        }
        if (!d_T_solver_db) d_T_solver_db = new MemoryDatabase("T_solver_db");

        if (input_db->keyExists("T_precond_type"))
        {
            d_T_precond_type = input_db->getString("T_precond_type");
            if (input_db->keyExists("T_precond_db")) d_T_precond_db = input_db->getDatabase("T_precond_db");
        }
        if (!d_T_precond_db) d_T_precond_db = new MemoryDatabase("T_precond_db");

        if (input_db->keyExists("latent_heat")) d_latent_heat = input_db->getDouble("latent_heat");
        if (input_db->keyExists("phase_change")) d_phase_change = input_db->getBool("phase_change");
        if (input_db->keyExists("rho_liquid")) d_rho_liquid = input_db->getDouble("rho_liquid");
        if (input_db->keyExists("T_ref")) d_T_ref = input_db->getDouble("T_ref");

        if (input_db->keyExists("num_interface_cells"))
            d_num_interface_cells = input_db->getInteger("num_interface_cells");

        if (input_db->keyExists("M_lf")) d_M_lf = input_db->getDouble("M_lf");
        if (input_db->keyExists("lambda_lf")) d_lambda_lf = input_db->getDouble("lambda_lf");
        if (input_db->keyExists("eta_lf")) d_eta_lf = input_db->getDouble("eta_lf");
        if (input_db->keyExists("solve_energy")) d_solve_energy = input_db->getBool("solve_energy");

        /*if (input_db->keyExists("k_convective_difference_form"))
            d_k_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("k_convective_difference_form"));
        else if (input_db->keyExists("k_convective_difference_type"))
            d_k_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("k_convective_difference_type"));
        else if (input_db->keyExists("default_k_convective_difference_form"))
            d_k_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("default_k_convective_difference_form"));
        else if (input_db->keyExists("default_k_convective_difference_type"))
            d_k_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("default_k_convective_difference_type"));
        if (input_db->keyExists("k_convective_op_type"))
            d_k_convective_op_type = input_db->getString("k_convective_op_type");
        else if (input_db->keyExists("k_convective_operator_type"))
            d_k_convective_op_type = input_db->getString("k_convective_operator_type");
        else if (input_db->keyExists("default_k_convective_op_type"))
            d_k_convective_op_type = input_db->getString("default_k_convective_op_type");
        else if (input_db->keyExists("default_k_convective_operator_type"))
            d_k_convective_op_type = input_db->getString("default_k_convective_operator_type");

        if (input_db->keyExists("k_convective_op_db"))
            d_k_convective_op_input_db = input_db->getDatabase("k_convective_op_db");
        else if (input_db->keyExists("default_k_convective_op_db"))
            d_k_convective_op_input_db = input_db->getDatabase("default_k_convective_op_db");

        if (input_db->keyExists("w_convective_op_type"))
            d_w_convective_op_type = input_db->getString("w_convective_op_type");
        else if (input_db->keyExists("w_convective_operator_type"))
            d_w_convective_op_type = input_db->getString("w_convective_operator_type");
        else if (input_db->keyExists("default_w_convective_op_type"))
            d_w_convective_op_type = input_db->getString("default_w_convective_op_type");
        else if (input_db->keyExists("default_w_convective_operator_type"))
            d_w_convective_op_type = input_db->getString("default_w_convective_operator_type");

        if (input_db->keyExists("w_convective_op_db"))
            d_w_convective_op_input_db = input_db->getDatabase("w_convective_op_db");
        else if (input_db->keyExists("default_w_convective_op_db"))
            d_w_convective_op_input_db = input_db->getDatabase("default_w_convective_op_db");

        if (input_db->keyExists("w_convective_difference_form"))
            d_w_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("w_convective_difference_form"));
        else if (input_db->keyExists("w_convective_difference_type"))
            d_w_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("w_convective_difference_type"));
        else if (input_db->keyExists("default_w_convective_difference_form"))
            d_w_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("default_w_convective_difference_form"));
        else if (input_db->keyExists("default_w_convective_difference_type"))
            d_w_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("default_w_convective_difference_type"));
        if (input_db->keyExists("k_convective_time_stepping_type"))
            d_k_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("k_convective_time_stepping_type"));
        if (input_db->keyExists("w_convective_time_stepping_type"))
            d_w_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("w_convective_time_stepping_type"));
        if (input_db->keyExists("rho_coarsen_type")) d_rho_coarsen_type = input_db->getString("rho_coarsen_type");
        if (input_db->keyExists("rho_refine_type")) d_rho_refine_type = input_db->getString("rho_refine_type");
        if (input_db->keyExists("rho_bdry_extrap_type"))
            d_rho_bdry_extrap_type = input_db->getString("rho_bdry_extrap_type");*/
    }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
