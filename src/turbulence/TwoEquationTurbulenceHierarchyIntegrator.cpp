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

#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/TwoEquationTurbulenceHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h"

#include "ibtk/CCLaplaceOperator.h"
#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/PoissonSolver.h"

#include "BasePatchHierarchy.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellDataFactory.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideDataFactory.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include <boost/math/tools/roots.hpp>

#include <algorithm>
#include <deque>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

// FORTRAN ROUTINES
#if (NDIM == 2)
#define PRODUCTION IBAMR_FC_FUNC(production_2d, PRODUCTION_2D)
#define SST_BLENDING_FCN IBAMR_FC_FUNC(sst_blending_fcn_2d, SST_BLENDING_FCN_2D)
#define SST_F2_FCN IBAMR_FC_FUNC(sst_f2_fcn_2d, SST_F2_FCN_2D)
#define SST_MU_T_FCN IBAMR_FC_FUNC(sst_mu_t_fcn_2d, SST_MU_T_FCN_2D)
#define WALL_SHEAR_STRESS IBAMR_FC_FUNC(wall_shear_stress_2d, WALL_SHEAR_STRESS_2D)
#endif
#if (NDIM == 3)
#define PRODUCTION IBAMR_FC_FUNC(production_3d, PRODUCTION_3D)
#define SST_BLENDING_FCN IBAMR_FC_FUNC(sst_blending_fcn_3d, SST_BLENDING_FCN_3D)
#define SST_F2_FCN IBAMR_FC_FUNC(sst_f2_fcn_3d, SST_F2_FCN_3D)
#define SST_MU_T_FCN IBAMR_FC_FUNC(sst_mu_t_fcn_3d, SST_MU_T_FCN_3D)
#define WALL_SHEAR_STRESS IBAMR_FC_FUNC(wall_shear_stress_3d, WALL_SHEAR_STRESS_3D)
#endif

extern "C"
{
    void PRODUCTION(const double*,
                    const int&,
                    const double*,
                    const int&,
                    const double*,
                    const double*,
#if (NDIM == 3)
                    const double*,
#endif
                    const int&,
                    const int&,
#if (NDIM == 3)
                    const int&,
#endif
                    const int&,
                    const int&,
                    const int&,
                    const int&,
#if (NDIM == 3)
                    const int&,
                    const int&,
#endif
                    const double*);

    void SST_BLENDING_FCN(const double*,
                          const int&,
                          const double*,
                          const int&,
                          const double*,
                          const int&,
                          const double*,
                          const int&,
                          const double*,
                          const int&,
                          const double*,
                          const int&,
                          const double&,
                          const double&,
                          const int&,
                          const int&,
                          const int&,
                          const int&,
#if (NDIM == 3)
                          const int&,
                          const int&,
#endif
                          const double*);

    void SST_F2_FCN(const double*,
                    const int&,
                    const double*,
                    const int&,
                    const double*,
                    const int&,
                    const double*,
                    const int&,
                    const double*,
                    const int&,
                    const double*,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
                    const int&,
#if (NDIM == 3)

                    const int&,
                    const int&,
#endif
                    const double&);

    void SST_MU_T_FCN(const double*,
                      const int&,
                      const double*,
                      const int&,
                      const double*,
                      const int&,
                      const double*,
                      const int&,
                      const double*,
                      const int&,
                      const double*,
                      const int&,
                      const int&,
                      const int&,
                      const int&,
                      const int&,
#if (NDIM == 3)
                      const int&,
                      const int&,
#endif
                      const double&);

    void WALL_SHEAR_STRESS(const double*,
                           const int&,
                           const double*,
                           const int&,
                           const double*,
                           const int&,
#if (NDIM == 3)
                           const double*,
                           const int&,
#endif
                           const double*,
                           const int&,
                           const double*,
                           const int&,
                           const double*,
                           const int&,
                           const double*,
                           const int&,
                           const double*,
                           const int&,
                           const double*,
                           const int&,
#if (NDIM == 3)
                           const double*,
                           const int&,
#endif
                           const double*,
                           const int&,
                           const double*,
                           const int&,
#if (NDIM == 3)
                           const double*,
                           const int&,
#endif
                           const double&,
                           const double&,
                           const double&,
                           const int&,
                           const int&,
                           const int&,
                           const int&,
                           const int&,
#if (NDIM == 3)
                           const int&,
                           const int&,
#endif
                           const double*);
}
/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Turbulence model constants.

// sigma_k.
static const double SIGMA_K1 = 0.85;
static const double SIGMA_K2 = 1.0;

// sigma_w.
static const double SIGMA_W1 = 0.5;
static const double SIGMA_W2 = 0.856;

// beta.
static const double BETA_1 = 0.075;
static const double BETA_2 = 0.0828;

// beta_star.
static const double BETA_STAR = 0.09;

static const double B = 5.25;

// a1.
static const double A1 = 0.31;

// Constants appear in wall function.
static const double KAPPA = 0.4187;
static const double E = 9.793;
static const double DPLUS = 11.0;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int FACEG = 1;

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

struct FrictionVelocityFunctor
{
    std::pair<double, double> operator()(const double& U_tau);
    double kappa, E, d_plus, delta, nu;
    static double s_newton_guess, s_newton_min, s_newton_max;
    double U_s;
};

double FrictionVelocityFunctor::s_newton_guess = 1.0;
double FrictionVelocityFunctor::s_newton_max = 100.0;
double FrictionVelocityFunctor::s_newton_min = 0.01;

std::pair<double, double>
FrictionVelocityFunctor::operator()(const double& U_tau)
{
    const double c = log(E / kappa) / kappa;
    const double b = 0.5 * (d_plus * kappa / c + 1.0 / d_plus);
    const double kappa_star = kappa * delta / nu;
    const double b_star = b * delta / nu;
    double fx =
        (U_s / U_tau) - (log(1.0 + kappa_star * U_tau) / kappa) -
        c * (1.0 - exp(-delta * U_tau / (nu * d_plus)) - (exp(-b_star * U_tau) * delta * U_tau / (nu * d_plus)));
    double dx =
        -U_s / (U_tau * U_tau) - kappa_star / (kappa * (1.0 + kappa_star * U_tau)) -
        c * delta *
            (exp(-delta * U_tau / (nu * d_plus)) - exp(-b_star * U_tau) + b_star * U_tau * exp(-b_star * U_tau)) /
            (nu * d_plus);
    return std::make_pair(fx, dx);
}
} // namespace
  //
/////////////////////////////// PUBLIC ///////////////////////////////////////
TwoEquationTurbulenceHierarchyIntegrator::TwoEquationTurbulenceHierarchyIntegrator(const std::string& object_name,
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
    // if (from_restart) getFromRestart();

    // Set some default values.
    d_k_convective_op_input_db = new MemoryDatabase(d_object_name + "::k_convective_op_input_db");
    d_w_convective_op_input_db = new MemoryDatabase(d_object_name + "::w_convective_op_input_db");

    d_wall_location_index = input_db->getIntegerArray("wall_location_index");
    d_distance_to_virtual_point = input_db->getDoubleWithDefault("distance_to_virtual_point", 2.0);

    getFromInput(input_db, from_restart);

    // Get plotting options from database
    if (input_db->keyExists("output_k")) d_output_k = input_db->getBool("output_k");
    if (input_db->keyExists("output_w")) d_output_w = input_db->getBool("output_w");

    return;
} // TwoEquationTurbulenceHierarchyIntegrator

void
TwoEquationTurbulenceHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                        Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffSemiImplicitHierarchyIntegrator.
    AdvDiffSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Setup hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<NodeVariable<NDIM, double> > nc_var = new NodeVariable<NDIM, double>("nc_var");
    d_hier_nc_data_ops = hier_ops_manager->getOperationsDouble(nc_var, d_hierarchy, true);

    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Operators and solvers are maintained for each variable registered with the
    // integrator.
    if (d_k_solver_type == CCPoissonSolverManager::UNDEFINED)
    {
        d_k_solver_type = CCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER;
    }
    if (d_k_precond_type == CCPoissonSolverManager::UNDEFINED)
    {
        const int max_levels = gridding_alg->getMaxLevels();
        if (max_levels == 1)
        {
            d_k_precond_type = CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
        }
        else
        {
            d_k_precond_type = CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
        }
        d_k_precond_db->putInteger("max_iterations", 1);
    }

    if (d_w_solver_type == CCPoissonSolverManager::UNDEFINED)
    {
        d_w_solver_type = CCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER;
    }
    if (d_w_precond_type == CCPoissonSolverManager::UNDEFINED)
    {
        const int max_levels = gridding_alg->getMaxLevels();
        if (max_levels == 1)
        {
            d_w_precond_type = CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER;
        }
        else
        {
            d_w_precond_type = CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
        }
        d_w_precond_db->putInteger("max_iterations", 1);
    }

    d_k_solver = getHelmholtzSolverKEquation(d_k_var);
    d_w_solver = getHelmholtzSolverWEquation(d_w_var);
    d_k_rhs_op = getHelmholtzRHSOperatorKEquation(d_k_var);
    d_w_rhs_op = getHelmholtzRHSOperatorWEquation(d_w_var);

    // Register state variables that are maintained by the
    // TwoEquationTurbulenceHierarchyIntegrator.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> face_ghosts = FACEG;
    const IntVector<NDIM> no_ghosts = 0;

    if (d_k_var)
        registerVariable(d_k_current_idx,
                         d_k_new_idx,
                         d_k_scratch_idx,
                         d_k_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_k_init);
    if (d_w_var)
        registerVariable(d_w_current_idx,
                         d_w_new_idx,
                         d_w_scratch_idx,
                         d_w_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE",
                         d_w_init);

    int k_F_current_idx, k_F_scratch_idx, k_F_new_idx;
    if (d_k_F_var)
        registerVariable(k_F_current_idx,
                         k_F_new_idx,
                         k_F_scratch_idx,
                         d_k_F_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int w_F_current_idx, w_F_scratch_idx, w_F_new_idx;
    if (d_w_F_var)
        registerVariable(w_F_current_idx,
                         w_F_new_idx,
                         w_F_scratch_idx,
                         d_w_F_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int k_diff_coef_current_idx, k_diff_coef_scratch_idx, k_diff_coef_new_idx;
    if (d_k_diffusion_coef_var)
        registerVariable(k_diff_coef_current_idx,
                         k_diff_coef_new_idx,
                         k_diff_coef_scratch_idx,
                         d_k_diffusion_coef_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int w_diff_coef_current_idx, w_diff_coef_scratch_idx, w_diff_coef_new_idx;
    if (d_w_diffusion_coef_var)
        registerVariable(w_diff_coef_current_idx,
                         w_diff_coef_new_idx,
                         w_diff_coef_scratch_idx,
                         d_w_diffusion_coef_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int k_u_current_idx, k_u_scratch_idx, k_u_new_idx;
    if (d_k_u_var)
        registerVariable(k_u_current_idx,
                         k_u_new_idx,
                         k_u_scratch_idx,
                         d_k_u_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int w_u_current_idx, w_u_scratch_idx, w_u_new_idx;
    if (d_w_u_var)
        registerVariable(w_u_current_idx,
                         w_u_new_idx,
                         w_u_scratch_idx,
                         d_w_u_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int k_rhs_scratch_idx, w_rhs_scratch_idx, k_diff_coef_rhs_scratch_idx, w_diff_coef_rhs_scratch_idx;

    registerVariable(k_rhs_scratch_idx, d_k_rhs_var, cell_ghosts, getScratchContext());
    registerVariable(w_rhs_scratch_idx, d_w_rhs_var, cell_ghosts, getScratchContext());

    Pointer<CellDataFactory<NDIM, double> > k_factory = d_k_var->getPatchDataFactory();
    const int k_depth = k_factory->getDefaultDepth();
    d_k_diffusion_coef_rhs_var = new SideVariable<NDIM, double>(d_k_var->getName() + "::k_diff_rhs", k_depth);
    registerVariable(k_diff_coef_rhs_scratch_idx, d_k_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    Pointer<CellDataFactory<NDIM, double> > w_factory = d_w_var->getPatchDataFactory();
    const int w_depth = w_factory->getDefaultDepth();
    d_w_diffusion_coef_rhs_var = new SideVariable<NDIM, double>(d_w_var->getName() + "::w_diff_rhs", w_depth);
    registerVariable(w_diff_coef_rhs_scratch_idx, d_w_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    d_k_N_var = new CellVariable<NDIM, double>(d_k_var->getName() + "::N", k_depth);
    int k_N_scratch_idx;
    registerVariable(k_N_scratch_idx, d_k_N_var, cell_ghosts, getScratchContext());

    d_w_N_var = new CellVariable<NDIM, double>(d_w_var->getName() + "::N", w_depth);
    int w_N_scratch_idx;
    registerVariable(w_N_scratch_idx, d_w_N_var, cell_ghosts, getScratchContext());

    d_k_N_old_var = new CellVariable<NDIM, double>(d_k_var->getName() + "::k_N_old", k_depth);
    int k_N_old_current_idx, k_N_old_new_idx, k_N_old_scratch_idx;
    registerVariable(k_N_old_current_idx,
                     k_N_old_new_idx,
                     k_N_old_scratch_idx,
                     d_k_N_old_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_w_N_old_var = new CellVariable<NDIM, double>(d_w_var->getName() + "::w_N_old", w_depth);
    int w_N_old_current_idx, w_N_old_new_idx, w_N_old_scratch_idx;
    registerVariable(w_N_old_current_idx,
                     w_N_old_new_idx,
                     w_N_old_scratch_idx,
                     d_w_N_old_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_F1_var = new CellVariable<NDIM, double>("::F1");
    registerVariable(d_F1_scratch_idx, d_F1_var, no_ghosts, getScratchContext());

    d_F2_var = new CellVariable<NDIM, double>("::F2");
    registerVariable(d_F2_scratch_idx, d_F2_var, no_ghosts, getScratchContext());

    d_p_var = new CellVariable<NDIM, double>("::Production");
    registerVariable(d_p_scratch_idx, d_p_var, no_ghosts, getScratchContext());

    d_distance_to_closest_surface_var = new CellVariable<NDIM, double>("::Distance_to_closest_surface");
    registerVariable(
        d_distance_to_closest_surface_scratch_idx, d_distance_to_closest_surface_var, no_ghosts, getScratchContext());

    d_mu_eff_var = new CellVariable<NDIM, double>("::mu_eff");
    registerVariable(d_mu_eff_scratch_idx, d_mu_eff_var, cell_ghosts, getScratchContext());

    // Initialize and register variables used to compute k and w equations coefficient.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_k_temp_var = new CellVariable<NDIM, double>(d_k_var->getName() + "::temp");
    d_k_temp_idx = var_db->registerVariableAndContext(d_k_temp_var, getCurrentContext(), no_ghosts);

    d_k_temp_rhs_var = new CellVariable<NDIM, double>(d_k_var->getName() + "::temp_rhs");
    d_k_temp_rhs_idx = var_db->registerVariableAndContext(d_k_temp_rhs_var, getCurrentContext(), no_ghosts);

    d_k_dissipation_var = new CellVariable<NDIM, double>(d_k_var->getName() + "::dissipation");
    d_k_dissipation_idx = var_db->registerVariableAndContext(d_k_dissipation_var, getCurrentContext(), no_ghosts);

    d_w_temp_var = new CellVariable<NDIM, double>(d_w_var->getName() + "::temp");
    d_w_temp_idx = var_db->registerVariableAndContext(d_w_temp_var, getCurrentContext(), no_ghosts);

    d_w_temp_rhs_var = new CellVariable<NDIM, double>(d_w_var->getName() + "::temp_rhs");
    d_w_temp_rhs_idx = var_db->registerVariableAndContext(d_w_temp_rhs_var, getCurrentContext(), no_ghosts);

    d_w_dissipation_var = new CellVariable<NDIM, double>(d_w_var->getName() + "::dissipation");
    d_w_dissipation_idx = var_db->registerVariableAndContext(d_w_dissipation_var, getCurrentContext(), no_ghosts);

    d_k_C_var = new CellVariable<NDIM, double>(d_k_var->getName() + "::C");
    d_k_C_idx = var_db->registerVariableAndContext(d_k_C_var, getCurrentContext(), no_ghosts);

    d_w_C_var = new CellVariable<NDIM, double>(d_w_var->getName() + "::C");
    d_w_C_idx = var_db->registerVariableAndContext(d_w_C_var, getCurrentContext(), no_ghosts);

    d_rho_cc_var = new CellVariable<NDIM, double>("rho::cc", /*depth*/ 1);
    registerVariable(d_rho_cc_current_idx,
                     d_rho_cc_new_idx,
                     d_rho_cc_scratch_idx,
                     d_rho_cc_var,
                     cell_ghosts,
                     d_rho_coarsen_type,
                     d_rho_refine_type);

    // Registering a temporary cell-centered vector variable to be used in the interpolation
    // function.
    d_rho_vec_cc_var = new CellVariable<NDIM, double>("rho_vec_cc", NDIM);
    d_rho_vec_cc_idx = var_db->registerVariableAndContext(d_rho_vec_cc_var, getCurrentContext(), cell_ghosts);

    d_yplus_cc_var = new CellVariable<NDIM, double>("y_plus", NDIM);
    d_yplus_cc_idx = var_db->registerVariableAndContext(d_yplus_cc_var, getCurrentContext(), cell_ghosts);

    d_U_tau_cc_var = new CellVariable<NDIM, double>("U_tau", NDIM);
    d_U_tau_cc_idx = var_db->registerVariableAndContext(d_U_tau_cc_var, getCurrentContext(), cell_ghosts);

    // Setup coarsening communications algorithms, used in synchronizing refined
    // regions of coarse data with the underlying fine data.
    Pointer<CoarsenOperator<NDIM> > coarsen_operator_k =
        grid_geom->lookupCoarsenOperator(d_k_var, "CONSERVATIVE_COARSEN");
    getCoarsenAlgorithm(SYNCH_CURRENT_DATA_ALG)->registerCoarsen(d_k_current_idx, d_k_current_idx, coarsen_operator_k);
    getCoarsenAlgorithm(SYNCH_NEW_DATA_ALG)->registerCoarsen(d_k_new_idx, d_k_new_idx, coarsen_operator_k);

    Pointer<CoarsenOperator<NDIM> > coarsen_operator_w =
        grid_geom->lookupCoarsenOperator(d_w_var, "CONSERVATIVE_COARSEN");
    getCoarsenAlgorithm(SYNCH_CURRENT_DATA_ALG)->registerCoarsen(d_w_current_idx, d_w_current_idx, coarsen_operator_w);
    getCoarsenAlgorithm(SYNCH_NEW_DATA_ALG)->registerCoarsen(d_w_new_idx, d_w_new_idx, coarsen_operator_w);

    // Setup the convective operator.
    d_k_convective_op = getConvectiveOperatorKEquation(d_k_var);
    d_w_convective_op = getConvectiveOperatorWEquation(d_w_var);

    // Register variables for plotting.
    if (d_visit_writer)
    {
        if (d_output_k)
        {
            d_visit_writer->registerPlotQuantity("k", "SCALAR", d_k_current_idx, 0, d_k_scale);
        }

        if (d_output_w)
        {
            d_visit_writer->registerPlotQuantity("w", "SCALAR", d_w_current_idx, 0, d_w_scale);
        }
    }

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}

void
TwoEquationTurbulenceHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
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
        d_k_solver_needs_init = true;
        d_k_rhs_op_needs_init = true;
        d_w_solver_needs_init = true;
        d_w_rhs_op_needs_init = true;
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }
    // Allocate the scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_k_temp_idx)) level->allocatePatchData(d_k_temp_idx, current_time);
        if (!level->checkAllocated(d_k_temp_rhs_idx)) level->allocatePatchData(d_k_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_k_dissipation_idx)) level->allocatePatchData(d_k_dissipation_idx, current_time);
        if (!level->checkAllocated(d_w_temp_idx)) level->allocatePatchData(d_w_temp_idx, current_time);
        if (!level->checkAllocated(d_w_temp_rhs_idx)) level->allocatePatchData(d_w_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_w_dissipation_idx)) level->allocatePatchData(d_w_dissipation_idx, current_time);
        if (!level->checkAllocated(d_k_C_idx)) level->allocatePatchData(d_k_C_idx, current_time);
        if (!level->checkAllocated(d_w_C_idx)) level->allocatePatchData(d_w_C_idx, current_time);
        if (!level->checkAllocated(d_rho_vec_cc_idx)) level->allocatePatchData(d_rho_vec_cc_idx, current_time);
        if (!level->checkAllocated(d_yplus_cc_idx)) level->allocatePatchData(d_yplus_cc_idx, current_time);
        if (!level->checkAllocated(d_U_tau_cc_idx)) level->allocatePatchData(d_U_tau_cc_idx, current_time);
    }

    // Update the advection velocity for k.
    const int k_u_current_idx = var_db->mapVariableAndContextToIndex(d_k_u_var, getCurrentContext());
    const int k_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_u_var, getScratchContext());
    const int k_u_new_idx = var_db->mapVariableAndContextToIndex(d_k_u_var, getNewContext());

    if (d_k_u_fcn)
    {
        d_k_u_fcn->setDataOnPatchHierarchy(k_u_current_idx, d_k_u_var, d_hierarchy, current_time);
        d_k_u_fcn->setDataOnPatchHierarchy(k_u_new_idx, d_k_u_var, d_hierarchy, new_time);
    }
    else
    {
        d_hier_fc_data_ops->copyData(k_u_new_idx, k_u_current_idx);
    }
    d_hier_fc_data_ops->linearSum(k_u_scratch_idx, 0.5, k_u_current_idx, 0.5, k_u_new_idx);

    // Setup the operators and solvers and compute the right-hand-side terms.
    const int k_current_idx = var_db->mapVariableAndContextToIndex(d_k_var, getCurrentContext());
    const int k_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_var, getScratchContext());
    const int k_new_idx = var_db->mapVariableAndContextToIndex(d_k_var, getNewContext());
    const int k_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_rhs_var, getScratchContext());
    const int k_diff_coef_current_idx =
        (var_db->mapVariableAndContextToIndex(d_k_diffusion_coef_var, getCurrentContext()));
    const int k_diff_coef_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_k_diffusion_coef_var, getScratchContext()));
    const int k_diff_coef_rhs_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_k_diffusion_coef_rhs_var, getScratchContext()));

    // setting k equation diffusion timestepping type
    double alpha = 0.0;
    switch (d_k_diffusion_time_stepping_type)
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
                                 << enum_to_string<TimeSteppingType>(d_k_diffusion_time_stepping_type) << " \n"
                                 << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }
    PoissonSpecifications k_solver_spec(d_object_name + "::solver_spec::" + d_k_var->getName());
    PoissonSpecifications k_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_k_var->getName());

    // Interpolate the side-centered density maintained at INSVCStaggeredHierarchyIntegrator
    // to cell-centered
    interpolateSCMassDensityToCC(getCurrentContext());
    const int rho_cc_current_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getCurrentContext());
    const int rho_cc_new_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getNewContext());

    // set rho/dt
    d_hier_cc_data_ops->scale(d_k_temp_rhs_idx, 1.0 / dt, rho_cc_current_idx);
    k_rhs_op_spec.setCPatchDataId(d_k_temp_rhs_idx);

    // set rho*beta_star*omega
    d_w_current_idx = var_db->mapVariableAndContextToIndex(d_w_var, getCurrentContext());
    d_hier_cc_data_ops->multiply(d_k_dissipation_idx, rho_cc_current_idx, d_w_current_idx);
    d_hier_cc_data_ops->scale(d_k_C_idx, BETA_STAR, d_k_dissipation_idx);
    d_hier_cc_data_ops->add(d_k_C_idx, d_k_C_idx, d_k_temp_rhs_idx);
    k_solver_spec.setCPatchDataId(d_k_C_idx);

    d_mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
    const int mu_current_idx =
        var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getCurrentContext());
    const int mu_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getScratchContext());
    d_hier_cc_data_ops->copyData(mu_scratch_idx, mu_current_idx);
    d_mu_bdry_bc_fill_op->fillData(current_time);

    // Calculate F1 and blend the coefficients
    calculateBlendingFunction(current_time, getCurrentContext());
    d_mu_t_var = d_ins_hierarchy_integrator->getTurbulentViscosityVariable();
    const int mu_t_current_idx =
        var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getCurrentContext());
    const int mu_t_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getScratchContext());
    d_hier_cc_data_ops->copyData(mu_t_scratch_idx, mu_t_current_idx);
    d_mu_t_bdry_bc_fill_op->fillData(current_time);

    d_mu_eff_scratch_idx = var_db->mapVariableAndContextToIndex(d_mu_eff_var, getScratchContext());
    d_hier_cc_data_ops->axpy(
        d_mu_eff_scratch_idx, d_sigma_k, mu_t_scratch_idx, mu_scratch_idx, false /*interior_only*/);

    // Interpolate the cell-centered mu_eff to side-centered
    d_hier_math_ops->interp(k_diff_coef_current_idx,
                            d_k_diffusion_coef_var,
                            true,
                            d_mu_eff_scratch_idx,
                            d_mu_eff_var,
                            d_no_fill_op,
                            d_integrator_time);

    d_hier_sc_data_ops->scale(k_diff_coef_rhs_scratch_idx, (1.0 - alpha), k_diff_coef_current_idx);
    k_rhs_op_spec.setDPatchDataId(k_diff_coef_rhs_scratch_idx);

    d_hier_sc_data_ops->scale(k_diff_coef_scratch_idx, -alpha, k_diff_coef_current_idx);
    k_solver_spec.setDPatchDataId(k_diff_coef_scratch_idx);

    // Initialize the RHS operator and compute the RHS vector for k equation.
    Pointer<LaplaceOperator> k_rhs_op = d_k_rhs_op;
    k_rhs_op->setPoissonSpecifications(k_rhs_op_spec);
    k_rhs_op->setPhysicalBcCoef(d_k_bc_coef);
    k_rhs_op->setHomogeneousBc(false);
    k_rhs_op->setSolutionTime(current_time);
    k_rhs_op->setTimeInterval(current_time, new_time);
    if (d_k_rhs_op_needs_init)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the RHS operator for" << d_k_var->getName() << "\n";
        }
        k_rhs_op->initializeOperatorState(*d_k_sol, *d_k_rhs);
        d_k_rhs_op_needs_init = false;
    }
    d_hier_cc_data_ops->copyData(k_scratch_idx, k_current_idx, false);
    k_rhs_op->apply(*d_k_sol, *d_k_rhs);

    // Initialize the linear solver for k equation.
    Pointer<PoissonSolver> k_solver = d_k_solver;
    k_solver->setPoissonSpecifications(k_solver_spec);
    k_solver->setPhysicalBcCoef(d_k_bc_coef);
    k_solver->setHomogeneousBc(false);
    k_solver->setSolutionTime(new_time);
    k_solver->setTimeInterval(current_time, new_time);
    if (d_k_solver_needs_init)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the solvers for" << d_k_var->getName() << "\n";
        }
        k_solver->initializeSolverState(*d_k_sol, *d_k_rhs);
        d_k_solver_needs_init = false;
    }

    // Account for the convective difference term.
    if (d_k_u_var)
    {
        if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(d_k_convective_time_stepping_type))
        {
            d_k_convective_time_stepping_type = d_k_init_convective_time_stepping_type;
        }
        if ((num_cycles == 1) && (d_k_convective_time_stepping_type == MIDPOINT_RULE ||
                                  d_k_convective_time_stepping_type == TRAPEZOIDAL_RULE))
        {
            TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                     << "  time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_k_convective_time_stepping_type)
                                     << " requires num_cycles > 1.\n"
                                     << "  at current time step, num_cycles = " << num_cycles << "\n");
        }
        if (d_k_convective_op_needs_init)
        {
            d_k_convective_op->initializeOperatorState(*d_k_sol, *d_k_rhs);
            d_k_convective_op_needs_init = false;
        }

        d_k_convective_op->setAdvectionVelocity(k_u_current_idx);
        const int k_current_idx = var_db->mapVariableAndContextToIndex(d_k_var, getCurrentContext());
        const int k_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_var, getScratchContext());
        const int k_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_N_var, getScratchContext());
        d_hier_cc_data_ops->copyData(k_scratch_idx, k_current_idx);
        d_k_convective_op->setSolutionTime(current_time);
        d_k_convective_op->applyConvectiveOperator(k_scratch_idx, k_N_scratch_idx);
        const int k_N_old_new_idx = var_db->mapVariableAndContextToIndex(d_k_N_old_var, getNewContext());
        d_hier_cc_data_ops->copyData(k_N_old_new_idx, k_N_scratch_idx); // this is redundant
        if (d_k_convective_time_stepping_type == FORWARD_EULER)
        {
            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, -1.0, k_N_scratch_idx, k_rhs_scratch_idx);
        }
        else if (d_k_convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, -0.5, k_N_scratch_idx, k_rhs_scratch_idx);
        }
    }

    // Update the advection velocity for w.
    const int w_u_current_idx = var_db->mapVariableAndContextToIndex(d_w_u_var, getCurrentContext());
    const int w_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_u_var, getScratchContext());
    const int w_u_new_idx = var_db->mapVariableAndContextToIndex(d_w_u_var, getNewContext());

    if (d_w_u_fcn)
    {
        d_w_u_fcn->setDataOnPatchHierarchy(w_u_current_idx, d_w_u_var, d_hierarchy, current_time);
        d_w_u_fcn->setDataOnPatchHierarchy(w_u_new_idx, d_w_u_var, d_hierarchy, new_time);
    }
    else
    {
        d_hier_fc_data_ops->copyData(w_u_new_idx, w_u_current_idx);
    }
    d_hier_fc_data_ops->linearSum(w_u_scratch_idx, 0.5, w_u_current_idx, 0.5, w_u_new_idx);

    const int w_current_idx = var_db->mapVariableAndContextToIndex(d_w_var, getCurrentContext());
    const int w_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_var, getScratchContext());
    const int w_new_idx = var_db->mapVariableAndContextToIndex(d_w_var, getNewContext());
    const int w_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_rhs_var, getScratchContext());
    const int w_diff_coef_current_idx =
        (var_db->mapVariableAndContextToIndex(d_w_diffusion_coef_var, getCurrentContext()));
    const int w_diff_coef_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_w_diffusion_coef_var, getScratchContext()));
    const int w_diff_coef_rhs_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_w_diffusion_coef_rhs_var, getScratchContext()));

    // setting w equation diffusion timestepping type.
    switch (d_w_diffusion_time_stepping_type)
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
                                 << enum_to_string<TimeSteppingType>(d_w_diffusion_time_stepping_type) << " \n"
                                 << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }
    PoissonSpecifications w_solver_spec(d_object_name + "::solver_spec::" + d_w_var->getName());
    PoissonSpecifications w_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_w_var->getName());

    // set rho/dt.
    d_hier_cc_data_ops->scale(d_w_temp_rhs_idx, 1.0 / dt, rho_cc_current_idx);
    w_rhs_op_spec.setCPatchDataId(d_w_temp_rhs_idx);

    // set rho*beta*omega
    d_hier_cc_data_ops->multiply(d_w_dissipation_idx, rho_cc_current_idx, w_current_idx);
    d_hier_cc_data_ops->axpy(d_w_C_idx, d_beta, d_w_dissipation_idx, d_w_temp_rhs_idx);
    w_solver_spec.setCPatchDataId(d_w_C_idx);

    d_hier_cc_data_ops->axpy(
        d_mu_eff_scratch_idx, d_sigma_w, mu_t_scratch_idx, mu_scratch_idx, false /*interior_only*/);

    // Interpolate the cell-centered mu_eff to side-centered
    d_hier_math_ops->interp(w_diff_coef_current_idx,
                            d_w_diffusion_coef_var,
                            true,
                            d_mu_eff_scratch_idx,
                            d_mu_eff_var,
                            d_no_fill_op,
                            d_integrator_time);

    d_hier_sc_data_ops->scale(w_diff_coef_rhs_scratch_idx, (1.0 - alpha), w_diff_coef_current_idx);
    // d_hier_sc_data_ops->setToScalar(w_diff_coef_rhs_scratch_idx, 0.0);
    w_rhs_op_spec.setDPatchDataId(w_diff_coef_rhs_scratch_idx);

    d_hier_sc_data_ops->scale(w_diff_coef_scratch_idx, -alpha, w_diff_coef_current_idx);
    // d_hier_sc_data_ops->setToScalar(w_diff_coef_scratch_idx, 0.0);
    w_solver_spec.setDPatchDataId(w_diff_coef_scratch_idx);

    // Initialize the RHS operator and compute the RHS vector for w equation.
    Pointer<LaplaceOperator> w_rhs_op = d_w_rhs_op;
    w_rhs_op->setPoissonSpecifications(w_rhs_op_spec);
    w_rhs_op->setPhysicalBcCoef(d_w_bc_coef);
    w_rhs_op->setHomogeneousBc(false);
    w_rhs_op->setSolutionTime(current_time);
    w_rhs_op->setTimeInterval(current_time, new_time);
    if (d_w_rhs_op_needs_init)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the RHS operator for" << d_w_var->getName() << "\n";
        }
        w_rhs_op->initializeOperatorState(*d_w_sol, *d_w_rhs);
        d_w_rhs_op_needs_init = false;
    }
    d_hier_cc_data_ops->copyData(w_scratch_idx, w_current_idx, false);
    w_rhs_op->apply(*d_w_sol, *d_w_rhs);

    // Initialize the linear solver for w equation.
    Pointer<PoissonSolver> w_solver = d_w_solver;
    w_solver->setPoissonSpecifications(w_solver_spec);
    w_solver->setPhysicalBcCoef(d_w_bc_coef);
    w_solver->setHomogeneousBc(false);
    w_solver->setSolutionTime(new_time);
    w_solver->setTimeInterval(current_time, new_time);
    if (d_w_solver_needs_init)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the solvers for" << d_w_var->getName() << "\n";
        }
        w_solver->initializeSolverState(*d_w_sol, *d_w_rhs);
        d_w_solver_needs_init = false;
    }

    // Account for the convective difference term.
    if (d_w_u_var)
    {
        if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(d_w_convective_time_stepping_type))
        {
            d_w_convective_time_stepping_type = d_w_init_convective_time_stepping_type;
        }
        if ((num_cycles == 1) && (d_w_convective_time_stepping_type == MIDPOINT_RULE ||
                                  d_w_convective_time_stepping_type == TRAPEZOIDAL_RULE))
        {
            TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                     << "  time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_w_convective_time_stepping_type)
                                     << " requires num_cycles > 1.\n"
                                     << "  at current time step, num_cycles = " << num_cycles << "\n");
        }
        if (d_w_convective_op_needs_init)
        {
            d_w_convective_op->initializeOperatorState(*d_w_sol, *d_w_rhs);
            d_w_convective_op_needs_init = false;
        }

        d_w_convective_op->setAdvectionVelocity(w_u_current_idx);
        const int w_current_idx = var_db->mapVariableAndContextToIndex(d_w_var, getCurrentContext());
        const int w_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_var, getScratchContext());
        const int w_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_N_var, getScratchContext());
        d_hier_cc_data_ops->copyData(w_scratch_idx, w_current_idx);
        d_w_convective_op->setSolutionTime(current_time);
        d_w_convective_op->applyConvectiveOperator(w_scratch_idx, w_N_scratch_idx);
        const int w_N_old_new_idx = var_db->mapVariableAndContextToIndex(d_w_N_old_var, getNewContext());
        d_hier_cc_data_ops->copyData(w_N_old_new_idx, w_N_scratch_idx); // redundant
        if (d_w_convective_time_stepping_type == FORWARD_EULER)
        {
            d_hier_cc_data_ops->axpy(w_rhs_scratch_idx, -1.0, w_N_scratch_idx, w_rhs_scratch_idx);
        }
        else if (d_w_convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(w_rhs_scratch_idx, -0.5, w_N_scratch_idx, w_rhs_scratch_idx);
        }
    }

    computeWallShearStressFromWallLaw(current_time, getCurrentContext());

    // Set the initial guess for k, w and rho.
    d_hier_cc_data_ops->copyData(k_new_idx, k_current_idx);
    d_hier_cc_data_ops->copyData(w_new_idx, w_current_idx);
    d_hier_cc_data_ops->copyData(rho_cc_new_idx, rho_cc_current_idx);

    // compute wall shear stress from wall law.

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
TwoEquationTurbulenceHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                             const double new_time,
                                                             const int cycle_num)
{
    AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int ins_expected_num_cycles = d_ins_hierarchy_integrator->getNumberOfCycles();
    if ((d_current_num_cycles != ins_expected_num_cycles))
    {
        IBAMR_DO_ONCE({
            pout << "TwoEquationTurbulenceHierarchyIntegrator::integrateHierarchy():\n"
                 << "  WARNING: adv diff num_cycles = " << d_current_num_cycles
                 << " but the ins num_cycles is = " << ins_expected_num_cycles << ".\n"
                 << "The current implementation requires both turbulence hierarchy and ins hierarchy integrators have "
                    "the same number of cycles";
        });
    }

    // Setup the operators and solvers and compute the right-hand-side terms.
    const int k_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_var, getScratchContext());
    const int k_new_idx = var_db->mapVariableAndContextToIndex(d_k_var, getNewContext());
    const int k_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_rhs_var, getScratchContext());
    const int k_diff_coef_new_idx = var_db->mapVariableAndContextToIndex(d_k_diffusion_coef_var, getNewContext());
    const int k_diff_coef_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_k_diffusion_coef_var, getScratchContext());
    const int w_new_idx = var_db->mapVariableAndContextToIndex(d_w_var, getNewContext());

    if (cycle_num > 0) interpolateSCMassDensityToCC(getNewContext());
    const int rho_cc_new_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getNewContext());
    d_mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
    const int mu_new_idx = var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getNewContext());
    const int mu_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getScratchContext());
    d_mu_t_var = d_ins_hierarchy_integrator->getTurbulentViscosityVariable();
    const int mu_t_new_idx =
        var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getNewContext());
    const int mu_t_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getScratchContext());
    d_mu_eff_scratch_idx = var_db->mapVariableAndContextToIndex(d_mu_eff_var, getScratchContext());
    if (cycle_num > 0)
    {
        // set rho/dt.
        d_hier_cc_data_ops->scale(d_k_temp_idx, 1.0 / dt, rho_cc_new_idx);

        // set rho*beta_star*omega.
        d_hier_cc_data_ops->multiply(d_k_dissipation_idx, rho_cc_new_idx, w_new_idx);
        d_hier_cc_data_ops->scale(d_k_C_idx, BETA_STAR, d_k_dissipation_idx);
        d_hier_cc_data_ops->add(d_k_C_idx, d_k_C_idx, d_k_temp_idx);

        double alpha = 0.0;
        switch (d_k_diffusion_time_stepping_type)
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
                                     << enum_to_string<TimeSteppingType>(d_k_diffusion_time_stepping_type) << " \n"
                                     << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
        }

        // Calculate F1 and blend the coefficients.
        calculateBlendingFunction(new_time, getNewContext());

        // Update the diffusion coefficient.
        d_hier_cc_data_ops->copyData(mu_t_scratch_idx, mu_t_new_idx);
        d_hier_cc_data_ops->copyData(mu_scratch_idx, mu_new_idx);
        d_mu_t_bdry_bc_fill_op->fillData(new_time);
        d_mu_bdry_bc_fill_op->fillData(new_time);
        d_hier_cc_data_ops->axpy(
            d_mu_eff_scratch_idx, d_sigma_k, mu_t_scratch_idx, mu_scratch_idx, /*interior_only*/ false);

        // Interpolate the cell-centered mu_eff to side-centered
        d_hier_math_ops->interp(k_diff_coef_new_idx,
                                d_k_diffusion_coef_var,
                                true,
                                d_mu_eff_scratch_idx,
                                d_mu_eff_var,
                                d_no_fill_op,
                                d_integrator_time);

        d_hier_sc_data_ops->scale(k_diff_coef_scratch_idx, -alpha, k_diff_coef_new_idx);

        // Update the advection velocity for k.
        const int k_u_current_idx = var_db->mapVariableAndContextToIndex(d_k_u_var, getCurrentContext());
        const int k_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_u_var, getScratchContext());
        const int k_u_new_idx = var_db->mapVariableAndContextToIndex(d_k_u_var, getNewContext());
        if (d_k_u_fcn)
        {
            d_k_u_fcn->setDataOnPatchHierarchy(k_u_new_idx, d_k_u_var, d_hierarchy, new_time);
        }
        d_hier_fc_data_ops->linearSum(k_u_scratch_idx, 0.5, k_u_current_idx, 0.5, k_u_new_idx);
    }

    // account for the convective term
    TimeSteppingType convective_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
    if (d_k_u_var)
    {
        convective_time_stepping_type = d_k_convective_time_stepping_type;
        if (is_multistep_time_stepping_type(convective_time_stepping_type))
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(convective_time_stepping_type == ADAMS_BASHFORTH);
#endif
            if (getIntegratorStep() == 0)
            {
                convective_time_stepping_type = d_k_init_convective_time_stepping_type;
            }
            else if (cycle_num > 0)
            {
                convective_time_stepping_type = MIDPOINT_RULE;
                IBAMR_DO_ONCE({
                    pout << "TwoEquationTurbulenceHierarchyIntegrator::"
                            "integrateHierarchy():"
                            "\n"
                         << "  WARNING: convective_time_stepping_type = "
                         << enum_to_string<TimeSteppingType>(d_k_convective_time_stepping_type)
                         << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                         << "           using " << enum_to_string<TimeSteppingType>(d_k_convective_time_stepping_type)
                         << " only for the first cycle in each time step;\n"
                         << "           using " << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                         << " for subsequent cycles.\n";
                });
            }
        }
        const int k_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_N_var, getScratchContext());
        if (cycle_num > 0)
        {
            if (convective_time_stepping_type == MIDPOINT_RULE)
            {
                const int k_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_u_var, getScratchContext());
                d_k_convective_op->setAdvectionVelocity(k_u_scratch_idx);
                const int k_current_idx = var_db->mapVariableAndContextToIndex(d_k_var, getCurrentContext());
                const int k_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_var, getScratchContext());
                const int k_new_idx = var_db->mapVariableAndContextToIndex(d_k_var, getNewContext());
                d_hier_cc_data_ops->linearSum(k_scratch_idx, 0.5, k_current_idx, 0.5, k_new_idx);
                d_k_convective_op->setSolutionTime(half_time);
                d_k_convective_op->applyConvectiveOperator(k_scratch_idx, k_N_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                const int k_u_new_idx = var_db->mapVariableAndContextToIndex(d_k_u_var, getNewContext());
                d_k_convective_op->setAdvectionVelocity(k_u_new_idx);
                const int k_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_var, getScratchContext());
                const int k_new_idx = var_db->mapVariableAndContextToIndex(d_k_var, getNewContext());
                d_hier_cc_data_ops->copyData(k_scratch_idx, k_new_idx);
                d_k_convective_op->setSolutionTime(new_time);
                d_k_convective_op->applyConvectiveOperator(k_scratch_idx, k_N_scratch_idx);
            }
        }
        if (convective_time_stepping_type == ADAMS_BASHFORTH)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(cycle_num == 0);
#endif
            const int k_N_old_current_idx = var_db->mapVariableAndContextToIndex(d_k_N_old_var, getCurrentContext());
            const double omega = dt / d_dt_previous[0];
            d_hier_cc_data_ops->linearSum(
                k_N_scratch_idx, 1.0 + 0.5 * omega, k_N_scratch_idx, -0.5 * omega, k_N_old_current_idx);
        }
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, -1.0, k_N_scratch_idx, k_rhs_scratch_idx);
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, -0.5, k_N_scratch_idx, k_rhs_scratch_idx);
        }
    }
    // Compute the production term.
    calculateTurbulentKEProduction(new_time);
    computeWallShearStressFromWallLaw(new_time, getNewContext());
    // Account for forcing terms.
    const int k_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_F_var, getScratchContext());
    const int k_F_new_idx = var_db->mapVariableAndContextToIndex(d_k_F_var, getNewContext());

    if (d_k_F_var)
    {
        d_k_F_fcn->setDataOnPatchHierarchy(k_F_scratch_idx, d_k_F_var, d_hierarchy, half_time);
        d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, +1.0, k_F_scratch_idx, k_rhs_scratch_idx);
    }

    // Solve for k(n+1).
    Pointer<PoissonSolver> k_solver = d_k_solver;
    k_solver->solveSystem(*d_k_sol, *d_k_rhs);
    d_hier_cc_data_ops->copyData(k_new_idx, k_scratch_idx);
    if (d_enable_logging && d_enable_logging_solver_iterations)
        plog << d_object_name << ":" << d_k_var->getName()
             << "::integrateHierarchy():diffusion solve number of iterations = " << k_solver->getNumIterations()
             << "\n";
    if (d_enable_logging)
        plog << d_object_name << ":" << d_k_var->getName()
             << "::integrateHierarchy():diffusion solve residual norm        = " << k_solver->getResidualNorm() << "\n";
    if (k_solver->getNumIterations() == k_solver->getMaxIterations())
    {
        pout << d_object_name << ":" << d_k_var->getName()
             << "::integrateHierarchy():WARNING: linear solver iterations == max iterations\n";
    }
    // Reset the right-hand side vector.
    if (d_k_u_var)
    {
        const int k_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_N_var, getScratchContext());
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, +1.0, k_N_scratch_idx, k_rhs_scratch_idx);
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, +0.5, k_N_scratch_idx, k_rhs_scratch_idx);
        }
    }
    if (d_k_F_fcn)
    {
        d_hier_cc_data_ops->axpy(k_rhs_scratch_idx, -1.0, k_F_scratch_idx, k_rhs_scratch_idx);
        d_hier_cc_data_ops->copyData(k_F_new_idx, k_F_scratch_idx);
    }

    const int w_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_var, getScratchContext());
    const int w_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_rhs_var, getScratchContext());
    const int w_diff_coef_new_idx = var_db->mapVariableAndContextToIndex(d_w_diffusion_coef_var, getNewContext());
    const int w_diff_coef_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_w_diffusion_coef_var, getScratchContext());

    if (cycle_num > 0)
    {
        // set rho/dt.
        d_hier_cc_data_ops->scale(d_w_temp_idx, 1.0 / dt, rho_cc_new_idx);

        // set rho*beta*omega.
        d_hier_cc_data_ops->multiply(d_w_dissipation_idx, rho_cc_new_idx, w_new_idx);
        d_hier_cc_data_ops->axpy(d_w_C_idx, d_beta, d_w_dissipation_idx, d_w_temp_idx);

        double alpha = 0.0;
        switch (d_w_diffusion_time_stepping_type)
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
                                     << enum_to_string<TimeSteppingType>(d_w_diffusion_time_stepping_type) << " \n"
                                     << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
        }

        // Update the diffusion coefficient.
        d_hier_cc_data_ops->axpy(
            d_mu_eff_scratch_idx, d_sigma_w, mu_t_scratch_idx, mu_scratch_idx, /*interior_only*/ false);

        // Interpolate the cell-centered mu_eff to side-centered.
        d_hier_math_ops->interp(w_diff_coef_new_idx,
                                d_w_diffusion_coef_var,
                                true,
                                d_mu_eff_scratch_idx,
                                d_mu_eff_var,
                                d_no_fill_op,
                                d_integrator_time);

        d_hier_sc_data_ops->scale(w_diff_coef_scratch_idx, -alpha, w_diff_coef_new_idx);
        // d_hier_sc_data_ops->setToScalar(w_diff_coef_scratch_idx, 0.0);

        // Update the advection velocity for w.
        const int w_u_current_idx = var_db->mapVariableAndContextToIndex(d_w_u_var, getCurrentContext());
        const int w_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_u_var, getScratchContext());
        const int w_u_new_idx = var_db->mapVariableAndContextToIndex(d_w_u_var, getNewContext());
        if (d_w_u_fcn)
        {
            d_w_u_fcn->setDataOnPatchHierarchy(w_u_new_idx, d_w_u_var, d_hierarchy, new_time);
        }
        d_hier_fc_data_ops->linearSum(w_u_scratch_idx, 0.5, w_u_current_idx, 0.5, w_u_new_idx);
    }
    // account for the convective term
    if (d_w_u_var)
    {
        convective_time_stepping_type = d_w_convective_time_stepping_type;
        if (is_multistep_time_stepping_type(convective_time_stepping_type))
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(convective_time_stepping_type == ADAMS_BASHFORTH);
#endif
            if (getIntegratorStep() == 0)
            {
                convective_time_stepping_type = d_w_init_convective_time_stepping_type;
            }
            else if (cycle_num > 0)
            {
                convective_time_stepping_type = MIDPOINT_RULE;
                IBAMR_DO_ONCE({
                    pout << "TwoEquationTurbulenceHierarchyIntegrator::"
                            "integrateHierarchy():"
                            "\n"
                         << "  WARNING: convective_time_stepping_type = "
                         << enum_to_string<TimeSteppingType>(d_w_convective_time_stepping_type)
                         << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                         << "           using " << enum_to_string<TimeSteppingType>(d_w_convective_time_stepping_type)
                         << " only for the first cycle in each time step;\n"
                         << "           using " << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                         << " for subsequent cycles.\n";
                });
            }
        }
        const int w_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_N_var, getScratchContext());
        if (cycle_num > 0)
        {
            if (convective_time_stepping_type == MIDPOINT_RULE)
            {
                const int w_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_u_var, getScratchContext());
                d_w_convective_op->setAdvectionVelocity(w_u_scratch_idx);
                const int w_current_idx = var_db->mapVariableAndContextToIndex(d_w_var, getCurrentContext());
                const int w_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_var, getScratchContext());
                const int w_new_idx = var_db->mapVariableAndContextToIndex(d_w_var, getNewContext());
                d_hier_cc_data_ops->linearSum(w_scratch_idx, 0.5, w_current_idx, 0.5, w_new_idx);
                d_w_convective_op->setSolutionTime(half_time);
                d_w_convective_op->applyConvectiveOperator(w_scratch_idx, w_N_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                const int w_u_new_idx = var_db->mapVariableAndContextToIndex(d_w_u_var, getNewContext());
                d_w_convective_op->setAdvectionVelocity(w_u_new_idx);
                const int w_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_var, getScratchContext());
                const int w_new_idx = var_db->mapVariableAndContextToIndex(d_w_var, getNewContext());
                d_hier_cc_data_ops->copyData(w_scratch_idx, w_new_idx);
                d_w_convective_op->setSolutionTime(new_time);
                d_w_convective_op->applyConvectiveOperator(w_scratch_idx, w_N_scratch_idx);
            }
        }
        if (convective_time_stepping_type == ADAMS_BASHFORTH)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(cycle_num == 0);
#endif
            const int w_N_old_current_idx = var_db->mapVariableAndContextToIndex(d_w_N_old_var, getCurrentContext());
            const double omega = dt / d_dt_previous[0];
            d_hier_cc_data_ops->linearSum(
                w_N_scratch_idx, 1.0 + 0.5 * omega, w_N_scratch_idx, -0.5 * omega, w_N_old_current_idx);
        }
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(w_rhs_scratch_idx, -1.0, w_N_scratch_idx, w_rhs_scratch_idx);
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(w_rhs_scratch_idx, -0.5, w_N_scratch_idx, w_rhs_scratch_idx);
        }
    }

    // Account for forcing terms.
    const int w_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_F_var, getScratchContext());
    const int w_F_new_idx = var_db->mapVariableAndContextToIndex(d_w_F_var, getNewContext());
    if (d_w_F_var)
    {
        d_w_F_fcn->setDataOnPatchHierarchy(w_F_scratch_idx, d_w_F_var, d_hierarchy, half_time);
        d_hier_cc_data_ops->axpy(w_rhs_scratch_idx, +1.0, w_F_scratch_idx, w_rhs_scratch_idx);
    }

    // Solve for w(n+1).
    Pointer<PoissonSolver> w_solver = d_w_solver;
    w_solver->solveSystem(*d_w_sol, *d_w_rhs);
    d_hier_cc_data_ops->copyData(w_new_idx, w_scratch_idx);

    if (d_enable_logging && d_enable_logging_solver_iterations)
        plog << d_object_name << ":" << d_w_var->getName()
             << "::integrateHierarchy(): diffusion solve number of iterations = " << w_solver->getNumIterations()
             << "\n";
    if (d_enable_logging)
        plog << d_object_name << ":" << d_w_var->getName()
             << "::integrateHierarchy(): diffusion solve residual norm        = " << w_solver->getResidualNorm()
             << "\n";
    if (w_solver->getNumIterations() == w_solver->getMaxIterations())
    {
        pout << d_object_name << ":" << d_w_var->getName() << "::integrateHierarchy():"
             << "  WARNING: linear solver iterations == max iterations\n";
    }
    // Reset the right-hand side vector.
    if (d_w_u_var)
    {
        const int w_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_N_var, getScratchContext());
        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(w_rhs_scratch_idx, +1.0, w_N_scratch_idx, w_rhs_scratch_idx);
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(w_rhs_scratch_idx, +0.5, w_N_scratch_idx, w_rhs_scratch_idx);
        }
    }
    if (d_w_F_fcn)
    {
        d_hier_cc_data_ops->axpy(w_rhs_scratch_idx, -1.0, w_F_scratch_idx, w_rhs_scratch_idx);
        d_hier_cc_data_ops->copyData(w_F_new_idx, w_F_scratch_idx);
    }

    // Find the U_tau and yplus from the wall law.
    // applyWallFunction(new_time);

    Pointer<SideVariable<NDIM, double> > yplus_sc_var = d_ins_hierarchy_integrator->getYplusVariable();
    Pointer<SideVariable<NDIM, double> > U_tau_sc_var = d_ins_hierarchy_integrator->getUtauVariable();
    const int yplus_sc_idx =
        var_db->mapVariableAndContextToIndex(yplus_sc_var, d_ins_hierarchy_integrator->getScratchContext());
    const int U_tau_sc_idx =
        var_db->mapVariableAndContextToIndex(U_tau_sc_var, d_ins_hierarchy_integrator->getScratchContext());

    const int yplus_cc_idx = var_db->mapVariableAndContextToIndex(d_yplus_cc_var, getCurrentContext());
    const int U_tau_cc_idx = var_db->mapVariableAndContextToIndex(d_U_tau_cc_var, getCurrentContext());
    static const bool synch_cf_interface = true;

    // Interpolate the side-centered yplus to cell.
    d_hier_math_ops->interp(
        yplus_cc_idx, d_yplus_cc_var, yplus_sc_idx, yplus_sc_var, d_no_fill_op, d_integrator_time, synch_cf_interface);

    // Interpolate the side-centered U_tau to cell.
    d_hier_math_ops->interp(
        U_tau_cc_idx, d_U_tau_cc_var, U_tau_sc_idx, U_tau_sc_var, d_no_fill_op, d_integrator_time, synch_cf_interface);

    // Based on the yplus, set the wall boundary conditions for k and w variables.
    // postprocessTurbulentVariablesBasedonYplus();

    // Postprocess w values in the near wall cells.
    postprocessTurbulentDissipationRate();

    // Update the turbulent viscosity.
    calculateF2();
    calculateTurbulentViscosity();

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);

    return;
} // integrateHierarchy

void
TwoEquationTurbulenceHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                        const double new_time,
                                                                        const bool skip_synchronize_new_state_data,
                                                                        const int num_cycles)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    AdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Deallocate any temporary data used to compute coefficients
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_k_temp_idx);
        level->deallocatePatchData(d_k_temp_rhs_idx);
        level->deallocatePatchData(d_k_dissipation_idx);
        level->deallocatePatchData(d_w_temp_idx);
        level->deallocatePatchData(d_w_temp_rhs_idx);
        level->deallocatePatchData(d_w_dissipation_idx);
        level->deallocatePatchData(d_k_C_idx);
        level->deallocatePatchData(d_w_C_idx);
        level->deallocatePatchData(d_rho_vec_cc_idx);
        // These two indexes are not deallocated because these two values are
        // required in example.cpp to plot the yplus and U_tau profile.
        // level->deallocatePatchData(d_yplus_cc_idx);
        // level->deallocatePatchData(d_U_tau_cc_idx);
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
}
void
TwoEquationTurbulenceHierarchyIntegrator::registerKVariable(Pointer<CellVariable<NDIM, double> > k_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(k_var);
#endif
    d_k_var = k_var;
    Pointer<CellDataFactory<NDIM, double> > k_factory = k_var->getPatchDataFactory();
    const int k_depth = k_factory->getDefaultDepth();
    Pointer<CellVariable<NDIM, double> > k_u_var = new CellVariable<NDIM, double>(k_var->getName() + "::u", k_depth);
    Pointer<CellVariable<NDIM, double> > k_rhs_var =
        new CellVariable<NDIM, double>(k_var->getName() + "::k_rhs", k_depth);
    Pointer<CellVariable<NDIM, double> > k_F_var = new CellVariable<NDIM, double>(k_var->getName() + "::F", k_depth);
    Pointer<SideVariable<NDIM, double> > k_diff_coef_var =
        new SideVariable<NDIM, double>(k_var->getName() + "::diff_coef", k_depth);

    // Set default values.
    d_k_u_var = nullptr;
    d_k_F_var = k_F_var;
    d_k_rhs_var = k_rhs_var;
    d_k_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    d_k_convective_difference_form = d_default_convective_difference_form;
    // if(!d_k_convective_time_stepping_type) d_k_convective_time_stepping_type =
    // d_default_convective_time_stepping_type;
    d_k_diffusion_coef_var = k_diff_coef_var;
    d_k_init = nullptr;
    d_k_F_fcn = nullptr;
    d_k_bc_coef = nullptr;
    d_k_reset_priority.push_back(std::numeric_limits<int>::max());
    return;
} // registerKVariable

void
TwoEquationTurbulenceHierarchyIntegrator::registerWVariable(Pointer<CellVariable<NDIM, double> > w_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(w_var);
#endif
    d_w_var = w_var;
    Pointer<CellDataFactory<NDIM, double> > w_factory = w_var->getPatchDataFactory();
    const int w_depth = w_factory->getDefaultDepth();
    Pointer<CellVariable<NDIM, double> > w_u_var = new CellVariable<NDIM, double>(w_var->getName() + "::u", w_depth);
    Pointer<CellVariable<NDIM, double> > w_rhs_var =
        new CellVariable<NDIM, double>(w_var->getName() + "::w_rhs", w_depth);
    Pointer<CellVariable<NDIM, double> > w_F_var = new CellVariable<NDIM, double>(w_var->getName() + "::F", w_depth);
    Pointer<SideVariable<NDIM, double> > w_diff_coef_var =
        new SideVariable<NDIM, double>(w_var->getName() + "::diff_coef", w_depth);

    // Set default values.
    d_w_u_var = nullptr;
    d_w_F_var = w_F_var;
    d_w_rhs_var = w_rhs_var;
    d_w_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    d_w_convective_difference_form = d_default_convective_difference_form;
    // d_w_convective_time_stepping_type = d_default_convective_time_stepping_type;
    d_w_diffusion_coef_var = w_diff_coef_var;
    d_w_is_diffusion_coef_variable = false;
    d_w_init = nullptr;
    d_w_bc_coef = nullptr;
    // d_w_bc_coef = std::vector<RobinBcCoefStrategy<NDIM>*>(w_depth, static_cast<RobinBcCoefStrategy<NDIM>*>(nullptr));
    d_w_reset_priority.push_back(std::numeric_limits<int>::max());
    return;
} // registerWVariable

Pointer<Variable<NDIM> >
TwoEquationTurbulenceHierarchyIntegrator::getKVariable() const
{
    return d_k_var;
} // getKVariable

Pointer<Variable<NDIM> >
TwoEquationTurbulenceHierarchyIntegrator::getWVariable() const
{
    return d_w_var;
} // getWVariable

void
TwoEquationTurbulenceHierarchyIntegrator::setAdvectionVelocityKEquation(Pointer<CellVariable<NDIM, double> > k_var,
                                                                        Pointer<FaceVariable<NDIM, double> > u_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
    TBOX_ASSERT(std::find(d_u_var.begin(), d_u_var.end(), u_var) != d_u_var.end());
#endif
    d_k_u_var = u_var;

    return;
} // setAdvectionVelocityKEquation
void
TwoEquationTurbulenceHierarchyIntegrator::setAdvectionVelocityWEquation(Pointer<CellVariable<NDIM, double> > w_var,
                                                                        Pointer<FaceVariable<NDIM, double> > u_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
    TBOX_ASSERT(std::find(d_u_var.begin(), d_u_var.end(), u_var) != d_u_var.end());
#endif
    d_w_u_var = u_var;

    return;
} // setAdvectionVelocityWEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setSourceTermFunctionKEquation(Pointer<CellVariable<NDIM, double> > k_var,
                                                                         Pointer<CartGridFunction> k_F_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
    TBOX_ASSERT(d_k_F_var);
#endif

    d_k_F_fcn = k_F_fcn;
    return;
} // setSourceTermFunctionKEquation
void
TwoEquationTurbulenceHierarchyIntegrator::setSourceTermFunctionWEquation(Pointer<CellVariable<NDIM, double> > w_F_var,
                                                                         Pointer<CartGridFunction> w_F_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
    TBOX_ASSERT(d_w_F_var);
#endif
    d_w_F_fcn = w_F_fcn;
    return;
} // setSourceTermFunctionWEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setDiffusionTimeSteppingTypeKEquation(
    Pointer<CellVariable<NDIM, double> > k_var,
    const TimeSteppingType diffusion_time_stepping_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
#endif
    d_k_diffusion_time_stepping_type = diffusion_time_stepping_type;
    return;
} // setDiffusionTimeSteppingTypeKEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setDiffusionTimeSteppingTypeWEquation(
    Pointer<CellVariable<NDIM, double> > w_var,
    const TimeSteppingType diffusion_time_stepping_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
#endif
    d_w_diffusion_time_stepping_type = diffusion_time_stepping_type;
    return;
} // setDiffusionTimeSteppingTypeWEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setConvectiveTimeSteppingTypeKEquation(
    Pointer<CellVariable<NDIM, double> > k_var,
    TimeSteppingType convective_time_stepping_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
#endif
    d_k_convective_time_stepping_type = convective_time_stepping_type;
    return;
} // setConvectiveTimeSteppingTypeKEquation
void
TwoEquationTurbulenceHierarchyIntegrator::setConvectiveTimeSteppingTypeWEquation(
    Pointer<CellVariable<NDIM, double> > w_var,
    TimeSteppingType convective_time_stepping_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
#endif
    d_w_convective_time_stepping_type = convective_time_stepping_type;
    return;
} // setConvectiveTimeSteppingTypeWEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setInitialConvectiveTimeSteppingTypeKEquation(
    Pointer<CellVariable<NDIM, double> > k_var,
    TimeSteppingType init_convective_time_stepping_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
#endif
    d_k_init_convective_time_stepping_type = init_convective_time_stepping_type;
    return;
} // setInitialConvectiveTimeSteppingTypeKEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setInitialConvectiveTimeSteppingTypeWEquation(
    Pointer<CellVariable<NDIM, double> > w_var,
    TimeSteppingType init_convective_time_stepping_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
#endif
    d_w_init_convective_time_stepping_type = init_convective_time_stepping_type;
    return;
} // setInitialConvectiveTimeSteppingTypeWEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setInitialConditionsKEquation(Pointer<CellVariable<NDIM, double> > k_var,
                                                                        Pointer<IBTK::CartGridFunction> k_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
#endif
    d_k_init = k_init;
    return;
} // setInitialConditionsKEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setInitialConditionsWEquation(Pointer<CellVariable<NDIM, double> > w_var,
                                                                        Pointer<IBTK::CartGridFunction> w_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
#endif
    d_w_init = w_init;
    return;
} // setInitialConditionsWEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setPhysicalBcCoefKEquation(Pointer<CellVariable<NDIM, double> > k_var,
                                                                     RobinBcCoefStrategy<NDIM>* k_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
#endif
    d_k_bc_coef = k_bc_coef;
    return;
} // setPhysicalBcCoefKEquation

RobinBcCoefStrategy<NDIM>*
TwoEquationTurbulenceHierarchyIntegrator::getPhysicalBcCoefKEquation()
{
    return d_k_bc_coef;
} // getPhysicalBcCoefKEquation

void
TwoEquationTurbulenceHierarchyIntegrator::setPhysicalBcCoefWEquation(Pointer<CellVariable<NDIM, double> > w_var,
                                                                     RobinBcCoefStrategy<NDIM>* w_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
#endif
    d_w_bc_coef = w_bc_coef;
    return;
} // setPhysicalBcCoefWEquation

RobinBcCoefStrategy<NDIM>*
TwoEquationTurbulenceHierarchyIntegrator::getPhysicalBcCoefWEquation()
{
    return d_w_bc_coef;
} // getPhysicalBcCoefWEquation

Pointer<PoissonSolver>
TwoEquationTurbulenceHierarchyIntegrator::getHelmholtzSolverKEquation(Pointer<CellVariable<NDIM, double> > k_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
#endif
    if (!d_k_solver)
    {
        const std::string& name = k_var->getName();
        d_k_solver =
            CCPoissonSolverManager::getManager()->allocateSolver(d_k_solver_type,
                                                                 d_object_name + "::helmholtz_solver::" + name,
                                                                 d_k_solver_db,
                                                                 "k_turbulence_",
                                                                 d_k_precond_type,
                                                                 d_object_name + "::helmholtz_precond::" + name,
                                                                 d_k_precond_db,
                                                                 "k_turbulence_pc_");
        d_k_solver_needs_init = true;
    }
    return d_k_solver;
} // getHelmholtzSolverKEquation

Pointer<PoissonSolver>
TwoEquationTurbulenceHierarchyIntegrator::getHelmholtzSolverWEquation(Pointer<CellVariable<NDIM, double> > w_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
#endif
    if (!d_w_solver)
    {
        const std::string& name = w_var->getName();
        d_w_solver =
            CCPoissonSolverManager::getManager()->allocateSolver(d_w_solver_type,
                                                                 d_object_name + "::helmholtz_solver::" + name,
                                                                 d_w_solver_db,
                                                                 "w_turbulence_",
                                                                 d_w_precond_type,
                                                                 d_object_name + "::helmholtz_precond::" + name,
                                                                 d_w_precond_db,
                                                                 "w_turbulence_pc_");
        d_w_solver_needs_init = true;
    }
    return d_w_solver;
} // getHelmholtzSolverWEquation

Pointer<LaplaceOperator>
TwoEquationTurbulenceHierarchyIntegrator::getHelmholtzRHSOperatorKEquation(Pointer<CellVariable<NDIM, double> > k_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
#endif
    const std::string& name = k_var->getName();
    if (!d_k_rhs_op)
    {
        d_k_rhs_op = new CCLaplaceOperator(d_object_name + "::helmholtz_rhs_op::" + name, /*homogeneous_bc*/ false);
        d_k_rhs_op_needs_init = true;
    }
    return d_k_rhs_op;
} // getHelmholtzRHSOperatorKEquation

Pointer<LaplaceOperator>
TwoEquationTurbulenceHierarchyIntegrator::getHelmholtzRHSOperatorWEquation(Pointer<CellVariable<NDIM, double> > w_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
#endif
    const std::string& name = w_var->getName();
    if (!d_w_rhs_op)
    {
        d_w_rhs_op = new CCLaplaceOperator(d_object_name + "::helmholtz_rhs_op::" + name, /*homogeneous_bc*/ false);
        d_w_rhs_op_needs_init = true;
    }
    return d_w_rhs_op;
} // getHelmholtzRHSOperatorWEquation

Pointer<ConvectiveOperator>
TwoEquationTurbulenceHierarchyIntegrator::getConvectiveOperatorKEquation(Pointer<CellVariable<NDIM, double> > k_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_k_var);
#endif
    if (!d_k_convective_op)
    {
        std::vector<RobinBcCoefStrategy<NDIM>*> k_bc_coefs(1, d_k_bc_coef);
        AdvDiffConvectiveOperatorManager* k_convective_op_manager = AdvDiffConvectiveOperatorManager::getManager();
        d_k_convective_op = k_convective_op_manager->allocateOperator(d_k_convective_op_type,
                                                                      d_object_name + "::KConvectiveOperator",
                                                                      d_k_var,
                                                                      d_k_convective_op_input_db,
                                                                      d_k_convective_difference_form,
                                                                      k_bc_coefs);
        d_k_convective_op_needs_init = true;
    }
    return d_k_convective_op;
} // getConvectiveOperatorKEquation

Pointer<ConvectiveOperator>
TwoEquationTurbulenceHierarchyIntegrator::getConvectiveOperatorWEquation(Pointer<CellVariable<NDIM, double> > w_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_w_var);
#endif
    if (!d_w_convective_op)
    {
        std::vector<RobinBcCoefStrategy<NDIM>*> w_bc_coefs(1, d_w_bc_coef);
        AdvDiffConvectiveOperatorManager* w_convective_op_manager = AdvDiffConvectiveOperatorManager::getManager();
        d_w_convective_op = w_convective_op_manager->allocateOperator(d_w_convective_op_type,
                                                                      d_object_name + "::WConvectiveOperator",
                                                                      d_w_var,
                                                                      d_w_convective_op_input_db,
                                                                      d_w_convective_difference_form,
                                                                      w_bc_coefs);
        d_w_convective_op_needs_init = true;
    }
    return d_w_convective_op;
} // getConvectiveOperatorWEquation

void
TwoEquationTurbulenceHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    AdvDiffSemiImplicitHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
        base_hierarchy, coarsest_level, finest_level);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_k_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_var, getScratchContext());
    d_w_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_var, getScratchContext());
    d_rho_cc_scratch_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getScratchContext());

    // Setup the patch boundary filling objects.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent k_bc_component(d_k_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     BDRY_EXTRAP_TYPE,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_k_bc_coef);
    d_k_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_k_bdry_bc_fill_op->initializeOperatorState(k_bc_component, d_hierarchy);

    InterpolationTransactionComponent w_bc_component(d_w_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     BDRY_EXTRAP_TYPE,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_w_bc_coef);
    d_w_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_w_bdry_bc_fill_op->initializeOperatorState(w_bc_component, d_hierarchy);

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* mu_bc_coef = d_ins_hierarchy_integrator->getViscosityBoundaryConditions();
    Pointer<CellVariable<NDIM, double> > mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
    const int mu_scratch_idx =
        var_db->mapVariableAndContextToIndex(mu_var, d_ins_hierarchy_integrator->getScratchContext());
    InterpolationTransactionComponent mu_bc_component(mu_scratch_idx,
                                                      DATA_REFINE_TYPE,
                                                      USE_CF_INTERPOLATION,
                                                      DATA_COARSEN_TYPE,
                                                      BDRY_EXTRAP_TYPE,
                                                      CONSISTENT_TYPE_2_BDRY,
                                                      mu_bc_coef);
    d_mu_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_mu_bdry_bc_fill_op->initializeOperatorState(mu_bc_component, d_hierarchy);

    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* mu_t_bc_coef =
        d_ins_hierarchy_integrator->getTurbulentViscosityBoundaryConditions();
    Pointer<CellVariable<NDIM, double> > mu_t_var = d_ins_hierarchy_integrator->getTurbulentViscosityVariable();
    const int mu_t_scratch_idx =
        var_db->mapVariableAndContextToIndex(mu_t_var, d_ins_hierarchy_integrator->getScratchContext());
    InterpolationTransactionComponent mu_t_bc_component(mu_t_scratch_idx,
                                                        DATA_REFINE_TYPE,
                                                        USE_CF_INTERPOLATION,
                                                        DATA_COARSEN_TYPE,
                                                        BDRY_EXTRAP_TYPE,
                                                        CONSISTENT_TYPE_2_BDRY,
                                                        mu_t_bc_coef);
    d_mu_t_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_mu_t_bdry_bc_fill_op->initializeOperatorState(mu_t_bc_component, d_hierarchy);

    std::vector<RobinBcCoefStrategy<NDIM>*> rho_bc_coefs =
        d_ins_hierarchy_integrator->getMassDensityBoundaryConditions();
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent rho_ghost_cc_interpolation(d_rho_cc_scratch_idx,
                                                                 DATA_REFINE_TYPE,
                                                                 USE_CF_INTERPOLATION,
                                                                 DATA_COARSEN_TYPE,
                                                                 BDRY_EXTRAP_TYPE,
                                                                 CONSISTENT_TYPE_2_BDRY,
                                                                 rho_bc_coefs[0]);
    d_rho_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_rho_bdry_bc_fill_op->initializeOperatorState(rho_ghost_cc_interpolation, d_hierarchy);

    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();

    // Reset the solution and rhs vectors.
    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    const int k_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_var, getScratchContext());
    d_k_sol = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::sol_vec::" + d_k_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_k_sol->addComponent(d_k_var, k_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    const int k_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_k_rhs_var, getScratchContext());
    d_k_rhs = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::rhs_vec::" + d_k_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_k_rhs->addComponent(d_k_rhs_var, k_rhs_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    const int w_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_var, getScratchContext());
    d_w_sol = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::sol_vec::" + d_w_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_w_sol->addComponent(d_w_var, w_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    const int w_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_w_rhs_var, getScratchContext());
    d_w_rhs = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::rhs_vec::" + d_w_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_w_rhs->addComponent(d_w_rhs_var, w_rhs_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    d_k_solver_needs_init = true;
    d_w_solver_needs_init = true;

    d_k_convective_op_needs_init = true;
    d_w_convective_op_needs_init = true;
    return;
}

Pointer<CellVariable<NDIM, double> >
TwoEquationTurbulenceHierarchyIntegrator::getCellCenteredMassDensityVariable()
{
#if (!NDEBUG)
    TBOX_ASSERT(d_rho_cc_var);
#endif
    return d_rho_cc_var;
} // getCellCenteredMassDensityVariable

Pointer<CellVariable<NDIM, double> >
TwoEquationTurbulenceHierarchyIntegrator::getCellCenteredUtauVariable()
{
    return d_U_tau_cc_var;
} // getCellCenteredUtauVariable

Pointer<CellVariable<NDIM, double> >
TwoEquationTurbulenceHierarchyIntegrator::getCellCenteredYplusVariable()
{
    return d_yplus_cc_var;
} // getCellCenteredYplusVariables

void
TwoEquationTurbulenceHierarchyIntegrator::registerINSVCStaggeredHierarchyIntegrator(
    Pointer<INSVCStaggeredHierarchyIntegrator> ins_cons_hier_integrator)
{
    d_ins_hierarchy_integrator = ins_cons_hier_integrator;
    return;
} // registerINSVCStaggeredHierarchyIntegrator

void
TwoEquationTurbulenceHierarchyIntegrator::calculateTurbulentViscosity()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getNewContext());

        d_mu_t_var = d_ins_hierarchy_integrator->getTurbulentViscosityVariable();
        const int mu_t_new_idx =
            var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getNewContext());

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > k_data = patch->getPatchData(d_k_new_idx);
            Pointer<CellData<NDIM, double> > w_data = patch->getPatchData(d_w_new_idx);
            Pointer<CellData<NDIM, double> > F2_data = patch->getPatchData(d_F2_scratch_idx);
            Pointer<CellData<NDIM, double> > P_data = patch->getPatchData(d_p_scratch_idx);
            Pointer<CellData<NDIM, double> > mu_t_data = patch->getPatchData(mu_t_new_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_new_idx);

            SST_MU_T_FCN(mu_t_data->getPointer(),
                         mu_t_data->getGhostCellWidth().max(),
                         k_data->getPointer(),
                         k_data->getGhostCellWidth().max(),
                         w_data->getPointer(),
                         w_data->getGhostCellWidth().max(),
                         rho_data->getPointer(),
                         rho_data->getGhostCellWidth().max(),
                         F2_data->getPointer(),
                         F2_data->getGhostCellWidth().max(),
                         P_data->getPointer(),
                         P_data->getGhostCellWidth().max(),

                         patch_box.lower(0),
                         patch_box.upper(0),
                         patch_box.lower(1),
                         patch_box.upper(1),
#if (NDIM == 3)
                         patch_box.lower(2),
                         patch_box.upper(2),
#endif
                         A1);
        }
    }

    return;
} // calculateTurbulentViscosity

void
TwoEquationTurbulenceHierarchyIntegrator::calculateTurbulentKEProduction(const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        d_p_scratch_idx = var_db->mapVariableAndContextToIndex(d_p_var, getScratchContext());
        Pointer<SideVariable<NDIM, double> > U_var = d_ins_hierarchy_integrator->getVelocityVariable();
        const int U_new_idx = var_db->mapVariableAndContextToIndex(U_var, d_ins_hierarchy_integrator->getNewContext());
        const int U_scratch_idx =
            var_db->mapVariableAndContextToIndex(U_var, d_ins_hierarchy_integrator->getScratchContext());
        d_hier_sc_data_ops->copyData(U_scratch_idx, U_new_idx);

        // filling ghost cells for velocity
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> U_bc_coefs =
            d_ins_hierarchy_integrator->getVelocityBoundaryConditions();
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        InterpolationTransactionComponent U_ghost_cc_interpolation(
            U_scratch_idx, "NONE", true, "CUBIC_COARSEN", "LINEAR", false, U_bc_coefs);
        HierarchyGhostCellInterpolation U_ghost_cell_fill_op;
        U_ghost_cell_fill_op.initializeOperatorState(U_ghost_cc_interpolation, d_hierarchy, coarsest_ln, finest_ln);
        U_ghost_cell_fill_op.fillData(data_time);

        // Get the data index of the turbulent viscosity associated with new index
        d_mu_t_var = d_ins_hierarchy_integrator->getTurbulentViscosityVariable();
        const int mu_t_new_idx =
            var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getNewContext());

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(U_scratch_idx);

            const IntVector<NDIM>& u_ghost_cells = u_data->getGhostCellWidth();
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            Pointer<CellData<NDIM, double> > mu_t_data = patch->getPatchData(mu_t_new_idx);
            Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(d_p_scratch_idx);

            // routine to calculate production term
            PRODUCTION(p_data->getPointer(),
                       p_data->getGhostCellWidth().max(),
                       mu_t_data->getPointer(),
                       mu_t_data->getGhostCellWidth().max(),
                       u_data->getPointer(0),
                       u_data->getPointer(1),
#if (NDIM == 3)
                       u_data->getPointer(2),
#endif
                       u_ghost_cells(0),
                       u_ghost_cells(1),
#if (NDIM == 3)
                       u_ghost_cells(2),
#endif
                       patch_box.lower(0),
                       patch_box.upper(0),
                       patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2),
                       patch_box.upper(2),
#endif
                       dx);
        }
    }

    return;
} // calculateTurbulentKEProduction

void
TwoEquationTurbulenceHierarchyIntegrator::calculateF2()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();

            // Get the data index of the density and viscosity associated with new index
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getNewContext());
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_new_idx);
            d_mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
            const int mu_new_idx =
                var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getNewContext());
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_new_idx);

            Pointer<CellData<NDIM, double> > k_data = patch->getPatchData(d_k_new_idx);
            Pointer<CellData<NDIM, double> > w_data = patch->getPatchData(d_w_new_idx);
            Pointer<CellData<NDIM, double> > F2_data = patch->getPatchData(d_F2_scratch_idx);
            Pointer<CellData<NDIM, double> > d_data = patch->getPatchData(d_distance_to_closest_surface_scratch_idx);

            // routine to calculate blending function
            SST_F2_FCN(F2_data->getPointer(),
                       F2_data->getGhostCellWidth().max(),
                       k_data->getPointer(),
                       k_data->getGhostCellWidth().max(),
                       w_data->getPointer(),
                       w_data->getGhostCellWidth().max(),
                       mu_data->getPointer(),
                       mu_data->getGhostCellWidth().max(),
                       rho_data->getPointer(),
                       rho_data->getGhostCellWidth().max(),
                       d_data->getPointer(),
                       d_data->getGhostCellWidth().max(),
                       patch_box.lower(0),
                       patch_box.upper(0),
                       patch_box.lower(1),
                       patch_box.upper(1),
#if (NDIM == 3)
                       patch_box.lower(2),
                       patch_box.upper(2),
#endif
                       BETA_STAR);
        }
    }
    return;
} // calculateF2

void
TwoEquationTurbulenceHierarchyIntegrator::calculateBlendingFunction(const double data_time,
                                                                    const Pointer<VariableContext> ctx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_F1_var);
#endif
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int level_number = coarsest_ln; level_number <= finest_ln; ++level_number)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_number);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const SAMRAI::hier::Index<NDIM>& patch_lower_index = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            const double* grid_lower = grid_geom->getXLower();
            const double* grid_upper = grid_geom->getXUpper();
            const double* geometry_coord[2] = { grid_lower, grid_upper };

            Pointer<CellData<NDIM, double> > d_data = patch->getPatchData(d_distance_to_closest_surface_scratch_idx);
            double distance = 0.0;
            // Loop over cells
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                IBTK::VectorNd cell_coord;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    cell_coord[d] = patch_X_lower[d] + patch_dx[d] * (ci(d) - patch_lower_index(d) + 0.5);
                }

                std::vector<double> distance_to_surface(2 * NDIM, std::numeric_limits<int>::max());

                for (int i = 0; i < d_wall_location_index.size(); i++)
                {
                    // since we need both modulo and division here, using div() function.
                    std::div_t coord = std::div(d_wall_location_index[i], 2);
                    distance = std::abs(geometry_coord[coord.rem][coord.quot] - cell_coord[coord.quot]);
                    distance_to_surface[i] = distance;
                }
                (*d_data)(ci) = *min_element(distance_to_surface.begin(), distance_to_surface.end());
            }

            // Calculating F1
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

            // filling ghost cells for turbulent kinetic energy
            const int k_idx = var_db->mapVariableAndContextToIndex(d_k_var, ctx);
            d_hier_cc_data_ops->copyData(d_k_scratch_idx, k_idx);
            d_k_bdry_bc_fill_op->fillData(data_time);
            Pointer<CellData<NDIM, double> > k_data = patch->getPatchData(d_k_scratch_idx);

            const int w_idx = var_db->mapVariableAndContextToIndex(d_w_var, ctx);
            d_hier_cc_data_ops->copyData(d_w_scratch_idx, w_idx);
            d_w_bdry_bc_fill_op->fillData(data_time);

            Pointer<CellData<NDIM, double> > w_data = patch->getPatchData(d_w_scratch_idx);
            Pointer<CellData<NDIM, double> > F1_data = patch->getPatchData(d_F1_scratch_idx);

            const int rho_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, ctx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

            d_mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
            int mu_idx = 0;
            if (ctx == getNewContext())
            {
                mu_idx = var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getNewContext());
            }
            else if (ctx == getCurrentContext())
            {
                mu_idx =
                    var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getCurrentContext());
            }
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

            // Routine to calculate blending function.
            SST_BLENDING_FCN(F1_data->getPointer(),
                             F1_data->getGhostCellWidth().max(),
                             k_data->getPointer(),
                             k_data->getGhostCellWidth().max(),
                             w_data->getPointer(),
                             w_data->getGhostCellWidth().max(),
                             mu_data->getPointer(),
                             mu_data->getGhostCellWidth().max(),
                             rho_data->getPointer(),
                             rho_data->getGhostCellWidth().max(),
                             d_data->getPointer(),
                             d_data->getGhostCellWidth().max(),
                             BETA_STAR,
                             SIGMA_W2,
                             patch_box.lower(0),
                             patch_box.upper(0),
                             patch_box.lower(1),
                             patch_box.upper(1),
#if (NDIM == 3)
                             patch_box.lower(2),
                             patch_box.upper(2),
#endif
                             patch_dx);

            // Perform the F1 weighting to turbulence coefficients.
            // Loop over cells.
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                d_sigma_k = (*F1_data)(ci)*SIGMA_K1 + (1.0 - (*F1_data)(ci)) * SIGMA_K2;
                d_sigma_w = (*F1_data)(ci)*SIGMA_W1 + (1.0 - (*F1_data)(ci)) * SIGMA_W2;
                d_beta = (*F1_data)(ci)*BETA_1 + (1.0 - (*F1_data)(ci)) * BETA_2;
            }
        }
    }

    return;
} // calculateBlendingFunction

void
TwoEquationTurbulenceHierarchyIntegrator::interpolateSCMassDensityToCC(Pointer<VariableContext> ctx)
{
#if (!NDEBUG)
    TBOX_ASSERT(d_rho_cc_var);
#endif
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double> > rho_sc_var = d_ins_hierarchy_integrator->getMassDensityVariable();
    int rho_sc_idx = 0;
    if (ctx == getNewContext())
    {
        rho_sc_idx = var_db->mapVariableAndContextToIndex(rho_sc_var, d_ins_hierarchy_integrator->getNewContext());
    }
    else if (ctx == getCurrentContext())
    {
        rho_sc_idx = var_db->mapVariableAndContextToIndex(rho_sc_var, d_ins_hierarchy_integrator->getCurrentContext());
    }

    const int rho_cc_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, ctx);
    const int rho_vec_cc_idx = var_db->mapVariableAndContextToIndex(d_rho_vec_cc_var, getCurrentContext());
    static const bool synch_cf_interface = true;

    d_hier_math_ops->interp(
        rho_vec_cc_idx, d_rho_vec_cc_var, rho_sc_idx, rho_sc_var, d_no_fill_op, d_integrator_time, synch_cf_interface);

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
    return;
} // interpolateSCMassDensityToCC

void
TwoEquationTurbulenceHierarchyIntegrator::applyWallFunction(const double data_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    FrictionVelocityFunctor friction_velocity;
    friction_velocity.d_plus = DPLUS;
    friction_velocity.kappa = KAPPA;
    friction_velocity.E = E;

    for (int ln = coarsest_ln; ln <= finest_ln; ln++)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<SideVariable<NDIM, double> > U_var = d_ins_hierarchy_integrator->getVelocityVariable();
        const int U_new_idx = var_db->mapVariableAndContextToIndex(U_var, d_ins_hierarchy_integrator->getNewContext());
        const int U_scratch_idx =
            var_db->mapVariableAndContextToIndex(U_var, d_ins_hierarchy_integrator->getScratchContext());
        d_hier_sc_data_ops->copyData(U_scratch_idx, U_new_idx);

        Pointer<SideVariable<NDIM, double> > U_tau_sc_var = d_ins_hierarchy_integrator->getUtauVariable();
        const int U_tau_sc_idx =
            var_db->mapVariableAndContextToIndex(U_tau_sc_var, d_ins_hierarchy_integrator->getScratchContext());
        d_hier_sc_data_ops->setToScalar(U_tau_sc_idx, 0.0);
        Pointer<SideVariable<NDIM, double> > yplus_sc_var = d_ins_hierarchy_integrator->getYplusVariable();
        const int yplus_sc_idx =
            var_db->mapVariableAndContextToIndex(yplus_sc_var, d_ins_hierarchy_integrator->getScratchContext());
        d_hier_sc_data_ops->setToScalar(yplus_sc_idx, 0.0);
        Pointer<CellVariable<NDIM, double> > mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
        const int mu_new_idx =
            var_db->mapVariableAndContextToIndex(mu_var, d_ins_hierarchy_integrator->getNewContext());
        const int mu_scratch_idx =
            var_db->mapVariableAndContextToIndex(mu_var, d_ins_hierarchy_integrator->getScratchContext());
        d_hier_cc_data_ops->copyData(mu_scratch_idx, mu_new_idx);
        d_mu_bdry_bc_fill_op->fillData(data_time);

        const int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getNewContext());
        const int rho_scratch_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getScratchContext());
        d_hier_cc_data_ops->copyData(rho_scratch_idx, rho_new_idx);
        d_rho_bdry_bc_fill_op->fillData(data_time);

        // filling ghost cells for velocity
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> U_bc_coefs =
            d_ins_hierarchy_integrator->getVelocityBoundaryConditions();
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        InterpolationTransactionComponent U_ghost_cc_interpolation(
            U_scratch_idx, "NONE", true, "CUBIC_COARSEN", "LINEAR", false, U_bc_coefs);
        HierarchyGhostCellInterpolation U_ghost_cell_fill_op;
        U_ghost_cell_fill_op.initializeOperatorState(U_ghost_cc_interpolation, d_hierarchy, coarsest_ln, finest_ln);
        U_ghost_cell_fill_op.fillData(data_time);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            const double* grid_lower = grid_geom->getXLower();
            const double* grid_upper = grid_geom->getXUpper();
            Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_new_idx);
            Pointer<SideData<NDIM, double> > U_tau_data = patch->getPatchData(U_tau_sc_idx);
            Pointer<SideData<NDIM, double> > yplus_sc_data = patch->getPatchData(yplus_sc_idx);

            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_scratch_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_scratch_idx);

            IBTK::VectorNd number_of_indices(NDIM);
            for (int d = 0; d < NDIM; d++)
            {
                number_of_indices(d) = (grid_upper[d] - grid_lower[d]) / patch_dx[d];
            }
            friction_velocity.delta = d_distance_to_virtual_point * patch_dx[1];

            // Set up Newton itertion
            const double guess = FrictionVelocityFunctor::s_newton_guess;
            const double min = FrictionVelocityFunctor::s_newton_min;
            const double max = FrictionVelocityFunctor::s_newton_max;

            // Accuracy doubles at each step, so stop when just over half of the
            // digits are correct, and rely on that step to polish off the remainder:
            auto get_digits = static_cast<int>(std::numeric_limits<double>::digits * 0.6);
            const boost::uintmax_t maxit = 50;
            boost::uintmax_t it = maxit;

            if (!patch_geom->getTouchesRegularBoundary()) continue;

            for (unsigned int i = 0; i < d_wall_location_index.size(); i++)
            {
                const unsigned int axis = d_wall_location_index[i] / 2;
                const unsigned int side = d_wall_location_index[i] % 2;
                const bool is_lower = side == 0;
                if (patch_geom->getTouchesRegularBoundary(axis, side))
                {
                    Box<NDIM> bdry_box = patch_box;
                    if (is_lower)
                    {
                        bdry_box.upper(axis) = patch_box.lower(axis) + d_distance_to_virtual_point;
                    }
                    else
                    {
                        bdry_box.lower(axis) = patch_box.upper(axis) - d_distance_to_virtual_point;
                    }

                    Box<NDIM> trim_box = patch_box * bdry_box;
                    for (unsigned int d = 0; d < NDIM; d++)
                    {
                        for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(trim_box, d)); b; b++)
                        {
                            SideIndex<NDIM> si(b(), d, SideIndex<NDIM>::Lower);
                            // only for horizontal boundaries.
                            IntVector<NDIM> index_1(0, 1);
                            IntVector<NDIM> index_2(0, 2);
                            SideIndex<NDIM> si_1, si_2;

                            // bottom boundary.
                            if (axis == 1 && side == 0 && si(1) == 0 && d == 0)
                            {
                                si_1 = si + index_1;
                                si_2 = si + index_2;

                                const double mu = 0.5 * ((*mu_data)(si.toCell(0)) + (*mu_data)(si.toCell(1)));
                                const double rho = 0.5 * ((*rho_data)(si.toCell(0)) + (*rho_data)(si.toCell(1)));
                                friction_velocity.nu = mu / rho;
                                friction_velocity.U_s = 0.5 * ((*U_data)(si_1) + (*U_data)(si_2));
                                (*U_tau_data)(si) = boost::math::tools::newton_raphson_iterate(
                                    friction_velocity, guess, min, max, get_digits, it);
                                (*yplus_sc_data)(si) = rho * (*U_tau_data)(si)*patch_dx[axis] / (2.0 * mu);
                            }
                            // top boundary
                            else if (axis == 1 && side == 1 && si(1) == number_of_indices(axis) - 1 && d == 0)
                            {
                                si_1 = si - index_1;
                                si_2 = si - index_2;

                                const double mu = 0.5 * ((*mu_data)(si.toCell(0)) + (*mu_data)(si.toCell(1)));
                                const double rho = 0.5 * ((*rho_data)(si.toCell(0)) + (*rho_data)(si.toCell(1)));

                                friction_velocity.nu = mu / rho;
                                friction_velocity.U_s = 0.5 * ((*U_data)(si_1) + (*U_data)(si_2));
                                (*U_tau_data)(si) = boost::math::tools::newton_raphson_iterate(
                                    friction_velocity, guess, min, max, get_digits, it);
                                (*yplus_sc_data)(si) = rho * (*U_tau_data)(si)*patch_dx[axis] / (2.0 * mu);
                            }
                        }
                    }
                }
            }
        }
    }
    return;
} // applyWallFunction

void
TwoEquationTurbulenceHierarchyIntegrator::postprocessTurbulentVariablesBasedonYplus()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ln++)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int yplus_cc_idx = var_db->mapVariableAndContextToIndex(d_yplus_cc_var, getCurrentContext());
        const int U_tau_cc_idx = var_db->mapVariableAndContextToIndex(d_U_tau_cc_var, getCurrentContext());
        const int k_new_idx = var_db->mapVariableAndContextToIndex(d_k_var, getNewContext());
        const int w_new_idx = var_db->mapVariableAndContextToIndex(d_w_var, getNewContext());
        const int distance_to_closest_surface_idx =
            var_db->mapVariableAndContextToIndex(d_distance_to_closest_surface_var, getScratchContext());

        Pointer<CellVariable<NDIM, double> > mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
        const int mu_new_idx =
            var_db->mapVariableAndContextToIndex(mu_var, d_ins_hierarchy_integrator->getNewContext());
        const int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getNewContext());

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            const double* grid_lower = grid_geom->getXLower();
            const double* grid_upper = grid_geom->getXUpper();

            Pointer<CellData<NDIM, double> > k_data = patch->getPatchData(k_new_idx);
            Pointer<CellData<NDIM, double> > w_data = patch->getPatchData(w_new_idx);
            Pointer<CellData<NDIM, double> > distance_to_closest_surface_data =
                patch->getPatchData(distance_to_closest_surface_idx);

            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_new_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_new_idx);

            Pointer<CellData<NDIM, double> > U_tau_data = patch->getPatchData(U_tau_cc_idx);
            Pointer<CellData<NDIM, double> > yplus_data = patch->getPatchData(yplus_cc_idx);

            IBTK::VectorNd number_of_indices(NDIM);
            for (int d = 0; d < NDIM; d++)
            {
                number_of_indices(d) = (grid_upper[d] - grid_lower[d]) / patch_dx[d];
            }
            if (!patch_geom->getTouchesRegularBoundary()) continue;

            for (int i = 0; i < d_wall_location_index.size(); i++)
            {
                const unsigned int axis = d_wall_location_index[i] / 2;
                const unsigned int side = d_wall_location_index[i] % 2;
                const bool is_lower = side == 0;
                if (patch_geom->getTouchesRegularBoundary(axis, side))
                {
                    Box<NDIM> bdry_box = patch_box;
                    if (is_lower)
                    {
                        bdry_box.upper(axis) = patch_box.lower(axis) + d_distance_to_virtual_point;
                    }
                    else
                    {
                        bdry_box.lower(axis) = patch_box.upper(axis) - d_distance_to_virtual_point;
                    }

                    Box<NDIM> trim_box = patch_box * bdry_box;

                    for (Box<NDIM>::Iterator b(trim_box); b; b++)
                    {
                        CellIndex<NDIM> ci(b());
                        IntVector<NDIM> index_1(0, 1);
                        CellIndex<NDIM> ci_1;
                        if (ci(1) == 0 || ci(1) == number_of_indices(1) - 1)
                        {
                            ci_1 = (ci(1) == 0) ? ci + index_1 : ci - index_1;
                            double sum = 0.0, sum_1 = 0.0;
                            for (unsigned int d = 0; d < NDIM; d++)
                            {
                                sum += (*yplus_data)(ci, d);
                                sum_1 += (*U_tau_data)(ci, d);
                            }
                            // True only for equidistant cells.
                            double yplus = 0.5 * sum;
                            double U_tau = 0.5 * sum_1;

                            if (yplus <= 11.0)
                            {
                                // lineear interpolation.
                                (*k_data)(ci) = (*distance_to_closest_surface_data)(ci) * (*k_data)(ci_1) /
                                                (*distance_to_closest_surface_data)(ci_1);

                                (*w_data)(ci) = 6.0 * (*mu_data)(ci) /
                                                ((*rho_data)(ci)*BETA_1 * (*distance_to_closest_surface_data)(ci) *
                                                 (*distance_to_closest_surface_data)(ci));
                            }
                            else
                            {
                                (*k_data)(ci) = U_tau * U_tau / sqrt(BETA_STAR);

                                (*w_data)(ci) =
                                    U_tau / sqrt(BETA_STAR * KAPPA * (*distance_to_closest_surface_data)(ci));
                            }
                        }
                    }
                }
            }
        }
    }
    return;
} // postProcessTurbulentVariablesBasedonYplus

void
TwoEquationTurbulenceHierarchyIntegrator::postprocessTurbulentDissipationRate()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ln++)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int k_new_idx = var_db->mapVariableAndContextToIndex(d_k_var, getNewContext());
        const int w_new_idx = var_db->mapVariableAndContextToIndex(d_w_var, getNewContext());
        Pointer<CellVariable<NDIM, double> > mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
        const int mu_new_idx =
            var_db->mapVariableAndContextToIndex(mu_var, d_ins_hierarchy_integrator->getNewContext());
        const int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getNewContext());

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            const double* grid_lower = grid_geom->getXLower();
            const double* grid_upper = grid_geom->getXUpper();

            Pointer<CellData<NDIM, double> > k_data = patch->getPatchData(k_new_idx);
            Pointer<CellData<NDIM, double> > w_data = patch->getPatchData(w_new_idx);
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_new_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_new_idx);

            IBTK::VectorNd number_of_indices(NDIM);
            for (int d = 0; d < NDIM; d++)
            {
                number_of_indices(d) = (grid_upper[d] - grid_lower[d]) / patch_dx[d];
            }
            if (!patch_geom->getTouchesRegularBoundary()) continue;
            for (int i = 0; i < d_wall_location_index.size(); i++)
            {
                const unsigned int axis = d_wall_location_index[i] / 2;
                const unsigned int side = d_wall_location_index[i] % 2;
                const bool is_lower = side == 0;
                if (patch_geom->getTouchesRegularBoundary(axis, side))
                {
                    Box<NDIM> bdry_box = patch_box;
                    if (is_lower)
                    {
                        bdry_box.upper(axis) = patch_box.lower(axis);
                    }
                    else
                    {
                        bdry_box.lower(axis) = patch_box.upper(axis);
                    }

                    Box<NDIM> trim_box = patch_box * bdry_box;
                    for (Box<NDIM>::Iterator b(trim_box); b; b++)
                    {
                        CellIndex<NDIM> ci(b());
                        // This is only for horizontal boundaries.
                        if (ci(1) == 0 || ci(1) == number_of_indices(1) - 1)
                        {
                            const double w_vis =
                                6.0 * 4.0 * (*mu_data)(ci) / (BETA_1 * patch_dx[1] * patch_dx[1] * (*rho_data)(ci));
                            const double U_star = std::pow(BETA_STAR, 0.25) * std::pow(((*k_data)(ci)), 0.5);
                            const double w_log = 2.0 * U_star / (KAPPA * patch_dx[1] * std::sqrt(BETA_STAR));
                            (*w_data)(ci) = std::sqrt((w_vis * w_vis) + (w_log * w_log));
                            // std::cout << "w data in top and bottom near wall cell is\t" << ci << "\t" <<
                            // (*w_data)(ci)
                            //<< std::endl;
                        }
                    }
                }
            }
        }
    }
    return;
} // postprocessTurbulentDissipationRate

void
TwoEquationTurbulenceHierarchyIntegrator::computeWallShearStressFromWallLaw(const double data_time,
                                                                            const Pointer<VariableContext> ctx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ln++)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        // filling ghost cells for k, rho, mu and mu_t. Because these cell-centered quantities will be interpolated to
        // node or side.
        int k_idx = var_db->mapVariableAndContextToIndex(d_k_var, ctx);
        d_hier_cc_data_ops->copyData(d_k_scratch_idx, k_idx);
        d_k_bdry_bc_fill_op->fillData(data_time);

        const int rho_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, ctx);
        const int rho_scratch_idx = var_db->mapVariableAndContextToIndex(d_rho_cc_var, getScratchContext());
        d_hier_cc_data_ops->copyData(rho_scratch_idx, rho_idx);
        d_rho_bdry_bc_fill_op->fillData(data_time);

        d_mu_var = d_ins_hierarchy_integrator->getViscosityVariable();
        int mu_idx = 0;
        if (ctx == getNewContext())
        {
            mu_idx = var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getNewContext());
        }
        else if (ctx == getCurrentContext())
        {
            mu_idx = var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getCurrentContext());
        }
        const int mu_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_mu_var, d_ins_hierarchy_integrator->getScratchContext());
        d_hier_cc_data_ops->copyData(mu_scratch_idx, mu_idx);
        d_mu_bdry_bc_fill_op->fillData(data_time);

        d_mu_t_var = d_ins_hierarchy_integrator->getTurbulentViscosityVariable();
        int mu_t_idx = 0;
        if (ctx == getNewContext())
        {
            mu_t_idx = var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getNewContext());
        }
        else if (ctx == getCurrentContext())
        {
            mu_t_idx =
                var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getCurrentContext());
        }
        const int mu_t_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_mu_t_var, d_ins_hierarchy_integrator->getScratchContext());
        d_hier_cc_data_ops->copyData(mu_t_scratch_idx, mu_t_idx);
        d_mu_t_bdry_bc_fill_op->fillData(data_time);

        Pointer<SideVariable<NDIM, double> > U_var = d_ins_hierarchy_integrator->getVelocityVariable();
        Pointer<NodeVariable<NDIM, double> > tau_w_var = d_ins_hierarchy_integrator->getTauwVariable();

        int U_idx = 0;
        if (ctx == getNewContext())
        {
            U_idx = var_db->mapVariableAndContextToIndex(U_var, d_ins_hierarchy_integrator->getNewContext());
        }
        else if (ctx == getCurrentContext())
        {
            U_idx = var_db->mapVariableAndContextToIndex(U_var, d_ins_hierarchy_integrator->getCurrentContext());
        }

        const int tau_w_idx =
            var_db->mapVariableAndContextToIndex(tau_w_var, d_ins_hierarchy_integrator->getScratchContext());
        // d_hier_nc_data_ops->setToScalar(tau_w_idx, 0.0);

        Pointer<SideVariable<NDIM, double> > yplus_var = d_ins_hierarchy_integrator->getYplusVariable();
        Pointer<SideVariable<NDIM, double> > U_tau_var = d_ins_hierarchy_integrator->getUtauVariable();
        const int yplus_idx =
            var_db->mapVariableAndContextToIndex(yplus_var, d_ins_hierarchy_integrator->getScratchContext());
        // d_hier_sc_data_ops->setToScalar(yplus_idx, 0.0);
        const int U_tau_idx =
            var_db->mapVariableAndContextToIndex(U_tau_var, d_ins_hierarchy_integrator->getScratchContext());
        // d_hier_sc_data_ops->setToScalar(U_tau_idx, 0.0);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            Pointer<CellData<NDIM, double> > k_data = patch->getPatchData(d_k_scratch_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_scratch_idx);
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_scratch_idx);
            Pointer<CellData<NDIM, double> > mu_t_data = patch->getPatchData(mu_t_scratch_idx);

            Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
            const IntVector<NDIM>& U_ghost_cells = U_data->getGhostCellWidth();
            Pointer<NodeData<NDIM, double> > tau_w_data = patch->getPatchData(tau_w_idx);
            Pointer<SideData<NDIM, double> > yplus_data = patch->getPatchData(yplus_idx);
            const IntVector<NDIM>& yplus_ghost_cells = yplus_data->getGhostCellWidth();
            Pointer<SideData<NDIM, double> > U_tau_data = patch->getPatchData(U_tau_idx);
            const IntVector<NDIM>& U_tau_ghost_cells = U_tau_data->getGhostCellWidth();

            for (unsigned int i = 0; i < d_wall_location_index.size(); i++)
            {
                const unsigned int axis = d_wall_location_index[i] / 2;
                const unsigned int side = d_wall_location_index[i] % 2;
                const bool is_lower = side == 0;
                if (patch_geom->getTouchesRegularBoundary(axis, side))
                {
                    Box<NDIM> bdry_box = patch_box;
                    if (is_lower)
                    {
                        bdry_box.upper(axis) = patch_box.lower(axis);
                    }
                    else
                    {
                        bdry_box.lower(axis) = patch_box.upper(axis);
                    }

                    Box<NDIM> trim_box = patch_box * bdry_box;
                    WALL_SHEAR_STRESS(tau_w_data->getPointer(),
                                      tau_w_data->getGhostCellWidth().max(),
                                      U_data->getPointer(0),
                                      U_ghost_cells(0),
                                      U_data->getPointer(1),
                                      U_ghost_cells(1),
#if (NDIM == 3)
                                      U_data->getPointer(2),
                                      U_ghost_cells(2),
#endif
                                      k_data->getPointer(),
                                      k_data->getGhostCellWidth().max(),
                                      mu_t_data->getPointer(),
                                      mu_t_data->getGhostCellWidth().max(),
                                      rho_data->getPointer(),
                                      rho_data->getGhostCellWidth().max(),
                                      mu_data->getPointer(),
                                      mu_data->getGhostCellWidth().max(),
                                      U_tau_data->getPointer(0),
                                      U_tau_ghost_cells(0),
                                      U_tau_data->getPointer(1),
                                      U_tau_ghost_cells(1),
#if (NDIM == 3)
                                      U_tau_data->getPointer(2),
                                      U_tau_ghost_cells(2),
#endif
                                      yplus_data->getPointer(0),
                                      yplus_ghost_cells(0),
                                      yplus_data->getPointer(1),
                                      yplus_ghost_cells(1),
#if (NDIM == 3)
                                      yplus_data->getPointer(2),
                                      yplus_ghost_cells(2),
#endif
                                      KAPPA,
                                      BETA_STAR,
                                      B,
                                      d_wall_location_index[i],
                                      trim_box.lower(0),
                                      trim_box.upper(0),
                                      trim_box.lower(1),
                                      trim_box.upper(1),
#if (NDIM == 3)
                                      trim_box.lower(2),
                                      trim_box.upper(2),
#endif

                                      patch_dx);
                }
            }
        }
    }
    return;
} // computeWallShearStressFromWallLaw

void
TwoEquationTurbulenceHierarchyIntegrator::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db);
#endif
    // Read in data members from input database.
    if (!is_from_restart)
    {
        d_k_solver_type = CCPoissonSolverManager::UNDEFINED;
        d_k_precond_type = CCPoissonSolverManager::UNDEFINED;
        if (input_db->keyExists("k_solver_type"))
        {
            d_k_solver_type = input_db->getString("k_solver_type");
            if (input_db->keyExists("k_solver_db")) d_k_solver_db = input_db->getDatabase("k_solver_db");
        }
        if (!d_k_solver_db) d_k_solver_db = new MemoryDatabase("k_solver_db");

        if (input_db->keyExists("k_precond_type"))
        {
            d_k_precond_type = input_db->getString("k_precond_type");
            if (input_db->keyExists("k_precond_db")) d_k_precond_db = input_db->getDatabase("k_precond_db");
        }
        if (!d_k_precond_db) d_k_precond_db = new MemoryDatabase("k_precond_db");

        d_w_solver_type = CCPoissonSolverManager::UNDEFINED;
        d_w_precond_type = CCPoissonSolverManager::UNDEFINED;
        if (input_db->keyExists("w_solver_type"))
        {
            d_w_solver_type = input_db->getString("w_solver_type");
            if (input_db->keyExists("w_solver_db")) d_w_solver_db = input_db->getDatabase("w_solver_db");
        }
        if (!d_w_solver_db) d_w_solver_db = new MemoryDatabase("w_solver_db");

        if (input_db->keyExists("w_precond_type"))
        {
            d_w_precond_type = input_db->getString("w_precond_type");
            if (input_db->keyExists("w_precond_db")) d_w_precond_db = input_db->getDatabase("w_precond_db");
        }
        if (!d_w_precond_db) d_w_precond_db = new MemoryDatabase("w_precond_db");

        if (input_db->keyExists("k_convective_difference_form"))
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
            d_rho_bdry_extrap_type = input_db->getString("rho_bdry_extrap_type");
    }
}

} // namespace IBAMR
