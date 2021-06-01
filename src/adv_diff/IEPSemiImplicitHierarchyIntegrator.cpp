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

    // Initialize conservative mass and transported quantity integrator.
    if (d_solve_mass_conservation)
        d_rho_p_integrator = new AdvDiffConservativeMassTransportQuantityIntegrator(
            "AdvDiffConservativeMassTransportQuantityIntegrator::MassMomentumIntegrator",
            input_db->getDatabase("mass_transport_integrator_db"));

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

    int lf_u_current_idx, lf_u_scratch_idx, lf_u_new_idx;
    if (d_lf_u_var)
        registerVariable(lf_u_current_idx,
                         lf_u_new_idx,
                         lf_u_scratch_idx,
                         d_lf_u_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int T_u_current_idx, T_u_scratch_idx, T_u_new_idx;
    if (d_T_u_var)
        registerVariable(T_u_current_idx,
                         T_u_new_idx,
                         T_u_scratch_idx,
                         d_T_u_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int lf_diffusion_coef_rhs_scratch_idx, T_diffusion_coef_rhs_scratch_idx;
    d_lf_diffusion_coef_rhs_var = new SideVariable<NDIM, double>(d_lf_var->getName() + "::Diff");
    registerVariable(lf_diffusion_coef_rhs_scratch_idx, d_lf_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    d_T_diffusion_coef_rhs_var = new SideVariable<NDIM, double>(d_T_var->getName() + "::Diff");
    registerVariable(T_diffusion_coef_rhs_scratch_idx, d_T_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    int lf_rhs_scratch_idx, T_rhs_scratch_idx;
    registerVariable(lf_rhs_scratch_idx, d_lf_rhs_var, cell_ghosts, getScratchContext());
    registerVariable(T_rhs_scratch_idx, d_T_rhs_var, cell_ghosts, getScratchContext());

    int lf_H_scratch_idx;
    d_lf_H_var = new CellVariable<NDIM, double>(d_object_name + "::lf_H");
    registerVariable(lf_H_scratch_idx, d_lf_H_var, cell_ghosts, getScratchContext());

    int div_u_scratch_idx;
    d_div_u_var = new CellVariable<NDIM, double>(d_object_name + "::div_u");
    registerVariable(div_u_scratch_idx, d_div_u_var, no_ghosts, getScratchContext());

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_lf_C_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::C");
    d_lf_C_idx = var_db->registerVariableAndContext(d_lf_C_var, getCurrentContext(), no_ghosts);

    d_T_C_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::C");
    d_T_C_idx = var_db->registerVariableAndContext(d_T_C_var, getCurrentContext(), no_ghosts);

    d_lf_temp_rhs_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::temp_rhs");
    d_lf_temp_rhs_idx = var_db->registerVariableAndContext(d_lf_temp_rhs_var, getCurrentContext(), no_ghosts);

    d_T_temp_rhs_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::temp_rhs");
    d_T_temp_rhs_idx = var_db->registerVariableAndContext(d_T_temp_rhs_var, getCurrentContext(), no_ghosts);

    d_lf_N_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::N");
    int lf_N_scratch_idx;
    registerVariable(lf_N_scratch_idx, d_lf_N_var, cell_ghosts, getScratchContext());

    d_T_N_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::N");
    int T_N_scratch_idx;
    registerVariable(T_N_scratch_idx, d_T_N_var, cell_ghosts, getScratchContext());

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

    d_T_N_old_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::T_N_old");
    int T_N_old_current_idx, T_N_old_new_idx, T_N_old_scratch_idx;
    registerVariable(T_N_old_current_idx,
                     T_N_old_new_idx,
                     T_N_old_scratch_idx,
                     d_T_N_old_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

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

    if (d_H_var)
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

    d_chemical_potential_var = new CellVariable<NDIM, double>("chemical_potential_var");
    d_chemical_potential_idx =
        var_db->registerVariableAndContext(d_chemical_potential_var, getCurrentContext(), cell_ghosts);

    d_lf_pre_var = new CellVariable<NDIM, double>("lf_pre_var");
    d_lf_pre_idx = var_db->registerVariableAndContext(d_lf_pre_var, getCurrentContext());

    d_grad_lf_var = new SideVariable<NDIM, double>(d_object_name + "::grad_lf");
    d_grad_lf_idx =
        var_db->registerVariableAndContext(d_grad_lf_var, var_db->getContext(d_object_name + "grad_lf::SCRATCH"));

    d_H_sc_idx = var_db->registerVariableAndContext(d_grad_lf_var, var_db->getContext(d_object_name + "H_sc::SCRATCH"));

    if (d_visit_writer)
    {
        d_visit_writer->registerPlotQuantity("D", "VECTOR", d_D_cc_current_idx, 0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d == 0) d_visit_writer->registerPlotQuantity("D_x", "SCALAR", d_D_cc_current_idx, d);
            if (d == 1) d_visit_writer->registerPlotQuantity("D_y", "SCALAR", d_D_cc_current_idx, d);
            if (d == 2) d_visit_writer->registerPlotQuantity("D_z", "SCALAR", d_D_cc_current_idx, d);
        }
    }

    // Setup the convective operator.
    d_lf_convective_op = getConvectiveOperatorLiquidFractionEquation(d_lf_var);
    d_T_convective_op = getConvectiveOperatorTemperatureEquation(d_T_var);

    d_U_old_var = new FaceVariable<NDIM, double>(d_object_name + "::U_old");
    registerVariable(d_U_old_current_idx,
                     d_U_old_new_idx,
                     d_U_old_scratch_idx,
                     d_U_old_var,
                     face_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_cp_old_var = new CellVariable<NDIM, double>(d_object_name + "::cp_old");
    registerVariable(d_cp_old_current_idx,
                     d_cp_old_new_idx,
                     d_cp_old_scratch_idx,
                     d_cp_old_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_T_old_var = new CellVariable<NDIM, double>(d_object_name + "::T_old");
    registerVariable(d_T_old_current_idx,
                     d_T_old_new_idx,
                     d_T_old_scratch_idx,
                     d_T_old_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // computing AC LHS error
    d_lf_lhs_var = new CellVariable<NDIM, double>(d_object_name + "::lf_LHS");
    d_lf_lhs_idx = var_db->registerVariableAndContext(d_lf_lhs_var, getCurrentContext());

    d_lf_lhs_N_var = new CellVariable<NDIM, double>(d_lf_var->getName() + "::lhs_N");
    registerVariable(d_lf_lhs_N_scratch_idx, d_lf_lhs_N_var, cell_ghosts, getScratchContext());

    // Register variables for plotting.
    if (d_visit_writer)
    {

        if (d_output_lf)
            d_visit_writer->registerPlotQuantity("liquid_fraction", "SCALAR", d_lf_current_idx, 0);

        if (d_output_T)
            d_visit_writer->registerPlotQuantity("Temperature", "SCALAR", d_T_current_idx, 0);

        if (d_output_H) d_visit_writer->registerPlotQuantity("Heaviside", "SCALAR", d_H_current_idx, 0);
    }
    if (d_solve_mass_conservation)
    {
        // Set various objects with conservative time integrator.
        d_rho_p_integrator->setCellCenteredDensityBoundaryConditions(d_rho_bc_coef);
        d_rho_p_integrator->setCellCenteredSpecificHeatBoundaryConditions(d_rho_bc_coef);
        d_rho_p_integrator->setCellCenteredTemperatureBoundaryConditions(d_T_bc_coef);
        if (d_S_fcn) d_rho_p_integrator->setMassDensitySourceTerm(d_S_fcn);
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
        if (!level->checkAllocated(d_C_rhs_scratch_idx)) level->allocatePatchData(d_C_rhs_scratch_idx, current_time);
        if (!level->checkAllocated(d_lf_C_idx)) level->allocatePatchData(d_lf_C_idx, current_time);
        if (!level->checkAllocated(d_T_C_idx)) level->allocatePatchData(d_T_C_idx, current_time);
        if (!level->checkAllocated(d_lf_temp_rhs_idx)) level->allocatePatchData(d_lf_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_T_temp_rhs_idx)) level->allocatePatchData(d_T_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_g_firstder_idx))
            level->allocatePatchData(d_g_firstder_idx, current_time); // should this be new time?
        if (!level->checkAllocated(d_g_secondder_idx)) level->allocatePatchData(d_g_secondder_idx, current_time);
        if (!level->checkAllocated(d_p_firstder_idx)) level->allocatePatchData(d_p_firstder_idx, current_time);
        if (!level->checkAllocated(d_chemical_potential_idx))
            level->allocatePatchData(d_chemical_potential_idx, current_time);
        if (!level->checkAllocated(d_lf_pre_idx)) level->allocatePatchData(d_lf_pre_idx, current_time);
        if (!level->checkAllocated(d_grad_lf_idx)) level->allocatePatchData(d_grad_lf_idx, current_time);
        if (!level->checkAllocated(d_H_sc_idx)) level->allocatePatchData(d_H_sc_idx, current_time);
        if (!level->checkAllocated(d_lf_lhs_idx)) level->allocatePatchData(d_lf_lhs_idx, current_time);
    }

    if (d_lf_u_var)
    {
        // Update the advection velocity.
        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
        const int lf_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getScratchContext());
        const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());

        if (d_lf_u_fcn)
        {
            d_lf_u_fcn->setDataOnPatchHierarchy(lf_u_current_idx, d_lf_u_var, d_hierarchy, current_time);
            d_lf_u_fcn->setDataOnPatchHierarchy(lf_u_new_idx, d_lf_u_var, d_hierarchy, new_time);
        }
        else
        {
            d_hier_fc_data_ops->copyData(lf_u_new_idx, lf_u_current_idx);
        }
        d_hier_fc_data_ops->linearSum(lf_u_scratch_idx, 0.5, lf_u_current_idx, 0.5, lf_u_new_idx);
    }
    // Setup the operators and solvers and compute the right-hand-side terms.
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
    const int lf_N_scratch_idx = (var_db->mapVariableAndContextToIndex(d_lf_N_var, getScratchContext()));

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
    PoissonSpecifications lf_solver_spec(d_object_name + "::solver_spec::" + d_lf_var->getName());
    PoissonSpecifications lf_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_lf_var->getName());

    // set C coefficients.
    // compute heaviside based on phi^n.
    int ls_current_idx;
    for (auto Q_var : d_Q_var)
    {
        ls_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
    }
    computeHeavisideFunction(d_H_current_idx, ls_current_idx);
    computeDoubleWellPotential(d_g_firstder_idx, d_g_secondder_idx, d_lf_current_idx);

    d_hier_cc_data_ops->scale(d_lf_C_idx, 1.0 / dt, d_H_current_idx);
    d_hier_cc_data_ops->multiply(d_C_rhs_scratch_idx, d_g_secondder_idx, d_H_current_idx);
    d_hier_cc_data_ops->scale(d_C_rhs_scratch_idx, d_M_lf * d_lambda_lf / std::pow(d_eta_lf, 2.0), d_C_rhs_scratch_idx);
    d_hier_cc_data_ops->add(d_lf_C_idx, d_C_rhs_scratch_idx, d_lf_C_idx);
    lf_solver_spec.setCPatchDataId(d_lf_C_idx);

    d_hier_cc_data_ops->scale(d_lf_temp_rhs_idx, 1.0 / dt, d_H_current_idx);
    lf_rhs_op_spec.setCPatchDataId(d_lf_temp_rhs_idx);
    //    std::cout << "L2 norm of d_lf_temp_rhs_idx\t" << d_hier_cc_data_ops->L2Norm(d_lf_temp_rhs_idx) << std::endl;
    // set D coefficients.
    // Interpolate the cell-centered Heaviside to side-centered.
    d_hier_cc_data_ops->copyData(d_H_scratch_idx, d_H_current_idx);
    d_H_bdry_bc_fill_op->fillData(current_time);

    interpolateCCHeaviside(lf_diff_coef_current_idx, d_H_scratch_idx);

    d_hier_sc_data_ops->scale(lf_diff_coef_current_idx, d_M_lf * d_lambda_lf, lf_diff_coef_current_idx);

    d_hier_sc_data_ops->scale(lf_diff_coef_scratch_idx, -alpha, lf_diff_coef_current_idx);
    lf_solver_spec.setDPatchDataId(lf_diff_coef_scratch_idx);

    d_hier_sc_data_ops->scale(lf_diff_coef_rhs_scratch_idx, (1.0 - alpha), lf_diff_coef_current_idx);
    lf_rhs_op_spec.setDPatchDataId(lf_diff_coef_rhs_scratch_idx);
    //    std::cout << "L2 norm of lf_diff_coef_rhs_scratch_idx\t" <<
    //    d_hier_sc_data_ops->L2Norm(lf_diff_coef_rhs_scratch_idx)
    //              << std::endl;
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

    if (d_lf_u_var)
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
        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
        d_lf_convective_op->setAdvectionVelocity(lf_u_current_idx);
        const int lf_H_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_H_var, getScratchContext());
        /// d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_current_idx);
        d_hier_cc_data_ops->multiply(lf_H_scratch_idx, lf_current_idx, d_H_current_idx);
        d_lf_convective_op->setSolutionTime(current_time);
        d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_N_scratch_idx);

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
    d_hier_cc_data_ops->copyData(d_H_new_idx, d_H_current_idx);
    if (d_T_var) d_hier_cc_data_ops->copyData(d_T_new_idx, d_T_current_idx);

    if (d_solve_energy)
    {
        const int rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());
        const int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());

        // Note that we always reset current context of state variables here, if
        // necessary.
        const double apply_time = current_time;
        for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
        {
            d_reset_rho_fcns[k](rho_current_idx,
                                d_rho_var,
                                d_hier_math_ops,
                                -1 /*cycle_num*/,
                                apply_time,
                                current_time,
                                new_time,
                                d_reset_rho_fcns_ctx[k]);
        }

        const int Cp_current_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getCurrentContext());
        const int Cp_scratch_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getScratchContext());
        const int Cp_new_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getNewContext());
        // std::cout << "L2 norm of Cp_current_idx at cycle 0\t" << d_hier_cc_data_ops->L2Norm(Cp_current_idx) <<
        // std::endl;
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

        if (d_solve_mass_conservation)
        {
            const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());

            // Keep track of the time-lagged velocity, specific heat and temperature.
            d_hier_fc_data_ops->copyData(d_U_old_new_idx, lf_u_current_idx);
            d_hier_cc_data_ops->copyData(d_cp_old_new_idx, Cp_current_idx);
            d_hier_cc_data_ops->copyData(d_T_old_new_idx, T_current_idx);

            d_rho_p_integrator->setSolutionTime(current_time);
            d_rho_p_integrator->setTimeInterval(current_time, new_time);

            // For conservative discretization, an approximation to rho^{n+1}
            // will be computed from rho^{n}, which requires additional options to be
            // set.

            // Set the rho^{n} density
            d_rho_p_integrator->setCellCenteredDensityPatchDataIndex(rho_current_idx);
            // d_rho_p_integrator->setCellCenteredSpecificHeatPatchDataIndex(Cp_current_idx);
            // d_rho_p_integrator->setCellCenteredTemperaturePatchDataIndex(T_current_idx);

            // Set the convective derivative patch data index.
            const int T_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_N_var, getScratchContext());
            d_rho_p_integrator->setCellCenteredConvectiveDerivativePatchDataIndex(T_N_scratch_idx);

            // Data for the conservative time integrator is for cycle 0
            const int cycle_num = 0;
            d_rho_p_integrator->setCycleNumber(cycle_num);

            // Set the velocities used to update the density and the previous time step
            // size
            if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ -1, /*current*/ lf_u_current_idx, /*new*/ -1);
                d_rho_p_integrator->setSpecificHeatPatchDataIndices(
                    /*old*/ -1, /*current*/ Cp_current_idx, /*new*/ -1);
                d_rho_p_integrator->setTemperaturePatchDataIndices(
                    /*old*/ -1, /*current*/ T_current_idx, /*new*/ -1);
            }
            else
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ d_U_old_current_idx, /*current*/ lf_u_current_idx, /*new*/ -1);
                d_rho_p_integrator->setSpecificHeatPatchDataIndices(
                    /*old*/ d_cp_old_current_idx, /*current*/ Cp_current_idx, /*new*/ -1);
                d_rho_p_integrator->setTemperaturePatchDataIndices(
                    /*old*/ d_T_old_current_idx, /*current*/ T_current_idx, /*new*/ -1);
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
        PoissonSpecifications T_solver_spec(d_object_name + "::solver_spec::" + d_lf_var->getName());
        PoissonSpecifications T_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_lf_var->getName());

        // set rho*Cp/dt + K*lambda.
        const double lambda = 0.0;
        d_hier_cc_data_ops->multiply(d_C_current_idx, rho_current_idx, Cp_current_idx);
        d_hier_cc_data_ops->scale(d_C_current_idx, 1.0 / dt, d_C_current_idx);
        d_hier_cc_data_ops->copyData(d_T_C_idx, d_C_current_idx);
        T_solver_spec.setCPatchDataId(d_T_C_idx);

        //    for (unsigned k = 0; k < d_reset_Cp_fcns.size(); ++k)
        //    {
        //        d_reset_Cp_fcns[k](Cp_current_idx,
        //                           d_Cp_var,
        //                           d_hier_math_ops,
        //                           -1 /*cycle_num*/,
        //                           apply_time,
        //                           current_time,
        //                           new_time,
        //                           d_reset_Cp_fcns_ctx[k]);
        //    }

        d_hier_cc_data_ops->copyData(d_T_temp_rhs_idx, d_C_current_idx);
        T_rhs_op_spec.setCPatchDataId(d_T_temp_rhs_idx);
        //    std::cout << "L2 norm of d_T_temp_rhs_idx\t" << d_hier_cc_data_ops->L2Norm(d_T_temp_rhs_idx) << std::endl;
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

            // for plotting purpose.
            static const bool synch_cf_interface = true;
            d_hier_math_ops->interp(d_D_cc_new_idx,
                                    d_D_cc_var,
                                    T_diff_coef_current_idx,
                                    d_T_diffusion_coef_var,
                                    d_no_fill_op,
                                    d_integrator_time,
                                    synch_cf_interface);
        }

        d_hier_sc_data_ops->scale(T_diff_coef_scratch_idx, -alpha, T_diff_coef_current_idx);
        T_solver_spec.setDPatchDataId(T_diff_coef_scratch_idx);

        d_hier_sc_data_ops->scale(T_diff_coef_rhs_scratch_idx, (1.0 - alpha), T_diff_coef_current_idx);
        T_rhs_op_spec.setDPatchDataId(T_diff_coef_rhs_scratch_idx);
        //    std::cout << "L2 norm of T_diff_coef_rhs_scratch_idx\t" <<
        //    d_hier_sc_data_ops->L2Norm(T_diff_coef_rhs_scratch_idx) << std::endl;
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

        if (d_T_convective_op_needs_init)
        {
            d_T_convective_op->initializeOperatorState(*d_T_sol, *d_T_rhs);
            d_T_convective_op_needs_init = false;
        }

        d_hier_cc_data_ops->copyData(rho_new_idx, rho_current_idx);
        d_hier_cc_data_ops->copyData(Cp_new_idx, Cp_current_idx);
        d_hier_cc_data_ops->copyData(T_new_idx, T_current_idx);

        // Add div (u H lf).
        d_hier_cc_data_ops->scale(d_T_lf_N_scratch_idx, d_rho_liquid * d_latent_heat, lf_N_scratch_idx);
        if (d_lf_convective_time_stepping_type == FORWARD_EULER)
        {
            d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -1.0, d_T_lf_N_scratch_idx, T_rhs_scratch_idx);
            //            std::cout << "L2 norm of lf_N_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(lf_N_scratch_idx)
            //            << std::endl;
        }
        else if (d_lf_convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -0.5, d_T_lf_N_scratch_idx, T_rhs_scratch_idx);
        }
    }
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
    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
    const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
    const int lf_current_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getCurrentContext());
    const int lf_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_rhs_var, getScratchContext());
    const int lf_diff_coef_new_idx = (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_var, getNewContext()));
    const int lf_diff_coef_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_var, getScratchContext()));
    const int lf_diff_coef_rhs_scratch_idx =
        (var_db->mapVariableAndContextToIndex(d_lf_diffusion_coef_rhs_var, getScratchContext()));
    const int lf_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_N_var, getScratchContext());

    // compute heaviside based on phi^n+1.
    int ls_new_idx;
    for (auto Q_var : d_Q_var)
    {
        ls_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
    }
    computeHeavisideFunction(d_H_new_idx, ls_new_idx);

    // update C coefficients.
    if (cycle_num > 0) computeDoubleWellPotential(d_g_firstder_idx, d_g_secondder_idx, d_lf_new_idx);

    d_hier_cc_data_ops->scale(d_lf_C_idx, 1.0 / dt, d_H_new_idx);
    d_hier_cc_data_ops->multiply(d_C_rhs_scratch_idx, d_g_secondder_idx, d_H_new_idx);
    d_hier_cc_data_ops->scale(d_C_rhs_scratch_idx, d_M_lf * d_lambda_lf / std::pow(d_eta_lf, 2.0), d_C_rhs_scratch_idx);
    d_hier_cc_data_ops->add(d_lf_C_idx, d_C_rhs_scratch_idx, d_lf_C_idx);
    //    std::cout << "L2 norm of d_lf_C_idx\t" << d_hier_cc_data_ops->L2Norm(d_lf_C_idx) << std::endl;
    // update D coefficients.
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

    // Interpolate the cell-centered Heaviside to side-centered.
    d_hier_cc_data_ops->copyData(d_H_scratch_idx, d_H_new_idx);
    d_H_bdry_bc_fill_op->fillData(new_time);

    interpolateCCHeaviside(lf_diff_coef_new_idx, d_H_scratch_idx);
    d_hier_sc_data_ops->copyData(d_H_sc_idx, lf_diff_coef_new_idx);
    //    d_hier_math_ops->interp(lf_diff_coef_new_idx,
    //                            d_lf_diffusion_coef_var,
    //                            true,
    //                            d_H_scratch_idx,
    //                            d_H_var,
    //                            d_no_fill_op,
    //                            d_integrator_time);

    d_hier_sc_data_ops->scale(lf_diff_coef_new_idx, d_M_lf * d_lambda_lf, lf_diff_coef_new_idx);
    d_hier_sc_data_ops->scale(lf_diff_coef_scratch_idx, -alpha, lf_diff_coef_new_idx);
    //    std::cout << "L2 norm of lf_diff_coef_scratch_idx\t" << d_hier_sc_data_ops->L2Norm(lf_diff_coef_scratch_idx)
    //              << std::endl;
    if (cycle_num > 0)
    {
        // Update the advection velocity for lf.
        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
        const int lf_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getScratchContext());
        const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
        if (d_lf_u_fcn)
        {
            d_lf_u_fcn->setDataOnPatchHierarchy(lf_u_new_idx, d_lf_u_var, d_hierarchy, new_time);
        }
        d_hier_fc_data_ops->linearSum(lf_u_scratch_idx, 0.5, lf_u_current_idx, 0.5, lf_u_new_idx);
    }
    //    std::cout << "L2 norm of lf_rhs_scratch_idx before convective term\t"
    //              << d_hier_cc_data_ops->L2Norm(lf_rhs_scratch_idx) << std::endl;

    // Account for the convective difference term.
    const int lf_H_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_H_var, getScratchContext());
    TimeSteppingType convective_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
    if (d_lf_u_var)
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
                    pout << "IEPSemiImplicitHierarchyIntegrator::"
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
                const int lf_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getScratchContext());
                d_lf_convective_op->setAdvectionVelocity(lf_u_scratch_idx);
                d_hier_cc_data_ops->linearSum(lf_scratch_idx, 0.5, lf_current_idx, 0.5, lf_new_idx);
                d_hier_cc_data_ops->linearSum(d_H_scratch_idx, 0.5, d_H_current_idx, 0.5, d_H_new_idx);
                d_hier_cc_data_ops->multiply(lf_H_scratch_idx, d_H_scratch_idx, lf_scratch_idx);
                d_lf_convective_op->setSolutionTime(half_time);
                d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_N_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
                d_lf_convective_op->setAdvectionVelocity(lf_u_new_idx);
                // d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_new_idx);
                d_hier_cc_data_ops->multiply(lf_H_scratch_idx, lf_new_idx, d_H_new_idx);
                d_lf_convective_op->setSolutionTime(new_time);
                d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_N_scratch_idx);
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
        //        std::cout << "L2 norm of lf_N_scratch_idx after convection term\t"
        //                  << d_hier_cc_data_ops->L2Norm(lf_rhs_scratch_idx) << std::endl;
    }

    // Account for forcing terms.
    const int lf_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_F_var, getScratchContext());
    const int lf_F_new_idx = var_db->mapVariableAndContextToIndex(d_lf_F_var, getNewContext());

    if (d_lf_F_var)
    {
        computeInterpolationFunction(d_p_firstder_idx, d_lf_new_idx, d_T_new_idx);
        computeLiquidFractionSourceTerm(lf_F_scratch_idx);

        if (d_lf_u_var)
        {
            int div_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_div_u_var, getScratchContext());
            const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
            d_hier_math_ops->div(div_u_scratch_idx,
                                 d_div_u_var,
                                 1.0,
                                 lf_u_new_idx,
                                 d_lf_u_var,
                                 d_no_fill_op,
                                 d_integrator_time,
                                 /*synch_cf_bdry*/ false);

            d_hier_cc_data_ops->multiply(lf_H_scratch_idx, lf_new_idx, d_H_new_idx);
            d_hier_cc_data_ops->multiply(div_u_scratch_idx, lf_H_scratch_idx, div_u_scratch_idx);
            d_hier_cc_data_ops->axpy(lf_F_scratch_idx, +1.0, div_u_scratch_idx, lf_F_scratch_idx);
        }
        //        std::cout << "L2 norm of lf_F_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(lf_F_scratch_idx) <<
        //        std::endl;
        d_hier_cc_data_ops->axpy(lf_rhs_scratch_idx, +1.0, lf_F_scratch_idx, lf_rhs_scratch_idx);
        //        std::cout << "L2 norm of lf_rhs_scratch_idx after all contribution\t"
        //                  << d_hier_cc_data_ops->L2Norm(lf_rhs_scratch_idx) << std::endl;
    }

    // This will be used while computing chemical potential.
    d_hier_cc_data_ops->copyData(d_lf_pre_idx, lf_new_idx);

    // Solve for lf(n+1).
    Pointer<PoissonSolver> lf_solver = d_lf_solver;
    lf_solver->solveSystem(*d_lf_sol, *d_lf_rhs);
    d_hier_cc_data_ops->copyData(lf_new_idx, lf_scratch_idx);
    //    std::cout << "L2 norm of lf_new_idx\t" << d_hier_cc_data_ops->L2Norm(lf_new_idx) << std::endl;
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
    if (d_lf_u_var)
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

    // compute LHS of AC equation using lf^n+1,m+1
    computeLHSOfLiquidFractionEquation(
        d_lf_lhs_idx, lf_N_scratch_idx, dt, cycle_num, new_time, current_time, half_time);
    // Compute chemical potential using lf^n+1.
    computeChemicalPotential(d_chemical_potential_idx, d_H_sc_idx, new_time);

    if (d_solve_energy)
    {
        int rho_new_idx, rho_scratch_idx, rho_current_idx;
        rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());
        rho_scratch_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getScratchContext());
        rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());

        int Cp_new_idx, Cp_scratch_idx, Cp_current_idx;
        Cp_new_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getNewContext());
        Cp_scratch_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getScratchContext());
        Cp_current_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getCurrentContext());

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
        // Update N_idx if necessary
        if (cycle_num > 0)
        {
            const double dt = new_time - current_time;
            const double half_time = current_time + 0.5 * dt;
            d_rho_p_integrator->setSolutionTime(half_time);

            // Set the cycle number
            d_rho_p_integrator->setCycleNumber(cycle_num);

            // Set the patch data index for convective derivative.
            d_rho_p_integrator->setCellCenteredConvectiveDerivativePatchDataIndex(T_N_scratch_idx);

            // Always set to current because we want to update rho^{n} to rho^{n+1}
            d_rho_p_integrator->setCellCenteredDensityPatchDataIndex(rho_current_idx);

            const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
            const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());

            // Set the velocities used to update the density
            if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ -1, /*current*/ lf_u_current_idx, /*new*/ lf_u_new_idx);
                d_rho_p_integrator->setSpecificHeatPatchDataIndices(
                    /*old*/ -1, /*current*/ Cp_current_idx, /*new*/ Cp_new_idx);
                d_rho_p_integrator->setTemperaturePatchDataIndices(
                    /*old*/ -1, /*current*/ T_current_idx, /*new*/ T_new_idx);
            }
            else
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ d_U_old_current_idx,
                    /*current*/ lf_u_current_idx,
                    /*new*/ lf_u_new_idx);
                d_rho_p_integrator->setSpecificHeatPatchDataIndices(
                    /*old*/ d_cp_old_current_idx,
                    /*current*/ Cp_current_idx,
                    /*new*/ Cp_new_idx);
                d_rho_p_integrator->setTemperaturePatchDataIndices(
                    /*old*/ d_T_old_current_idx,
                    /*current*/ T_current_idx,
                    /*new*/ T_new_idx);

                d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
            }

            d_rho_p_integrator->integrate(dt);
        }

        const int rho_cc_new_idx = d_rho_p_integrator->getUpdatedCellCenteredDensityPatchDataIndex();
        d_hier_cc_data_ops->copyData(rho_new_idx,
                                     rho_cc_new_idx,
                                     /*interior_only*/ true);

        // set rho*Cp/dt + K*lambda.
        const double lambda = 0.0;
        d_hier_cc_data_ops->multiply(d_C_new_idx, rho_new_idx, Cp_new_idx);
        d_hier_cc_data_ops->scale(d_C_new_idx, 1.0 / dt, d_C_new_idx);
        d_hier_cc_data_ops->copyData(d_T_C_idx, d_C_new_idx);
        //       std::cout << "L2 norm of d_T_C_idx\t" << d_hier_cc_data_ops->L2Norm(d_T_C_idx) << std::endl;

        if (cycle_num > 0)
        {
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
                d_reset_kappa_fcns[k](T_diff_coef_new_idx,
                                      d_T_diffusion_coef_var,
                                      d_hier_math_ops,
                                      -1 /*cycle_num*/,
                                      apply_time,
                                      current_time,
                                      new_time,
                                      d_reset_kappa_fcns_ctx[k]);

                // for plotting purpose.
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
            //            std::cout << "L2 norm of T_diff_coef_scratch_idx\t" <<
            //            d_hier_sc_data_ops->L2Norm(T_diff_coef_scratch_idx) << std::endl;
        }

        //        std::cout << "L2 norm of T_rhs_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(T_rhs_scratch_idx) <<
        //        std::endl; std::cout << "L2 norm of T_N_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(T_N_scratch_idx)
        //        << std::endl;
        // Account for the convective acceleration term N_full.
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

        // Account for forcing terms.
        const int T_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_F_var, getScratchContext());
        const int T_F_new_idx = var_db->mapVariableAndContextToIndex(d_T_F_var, getNewContext());

        if (d_T_F_var)
        {
            computeTemperatureSourceTerm(T_F_scratch_idx, dt);
            d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, +1.0, T_F_scratch_idx, T_rhs_scratch_idx);
            //            std::cout << "L2 norm of T_rhs_scratch_idx after F\t" <<
            //            d_hier_cc_data_ops->L2Norm(T_rhs_scratch_idx)
            //                      << std::endl;
        }

        // Solve for T(n+1).
        Pointer<PoissonSolver> T_solver = d_T_solver;
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
        level->deallocatePatchData(d_g_firstder_idx);
        level->deallocatePatchData(d_g_secondder_idx);
        level->deallocatePatchData(d_p_firstder_idx);
        level->deallocatePatchData(d_chemical_potential_idx);
        level->deallocatePatchData(d_lf_pre_idx);
        level->deallocatePatchData(d_grad_lf_idx);
        level->deallocatePatchData(d_H_sc_idx);
        level->deallocatePatchData(d_lf_lhs_idx);
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

    d_lf_convective_op_needs_init = true;
    d_T_convective_op_needs_init = true;
    if (d_solve_mass_conservation)
    {
        d_rho_p_integrator->setHierarchyMathOps(d_hier_math_ops);
        d_rho_p_integrator->initializeTimeIntegrator(base_hierarchy);
    }
    return;
}

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
    d_lf_u_var = nullptr;
    d_lf_F_var = lf_F_var;
    d_lf_rhs_var = lf_rhs_var;
    d_lf_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    d_lf_convective_difference_form = d_default_convective_difference_form;
    if (!d_lf_convective_time_stepping_type)
        d_lf_convective_time_stepping_type = d_default_convective_time_stepping_type;
    d_lf_convective_op_type = d_default_convective_op_type;
    d_lf_convective_op_input_db = d_default_convective_op_input_db;
    d_lf_diffusion_coef_var = lf_diff_coef_var;
    d_lf_init = nullptr;
    d_lf_F_fcn = nullptr;
    d_lf_bc_coef = nullptr;
    return;
} // registerLiquidFractionVariable

void
IEPSemiImplicitHierarchyIntegrator::registerHeavisideVariable(Pointer<CellVariable<NDIM, double> > H_var,
                                                              const bool output_H_var)
{
    d_H_var = H_var;
    d_output_H = output_H_var;
    return;
} // registerHeavisideVariable

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
    d_T_u_var = nullptr;
    d_T_F_var = T_F_var;
    d_T_rhs_var = T_rhs_var;
    d_T_diffusion_time_stepping_type = d_default_diffusion_time_stepping_type;
    d_T_convective_difference_form = d_default_convective_difference_form;
    if (!d_T_convective_time_stepping_type) d_T_convective_time_stepping_type = d_default_convective_time_stepping_type;
    d_T_convective_op_type = d_default_convective_op_type;
    d_T_convective_op_input_db = d_default_convective_op_input_db;
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

Pointer<ConvectiveOperator>
IEPSemiImplicitHierarchyIntegrator::getConvectiveOperatorLiquidFractionEquation(
    Pointer<CellVariable<NDIM, double> > lf_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_lf_var);
#endif
    if (!d_lf_convective_op)
    {
        std::vector<RobinBcCoefStrategy<NDIM>*> lf_bc_coefs(1, d_lf_bc_coef);
        AdvDiffConvectiveOperatorManager* lf_convective_op_manager = AdvDiffConvectiveOperatorManager::getManager();
        d_lf_convective_op = lf_convective_op_manager->allocateOperator(d_lf_convective_op_type,
                                                                        d_object_name + "::lfConvectiveOperator",
                                                                        d_lf_var,
                                                                        d_lf_convective_op_input_db,
                                                                        d_lf_convective_difference_form,
                                                                        lf_bc_coefs);
        d_lf_convective_op_needs_init = true;
    }
    return d_lf_convective_op;
} // getConvectiveOperatorLiquidFractionEquation

Pointer<ConvectiveOperator>
IEPSemiImplicitHierarchyIntegrator::getConvectiveOperatorTemperatureEquation(Pointer<CellVariable<NDIM, double> > T_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_T_var);
#endif
    if (!d_T_convective_op)
    {
        std::vector<RobinBcCoefStrategy<NDIM>*> T_bc_coefs(1, d_T_bc_coef);
        AdvDiffConvectiveOperatorManager* T_convective_op_manager = AdvDiffConvectiveOperatorManager::getManager();
        d_T_convective_op = T_convective_op_manager->allocateOperator(d_T_convective_op_type,
                                                                      d_object_name + "::WConvectiveOperator",
                                                                      d_T_var,
                                                                      d_T_convective_op_input_db,
                                                                      d_T_convective_difference_form,
                                                                      T_bc_coefs);
        d_T_convective_op_needs_init = true;
    }
    return d_T_convective_op;
} // getConvectiveOperatorTemperatureEquation

int
IEPSemiImplicitHierarchyIntegrator::getChemicalPotentialIndex()
{
    return d_chemical_potential_idx;
} // getChemicalPotentialIndex

void
IEPSemiImplicitHierarchyIntegrator::registerMassDensityBoundaryConditions(RobinBcCoefStrategy<NDIM>*& rho_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_bc_coef = rho_bc_coef;
    return;
} // registerMassDensityBoundaryConditions

void
IEPSemiImplicitHierarchyIntegrator::registerMassDensitySourceTerm(Pointer<CartGridFunction> S_fcn)
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
        TBOX_ERROR(d_object_name << "::IEPSemiImplicitHierarchyIntegrator():\n"
                                 << " present implementation allows for only one mass density source\n"
                                 << " term to be set. Consider combining source terms into single "
                                    "CartGridFunction.\n");
    }
    return;
} // registerMassDensitySourceTerm

void
IEPSemiImplicitHierarchyIntegrator::setAdvectionVelocityLiquidFractionEquation(
    Pointer<CellVariable<NDIM, double> > lf_var,
    Pointer<FaceVariable<NDIM, double> > u_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_lf_var);
#endif
    d_lf_u_var = u_var;

    return;
} // setAdvectionVelocityLiquidFractionEquation
void
IEPSemiImplicitHierarchyIntegrator::setAdvectionVelocityTemperatureEquation(Pointer<CellVariable<NDIM, double> > T_var,
                                                                            Pointer<FaceVariable<NDIM, double> > u_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_T_var);
#endif
    d_T_u_var = u_var;

    return;
} // setAdvectionVelocityTemperatureEquation

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
                }
                // std::cout << "p' value is\t" << (*p_firstder_data)(ci) <<  std::endl;
            }
        }
    }
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
IEPSemiImplicitHierarchyIntegrator::computeChemicalPotential(int chemical_potential_idx,
                                                             const int H_sc_idx,
                                                             const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
    int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
    int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());
    int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());

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
    d_hier_sc_data_ops->multiply(d_grad_lf_idx, d_grad_lf_idx, H_sc_idx);

    // compute div(H*grad_lf).
    d_hier_math_ops->div(
        chemical_potential_idx, d_chemical_potential_var, 1.0, d_grad_lf_idx, d_grad_lf_var, nullptr, new_time, false);

    // update p' and g'.
    // Directly using p' calculated using var_phi^n+1,m and g'.
    // computeInterpolationFunction(d_p_firstder_idx, lf_new_idx, T_new_idx);
    // computeDoubleWellPotential(d_g_firstder_idx, d_g_secondder_idx, lf_new_idx);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_new_idx);
            Pointer<CellData<NDIM, double> > lf_pre_data = patch->getPatchData(d_lf_pre_idx);
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_new_idx);
            Pointer<CellData<NDIM, double> > p_firstder_data = patch->getPatchData(d_p_firstder_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double> > g_firstder_data = patch->getPatchData(d_g_firstder_idx);
            Pointer<CellData<NDIM, double> > g_secondder_data = patch->getPatchData(d_g_secondder_idx);
            Pointer<CellData<NDIM, double> > chemical_potential_data = patch->getPatchData(chemical_potential_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                double F = d_rho_liquid * d_latent_heat * (*H_data)(ci) * (*p_firstder_data)(ci) *
                           (d_T_ref - (*T_data)(ci)) / d_T_ref;

                double g_first_der_linear =
                    (*g_firstder_data)(ci) + (*g_secondder_data)(ci) * ((*lf_data)(ci) - (*lf_pre_data)(ci));

                (*chemical_potential_data)(ci) =
                    -d_lambda_lf * (*chemical_potential_data)(ci) +
                    (d_lambda_lf * (*H_data)(ci) / std::pow(d_eta_lf, 2.0) * g_first_der_linear) + F;
            }
        }
    }

    return;
} // computeChemicalPotential

void
IEPSemiImplicitHierarchyIntegrator::computeLHSOfLiquidFractionEquation(int lf_lhs_idx,
                                                                       const int lf_N_scratch_idx,
                                                                       const double dt,
                                                                       const int cycle_num,
                                                                       const double new_time,
                                                                       const double current_time,
                                                                       const double half_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

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
            Pointer<CellData<NDIM, double> > lf_lhs_data = patch->getPatchData(lf_lhs_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*lf_lhs_data)(ci) =
                    (((*H_new_data)(ci) * (*lf_new_data)(ci)) - ((*H_current_data)(ci) * (*lf_current_data)(ci))) / dt;
            }
        }
    }
    d_hier_cc_data_ops->axpy(lf_N_scratch_idx, +1.0, lf_lhs_idx, lf_N_scratch_idx);

    const int lf_lhs_N_scratch_idx = (var_db->mapVariableAndContextToIndex(d_lf_lhs_N_var, getScratchContext()));
    if (d_lf_u_var)
    {
        if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(d_lf_convective_time_stepping_type))
        {
            d_lf_convective_time_stepping_type = d_lf_init_convective_time_stepping_type;
        }

        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
        d_lf_convective_op->setAdvectionVelocity(lf_u_current_idx);
        const int lf_H_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_H_var, getScratchContext());
        /// d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_current_idx);
        d_hier_cc_data_ops->multiply(lf_H_scratch_idx, d_lf_current_idx, d_H_current_idx);
        d_lf_convective_op->setSolutionTime(current_time);
        d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_lhs_N_scratch_idx);

        const int lf_N_old_new_idx = var_db->mapVariableAndContextToIndex(d_lf_N_old_var, getNewContext());
        d_hier_cc_data_ops->copyData(lf_N_old_new_idx, lf_lhs_N_scratch_idx);

        if (d_lf_convective_time_stepping_type == FORWARD_EULER)
        {
            d_hier_cc_data_ops->axpy(lf_lhs_idx, 1.0, lf_lhs_N_scratch_idx, lf_lhs_idx);
        }
        else if (d_lf_convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(lf_lhs_idx, 0.5, lf_lhs_N_scratch_idx, lf_lhs_idx);
        }
    }

    if (cycle_num > 0)
    {
        // Update the advection velocity for lf.
        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
        const int lf_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getScratchContext());
        const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
        if (d_lf_u_fcn)
        {
            d_lf_u_fcn->setDataOnPatchHierarchy(lf_u_new_idx, d_lf_u_var, d_hierarchy, new_time);
        }
        d_hier_fc_data_ops->linearSum(lf_u_scratch_idx, 0.5, lf_u_current_idx, 0.5, lf_u_new_idx);
    }

    // Account for the convective difference term.
    const int lf_H_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_H_var, getScratchContext());
    TimeSteppingType convective_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
    if (d_lf_u_var)
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
                    pout << "IEPSemiImplicitHierarchyIntegrator::"
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
                const int lf_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getScratchContext());
                d_lf_convective_op->setAdvectionVelocity(lf_u_scratch_idx);
                d_hier_cc_data_ops->linearSum(d_lf_scratch_idx, 0.5, d_lf_current_idx, 0.5, d_lf_new_idx);
                d_hier_cc_data_ops->linearSum(d_H_scratch_idx, 0.5, d_H_current_idx, 0.5, d_H_new_idx);
                d_hier_cc_data_ops->multiply(lf_H_scratch_idx, d_H_scratch_idx, d_lf_scratch_idx);
                d_lf_convective_op->setSolutionTime(half_time);
                d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_lhs_N_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
                d_lf_convective_op->setAdvectionVelocity(lf_u_new_idx);
                // d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_new_idx);
                d_hier_cc_data_ops->multiply(lf_H_scratch_idx, d_lf_new_idx, d_H_new_idx);
                d_lf_convective_op->setSolutionTime(new_time);
                d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_lhs_N_scratch_idx);
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
                lf_lhs_N_scratch_idx, 1.0 + 0.5 * omega, lf_lhs_N_scratch_idx, -0.5 * omega, lf_N_old_current_idx);
        }

        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
        {
            d_hier_cc_data_ops->axpy(lf_lhs_idx, 1.0, lf_lhs_N_scratch_idx, lf_lhs_idx);
        }
        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
        {
            d_hier_cc_data_ops->axpy(lf_lhs_idx, 0.5, lf_lhs_N_scratch_idx, lf_lhs_idx);
        }
        std::cout << "L2 norm of lf_lhs_idx\t" << d_hier_cc_data_ops->L2Norm(lf_lhs_idx) << std::endl;
        std::cout << "L2 norm of lf_N_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(lf_N_scratch_idx) << std::endl;
    }

    return;
} // computeLHSOfLiquidFractionEquation

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
        if (input_db->keyExists("rho_liquid")) d_rho_liquid = input_db->getDouble("rho_liquid");
        if (input_db->keyExists("T_ref")) d_T_ref = input_db->getDouble("T_ref");

        if (input_db->keyExists("num_interface_cells"))
            d_num_interface_cells = input_db->getInteger("num_interface_cells");

        if (input_db->keyExists("M_lf")) d_M_lf = input_db->getDouble("M_lf");
        if (input_db->keyExists("lambda_lf")) d_lambda_lf = input_db->getDouble("lambda_lf");
        if (input_db->keyExists("eta_lf")) d_eta_lf = input_db->getDouble("eta_lf");
        if (input_db->keyExists("solve_energy")) d_solve_energy = input_db->getBool("solve_energy");
        if (input_db->keyExists("solve_mass_conservation"))
            d_solve_mass_conservation = input_db->getBool("solve_mass_conservation");

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

        if (input_db->keyExists("T_convective_op_type"))
            d_T_convective_op_type = input_db->getString("T_convective_op_type");
        else if (input_db->keyExists("T_convective_operator_type"))
            d_T_convective_op_type = input_db->getString("T_convective_operator_type");
        else if (input_db->keyExists("default_T_convective_op_type"))
            d_T_convective_op_type = input_db->getString("default_T_convective_op_type");
        else if (input_db->keyExists("default_T_convective_operator_type"))
            d_T_convective_op_type = input_db->getString("default_T_convective_operator_type");

        if (input_db->keyExists("T_convective_op_db"))
            d_T_convective_op_input_db = input_db->getDatabase("T_convective_op_db");
        else if (input_db->keyExists("default_T_convective_op_db"))
            d_T_convective_op_input_db = input_db->getDatabase("default_T_convective_op_db");

        if (input_db->keyExists("T_convective_difference_form"))
            d_T_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("T_convective_difference_form"));
        else if (input_db->keyExists("T_convective_difference_type"))
            d_T_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("T_convective_difference_type"));
        else if (input_db->keyExists("default_T_convective_difference_form"))
            d_T_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("default_T_convective_difference_form"));
        else if (input_db->keyExists("default_T_convective_difference_type"))
            d_T_convective_difference_form =
                string_to_enum<ConvectiveDifferencingType>(input_db->getString("default_T_convective_difference_type"));
        if (input_db->keyExists("lf_convective_time_stepping_type"))
            d_lf_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("lf_convective_time_stepping_type"));
        if (input_db->keyExists("T_convective_time_stepping_type"))
            d_T_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("T_convective_time_stepping_type"));
    }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
