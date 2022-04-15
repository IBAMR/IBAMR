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

#include "ibamr/AdvDiffConservativeCUIConvectiveOperator.h"
#include "ibamr/AdvDiffConservativeMassTransportQuantityIntegrator.h"
#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/IEPSemiImplicitHierarchyIntegrator.h"
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

// FORTRAN ROUTINES
#if (NDIM == 2)
#define C_TO_S_CWISE_INTERP_FC IBTK_FC_FUNC(ctoscwiseinterp2nd2d, CTOSCWISEINTERP2ND2D)
#define C_TO_S_CWISE_HARMONIC_INTERP_FC IBTK_FC_FUNC(ctoscwiseharmonicinterp2nd2d, CTOSCWISEARMONICINTERP2ND2D)
#endif // if (NDIM == 2)

#if (NDIM == 3)
#define C_TO_S_CWISE_INTERP_FC IBTK_FC_FUNC(ctoscwiseinterp2nd3d, CTOSCWISEINTERP2ND3D)
#define C_TO_S_CWISE_HARMONIC_INTERP_FC IBTK_FC_FUNC(ctoscwiseharmonicinterp2nd3d, CTOSCWISEHARMONICINTERP2ND3D)
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

    void C_TO_S_CWISE_HARMONIC_INTERP_FC(double* u0,
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
// Version of INSHierarchyIntegrator restart file data.
static const int IEP_HIERARCHY_INTEGRATOR_VERSION = 4;

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

    if (from_restart) getFromRestart();
    getFromInput(input_db, from_restart);

    if (!(d_lf_convective_difference_form == CONSERVATIVE && d_T_convective_difference_form == CONSERVATIVE))
    {
        TBOX_ERROR(d_object_name << "::IEPSemiImplicitHierarchyIntegrator():\n"
                                 << " variable coefficient discretization\n"
                                 << " requires CONSERVATIVE convective difference form\n");
    }

    // Initialize conservative mass and transported quantity integrator.
    if (d_solve_mass_conservation)
        d_rho_p_integrator = new AdvDiffConservativeMassTransportQuantityIntegrator(
            "AdvDiffConservativeMassTransportQuantityIntegrator::MassTransportIntegrator",
            input_db->getDatabase("mass_transport_integrator_db"));

    if (!(d_lf_convective_op_type == "CONSERVATIVE_CUI"))
    {
        TBOX_ERROR(d_object_name << "::IEPSemiImplicitHierarchyIntegrator():\n"
                                 << " current implementation supports only\n"
                                 << " CONSERVATIVE CUI limiter for Allen-Cahn equation\n");
    }

    if (!(d_T_convective_op_type == "CUI" || d_T_convective_op_type == "PPM"))
    {
        TBOX_ERROR(d_object_name << "::IEPSemiImplicitHierarchyIntegrator():\n"
                                 << " current implementation supports only\n"
                                 << " CUI and PPM convective limiters for energy equation\n");
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

    // Get the interpolation type for the material properties
    if (input_db->keyExists("vc_interpolation_type"))
    {
        d_k_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("vc_interpolation_type"));
    }
    if (input_db->keyExists("k_vc_interpolation_type"))
    {
        d_k_vc_interp_type = IBTK::string_to_enum<VCInterpType>(input_db->getString("k_vc_interpolation_type"));
    }
    switch (d_k_vc_interp_type)
    {
    case VC_HARMONIC_INTERP:
    case VC_AVERAGE_INTERP:
        break;
    default:
        TBOX_ERROR(d_object_name << "::IEPSemiImplicitHierarchyIntegrator():\n"
                                 << "  unsupported conductivity interpolation type: "
                                 << IBTK::enum_to_string<VCInterpType>(d_k_vc_interp_type) << " \n"
                                 << "  valid choices are: VC_HARMONIC_INTERP, VC_AVERAGE_INTERP\n");
    }

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

    registerVariable(d_h_current_idx,
                     d_h_new_idx,
                     d_h_scratch_idx,
                     d_h_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");
    if (d_visit_writer) d_visit_writer->registerPlotQuantity(d_h_var->getName(), "SCALAR", d_h_current_idx);

    int lf_F_current_idx, lf_F_scratch_idx, lf_F_new_idx;
    registerVariable(lf_F_current_idx,
                     lf_F_new_idx,
                     lf_F_scratch_idx,
                     d_lf_F_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    int T_F_current_idx, T_F_scratch_idx, T_F_new_idx;
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

    if (d_T_diffusion_coef_cc_var)
        registerVariable(d_T_diff_coef_cc_current_idx,
                         d_T_diff_coef_cc_new_idx,
                         d_T_diff_coef_cc_scratch_idx,
                         d_T_diffusion_coef_cc_var,
                         cell_ghosts,
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

    int T_diffusion_coef_rhs_scratch_idx;
    d_T_diffusion_coef_rhs_var = new SideVariable<NDIM, double>(d_T_var->getName() + "::Diff");
    registerVariable(T_diffusion_coef_rhs_scratch_idx, d_T_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    int T_rhs_scratch_idx, lf_rhs_scratch_idx;
    registerVariable(T_rhs_scratch_idx, d_T_rhs_var, cell_ghosts, getScratchContext());
    registerVariable(lf_rhs_scratch_idx, d_lf_rhs_var, cell_ghosts, getScratchContext());

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_T_C_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::C");
    d_T_C_idx = var_db->registerVariableAndContext(d_T_C_var, getCurrentContext(), no_ghosts);

    d_H_pre_var = new CellVariable<NDIM, double>(d_object_name + "::H_pre");
    d_H_pre_idx = var_db->registerVariableAndContext(d_H_pre_var, getCurrentContext(), no_ghosts);

    d_T_temp_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::temp");
    d_T_temp_idx = var_db->registerVariableAndContext(d_T_temp_var, getCurrentContext(), no_ghosts);

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

    d_D_cc_var = new CellVariable<NDIM, double>("D_cc", NDIM);
    registerVariable(d_D_cc_current_idx,
                     d_D_cc_new_idx,
                     d_D_cc_scratch_idx,
                     d_D_cc_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_C_var = new CellVariable<NDIM, double>("C_var");
    registerVariable(d_C_current_idx,
                     d_C_new_idx,
                     d_C_scratch_idx,
                     d_C_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_lf_material_derivative_var = new CellVariable<NDIM, double>(d_object_name + "::lf_material_derivative_var");
    registerVariable(d_lf_material_derivative_idx, d_lf_material_derivative_var, no_ghosts, getCurrentContext());

    d_div_U_F_var = new CellVariable<NDIM, double>(d_object_name + "::div_U_F_var");
    registerVariable(d_div_U_F_idx, d_div_U_F_var, no_ghosts, getCurrentContext());

    d_div_U_F_deno_var = new CellVariable<NDIM, double>(d_object_name + "::div_U_F_deno_var");
    registerVariable(d_div_U_F_deno_idx, d_div_U_F_deno_var, no_ghosts, getCurrentContext());

    d_div_U_F_diff_var = new CellVariable<NDIM, double>(d_object_name + "::div_U_F_diff_var");
    registerVariable(d_div_U_F_diff_idx, d_div_U_F_diff_var, no_ghosts, getCurrentContext());
    //
    d_lf_diffusion_var = new CellVariable<NDIM, double>(d_object_name + "::liquid_fraction_diffusion_var");
    registerVariable(d_lf_diffusion_idx, d_lf_diffusion_var, no_ghosts, getCurrentContext());

    if (d_visit_writer && d_solve_mass_conservation)
    {
        d_visit_writer->registerPlotQuantity("zeta", "SCALAR", d_div_U_F_idx); // In enthalpy porosity
                                                                               // zeta = H* D\varphi/Dt
        d_visit_writer->registerPlotQuantity("div_U_F_deno", "SCALAR", d_div_U_F_deno_idx);
        d_visit_writer->registerPlotQuantity("div_U_F_diff", "SCALAR", d_div_U_F_diff_idx);
    }

    d_lf_pre_var = new CellVariable<NDIM, double>(d_object_name + "::lf_pre_var");
    d_lf_pre_idx = var_db->registerVariableAndContext(d_lf_pre_var, getCurrentContext());

    d_T_pre_var = new CellVariable<NDIM, double>(d_object_name + "::T_pre_var");
    d_T_pre_idx = var_db->registerVariableAndContext(d_T_pre_var, getCurrentContext());

    d_dh_dT_var = new CellVariable<NDIM, double>(d_object_name + "::dh_dT_var");
    d_dh_dT_scratch_idx = var_db->registerVariableAndContext(d_dh_dT_var, getCurrentContext(), no_ghosts);

    d_drho_dT_var = new CellVariable<NDIM, double>(d_object_name + "::drho_dT_var");
    d_drho_dT_scratch_idx = var_db->registerVariableAndContext(d_drho_dT_var, getCurrentContext(), no_ghosts);

    d_grad_T_var = new SideVariable<NDIM, double>(d_object_name + "::grad_T");
    d_grad_T_idx =
        var_db->registerVariableAndContext(d_grad_T_var, var_db->getContext(d_object_name + "grad_T::SCRATCH"));

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
    d_nonconser_convective_op = getNonConservativeConvectiveOperator(d_lf_var, d_lf_bc_coef);
    Pointer<AdvDiffConservativeCUIConvectiveOperator> lf_conser_convective_op = d_lf_convective_op;
    lf_conser_convective_op->setHeavisideVariable(d_H_var);
    // std::vector<RobinBcCoefStrategy<NDIM>*> H_bc_coef = getPhysicalBcCoefs(d_H_var);
    // lf_conser_convective_op->setHeavisideBoundaryConditions(d_H_bc_coef);
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

    d_h_old_var = new CellVariable<NDIM, double>(d_object_name + "::h_old");
    registerVariable(d_h_old_current_idx,
                     d_h_old_new_idx,
                     d_h_old_scratch_idx,
                     d_h_old_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // computing magnitude of mass conservation equation pointwise
    d_M_var = new CellVariable<NDIM, double>(d_object_name + "::mass_conservation");
    registerVariable(d_M_idx, d_M_var, no_ghosts, getCurrentContext());

    // computing magnitude of mass conservation equation pointwise
    d_updated_rho_var = new CellVariable<NDIM, double>(d_object_name + "::updated_rho");
    registerVariable(d_updated_rho_cc_idx, d_updated_rho_var, no_ghosts, getCurrentContext());

    // Register variables for plotting.
    if (d_visit_writer)
    {

        if (d_output_lf)
            d_visit_writer->registerPlotQuantity("liquid_fraction", "SCALAR", d_lf_current_idx, 0);

        if (d_output_T)
            d_visit_writer->registerPlotQuantity("Temperature", "SCALAR", d_T_current_idx, 0);
    }
    if (d_solve_mass_conservation)
    {
        // Set various objects with conservative time integrator.
        Pointer<AdvDiffConservativeMassTransportQuantityIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
        rho_p_cc_integrator->setCellCenteredDensityBoundaryConditions(d_rho_bc_coef);
        //        rho_p_cc_integrator->setCellCenteredMaterialPropertyBoundaryConditions(d_rho_bc_coef);
        rho_p_cc_integrator->setCellCenteredTransportQuantityBoundaryConditions(d_h_bc_coef);
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
        d_T_solver_needs_init = true;
        d_T_rhs_op_needs_init = true;
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_T_C_idx)) level->allocatePatchData(d_T_C_idx, current_time);
        if (!level->checkAllocated(d_T_temp_idx)) level->allocatePatchData(d_T_temp_idx, current_time);
        if (!level->checkAllocated(d_T_temp_rhs_idx)) level->allocatePatchData(d_T_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_lf_pre_idx)) level->allocatePatchData(d_lf_pre_idx, current_time);
        //        if (!level->checkAllocated(d_lf_sc_idx)) level->allocatePatchData(d_lf_sc_idx, current_time);
        if (!level->checkAllocated(d_H_pre_idx)) level->allocatePatchData(d_H_pre_idx, current_time);
        if (!level->checkAllocated(d_dh_dT_scratch_idx)) level->allocatePatchData(d_dh_dT_scratch_idx, current_time);
        if (!level->checkAllocated(d_drho_dT_scratch_idx))
            level->allocatePatchData(d_drho_dT_scratch_idx, current_time);
        if (!level->checkAllocated(d_T_pre_idx)) level->allocatePatchData(d_T_pre_idx, current_time);
        if (!level->checkAllocated(d_grad_T_idx)) level->allocatePatchData(d_grad_T_idx, current_time);
    }

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

        const int h_current_idx = var_db->mapVariableAndContextToIndex(d_h_var, getCurrentContext());
        const int h_scratch_idx = var_db->mapVariableAndContextToIndex(d_h_var, getScratchContext());
        const int h_new_idx = var_db->mapVariableAndContextToIndex(d_h_var, getNewContext());

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

        // Setup the operators and solvers and compute the right-hand-side terms.
        const int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());
        const int H_current_idx = var_db->mapVariableAndContextToIndex(d_H_var, getCurrentContext());
        const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
        const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
        const int lf_current_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getCurrentContext());
        const int lf_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_N_var, getScratchContext());
        d_hier_cc_data_ops->copyData(d_H_pre_idx, H_current_idx);

        // Initialize h only at the start of the simulation.
        if (initial_time)
            computeEnthalpyBasedOnNonLinearTemperature(
                h_current_idx, T_current_idx, rho_current_idx, lf_current_idx, H_current_idx);

        //        const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
        //        std::cout << "L2 norm of h_current_idx\t" << d_hier_cc_data_ops->L2Norm(h_current_idx, wgt_cc_idx) <<
        //        std::endl;

        if (d_solve_mass_conservation)
        {
            const int T_u_current_idx = var_db->mapVariableAndContextToIndex(d_T_u_var, getCurrentContext());

            // Keep track of the time-lagged velocity, specific heat and temperature.
            d_hier_fc_data_ops->copyData(d_U_old_new_idx, T_u_current_idx);
            //            d_hier_cc_data_ops->copyData(d_cp_old_new_idx, Cp_current_idx);
            d_hier_cc_data_ops->copyData(d_T_old_new_idx, T_current_idx);

            d_rho_p_integrator->setSolutionTime(current_time);
            d_rho_p_integrator->setTimeInterval(current_time, new_time);

            // For conservative discretization, an approximation to rho^{n+1}
            // will be computed from rho^{n}, which requires additional options to be
            // set.

            // Set the rho^{n} density
            Pointer<AdvDiffConservativeMassTransportQuantityIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
            d_rho_p_integrator->setDensityPatchDataIndex(rho_current_idx);
            // d_rho_p_integrator->setCellCenteredMaterialPropertyPatchDataIndex(Cp_current_idx);
            // d_rho_p_integrator->setCellCenteredTemperaturePatchDataIndex(T_current_idx);

            // Set the convective derivative patch data index.
            const int T_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_N_var, getScratchContext());
            d_rho_p_integrator->setConvectiveDerivativePatchDataIndex(T_N_scratch_idx);
            rho_p_cc_integrator->setMassConservationPatchDataIndex(d_M_idx);

            // Data for the conservative time integrator is for cycle 0
            const int cycle_num = 0;
            d_rho_p_integrator->setCycleNumber(cycle_num);

            // Set the velocities used to update the density and the previous time step
            // size
            if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ -1, /*current*/ T_u_current_idx, /*new*/ -1);
                rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                    /*old*/ -1, /*current*/ h_current_idx, /*new*/ -1);
            }
            else
            {
                d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                    /*old*/ d_U_old_current_idx, /*current*/ T_u_current_idx, /*new*/ -1);
                rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                    /*old*/ d_h_old_current_idx, /*current*/ h_current_idx, /*new*/ -1);
                d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
            }

            // Integrate density and convective term of energy equation.
            d_rho_p_integrator->integrate(dt);
            //            std::cout << "max norm of d_M_idx \t" << d_hier_cc_data_ops->maxNorm(d_M_idx, wgt_cc_idx) <<
            //            std::endl;
        }
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
        PoissonSpecifications T_rhs_op_spec(d_object_name + "::rhs_op_spec::" + d_T_var->getName());

        // There is no coefficient in the RHS for the enthalpy formulation.
        T_rhs_op_spec.setCZero();
        for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
        {
            d_reset_kappa_fcns[k](d_T_diff_coef_cc_current_idx,
                                  d_T_diffusion_coef_cc_var,
                                  d_hier_math_ops,
                                  -1 /*cycle_num*/,
                                  apply_time,
                                  current_time,
                                  new_time,
                                  d_reset_kappa_fcns_ctx[k]);
        }

        // Interpolate the cell-centered diffusion coef to side-centered.
        d_hier_cc_data_ops->copyData(d_T_diff_coef_cc_scratch_idx, d_T_diff_coef_cc_current_idx);
        d_k_bdry_bc_fill_op->fillData(current_time);

        //        std::cout << "L2 norm of T_diff_coef_cc_scratch_idx\t"
        //                  << d_hier_cc_data_ops->L2Norm(d_T_diff_coef_cc_scratch_idx) << std::endl;

        if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
        {
            interpolateCCToSC(T_diff_coef_current_idx, d_T_diff_coef_cc_scratch_idx);
            //            std::cout << "L2 norm of T_diff_coef_current_idx just after interpolating\t"
            //                      << d_hier_sc_data_ops->L2Norm(T_diff_coef_current_idx) << std::endl;
        }
        else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
        {
            interpolateCCToSCHarmonic(T_diff_coef_current_idx, d_T_diff_coef_cc_scratch_idx);
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

        //        d_hier_sc_data_ops->scale(T_diff_coef_scratch_idx, -alpha, T_diff_coef_current_idx);
        //        d_hier_sc_data_ops->scale(T_diff_coef_scratch_idx, dt, T_diff_coef_current_idx);
        //        T_solver_spec.setDPatchDataId(T_diff_coef_scratch_idx);

        d_hier_sc_data_ops->scale(T_diff_coef_rhs_scratch_idx, (1.0 - alpha), T_diff_coef_current_idx);
        T_rhs_op_spec.setDPatchDataId(T_diff_coef_rhs_scratch_idx);
        //        std::cout << "L2 norm of T_diff_coef_scratch_idx\t" <<
        //        d_hier_sc_data_ops->L2Norm(T_diff_coef_scratch_idx)
        //                  << std::endl;

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

        if (d_T_convective_op_needs_init)
        {
            d_T_convective_op->initializeOperatorState(*d_T_sol, *d_T_rhs);
            d_T_convective_op_needs_init = false;
        }

        if (d_nonconser_convective_op_needs_init)
        {
            d_nonconser_convective_op->initializeOperatorState(*d_lf_sol, *d_lf_rhs);
            d_nonconser_convective_op_needs_init = false;
        }

        d_hier_cc_data_ops->copyData(rho_new_idx, rho_current_idx);
        d_hier_cc_data_ops->copyData(h_new_idx, h_current_idx);
        d_hier_cc_data_ops->copyData(T_new_idx, T_current_idx);
        d_hier_cc_data_ops->copyData(d_T_diff_coef_cc_new_idx, d_T_diff_coef_cc_current_idx);
        d_hier_cc_data_ops->copyData(lf_new_idx, lf_current_idx);

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
    const int lf_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_N_var, getScratchContext());
    const int H_current_idx = var_db->mapVariableAndContextToIndex(d_H_var, getCurrentContext());
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    const int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());

    //    if (cycle_num > 0)
    //    {
    //        // Update the advection velocity for lf.
    //        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
    //        const int lf_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getScratchContext());
    //        const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
    //        if (d_lf_u_fcn)
    //        {
    //            d_lf_u_fcn->setDataOnPatchHierarchy(lf_u_new_idx, d_lf_u_var, d_hierarchy, new_time);
    //        }
    //        d_hier_fc_data_ops->linearSum(lf_u_scratch_idx, 0.5, lf_u_current_idx, 0.5, lf_u_new_idx);
    //    }

    int rho_new_idx, rho_scratch_idx, rho_current_idx;
    rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());
    rho_scratch_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getScratchContext());
    rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());

    int Cp_new_idx, Cp_scratch_idx, Cp_current_idx;
    Cp_new_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getNewContext());
    Cp_scratch_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getScratchContext());
    Cp_current_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getCurrentContext());

    int h_new_idx, h_scratch_idx, h_current_idx;
    h_new_idx = var_db->mapVariableAndContextToIndex(d_h_var, getNewContext());
    h_scratch_idx = var_db->mapVariableAndContextToIndex(d_h_var, getScratchContext());
    h_current_idx = var_db->mapVariableAndContextToIndex(d_h_var, getCurrentContext());

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

    // update density.
    double apply_time = new_time;
    //    for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
    //    {
    //        d_reset_rho_fcns[k](rho_new_idx,
    //                            d_rho_var,
    //                            d_hier_math_ops,
    //                            -1 /*cycle_num*/,
    //                            apply_time,
    //                            current_time,
    //                            new_time,
    //                            d_reset_rho_fcns_ctx[k]);
    //    }

    // In the special case of a conservative discretization form, the updated
    // density is calculated by application of the mass and convective
    // momentum integrator.
    const int T_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_N_var, getScratchContext());

    Pointer<AdvDiffConservativeMassTransportQuantityIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
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

        const int T_u_current_idx = var_db->mapVariableAndContextToIndex(d_T_u_var, getCurrentContext());
        const int T_u_new_idx = var_db->mapVariableAndContextToIndex(d_T_u_var, getNewContext());

        // Set the velocities used to update the density
        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ -1, /*current*/ T_u_current_idx, /*new*/ T_u_new_idx);
            //            rho_p_cc_integrator->setMaterialPropertyPatchDataIndices(
            //                /*old*/ -1, /*current*/ Cp_current_idx, /*new*/ Cp_new_idx);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*old*/ -1, /*current*/ h_current_idx, /*new*/ h_new_idx);
        }
        else
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ d_U_old_current_idx,
                /*current*/ T_u_current_idx,
                /*new*/ T_u_new_idx);
            //            rho_p_cc_integrator->setMaterialPropertyPatchDataIndices(
            //                /*old*/ d_cp_old_current_idx,
            //                /*current*/ Cp_current_idx,
            //                /*new*/ Cp_new_idx);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*old*/ d_h_old_current_idx,
                /*current*/ h_current_idx,
                /*new*/ h_new_idx);

            d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        d_rho_p_integrator->integrate(dt);
    }

        d_updated_rho_cc_idx = d_rho_p_integrator ? d_rho_p_integrator->getUpdatedDensityPatchDataIndex() : rho_new_idx;
        d_hier_cc_data_ops->copyData(rho_new_idx,
                                     d_updated_rho_cc_idx,
                                     /*interior_only*/ true);

        const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
        //        std::cout << "L2 norm of T_rhs_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(T_rhs_scratch_idx,
        //        wgt_cc_idx)
        //                  << "for cycle\t" << cycle_num << "\n";

        // Account for the convective acceleration term N_full.
        if (d_T_u_var) d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -1.0, T_N_scratch_idx, T_rhs_scratch_idx);

        PoissonSpecifications T_solver_spec(d_object_name + "::solver_spec::" + d_T_var->getName());

        double lf_relative_iteration_error = 1.0;
        double inner_iterations = 1.0;

        // inner iteration for Newton scheme.
        while (lf_relative_iteration_error >= d_lf_iteration_error_tolerance &&
               inner_iterations <= d_max_inner_iterations)
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

            //        apply_time = new_time;
            //        for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
            //        {
            //            d_reset_kappa_fcns[k](d_T_diff_coef_cc_new_idx,
            //                                  d_T_diffusion_coef_cc_var,
            //                                  d_hier_math_ops,
            //                                  -1 /*cycle_num*/,
            //                                  apply_time,
            //                                  current_time,
            //                                  new_time,
            //                                  d_reset_kappa_fcns_ctx[k]);
            //
            //        }

            // Interpolate the cell-centered diffusion coef to side-centered.
            d_hier_cc_data_ops->copyData(d_T_diff_coef_cc_scratch_idx, d_T_diff_coef_cc_new_idx);
            d_k_bdry_bc_fill_op->fillData(new_time);

            if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
            {
                interpolateCCToSC(T_diff_coef_new_idx, d_T_diff_coef_cc_scratch_idx);
            }
            else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
            {
                interpolateCCToSCHarmonic(T_diff_coef_new_idx, d_T_diff_coef_cc_scratch_idx);
            }
            else
            {
                TBOX_ERROR("this statement should not be reached");
            }
            //        std::cout << "L2 norm of T_diff_coef_new_idx\t" << d_hier_sc_data_ops->L2Norm(T_diff_coef_new_idx,
            //        wgt_cc_idx) << "for cycle\t" << cycle_num << "\n";
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

            // To use same H used in advection of H in the convective term.
            d_hier_cc_data_ops->copyData(d_H_pre_idx, H_new_idx); //

            computeEnthalpyDerivative(d_dh_dT_scratch_idx, T_new_idx, H_new_idx);
            //            std::cout << "L2 norm of d_dh_dT_scratch_idx\t"
            //                      << d_hier_cc_data_ops->L2Norm(d_dh_dT_scratch_idx, wgt_cc_idx) << "for cycle\t" <<
            //                      cycle_num
            //                      << "\n";

            // set rho*Cp/dt + K*lambda.
            const double lambda = 0.0;
            d_hier_cc_data_ops->multiply(d_C_new_idx, rho_new_idx, d_dh_dT_scratch_idx);
            d_hier_cc_data_ops->copyData(d_T_C_idx, d_C_new_idx);
            d_hier_cc_data_ops->scale(d_T_C_idx, 1.0 / dt, d_T_C_idx);

            //        std::cout << "L2 norm of d_T_C_idx\t" << d_hier_cc_data_ops->L2Norm(d_T_C_idx, wgt_cc_idx) << "for
            //        cycle\t" << cycle_num << "\n";
            T_solver_spec.setCPatchDataId(d_T_C_idx);

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

            // Account for forcing terms.
            const int T_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_F_var, getScratchContext());
            const int T_F_new_idx = var_db->mapVariableAndContextToIndex(d_T_F_var, getNewContext());
            if (d_T_F_fcn)
            {
                d_T_F_fcn->setDataOnPatchHierarchy(T_F_scratch_idx, d_T_F_var, d_hierarchy, half_time);
            }
            else
                d_hier_cc_data_ops->setToScalar(T_F_scratch_idx, 0.0);

            //        std::cout << "L2 norm of T_rhs_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(T_rhs_scratch_idx,
            //        wgt_cc_idx) << "for cycle\t" << cycle_num << "\n";
            // Add Allen-Cahn advection term.
            computeTemperatureSourceTerm(T_F_scratch_idx, dt);
            d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, +1.0, T_F_scratch_idx, T_rhs_scratch_idx);

            // Storing T^n+1,m.
            d_hier_cc_data_ops->copyData(d_T_pre_idx, T_new_idx);

            // Solve for T(n+1, m+1).
            T_solver->solveSystem(*d_T_sol, *d_T_rhs);
            d_hier_cc_data_ops->copyData(T_new_idx, T_scratch_idx);

            // Update Temperature. T^n+1, m+1 = T^n+1,m + \Delta T.
            //    std::cout << "L2 norm of T and delta T after solving\t" << d_hier_cc_data_ops->L2Norm(T_new_idx,
            //               wgt_cc_idx) << "\t" << d_hier_cc_data_ops->L2Norm(T_scratch_idx,
            //               wgt_cc_idx) << "for cycle\t" << cycle_num << "\n";
            //    d_hier_cc_data_ops->add(T_new_idx, T_scratch_idx, T_new_idx);
            //    std::cout << "L2 norm of Final Tafter solving\t" << d_hier_cc_data_ops->L2Norm(T_new_idx,
            //    wgt_cc_idx) << "\t" << d_hier_cc_data_ops->L2Norm(T_scratch_idx,
            //    wgt_cc_idx) << "for cycle\t" << cycle_num << "\n";

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
            updateEnthalpy(h_new_idx, T_new_idx, d_T_pre_idx);
            //            std::cout << "L2 norm of h_new_idx\t" << d_hier_cc_data_ops->L2Norm(h_new_idx, wgt_cc_idx) <<
            //            "\t"
            //                      << "for cycle\t" << cycle_num << "\n";

            // Find T^n+1,m+1 based on h^n+1, m+1.
            computeTemperatureBasedOnNonLinearEnthalpy(T_new_idx, h_new_idx, H_new_idx);

            // Find lf^n+1, m+1 based on h^n+1, m+1.
            d_hier_cc_data_ops->copyData(d_lf_pre_idx, lf_new_idx);
            computeLiquidFraction(lf_new_idx, h_new_idx, H_new_idx);
            //            std::cout << "L2 norm of lf_new_idx\t" << d_hier_cc_data_ops->L2Norm(lf_new_idx, wgt_cc_idx)
            //            << "\t"
            //                      << "for cycle\t" << cycle_num << "\n";

            //            // update density. Do not update density. We want to use from mass integrator.
            //            double apply_time = new_time;
            //            for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
            //            {
            //                d_reset_rho_fcns[k](rho_new_idx,
            //                                    d_rho_var,
            //                                    d_hier_math_ops,
            //                                    -1 /*cycle_num*/,
            //                                    apply_time,
            //                                    current_time,
            //                                    new_time,
            //                                    d_reset_rho_fcns_ctx[k]);
            //            }

            //             This will be used in continuity source term
            double apply_time = new_time;
            for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
            {
                d_reset_rho_fcns[k](d_updated_rho_cc_idx,
                                    d_updated_rho_var,
                                    d_hier_math_ops,
                                    -1 /*cycle_num*/,
                                    apply_time,
                                    current_time,
                                    new_time,
                                    d_reset_rho_fcns_ctx[k]);
            }

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

            // update conductivity.
            for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
            {
                d_reset_kappa_fcns[k](d_T_diff_coef_cc_new_idx,
                                      d_T_diffusion_coef_cc_var,
                                      d_hier_math_ops,
                                      -1 /*cycle_num*/,
                                      apply_time,
                                      current_time,
                                      new_time,
                                      d_reset_kappa_fcns_ctx[k]);
            }

            // Finding L2 norm of lf^m-1 iteration.
            lf_relative_iteration_error = d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx);

            // Finding lf^m - lf^m-1.
            d_hier_cc_data_ops->subtract(d_lf_pre_idx, lf_new_idx, d_lf_pre_idx);

            //            computeLiquidFractionRelativeError(lf_new_idx, d_lf_pre_idx); // you can remove this function.
            pout << "liquid fraction error: L2 norm:\t"
                 << d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error) << "\n"
                 << "L1 norm: "
                 << d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error) << "\n"
                 << "max norm: "
                 << d_hier_cc_data_ops->maxNorm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error) << "\n";

            lf_relative_iteration_error =
                d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error);
            inner_iterations++;

            if (d_T_F_var)
            {
                d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, -1.0, T_F_scratch_idx, T_rhs_scratch_idx);
                d_hier_cc_data_ops->copyData(T_F_new_idx, T_F_scratch_idx);
            }
        }


        // Reset the right-hand side vector.
        if (d_T_u_var) d_hier_cc_data_ops->axpy(T_rhs_scratch_idx, +1.0, T_N_scratch_idx, T_rhs_scratch_idx);

        // To be used in continuity equation.
        computeContinuitySourceTerm(
            d_div_U_F_idx, T_new_idx, h_new_idx, rho_new_idx, T_diff_coef_new_idx, H_new_idx, new_time);

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
        level->deallocatePatchData(d_T_C_idx);
        level->deallocatePatchData(d_T_temp_idx);
        level->deallocatePatchData(d_T_temp_rhs_idx);
        level->deallocatePatchData(d_lf_pre_idx);
        level->deallocatePatchData(d_T_pre_idx);
        level->deallocatePatchData(d_H_pre_idx);
        level->deallocatePatchData(d_dh_dT_scratch_idx);
        level->deallocatePatchData(d_drho_dT_scratch_idx);
        level->deallocatePatchData(d_grad_T_idx);
    }

    AdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

// Pointer<CellVariable<NDIM, double> >
// IEPSemiImplicitHierarchyIntegrator::getHeavisideVariable() const
//{
//    return d_H_var;
//} // getHeavisideVariable

void
IEPSemiImplicitHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    AdvDiffSemiImplicitHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
        base_hierarchy, coarsest_level, finest_level);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());
    std::vector<RobinBcCoefStrategy<NDIM>*> H_bc_coef = getPhysicalBcCoefs(d_H_var);

    // Setup the patch boundary filling objects.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent H_bc_component(H_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     BDRY_EXTRAP_TYPE,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     H_bc_coef);
    d_H_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_H_bdry_bc_fill_op->initializeOperatorState(H_bc_component, d_hierarchy);

    d_T_diff_coef_cc_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_cc_var, getScratchContext());
    InterpolationTransactionComponent k_bc_component(d_T_diff_coef_cc_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     BDRY_EXTRAP_TYPE,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_lf_bc_coef); // liquid fraction boundary condition is used.
    d_k_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_k_bdry_bc_fill_op->initializeOperatorState(k_bc_component, d_hierarchy);

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

    d_T_solver_needs_init = true;
    d_nonconser_convective_op_needs_init = true;
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
    d_lf_diffusion_coef_var = lf_diff_coef_var;
    d_lf_init = nullptr;
    d_lf_F_fcn = nullptr;
    d_lf_bc_coef = nullptr;
    return;
} // registerLiquidFractionVariable

void
IEPSemiImplicitHierarchyIntegrator::setHeavisideVariable(Pointer<CellVariable<NDIM, double> > H_var)
{
    d_H_var = H_var;
    return;
} // setHeavisideVariable

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
    Pointer<CellVariable<NDIM, double> > T_diff_coef_cc_var =
        new CellVariable<NDIM, double>(T_var->getName() + "::diff_coef_cc", T_depth);
    Pointer<CellVariable<NDIM, double> > h_var = new CellVariable<NDIM, double>("enthalpy", T_depth);

    // Set default values.
    d_h_var = h_var;
    d_T_u_var = nullptr;
    d_T_F_var = T_F_var;
    d_T_rhs_var = T_rhs_var;
    d_T_diffusion_coef_var = T_diff_coef_var;
    d_T_diffusion_coef_cc_var = T_diff_coef_cc_var;
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

// void
// IEPSemiImplicitHierarchyIntegrator::setPhysicalBcCoefHeavisideEquation(Pointer<CellVariable<NDIM, double> > H_var,
//                                                                       RobinBcCoefStrategy<NDIM>* H_bc_coef)
//{
//#if !defined(NDEBUG)
//    TBOX_ASSERT(d_H_var);
//#endif
//    d_H_bc_coef = H_bc_coef;
//    return;
//} // setPhysicalBcCoefHeavisideEquation

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

void
IEPSemiImplicitHierarchyIntegrator::setEnthalpyBcCoef(RobinBcCoefStrategy<NDIM>* h_bc_coef)
{
    d_h_bc_coef = h_bc_coef;
    return;
} // setEnthalpyBcCoef

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
IEPSemiImplicitHierarchyIntegrator::getNonConservativeConvectiveOperator(
    Pointer<CellVariable<NDIM, double> > Q_var,
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* Q_bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_var);
#endif
    if (!d_nonconser_convective_op)
    {
        std::vector<RobinBcCoefStrategy<NDIM>*> Q_bc_coef(1, Q_bc_coefs);
        AdvDiffConvectiveOperatorManager* nonconser_convective_op_manager =
            AdvDiffConvectiveOperatorManager::getManager();
        d_nonconser_convective_op =
            nonconser_convective_op_manager->allocateOperator("CUI",
                                                              d_object_name + "::lfConvectiveOperator",
                                                              Q_var,
                                                              d_nonconser_convective_op_input_db,
                                                              d_nonconser_convective_difference_form,
                                                              Q_bc_coef);
        d_nonconser_convective_op_needs_init = true;
    }
    return d_nonconser_convective_op;
} // getNonConservativeConvectiveOperator

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
                                                                      d_object_name + "::TConvectiveOperator",
                                                                      d_T_var,
                                                                      d_T_convective_op_input_db,
                                                                      d_T_convective_difference_form,
                                                                      T_bc_coefs);
        d_T_convective_op_needs_init = true;
    }
    return d_T_convective_op;
} // getConvectiveOperatorTemperatureEquation

void
IEPSemiImplicitHierarchyIntegrator::setTemperatureSourceTermFunction(Pointer<IBTK::CartGridFunction> T_F_fcn)
{
    d_T_F_fcn = T_F_fcn;
    return;
} // setTemperatureSourceTermFunction

int
IEPSemiImplicitHierarchyIntegrator::getLiquidFractionMaterialDerivativeIndex()
{
    return d_lf_material_derivative_idx;
} // getLiquidFractionMaterialDerivativeIndex

int
IEPSemiImplicitHierarchyIntegrator::getContinuitySourceTermIndex()
{
    return d_div_U_F_idx;
} // getContinuitySourceTermIndex

int
IEPSemiImplicitHierarchyIntegrator::getUpdatedDensityIndex()
{
    return d_updated_rho_cc_idx;
} // getUpdatedDensityIndex

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

void
IEPSemiImplicitHierarchyIntegrator::registerResetLiquidFractionFcn(ResetLiquidFractionFcnPtr callback, void* ctx)
{
    d_reset_liquid_fraction_fcns.push_back(callback);
    d_reset_liquid_fraction_fcns_ctx.push_back(ctx);
    return;
} // registerResetLiquidFractionFcn

void
IEPSemiImplicitHierarchyIntegrator::registerResetLiquidFractionDerivativeFcn(
    ResetLiquidFractionDerivativeFcnPtr callback,
    void* ctx)
{
    d_reset_dlf_dT_fcns.push_back(callback);
    d_reset_dlf_dT_fcns_ctx.push_back(ctx);
    return;
} // registerResetLiquidFractionDerivativeFcn

void
IEPSemiImplicitHierarchyIntegrator::registerResetLiquidFractionInverseFcn(ResetLiquidFractionInverseFcnPtr callback,
                                                                          void* ctx)
{
    d_reset_lf_inverse_fcns.push_back(callback);
    d_reset_lf_inverse_fcns_ctx.push_back(ctx);
    return;
} // registerResetLiquidFractionInverseFcn

void
IEPSemiImplicitHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
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
    d_lf_convective_op_input_db = db->putDatabase("d_lf_convective_op_db"); // I am not sure whether this is correct.

    db->putString("d_T_diffusion_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_T_diffusion_time_stepping_type));
    db->putString("d_T_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_T_convective_time_stepping_type));
    db->putString("d_T_init_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_T_init_convective_time_stepping_type));
    db->putString("d_T_convective_difference_form",
                  enum_to_string<ConvectiveDifferencingType>(d_T_convective_difference_form));
    db->putString("d_T_convective_op_type", d_T_convective_op_type);
    d_T_convective_op_input_db = db->putDatabase("d_T_convective_op_db"); // I am not sure whether this is correct.

    db->putDouble("d_latent_heat", d_latent_heat);
    db->putDouble("d_latent_heat_temp", d_latent_heat_temp);
    db->putDouble("d_rho_liquid", d_rho_liquid);
    db->putDouble("d_rho_solid", d_rho_solid);
    db->putDouble("d_cp_liquid", d_cp_liquid);
    db->putDouble("d_cp_solid", d_cp_solid);
    db->putDouble("d_cp_gas", d_cp_gas);
    db->putDouble("d_T_melt", d_T_melt);

    db->putDouble("d_M_lf", d_M_lf);
    db->putDouble("d_lambda_lf", d_lambda_lf);
    db->putDouble("d_eta_lf", d_eta_lf);
    db->putBool("d_solve_energy", d_solve_energy);
    db->putBool("d_solve_mass_conservation", d_solve_mass_conservation);

    db->putDouble("d_free_parameter", d_free_parameter);
    db->putDouble("d_eps", d_eps);
    db->putDouble("d_beta", d_beta);
    db->putDouble("d_H_diffusion_coefficient", d_H_diffusion_coefficient);
    db->putBool("d_apply_brinkman", d_apply_brinkman);
    if (d_apply_brinkman)
    {
        d_lf_brinkman_db = db->putDatabase("d_lf_brinkman_db"); // I am not sure whether this is correct.
        d_lf_brinkman_db->putDouble("d_lf_b", d_lf_b);
    }
    db->putBool("d_add_diffusion", d_add_diffusion);

    db->putDouble("d_liquidus_temperature", d_liquidus_temperature);
    db->putDouble("d_solidus_temperature", d_solidus_temperature);
    db->putBool("d_compute_T_material_derivative", d_compute_T_material_derivative);
    db->putInteger("d_max_inner_iterations", d_max_inner_iterations);
    db->putDouble("d_lf_iteration_error_tolerance", d_lf_iteration_error_tolerance);

    AdvDiffSemiImplicitHierarchyIntegrator::putToDatabaseSpecialized(db);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IEPSemiImplicitHierarchyIntegrator::computeDoubleWellPotential(int g_firstder_idx,
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

                if (MathUtilities<double>::equalEps(lf, 1.0) && T <= d_T_melt)
                {
                    (*p_firstder_data)(ci) = 1.0;
                }
                else if (MathUtilities<double>::equalEps(lf, 0.0) && T >= d_T_melt)
                {
                    (*p_firstder_data)(ci) = 1.0;
                }
                else
                {
                    (*p_firstder_data)(ci) = 30.0 * std::pow(lf, 4.0) - 60.0 * std::pow(lf, 3.0) + 30.0 * lf * lf;
                }
                // (*p_firstder_data)(ci) = 0.0;
                // std::cout << "p' value is\t" << (*p_firstder_data)(ci) <<  std::endl;
            }
        }
    }
    return;
} // computeInterpolationFunction

void
IEPSemiImplicitHierarchyIntegrator::interpolateCCToSC(int sc_idx, const int cc_idx)
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

            Pointer<SideData<NDIM, double> > sc_data = patch->getPatchData(sc_idx);
            Pointer<CellData<NDIM, double> > cc_data = patch->getPatchData(cc_idx);

            C_TO_S_CWISE_INTERP_FC(sc_data->getPointer(0),
                                   sc_data->getPointer(1),
#if (NDIM == 3)
                                   sc_data->getPointer(2),
#endif
                                   sc_data->getGhostCellWidth().max(),
                                   cc_data->getPointer(),
                                   cc_data->getGhostCellWidth().max(),
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
} // interpolateCCTOSC

void
IEPSemiImplicitHierarchyIntegrator::interpolateCCToSCHarmonic(int sc_idx, const int cc_idx)
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

            Pointer<SideData<NDIM, double> > sc_data = patch->getPatchData(sc_idx);
            Pointer<CellData<NDIM, double> > cc_data = patch->getPatchData(cc_idx);

            // Use this only k var
            C_TO_S_CWISE_HARMONIC_INTERP_FC(sc_data->getPointer(0),
                                            sc_data->getPointer(1),
#if (NDIM == 3)
                                            sc_data->getPointer(2),
#endif
                                            sc_data->getGhostCellWidth().max(),
                                            cc_data->getPointer(),
                                            cc_data->getGhostCellWidth().max(),
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
} // interpolateCCTOSCHarmonic

// void
// IEPSemiImplicitHierarchyIntegrator::computeLiquidFractionSourceTerm(int F_scratch_idx,
//                                                                    const double dt,
//                                                                    const double new_time)
//{
//    const int coarsest_ln = 0;
//    const int finest_ln = d_hierarchy->getFinestLevelNumber();
//    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
//
//    int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
//    int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
//    int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
//    int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());
//
//    // Filling ghost cells for heaviside and liquid fraction.
//    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
//    std::vector<InterpolationTransactionComponent> transaction_comps(2);
//    transaction_comps[0] = InterpolationTransactionComponent(lf_scratch_idx,
//                                                             lf_new_idx,
//                                                             "CONSERVATIVE_LINEAR_REFINE",
//                                                             false,
//                                                             "CONSERVATIVE_COARSEN",
//                                                             "LINEAR",
//                                                             false,
//                                                             d_lf_bc_coef);
//
//    transaction_comps[1] = InterpolationTransactionComponent(H_scratch_idx,
//                                                             H_new_idx,
//                                                             "CONSERVATIVE_LINEAR_REFINE",
//                                                             false,
//                                                             "CONSERVATIVE_COARSEN",
//                                                             "LINEAR",
//                                                             false,
//                                                             d_lf_bc_coef);
//
//    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
//    hier_bdry_fill->initializeOperatorState(transaction_comps, d_hierarchy);
//    hier_bdry_fill->fillData(new_time);
//
//    // perform gradient of heaviside.
//    d_hier_math_ops->grad(d_grad_H_idx, d_grad_H_var, true, 1.0, H_scratch_idx, d_H_var, nullptr, new_time);
//
//    // compute lf*D*grad_H.
//    d_hier_sc_data_ops->scale(d_grad_H_idx, d_H_diffusion_coefficient, d_grad_H_idx);
//
//    interpolateCCToSC(d_lf_sc_idx, lf_scratch_idx);
//    d_hier_sc_data_ops->multiply(d_grad_H_idx, d_lf_sc_idx, d_grad_H_idx);
//
//    // compute div(lf*D*grad_H).
//    d_hier_math_ops->div(
//        d_lf_diffusion_idx, d_lf_diffusion_var, 1.0, d_grad_H_idx, d_grad_H_var, nullptr, new_time, false);
//
//    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
//    {
//        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
//        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
//        {
//            Pointer<Patch<NDIM> > patch = level->getPatch(p());
//            const Box<NDIM>& patch_box = patch->getBox();
//
//            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(d_lf_new_idx);
//            Pointer<CellData<NDIM, double> > lf_cur_data = patch->getPatchData(d_lf_current_idx);
//            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(d_T_new_idx);
//            Pointer<CellData<NDIM, double> > p_firstder_data = patch->getPatchData(d_p_firstder_idx);
//            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);
//            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_scratch_idx);
//            Pointer<CellData<NDIM, double> > g_firstder_data = patch->getPatchData(d_g_firstder_idx);
//            Pointer<CellData<NDIM, double> > g_secondder_data = patch->getPatchData(d_g_secondder_idx);
//            Pointer<CellData<NDIM, double> > lf_diffusion_data = patch->getPatchData(d_lf_diffusion_idx);
//
//            for (Box<NDIM>::Iterator it(patch_box); it; it++)
//            {
//                CellIndex<NDIM> ci(it());
//
//                double F = -d_M_lf * d_rho_liquid * d_latent_heat * (*H_data)(ci) * (*p_firstder_data)(ci) *
//                           (d_T_melt - (*T_data)(ci)) / d_T_melt;
//                (*F_data)(ci) = F -
//                                (d_M_lf * d_lambda_lf * (*H_data)(ci) / std::pow(d_eta_lf, 2.0) *
//                                 ((*g_firstder_data)(ci) - ((*g_secondder_data)(ci) * (*lf_data)(ci)))) +
//                                (*lf_diffusion_data)(ci);
//
//                // Brinkman force
//                if (d_apply_brinkman) (*F_data)(ci) += d_beta / dt * (1.0 - (*H_data)(ci)) * d_lf_b;
//                //                if (d_apply_brinkman) (*F_data)(ci) += d_beta / dt * (1.0 - (*H_data)(ci)) *
//                //                (*lf_cur_data)(ci);
//            }
//        }
//    }
//    return;
//} // computeLiquidFractionSourceTerm

void
IEPSemiImplicitHierarchyIntegrator::computeTemperatureSourceTerm(int F_scratch_idx, const double dt)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int h_new_idx = var_db->mapVariableAndContextToIndex(d_h_var, getNewContext());
    int h_current_idx = var_db->mapVariableAndContextToIndex(d_h_var, getCurrentContext());
    int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());
    int dh_dT_scratch_idx = var_db->mapVariableAndContextToIndex(d_dh_dT_var, getCurrentContext());
    int rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());
    int rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_new_data = patch->getPatchData(T_new_idx);
            Pointer<CellData<NDIM, double> > h_new_data = patch->getPatchData(h_new_idx);
            Pointer<CellData<NDIM, double> > h_current_data = patch->getPatchData(h_current_idx);
            Pointer<CellData<NDIM, double> > rho_new_data = patch->getPatchData(rho_new_idx);
            Pointer<CellData<NDIM, double> > rho_current_data = patch->getPatchData(rho_current_idx);
            Pointer<CellData<NDIM, double> > dh_dT_data = patch->getPatchData(dh_dT_scratch_idx);
            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_scratch_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*F_data)(ci) += -1.0 / dt *
                                 ((*rho_new_data)(ci) * ((*h_new_data)(ci) - (*dh_dT_data)(ci) * (*T_new_data)(ci)) -
                                  (*rho_current_data)(ci) * (*h_current_data)(ci));

                //                (*F_data)(ci) +=
                //                d_rho_liquid * d_latent_heat_temp / dt * ((*H_current_data)(ci) *
                //                (*lf_current_data)(ci));
            }
        }
    }
    return;
} // computeTemperatureSourceTerm

void
IEPSemiImplicitHierarchyIntegrator::computeContinuitySourceTerm(int div_U_F_idx,
                                                                const int T_new_idx,
                                                                const int h_new_idx,
                                                                const int rho_new_idx,
                                                                const int T_diff_coef_idx,
                                                                const int H_new_idx,
                                                                const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, getScratchContext());

    // Filling ghost cells for temperature.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> T_transaction_comps(1);
    T_transaction_comps[0] = InterpolationTransactionComponent(T_scratch_idx,
                                                               T_new_idx,
                                                               "CONSERVATIVE_LINEAR_REFINE",
                                                               false,
                                                               "CONSERVATIVE_COARSEN",
                                                               "LINEAR",
                                                               false,
                                                               d_T_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(T_transaction_comps, d_hierarchy);
    hier_bdry_fill->fillData(new_time);

    // perform gradient of temperature.
    d_hier_math_ops->grad(d_grad_T_idx, d_grad_T_var, true, 1.0, T_scratch_idx, d_T_var, nullptr, new_time);

    // compute k*grad_T.
    d_hier_sc_data_ops->multiply(d_grad_T_idx, d_grad_T_idx, T_diff_coef_idx);

    // compute div(k*grad_T).
    d_hier_math_ops->div(div_U_F_idx, d_div_U_F_var, 1.0, d_grad_T_idx, d_grad_T_var, nullptr, new_time, false);

    // copying div(k*grad_T) for plotting purpose.
    d_hier_cc_data_ops->copyData(d_div_U_F_diff_idx, div_U_F_idx);

    const double h_s = d_cp_solid * d_solidus_temperature;
    const double cp_avg = (d_cp_solid + d_cp_liquid) / 2.0;
    const double h_l = cp_avg * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat_temp;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(h_new_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_new_idx);
            Pointer<CellData<NDIM, double> > div_U_F_data = patch->getPatchData(div_U_F_idx);
            Pointer<CellData<NDIM, double> > div_U_F_deno_data = patch->getPatchData(d_div_U_F_deno_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double> > lf_material_derivative_data =
                patch->getPatchData(d_lf_material_derivative_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*div_U_F_deno_data)(ci) = 0.0;

                if ((*h_data)(ci) >= h_s && (*h_data)(ci) <= h_l)
                {
                    // Plotting one over denominator.
                    (*div_U_F_deno_data)(ci) =
                        1.0 / ((*h_data)(ci) * (d_rho_liquid - d_rho_solid) - d_rho_liquid * h_l + d_rho_solid * h_s);

                    const double denominator = (*rho_data)(ci)*std::pow(
                        (*h_data)(ci) * (d_rho_liquid - d_rho_solid) - d_rho_liquid * h_l + d_rho_solid * h_s, 2.0);

                    (*lf_material_derivative_data)(ci) =
                        (*div_U_F_data)(ci)*d_rho_solid * d_rho_liquid * (h_l - h_s) /
                        denominator; // div k grad T rho_s*rho_l (h_l - h_s) / denominator
                }
                else
                    (*lf_material_derivative_data)(ci) = 0.0;

                (*div_U_F_data)(ci) = -(d_rho_liquid - d_rho_solid) * (*lf_material_derivative_data)(ci) *
                                      (*H_data)(ci) / (*rho_data)(ci);
            }
        }
    }

    return;
} // computeContinuitySourceTerm

void
IEPSemiImplicitHierarchyIntegrator::computeMaterialDerivativeOfLiquidFraction(int lf_material_derivative_idx,
                                                                              const double dt,
                                                                              const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
    int lf_current_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getCurrentContext());
    int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
    int lf_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_N_var, getScratchContext());
    int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());
    int T_current_idx = var_db->mapVariableAndContextToIndex(d_T_var, getCurrentContext());
    int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, getScratchContext());
    int T_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_N_var, getScratchContext());

    if (d_compute_T_material_derivative)
    {
        // compute (T^n+1 - T^n) / dt.
        d_hier_cc_data_ops->subtract(lf_material_derivative_idx, T_new_idx, T_current_idx);
        d_hier_cc_data_ops->scale(lf_material_derivative_idx, dt, lf_material_derivative_idx);

        // new time level T is used. This is different than how energy equation convective term is treated. Also lf
        //        and T convective op are conservative. We need non-conservative.
        if (d_T_u_var) // default the code uses
        {
            const int T_u_new_idx = var_db->mapVariableAndContextToIndex(d_T_u_var, getNewContext());
            d_nonconser_convective_op->setAdvectionVelocity(T_u_new_idx);
            d_hier_cc_data_ops->copyData(T_scratch_idx, T_new_idx);
            d_nonconser_convective_op->setSolutionTime(new_time);
            d_nonconser_convective_op->applyConvectiveOperator(T_scratch_idx, T_N_scratch_idx);

            // add u.grad T.
            d_hier_cc_data_ops->add(lf_material_derivative_idx, lf_material_derivative_idx, T_N_scratch_idx);
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();

                Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_new_idx);
                Pointer<CellData<NDIM, double> > lf_material_derivative_data =
                    patch->getPatchData(lf_material_derivative_idx);
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double T = (*T_data)(ci);

                    if (T <= d_liquidus_temperature && T >= d_solidus_temperature)
                        (*lf_material_derivative_data)(ci) *= 1.0 / (d_liquidus_temperature - d_solidus_temperature);
                    else
                        (*lf_material_derivative_data)(ci) = 0.0;
                }
            }
        }
    }
    else
    {
        // compute (varphi^n+1 - varphi^n) / dt.
        d_hier_cc_data_ops->subtract(lf_material_derivative_idx, lf_new_idx, lf_current_idx);
        d_hier_cc_data_ops->scale(lf_material_derivative_idx, dt, lf_material_derivative_idx);

        // new time level lf is used. This is different than how energy equation convective term is treated.
        if (d_lf_u_var)
        {
            const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
            d_nonconser_convective_op->setAdvectionVelocity(lf_u_new_idx);
            d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_new_idx);
            d_nonconser_convective_op->setSolutionTime(new_time);
            d_nonconser_convective_op->applyConvectiveOperator(lf_scratch_idx, lf_N_scratch_idx);

            // add u.grad varphi.
            d_hier_cc_data_ops->add(lf_material_derivative_idx, lf_material_derivative_idx, lf_N_scratch_idx);
        }
    }

    // Multiply H with Dvarphi/Dt
    d_hier_cc_data_ops->multiply(lf_material_derivative_idx, lf_material_derivative_idx, H_new_idx);
    return;
} // computeMaterialDerivativeOfLiquidFraction

// void
// IEPSemiImplicitHierarchyIntegrator::computeLHSOfLiquidFractionEquation(int lf_lhs_idx,
//                                                                       const int lf_N_scratch_idx,
//                                                                       const double dt,
//                                                                       const int cycle_num,
//                                                                       const double new_time,
//                                                                       const double current_time,
//                                                                       const double half_time)
//{
//    const int coarsest_ln = 0;
//    const int finest_ln = d_hierarchy->getFinestLevelNumber();
//    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
//
//    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
//    {
//        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
//        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
//        {
//            Pointer<Patch<NDIM> > patch = level->getPatch(p());
//            const Box<NDIM>& patch_box = patch->getBox();
//            Pointer<CellData<NDIM, double> > lf_new_data = patch->getPatchData(d_lf_new_idx);
//            Pointer<CellData<NDIM, double> > H_new_data = patch->getPatchData(d_H_new_idx);
//            Pointer<CellData<NDIM, double> > lf_current_data = patch->getPatchData(d_lf_current_idx);
//            Pointer<CellData<NDIM, double> > H_current_data = patch->getPatchData(d_H_current_idx);
//            Pointer<CellData<NDIM, double> > lf_lhs_data = patch->getPatchData(lf_lhs_idx);
//
//            for (Box<NDIM>::Iterator it(patch_box); it; it++)
//            {
//                CellIndex<NDIM> ci(it());
//
//                (*lf_lhs_data)(ci) =
//                    (((*H_new_data)(ci) * (*lf_new_data)(ci)) - ((*H_current_data)(ci) * (*lf_current_data)(ci))) /
//                    dt;
//            }
//        }
//    }
//    d_hier_cc_data_ops->axpy(lf_N_scratch_idx, +1.0, lf_lhs_idx, lf_N_scratch_idx);
//
//    const int lf_lhs_N_scratch_idx = (var_db->mapVariableAndContextToIndex(d_lf_lhs_N_var, getScratchContext()));
//    if (d_lf_u_var)
//    {
//        if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(d_lf_convective_time_stepping_type))
//        {
//            d_lf_convective_time_stepping_type = d_lf_init_convective_time_stepping_type;
//        }
//
//        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
//        d_lf_convective_op->setAdvectionVelocity(lf_u_current_idx);
//        const int lf_H_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_H_var, getScratchContext());
//        /// d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_current_idx);
//        d_hier_cc_data_ops->multiply(lf_H_scratch_idx, d_lf_current_idx, d_H_current_idx);
//        d_lf_convective_op->setSolutionTime(current_time);
//        d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_lhs_N_scratch_idx);
//
//        const int lf_N_old_new_idx = var_db->mapVariableAndContextToIndex(d_lf_N_old_var, getNewContext());
//        d_hier_cc_data_ops->copyData(lf_N_old_new_idx, lf_lhs_N_scratch_idx);
//
//        if (d_lf_convective_time_stepping_type == FORWARD_EULER)
//        {
//            d_hier_cc_data_ops->axpy(lf_lhs_idx, 1.0, lf_lhs_N_scratch_idx, lf_lhs_idx);
//        }
//        else if (d_lf_convective_time_stepping_type == TRAPEZOIDAL_RULE)
//        {
//            d_hier_cc_data_ops->axpy(lf_lhs_idx, 0.5, lf_lhs_N_scratch_idx, lf_lhs_idx);
//        }
//    }
//
//    if (cycle_num > 0)
//    {
//        // Update the advection velocity for lf.
//        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
//        const int lf_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getScratchContext());
//        const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
//        if (d_lf_u_fcn)
//        {
//            d_lf_u_fcn->setDataOnPatchHierarchy(lf_u_new_idx, d_lf_u_var, d_hierarchy, new_time);
//        }
//        d_hier_fc_data_ops->linearSum(lf_u_scratch_idx, 0.5, lf_u_current_idx, 0.5, lf_u_new_idx);
//    }
//
//    // Account for the convective difference term.
//    const int lf_H_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_H_var, getScratchContext());
//    TimeSteppingType convective_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
//    if (d_lf_u_var)
//    {
//        convective_time_stepping_type = d_lf_convective_time_stepping_type;
//        if (is_multistep_time_stepping_type(convective_time_stepping_type))
//        {
//#if !defined(NDEBUG)
//            TBOX_ASSERT(convective_time_stepping_type == ADAMS_BASHFORTH);
//#endif
//            if (getIntegratorStep() == 0)
//            {
//                convective_time_stepping_type = d_lf_init_convective_time_stepping_type;
//            }
//            else if (cycle_num > 0)
//            {
//                convective_time_stepping_type = MIDPOINT_RULE;
//                IBAMR_DO_ONCE({
//                    pout << "IEPSemiImplicitHierarchyIntegrator::"
//                            "integrateHierarchy():"
//                            "\n"
//                         << "  WARNING: convective_time_stepping_type = "
//                         << enum_to_string<TimeSteppingType>(d_lf_convective_time_stepping_type)
//                         << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
//                         << "           using " <<
//                         enum_to_string<TimeSteppingType>(d_lf_convective_time_stepping_type)
//                         << " only for the first cycle in each time step;\n"
//                         << "           using " << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
//                         << " for subsequent cycles.\n";
//                });
//            }
//        }
//
//        if (cycle_num > 0)
//        {
//            if (convective_time_stepping_type == MIDPOINT_RULE)
//            {
//                const int lf_u_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getScratchContext());
//                d_lf_convective_op->setAdvectionVelocity(lf_u_scratch_idx);
//                d_hier_cc_data_ops->linearSum(d_lf_scratch_idx, 0.5, d_lf_current_idx, 0.5, d_lf_new_idx);
//                d_hier_cc_data_ops->linearSum(d_H_scratch_idx, 0.5, d_H_current_idx, 0.5, d_H_new_idx);
//                d_hier_cc_data_ops->multiply(lf_H_scratch_idx, d_H_scratch_idx, d_lf_scratch_idx);
//                d_lf_convective_op->setSolutionTime(half_time);
//                d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_lhs_N_scratch_idx);
//            }
//            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
//            {
//                const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());
//                d_lf_convective_op->setAdvectionVelocity(lf_u_new_idx);
//                // d_hier_cc_data_ops->copyData(lf_scratch_idx, lf_new_idx);
//                d_hier_cc_data_ops->multiply(lf_H_scratch_idx, d_lf_new_idx, d_H_new_idx);
//                d_lf_convective_op->setSolutionTime(new_time);
//                d_lf_convective_op->applyConvectiveOperator(lf_H_scratch_idx, lf_lhs_N_scratch_idx);
//            }
//        }
//        if (convective_time_stepping_type == ADAMS_BASHFORTH)
//        {
//#if !defined(NDEBUG)
//            TBOX_ASSERT(cycle_num == 0);
//#endif
//            const int lf_N_old_current_idx = var_db->mapVariableAndContextToIndex(d_lf_N_old_var,
//            getCurrentContext()); const double omega = dt / d_dt_previous[0]; d_hier_cc_data_ops->linearSum(
//                lf_lhs_N_scratch_idx, 1.0 + 0.5 * omega, lf_lhs_N_scratch_idx, -0.5 * omega, lf_N_old_current_idx);
//        }
//
//        if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
//        {
//            d_hier_cc_data_ops->axpy(lf_lhs_idx, 1.0, lf_lhs_N_scratch_idx, lf_lhs_idx);
//        }
//        else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
//        {
//            d_hier_cc_data_ops->axpy(lf_lhs_idx, 0.5, lf_lhs_N_scratch_idx, lf_lhs_idx);
//        }
//        //        std::cout << "L2 norm of lf_lhs_idx\t" << d_hier_cc_data_ops->L2Norm(lf_lhs_idx) << std::endl;
//        //        std::cout << "L2 norm of lf_N_scratch_idx\t" << d_hier_cc_data_ops->L2Norm(lf_N_scratch_idx) <<
//        //        std::endl;
//    }
//
//    return;
//} // computeLHSOfLiquidFractionEquation

void
IEPSemiImplicitHierarchyIntegrator::boundLiquidFraction(int lf_new_idx)
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
            Pointer<CellData<NDIM, double> > lf_new_data = patch->getPatchData(lf_new_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if ((*lf_new_data)(ci) > 1.0)
                    (*lf_new_data)(ci) = 1.0;
                else if ((*lf_new_data)(ci) < 0.0)
                    (*lf_new_data)(ci) = 0.0;
            }
        }
    }
    return;
} // boundLiquidFraction

void
IEPSemiImplicitHierarchyIntegrator::computeLiquidFractionRelativeError(const int lf_new_idx, int lf_pre_idx)
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
            Pointer<CellData<NDIM, double> > lf_new_data = patch->getPatchData(lf_new_idx);
            Pointer<CellData<NDIM, double> > lf_pre_data = patch->getPatchData(lf_pre_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*lf_pre_data)(ci) = std::abs(((*lf_new_data)(ci) - (*lf_pre_data)(ci)));
            }
        }
    }
    return;
} // computeLiquidFractionRelativeError

void
IEPSemiImplicitHierarchyIntegrator::computeEnthalpyBasedOnNonLinearTemperature(int h_idx,
                                                                               const int T_idx,
                                                                               const int rho_idx,
                                                                               const int lf_idx,
                                                                               const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const double h_s = d_cp_solid * d_solidus_temperature;
    const double cp_avg = (d_cp_solid + d_cp_liquid) / 2.0;
    const double h_l = cp_avg * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat_temp;
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

                if ((*H_data)(ci) >= 0.5)
                {
                    if ((*T_data)(ci) < d_solidus_temperature)
                        (*h_data)(ci) = d_cp_solid * (*T_data)(ci);
                    else if ((*T_data)(ci) >= d_solidus_temperature && (*T_data)(ci) <= d_liquidus_temperature)
                        (*h_data)(ci) = cp_avg * ((*T_data)(ci)-d_solidus_temperature) + h_s +
                                        (*lf_data)(ci)*d_rho_liquid * d_latent_heat_temp / (*rho_data)(ci);
                    else
                        (*h_data)(ci) = d_cp_liquid * ((*T_data)(ci)-d_liquidus_temperature) + h_l;
                }
                else
                    (*h_data)(ci) = d_cp_gas * (*T_data)(ci);
            }
        }
    }
    return;
} // computeEnthalpyBasedOnNonLinearTemperature

void
IEPSemiImplicitHierarchyIntegrator::computeTemperatureBasedOnNonLinearEnthalpy(int T_idx,
                                                                               const int h_idx,
                                                                               const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const double h_s = d_cp_solid * d_solidus_temperature;
    const double cp_avg = (d_cp_solid + d_cp_liquid) / 2.0;
    const double h_l = cp_avg * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat_temp;
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

                if ((*H_data)(ci) >= 0.5)
                {
                    if ((*h_data)(ci) < h_s)
                        (*T_data)(ci) = (*h_data)(ci) / d_cp_solid;
                    else if ((*h_data)(ci) >= h_s && (*h_data)(ci) <= h_l)
                        (*T_data)(ci) = d_solidus_temperature + ((*h_data)(ci)-h_s) / (h_l - h_s) *
                                                                    (d_liquidus_temperature - d_solidus_temperature);
                    else
                        (*T_data)(ci) = d_liquidus_temperature + ((*h_data)(ci)-h_l) / d_cp_liquid;
                }
                else
                    (*T_data)(ci) = (*h_data)(ci) / d_cp_gas;
            }
        }
    }
    return;
} // computeTemperatureBasedOnNonLinearEnthalpy

void
IEPSemiImplicitHierarchyIntegrator::updateEnthalpy(int h_new_idx, const int T_new_idx, const int T_pre_idx)
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
IEPSemiImplicitHierarchyIntegrator::computeEnthalpyDerivative(int dh_dT_idx, const int T_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const double cp_avg = (d_cp_liquid + d_cp_solid) / 2.0;
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

                if ((*H_data)(ci) >= 0.5)
                {
                    if ((*T_data)(ci) < d_solidus_temperature)
                        (*dh_dT_data)(ci) = d_cp_solid;
                    else if ((*T_data)(ci) >= d_solidus_temperature && (*T_data)(ci) <= d_liquidus_temperature)
                        (*dh_dT_data)(ci) =
                            cp_avg + d_latent_heat_temp / (d_liquidus_temperature - d_solidus_temperature);
                    else
                        (*dh_dT_data)(ci) = d_cp_liquid;
                }
                else
                    (*dh_dT_data)(ci) = d_cp_gas;
            }
        }
    }
    return;
} // computeEnthalpyDerivative

void
IEPSemiImplicitHierarchyIntegrator::computeDensityDerivative(int drho_dT_idx, const int T_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const double mushy_interval = d_liquidus_temperature - d_solidus_temperature;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            Pointer<CellData<NDIM, double> > drho_dT_data = patch->getPatchData(drho_dT_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if ((*H_data)(ci) >= 0.5)
                {
                    if ((*T_data)(ci) < d_solidus_temperature)
                        (*drho_dT_data)(ci) = 0.0;
                    else if ((*T_data)(ci) >= d_solidus_temperature && (*T_data)(ci) <= d_liquidus_temperature)
                    {
                        const double denominator =
                            std::pow(1.0 + ((*T_data)(ci)-d_solidus_temperature) / mushy_interval *
                                               (d_rho_solid / d_rho_liquid - 1.0),
                                     2.0);
                        (*drho_dT_data)(ci) =
                            d_rho_solid / d_rho_liquid * (d_rho_liquid - d_rho_solid) / (denominator * mushy_interval);
                    }
                    else
                        (*drho_dT_data)(ci) = 0.0;
                }
                else
                    (*drho_dT_data)(ci) = 0.0;
            }
        }
    }
    return;
} // computeDensityDerivative

void
IEPSemiImplicitHierarchyIntegrator::computeLiquidFraction(int lf_idx, const int h_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const double h_s = d_cp_solid * d_solidus_temperature;
    const double cp_avg = (d_cp_solid + d_cp_liquid) / 2.0;
    const double h_l = cp_avg * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat_temp;
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

                if ((*H_data)(ci) >= 0.5)
                {
                    if ((*h_data)(ci) < h_s)
                        (*lf_data)(ci) = 0.0;
                    else if ((*h_data)(ci) > h_l)
                        (*lf_data)(ci) = 1.0;
                    else
                        (*lf_data)(ci) =
                            d_rho_solid * (h_s - (*h_data)(ci)) /
                            ((d_rho_liquid - d_rho_solid) * (*h_data)(ci)-d_rho_liquid * h_l + d_rho_solid * h_s);
                }
                else
                {
                    (*lf_data)(ci) = (*H_data)(ci);
                }
            }
        }
    }
    return;
} // computeLiquidFraction

void
IEPSemiImplicitHierarchyIntegrator::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db);
#endif
    // Read in data members from input database.
    if (!is_from_restart)
    {
        d_latent_heat = input_db->getDouble("latent_heat");
        d_latent_heat_temp = input_db->getDouble("latent_heat_temp");
        d_rho_liquid = input_db->getDouble("rho_liquid");
        d_rho_solid = input_db->getDouble("rho_solid");
        d_cp_liquid = input_db->getDouble("cp_liquid");
        d_cp_solid = input_db->getDouble("cp_solid");
        d_cp_gas = input_db->getDouble("cp_gas");
        d_T_melt = input_db->getDouble("T_melt");

        d_M_lf = input_db->getDouble("M_lf");
        d_lambda_lf = input_db->getDouble("lambda_lf");
        d_eta_lf = input_db->getDouble("eta_lf");

        if (input_db->keyExists("solve_energy")) d_solve_energy = input_db->getBool("solve_energy");
        if (input_db->keyExists("solve_mass_conservation"))
            d_solve_mass_conservation = input_db->getBool("solve_mass_conservation");

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

        if (input_db->keyExists("T_diffusion_time_stepping_type"))
            d_T_diffusion_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("T_diffusion_time_stepping_type"));

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

        if (input_db->keyExists("free_parameter")) d_free_parameter = input_db->getDouble("free_parameter");
        if (input_db->keyExists("eps")) d_eps = input_db->getDouble("eps");
        if (input_db->keyExists("beta")) d_beta = input_db->getDouble("beta");
        if (input_db->keyExists("H_diffusion_coefficient"))
            d_H_diffusion_coefficient = input_db->getDouble("H_diffusion_coefficient");
        if (input_db->keyExists("apply_brinkman")) d_apply_brinkman = input_db->getBool("apply_brinkman");
        if (input_db->keyExists("lf_brinkman_db"))
        {
            d_lf_brinkman_db = input_db->getDatabase("lf_brinkman_db");
            //        if (d_brinkman_db) d_brinkman_type = d_brinkman_db->getString("brinkman_type");
            d_lf_b = d_lf_brinkman_db->getDouble("lf_b");
        }
        if (input_db->keyExists("add_diffusion")) d_add_diffusion = input_db->getBool("add_diffusion");

        d_liquidus_temperature = input_db->getDouble("liquidus_temperature");
        d_solidus_temperature = input_db->getDouble("solidus_temperature");
        if (input_db->keyExists("compute_T_material_derivative"))
            d_compute_T_material_derivative = input_db->getBool("compute_T_material_derivative");
        //        if (input_db->keyExists("num_cycles"))
        //                   d_num_cycles = input_db->getInteger("num_cycles");
        if (input_db->keyExists("max_inner_iterations"))
            d_max_inner_iterations = input_db->getInteger("max_inner_iterations");

        if (input_db->keyExists("lf_iteration_error_tolerance"))
            d_lf_iteration_error_tolerance = input_db->getDouble("lf_iteration_error_tolerance");
    }
}

void
IEPSemiImplicitHierarchyIntegrator::getFromRestart()
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

    d_T_diffusion_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_T_diffusion_time_stepping_type"));
    d_T_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_T_convective_time_stepping_type"));
    d_T_init_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_T_init_convective_time_stepping_type"));
    d_T_convective_difference_form =
        string_to_enum<ConvectiveDifferencingType>(db->getString("d_T_convective_difference_form"));
    d_T_convective_op_type = db->getString("d_T_convective_op_type");
    d_T_convective_op_input_db = db->getDatabase("d_T_convective_op_db");

    d_latent_heat = db->getDouble("d_latent_heat");
    d_latent_heat_temp = db->getDouble("d_latent_heat_temp");
    d_rho_liquid = db->getDouble("d_rho_liquid");
    d_rho_solid = db->getDouble("d_rho_solid");
    d_cp_liquid = db->getDouble("d_cp_liquid");
    d_cp_solid = db->getDouble("d_cp_solid");
    d_cp_gas = db->getDouble("d_cp_gas");
    d_T_melt = db->getDouble("d_T_melt");

    d_M_lf = db->getDouble("d_M_lf");
    d_lambda_lf = db->getDouble("d_lambda_lf");
    d_eta_lf = db->getDouble("d_eta_lf");
    d_solve_energy = db->getBool("d_solve_energy");
    d_solve_mass_conservation = db->getBool("d_solve_mass_conservation");

    d_free_parameter = db->getDouble("d_free_parameter");
    d_eps = db->getDouble("d_eps");
    d_beta = db->getDouble("d_beta");
    d_H_diffusion_coefficient = db->getDouble("d_H_diffusion_coefficient");
    d_apply_brinkman = db->getBool("d_apply_brinkman");
    if (d_apply_brinkman)
    {
        d_lf_brinkman_db = db->getDatabase("d_lf_brinkman_db");
        d_lf_b = d_lf_brinkman_db->getDouble("d_lf_b");
    }
    d_add_diffusion = db->getBool("d_add_diffusion");

    d_liquidus_temperature = db->getDouble("d_liquidus_temperature");
    d_solidus_temperature = db->getDouble("d_solidus_temperature");
    d_compute_T_material_derivative = db->getBool("d_compute_T_material_derivative");
    d_max_inner_iterations = db->getInteger("d_max_inner_iterations");
    d_lf_iteration_error_tolerance = db->getDouble("d_lf_iteration_error_tolerance");
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
