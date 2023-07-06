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
#include "ibamr/AdvDiffConservativeMassScalarTransportRKIntegrator.h"
#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/PhaseChangeHierarchyIntegrator.h"
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
// Version of PhaseChangeHierarchyIntegrator restart file data.
static const int PC_HIERARCHY_INTEGRATOR_VERSION = 4;

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

PhaseChangeHierarchyIntegrator::PhaseChangeHierarchyIntegrator(const std::string& object_name,
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

    // Initialize conservative mass and transported quantity integrator.
    if (d_solve_mass_conservation)
        d_rho_p_integrator = new AdvDiffConservativeMassScalarTransportRKIntegrator(
            "AdvDiffConservativeMassScalarTransportRKIntegrator::MassTransportIntegrator",
            input_db->getDatabase("mass_scalar_transport_integrator_db"));

    if (!(d_T_convective_op_type == "CUI" || d_T_convective_op_type == "PPM"))
    {
        TBOX_ERROR(d_object_name << "::PhaseChangeHierarchyIntegrator():\n"
                                 << " current implementation supports only\n"
                                 << " CUI and PPM convective limiters for energy equation\n");
    }

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
        TBOX_ERROR(d_object_name << "::PhaseChangeHierarchyIntegrator():\n"
                                 << "  unsupported conductivity interpolation type: "
                                 << IBTK::enum_to_string<VCInterpType>(d_k_vc_interp_type) << " \n"
                                 << "  valid choices are: VC_HARMONIC_INTERP, VC_AVERAGE_INTERP\n");
    }

    return;
} // PhaseChangeHierarchyIntegrator

void
PhaseChangeHierarchyIntegrator::registerSpecificHeatVariable(Pointer<CellVariable<NDIM, double> > Cp_var,
                                                             const bool output_Cp)
{
    d_Cp_var = Cp_var;
    d_output_Cp = output_Cp;

    return;
} // registerSpecificHeatVariable

void
PhaseChangeHierarchyIntegrator::registerDensityVariable(Pointer<CellVariable<NDIM, double> > rho_var,
                                                        const bool output_rho)
{
    d_rho_var = rho_var;
    d_output_rho = output_rho;

    return;
} // registerDensityVariable

void
PhaseChangeHierarchyIntegrator::registerResetDensityFcn(ResetPhasePropertiesFcnPtr callback, void* ctx)
{
    d_reset_rho_fcns.push_back(callback);
    d_reset_rho_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidDensityFcn

void
PhaseChangeHierarchyIntegrator::registerResetSpecificHeatFcn(ResetPhasePropertiesFcnPtr callback, void* ctx)
{
    d_reset_Cp_fcns.push_back(callback);
    d_reset_Cp_fcns_ctx.push_back(ctx);
    return;
} // registerResetSpecificHeatFcn

void
PhaseChangeHierarchyIntegrator::registerResetDiffusionCoefficientFcn(ResetPhasePropertiesFcnPtr callback, void* ctx)
{
    d_reset_kappa_fcns.push_back(callback);
    d_reset_kappa_fcns_ctx.push_back(ctx);
    return;
} // registerResetDiffusionCoefficientFcn

void
PhaseChangeHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                              Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffSemiImplicitHierarchyIntegrator.
    AdvDiffSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Operators and solvers are maintained for each variable registered with the
    // integrator.
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

    d_T_solver = getEnergyEquationHelmholtzSolver(d_T_var);
    d_T_rhs_op = getEnergyEquationHelmholtzRHSOperator(d_T_var);

    // Register additional variables required for present time stepping algorithm.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> face_ghosts = FACEG;
    const IntVector<NDIM> no_ghosts = NOGHOSTS;

    registerVariable(d_lf_current_idx,
                     d_lf_new_idx,
                     d_lf_scratch_idx,
                     d_lf_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_lf_init);

    if (d_lf_gradient_var)
        registerVariable(d_lf_gradient_current_idx,
                         d_lf_gradient_new_idx,
                         d_lf_gradient_scratch_idx,
                         d_lf_gradient_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    registerVariable(d_T_current_idx,
                     d_T_new_idx,
                     d_T_scratch_idx,
                     d_T_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_T_init);

    // T_F contains forcing term of the temperature equation
    int T_F_current_idx, T_F_scratch_idx, T_F_new_idx;
    registerVariable(T_F_current_idx,
                     T_F_new_idx,
                     T_F_scratch_idx,
                     d_T_F_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    // T_diff contains the diffusion coefficient of the temperature equation
    int T_diff_coef_current_idx, T_diff_coef_new_idx, T_diff_coef_scratch_idx;
    registerVariable(T_diff_coef_current_idx,
                     T_diff_coef_new_idx,
                     T_diff_coef_scratch_idx,
                     d_T_diffusion_coef_var,
                     face_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    int T_diff_coef_cc_current_idx, T_diff_coef_cc_new_idx, T_diff_coef_cc_scratch_idx;
    registerVariable(T_diff_coef_cc_current_idx,
                     T_diff_coef_cc_new_idx,
                     T_diff_coef_cc_scratch_idx,
                     d_T_diffusion_coef_cc_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    int u_current_idx, u_scratch_idx, u_new_idx;
    if (d_u_adv_var)
        registerVariable(u_current_idx,
                         u_new_idx,
                         u_scratch_idx,
                         d_u_adv_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

    int T_diffusion_coef_rhs_scratch_idx;
    d_T_diffusion_coef_rhs_var = new SideVariable<NDIM, double>(d_T_var->getName() + "::Diff");
    registerVariable(T_diffusion_coef_rhs_scratch_idx, d_T_diffusion_coef_rhs_var, cell_ghosts, getScratchContext());

    int T_rhs_scratch_idx;
    registerVariable(T_rhs_scratch_idx, d_T_rhs_var, cell_ghosts, getScratchContext());

    // T_C contains the C coefficient of the temperature equation
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_T_C_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::C");
    d_T_C_idx = var_db->registerVariableAndContext(d_T_C_var, getCurrentContext(), no_ghosts);

    d_H_pre_var = new CellVariable<NDIM, double>("H::pre");
    d_H_pre_idx = var_db->registerVariableAndContext(d_H_pre_var, getCurrentContext(), no_ghosts);

    d_lf_pre_var = new CellVariable<NDIM, double>("lf_pre_var");
    d_lf_pre_idx = var_db->registerVariableAndContext(d_lf_pre_var, getCurrentContext());

    d_T_temp_rhs_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::temp_rhs");
    d_T_temp_rhs_idx = var_db->registerVariableAndContext(d_T_temp_rhs_var, getCurrentContext(), no_ghosts);

    d_T_N_var = new CellVariable<NDIM, double>(d_T_var->getName() + "::N");
    int T_N_scratch_idx;
    registerVariable(T_N_scratch_idx, d_T_N_var, cell_ghosts, getScratchContext());

    registerVariable(d_rho_current_idx,
                     d_rho_new_idx,
                     d_rho_scratch_idx,
                     d_rho_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE",
                     d_rho_init);

    registerVariable(d_Cp_current_idx,
                     d_Cp_new_idx,
                     d_Cp_scratch_idx,
                     d_Cp_var,
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

    d_C_rhs_scratch_idx = var_db->registerVariableAndContext(d_C_var, var_db->getContext("C_rhs"));

    d_D_cc_var = new CellVariable<NDIM, double>("D_cc", NDIM);
    registerVariable(d_D_cc_current_idx,
                     d_D_cc_new_idx,
                     d_D_cc_scratch_idx,
                     d_D_cc_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    if (d_visit_writer && d_output_temp_k)
    {
        d_visit_writer->registerPlotQuantity("temperature_kappa", "VECTOR", d_D_cc_current_idx, 0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d == 0) d_visit_writer->registerPlotQuantity("temp_kappa_x", "SCALAR", d_D_cc_current_idx, d);
            if (d == 1) d_visit_writer->registerPlotQuantity("temp_kappa_y", "SCALAR", d_D_cc_current_idx, d);
            if (d == 2) d_visit_writer->registerPlotQuantity("temp_kappa_z", "SCALAR", d_D_cc_current_idx, d);
        }
    }

    d_U_old_var = new FaceVariable<NDIM, double>(d_object_name + "::U_old");
    registerVariable(d_U_old_current_idx,
                     d_U_old_new_idx,
                     d_U_old_scratch_idx,
                     d_U_old_var,
                     face_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    d_updated_rho_var = new CellVariable<NDIM, double>(d_object_name + "::updated_rho");
    registerVariable(d_updated_rho_idx, d_updated_rho_var, no_ghosts, getCurrentContext());

    d_Div_U_F_var = new CellVariable<NDIM, double>(d_object_name + "::Div_U_F_var");
    registerVariable(d_Div_U_F_idx, d_Div_U_F_var, no_ghosts, getCurrentContext());

    // Register variables for plotting.
    if (d_visit_writer)
    {
        if (d_output_lf) d_visit_writer->registerPlotQuantity("liquid_fraction", "SCALAR", d_lf_current_idx, 0);

        if (d_output_lf_gradient)
            d_visit_writer->registerPlotQuantity("liquid_fraction_gradient", "SCALAR", d_lf_gradient_current_idx, 0);

        if (d_output_T) d_visit_writer->registerPlotQuantity("temperature", "SCALAR", d_T_current_idx, 0);

        if (d_output_rho) d_visit_writer->registerPlotQuantity("rho_cc", "SCALAR", d_rho_current_idx, 0);

        if (d_output_Cp) d_visit_writer->registerPlotQuantity("specific_heat", "SCALAR", d_Cp_current_idx, 0);

        if (d_output_Div_U_F) d_visit_writer->registerPlotQuantity("div_U_F", "SCALAR", d_Div_U_F_idx, 0);
    }
    if (d_solve_mass_conservation)
    {
        // Set various objects with conservative time integrator.
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
        rho_p_cc_integrator->setCellCenteredDensityBoundaryConditions(d_rho_bc_coef);
    }

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
PhaseChangeHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
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
        if (!level->checkAllocated(d_C_rhs_scratch_idx)) level->allocatePatchData(d_C_rhs_scratch_idx, current_time);
        if (!level->checkAllocated(d_T_C_idx)) level->allocatePatchData(d_T_C_idx, current_time);
        if (!level->checkAllocated(d_T_temp_rhs_idx)) level->allocatePatchData(d_T_temp_rhs_idx, current_time);
        if (!level->checkAllocated(d_H_pre_idx)) level->allocatePatchData(d_H_pre_idx, current_time);
        if (!level->checkAllocated(d_lf_pre_idx)) level->allocatePatchData(d_lf_pre_idx, current_time);
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

    const int Cp_current_idx = var_db->mapVariableAndContextToIndex(d_Cp_var, getCurrentContext());
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

    if (d_solve_mass_conservation)
    {
        const int u_current_idx = var_db->mapVariableAndContextToIndex(d_u_adv_var, getCurrentContext());

        // Keep track of the time-lagged velocity.
        d_hier_fc_data_ops->copyData(d_U_old_new_idx, u_current_idx);

        d_rho_p_integrator->setSolutionTime(current_time);
        d_rho_p_integrator->setTimeInterval(current_time, new_time);

        // For conservative discretization, an approximation to rho^{n+1}
        // will be computed from rho^{n}, which requires additional options to be
        // set.

        // Set the rho^{n} density
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
        d_rho_p_integrator->setDensityPatchDataIndex(rho_current_idx);

        // Set the convective derivative patch data index.
        const int T_N_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_N_var, getScratchContext());
        d_rho_p_integrator->setConvectiveDerivativePatchDataIndex(T_N_scratch_idx);

        // Data for the conservative time integrator is for cycle 0
        const int cycle_num = 0;
        d_rho_p_integrator->setCycleNumber(cycle_num);
    }

    d_hier_cc_data_ops->copyData(rho_new_idx, rho_current_idx);
    d_hier_cc_data_ops->copyData(Cp_new_idx, Cp_current_idx);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
PhaseChangeHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
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
        level->deallocatePatchData(d_C_rhs_scratch_idx);
        level->deallocatePatchData(d_T_C_idx);
        level->deallocatePatchData(d_T_temp_rhs_idx);
        level->deallocatePatchData(d_H_pre_idx);
        level->deallocatePatchData(d_lf_pre_idx);
    }

    if (d_lf_gradient_var)
    {
        // Ghost cell filling for liquid fraction
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getScratchContext());
        const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
        std::vector<RobinBcCoefStrategy<NDIM>*> H_bc_coef =
            getPhysicalBcCoefs(d_H_var); // Using H bc for now since I dont have lf_bc
        // in this class.
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        InterpolationTransactionComponent lf_bc_component(lf_scratch_idx,
                                                          lf_new_idx,
                                                          DATA_REFINE_TYPE,
                                                          USE_CF_INTERPOLATION,
                                                          DATA_COARSEN_TYPE,
                                                          d_bdry_extrap_type, // TODO: update variable name
                                                          CONSISTENT_TYPE_2_BDRY,
                                                          H_bc_coef);
        HierarchyGhostCellInterpolation lf_bdry_bc_fill_op;
        lf_bdry_bc_fill_op.initializeOperatorState(lf_bc_component, d_hierarchy);
        lf_bdry_bc_fill_op.fillData(new_time);

        // Find gradient of liquid fraction.
        d_hier_math_ops->grad(
            d_lf_gradient_new_idx, d_lf_gradient_var, 1.0, lf_scratch_idx, d_lf_var, nullptr, new_time);
    }

    AdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
PhaseChangeHierarchyIntegrator::registerLiquidFractionVariable(Pointer<CellVariable<NDIM, double> > lf_var,
                                                               const bool output_lf_var)
{
    d_lf_var = lf_var;
    d_output_lf = output_lf_var;

    // Set default values.
    d_lf_init = nullptr;
    return;
} // registerLiquidFractionVariable

void
PhaseChangeHierarchyIntegrator::registerLiquidFractionGradientVariable(
    Pointer<CellVariable<NDIM, double> > lf_gradient_var,
    const bool output_lf_gradient_var)
{
    d_lf_gradient_var = lf_gradient_var;
    d_output_lf_gradient = output_lf_gradient_var;

    return;
} // registerLiquidFractionGradientVariable

void
PhaseChangeHierarchyIntegrator::registerHeavisideVariable(Pointer<CellVariable<NDIM, double> > H_var)
{
    d_H_var = H_var;
    return;
} // registerHeavisideVariable

void
PhaseChangeHierarchyIntegrator::registerTemperatureVariable(Pointer<CellVariable<NDIM, double> > T_var,
                                                            const bool output_T_var)
{
    d_T_var = T_var;
    d_output_T = output_T_var;

    Pointer<CellDataFactory<NDIM, double> > T_factory = T_var->getPatchDataFactory();
    const int T_depth = T_factory->getDefaultDepth();
    Pointer<CellVariable<NDIM, double> > T_rhs_var =
        new CellVariable<NDIM, double>(T_var->getName() + "::T_rhs", T_depth);
    Pointer<CellVariable<NDIM, double> > T_F_var = new CellVariable<NDIM, double>(T_var->getName() + "::F", T_depth);
    Pointer<SideVariable<NDIM, double> > T_diff_coef_var =
        new SideVariable<NDIM, double>(T_var->getName() + "::diff_coef", T_depth);
    Pointer<CellVariable<NDIM, double> > T_diff_coef_cc_var =
        new CellVariable<NDIM, double>(T_var->getName() + "::diff_coef_cc", T_depth);

    // Set default values.
    d_u_adv_var = nullptr;
    d_T_F_var = T_F_var;
    d_T_rhs_var = T_rhs_var;
    d_T_diffusion_coef_var = T_diff_coef_var;
    d_T_diffusion_coef_cc_var = T_diff_coef_cc_var;
    d_T_init = nullptr;
    d_T_F_fcn = nullptr;
    d_T_bc_coef = nullptr;
    return;
} // registerTemperatureVariable

void
PhaseChangeHierarchyIntegrator::setLiquidFractionInitialCondition(Pointer<CellVariable<NDIM, double> > lf_var,
                                                                  Pointer<IBTK::CartGridFunction> lf_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(lf_var);
#endif
    d_lf_init = lf_init;
    return;
} // setLiquidFractionInitialCondition

void
PhaseChangeHierarchyIntegrator::setTemperatureInitialCondition(Pointer<CellVariable<NDIM, double> > T_var,
                                                               Pointer<IBTK::CartGridFunction> T_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(T_var);
#endif
    d_T_init = T_init;
    return;
} // setTemperatureInitialCondition

void
PhaseChangeHierarchyIntegrator::setDensityInitialCondition(Pointer<CellVariable<NDIM, double> > rho_var,
                                                           Pointer<IBTK::CartGridFunction> rho_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(rho_var);
#endif
    d_rho_init = rho_init;
    return;
} // setDensityInitialCondition

void
PhaseChangeHierarchyIntegrator::setTemperaturePhysicalBcCoef(Pointer<CellVariable<NDIM, double> > T_var,
                                                             RobinBcCoefStrategy<NDIM>* T_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(T_var);
#endif
    d_T_bc_coef = T_bc_coef;
    return;
} // setTemperaturePhysicalBcCoef

RobinBcCoefStrategy<NDIM>*
PhaseChangeHierarchyIntegrator::getTemperaturePhysicalBcCoef()
{
    return d_T_bc_coef;
} // getTemperaturePhysicalBcCoef

void
PhaseChangeHierarchyIntegrator::setEnergyEquationSourceTermFunction(Pointer<IBTK::CartGridFunction> T_F_fcn)
{
    d_T_F_fcn = T_F_fcn;
    return;
} // setTemperatureSourceTermFunction

int
PhaseChangeHierarchyIntegrator::getVelocityDivergencePatchDataIndex()
{
    return d_Div_U_F_idx;
} // getVelocityDivergencePatchDataIndex

void
PhaseChangeHierarchyIntegrator::registerMassDensityBoundaryConditions(RobinBcCoefStrategy<NDIM>*& rho_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_bc_coef = rho_bc_coef;
    return;
} // registerMassDensityBoundaryConditions

void
PhaseChangeHierarchyIntegrator::registerMassDensitySourceTerm(Pointer<CartGridFunction> S_fcn)
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
        TBOX_ERROR(d_object_name << "::PhaseChangeHierarchyIntegrator():\n"
                                 << " present implementation allows for only one mass density source\n"
                                 << " term to be set. Consider combining source terms into single "
                                    "CartGridFunction.\n");
    }
    return;
} // registerMassDensitySourceTerm

void
PhaseChangeHierarchyIntegrator::registerSpecificHeatBoundaryConditions(RobinBcCoefStrategy<NDIM>*& Cp_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_Cp_bc_coef = Cp_bc_coef;
    return;
} // registerSpecificHeatBoundaryConditions

void
PhaseChangeHierarchyIntegrator::registerThermalConductivityBoundaryConditions(RobinBcCoefStrategy<NDIM>*& k_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_k_bc_coef = k_bc_coef;
    return;
} // registerThermalConductivityBoundaryConditions

void
PhaseChangeHierarchyIntegrator::setAdvectionVelocity(Pointer<FaceVariable<NDIM, double> > u_var)
{
    d_u_adv_var = u_var;

    return;
} // setAdvectionVelocityTemperatureEquation

void
PhaseChangeHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    db->putInteger("PC_HIERARCHY_INTEGRATOR_VERSION", PC_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_T_diffusion_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_T_diffusion_time_stepping_type));
    db->putString("d_T_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_T_convective_time_stepping_type));
    db->putString("d_T_init_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_T_init_convective_time_stepping_type));
    db->putString("d_T_convective_difference_form",
                  enum_to_string<ConvectiveDifferencingType>(d_T_convective_difference_form));
    db->putString("d_T_convective_op_type", d_T_convective_op_type);
    d_T_convective_op_input_db = db->putDatabase("d_T_convective_op_db");

    db->putDouble("d_latent_heat", d_latent_heat);
    db->putDouble("d_rho_liquid", d_rho_liquid);
    db->putDouble("d_rho_solid", d_rho_solid);
    db->putDouble("d_T_melt", d_T_melt);
    db->putBool("d_solve_mass_conservation", d_solve_mass_conservation);
    db->putBool("d_output_Div_U_F", d_output_Div_U_F);
    db->putBool("d_output_temp_k", d_output_temp_k);

    AdvDiffSemiImplicitHierarchyIntegrator::putToDatabaseSpecialized(db);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PROTECTED ////////////////////////////////////

void
PhaseChangeHierarchyIntegrator::regridHierarchyBeginSpecialized()
{
    AdvDiffSemiImplicitHierarchyIntegrator::regridHierarchyBeginSpecialized();

    d_T_rhs_op->deallocateOperatorState();
    d_T_solver->deallocateSolverState();

    d_T_solver_needs_init = true;
    d_T_rhs_op_needs_init = true;
    d_T_convective_op_needs_init = true;

    return;
} // regridHierarchyBeginSpecialized

void
PhaseChangeHierarchyIntegrator::regridHierarchyEndSpecialized()
{
    AdvDiffSemiImplicitHierarchyIntegrator::regridHierarchyEndSpecialized();

    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    const int coarsest_hier_level = 0;
    d_hier_math_ops->setPatchHierarchy(d_hierarchy);
    d_hier_math_ops->resetLevels(coarsest_hier_level, finest_hier_level);

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

    int T_diff_coef_cc_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_T_diffusion_coef_cc_var, getScratchContext());
    InterpolationTransactionComponent k_bc_component(T_diff_coef_cc_scratch_idx,
                                                     DATA_REFINE_TYPE,
                                                     USE_CF_INTERPOLATION,
                                                     DATA_COARSEN_TYPE,
                                                     BDRY_EXTRAP_TYPE,
                                                     CONSISTENT_TYPE_2_BDRY,
                                                     d_k_bc_coef);
    d_k_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    d_k_bdry_bc_fill_op->initializeOperatorState(k_bc_component, d_hierarchy);

    // Reset the solution and rhs vectors.
    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, getScratchContext());
    d_T_sol = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::sol_vec::" + d_T_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_T_sol->addComponent(d_T_var, T_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    const int T_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_rhs_var, getScratchContext());
    d_T_rhs = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::rhs_vec::" + d_T_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_T_rhs->addComponent(d_T_rhs_var, T_rhs_scratch_idx, wgt_idx, d_hier_cc_data_ops);

    d_T_solver_needs_init = true;
    d_T_rhs_op_needs_init = true;
    d_T_convective_op_needs_init = true;

    if (d_solve_mass_conservation)
    {
        d_rho_p_integrator->setHierarchyMathOps(d_hier_math_ops);
        d_rho_p_integrator->initializeSTSIntegrator(d_hierarchy);
    }
    return;
} // regridHierarchyEndSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////
void
PhaseChangeHierarchyIntegrator::interpolateCCToSCSimpleAveraging(int sc_idx, const int cc_idx)
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
} // interpolateCCTOSCSimpleAveraging

void
PhaseChangeHierarchyIntegrator::interpolateCCToSCHarmonicAveraging(int sc_idx, const int cc_idx)
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
} // interpolateCCTOSCHarmonicAveraging

Pointer<PoissonSolver>
PhaseChangeHierarchyIntegrator::getEnergyEquationHelmholtzSolver(Pointer<CellVariable<NDIM, double> > T_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(T_var);
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
} // getEnergyEquationHelmholtzSolver

Pointer<LaplaceOperator>
PhaseChangeHierarchyIntegrator::getEnergyEquationHelmholtzRHSOperator(Pointer<CellVariable<NDIM, double> > T_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(T_var);
#endif
    const std::string& name = T_var->getName();
    if (!d_T_rhs_op)
    {
        d_T_rhs_op = new CCLaplaceOperator(d_object_name + "::helmholtz_rhs_op::" + name, /*homogeneous_bc*/ false);
        d_T_rhs_op_needs_init = true;
    }
    return d_T_rhs_op;
} // getEnergyEquationHelmholtzRHSOperator

void
PhaseChangeHierarchyIntegrator::boundLiquidFraction(int lf_new_idx)
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
PhaseChangeHierarchyIntegrator::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db);
#endif
    // Read in data members from input database.
    if (!is_from_restart)
    {
        d_latent_heat = input_db->getDouble("latent_heat");
        d_rho_liquid = input_db->getDouble("rho_liquid");
        d_rho_solid = input_db->getDouble("rho_solid");
        d_T_melt = input_db->getDouble("T_melt");

        if (input_db->keyExists("solve_mass_conservation"))
            d_solve_mass_conservation = input_db->getBool("solve_mass_conservation");

        if (input_db->keyExists("output_Div_U_F")) d_output_Div_U_F = input_db->getBool("output_Div_U_F");

        if (input_db->keyExists("output_temp_k")) d_output_temp_k = input_db->getBool("output_temp_k");

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
        if (input_db->keyExists("T_convective_time_stepping_type"))
            d_T_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(input_db->getString("T_convective_time_stepping_type"));
    }
}

void
PhaseChangeHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("PC_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != PC_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }

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
    d_rho_liquid = db->getDouble("d_rho_liquid");
    d_rho_solid = db->getDouble("d_rho_solid");
    d_T_melt = db->getDouble("d_T_melt");
    d_solve_mass_conservation = db->getBool("d_solve_mass_conservation");
    d_output_Div_U_F = db->getBool("d_output_Div_U_F");
    d_output_temp_k = db->getBool("d_output_temp_k");
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
