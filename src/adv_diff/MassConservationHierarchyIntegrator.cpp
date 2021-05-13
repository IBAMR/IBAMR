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
#include "ibamr/MassConservationHierarchyIntegrator.h"
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

MassConservationHierarchyIntegrator::MassConservationHierarchyIntegrator(const std::string& object_name,
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

    d_rho_p_integrator = new MassIntegrator("MassIntegrator", input_db->getDatabase("mass_transport_integrator_db"));

    return;
} // MassConservationHierarchyIntegrator

void
MassConservationHierarchyIntegrator::registerDensityVariable(Pointer<CellVariable<NDIM, double> > rho_var,
                                                             const bool output_rho)
{
    d_rho_var = rho_var;
    d_rho_output = output_rho;

    return;
} // registerDensityVariable

void
MassConservationHierarchyIntegrator::registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx)
{
    d_reset_rho_fcns.push_back(callback);
    d_reset_rho_fcns_ctx.push_back(ctx);
    return;
} // registerResetFluidDensityFcn

// void
// IEPSemiImplicitHierarchyIntegrator::registerINSVCStaggeredHierarchyIntegrator(
//    Pointer<INSVCStaggeredHierarchyIntegrator> ins_cons_hier_integrator)
//{
//    d_ins_hierarchy_integrator = ins_cons_hier_integrator;
//    return;
//} // registerINSVCStaggeredHierarchyIntegrator

void
MassConservationHierarchyIntegrator::setAdvectionVelocityLiquidFractionEquation(
    Pointer<FaceVariable<NDIM, double> > u_var)
{
    d_lf_u_var = u_var;

    return;
} // setAdvectionVelocityLiquidFractionEquation

void
MassConservationHierarchyIntegrator::setDensityInitialCondition(Pointer<CellVariable<NDIM, double> > rho_var,
                                                                Pointer<IBTK::CartGridFunction> rho_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_var);
#endif
    d_rho_init = rho_init;
    return;
} // setDensityInitialCondition

void
MassConservationHierarchyIntegrator::registerMassDensityBoundaryConditions(RobinBcCoefStrategy<NDIM>*& rho_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_bc_coef = rho_bc_coef;
    return;
} // registerMassDensityBoundaryConditions

void
MassConservationHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
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
    const IntVector<NDIM> face_ghosts = FACEG;
    const IntVector<NDIM> no_ghosts = NOGHOSTS;

    int lf_u_current_idx, lf_u_scratch_idx, lf_u_new_idx;
    if (d_lf_u_var)
        registerVariable(lf_u_current_idx,
                         lf_u_new_idx,
                         lf_u_scratch_idx,
                         d_lf_u_var,
                         face_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");

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

    d_U_old_var = new FaceVariable<NDIM, double>(d_object_name + "::U_old");
    registerVariable(d_U_old_current_idx,
                     d_U_old_new_idx,
                     d_U_old_scratch_idx,
                     d_U_old_var,
                     face_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    if (d_solve_mass_conservation)
    {
        // Set various objects with conservative time integrator.
        d_rho_p_integrator->setCellCenteredDensityBoundaryConditions(d_rho_bc_coef);
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
MassConservationHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
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
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
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

    if (d_solve_mass_conservation)
    {
        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());

        // Keep track of the time-lagged velocity, specific heat and temperature.
        d_hier_fc_data_ops->copyData(d_U_old_new_idx, lf_u_current_idx);

        d_rho_p_integrator->setSolutionTime(current_time);
        d_rho_p_integrator->setTimeInterval(current_time, new_time);

        // For conservative discretization, an approximation to rho^{n+1}
        // will be computed from rho^{n}, which requires additional options to be
        // set.

        // Set the rho^{n} density
        d_rho_p_integrator->setCellCenteredDensityPatchDataIndex(rho_current_idx);

        // Data for the conservative time integrator is for cycle 0
        const int cycle_num = 0;
        d_rho_p_integrator->setCycleNumber(cycle_num);

        // Set the velocities used to update the density and the previous time step
        // size
        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ -1, /*current*/ lf_u_current_idx, /*new*/ -1);
        }
        else
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ d_U_old_current_idx, /*current*/ lf_u_current_idx, /*new*/ -1);
            d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        // Integrate density and convective term of energy equation.
        d_rho_p_integrator->integrate(dt);
    }

    d_hier_cc_data_ops->copyData(rho_new_idx, rho_current_idx);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
MassConservationHierarchyIntegrator::integrateHierarchy(const double current_time,
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

    int rho_new_idx, rho_scratch_idx, rho_current_idx;
    rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getNewContext());
    rho_scratch_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getScratchContext());
    rho_current_idx = var_db->mapVariableAndContextToIndex(d_rho_var, getCurrentContext());

    // In the special case of a conservative discretization form, the updated
    // density is calculated by application of the mass and convective
    // momentum integrator.
    // Update N_idx if necessary
    if (cycle_num > 0)
    {
        const double dt = new_time - current_time;
        const double half_time = current_time + 0.5 * dt;
        d_rho_p_integrator->setSolutionTime(half_time);

        // Set the cycle number
        d_rho_p_integrator->setCycleNumber(cycle_num);

        // Always set to current because we want to update rho^{n} to rho^{n+1}
        d_rho_p_integrator->setCellCenteredDensityPatchDataIndex(rho_current_idx);

        const int lf_u_current_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getCurrentContext());
        const int lf_u_new_idx = var_db->mapVariableAndContextToIndex(d_lf_u_var, getNewContext());

        // Set the velocities used to update the density
        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ -1, /*current*/ lf_u_current_idx, /*new*/ lf_u_new_idx);
        }
        else
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ d_U_old_current_idx,
                /*current*/ lf_u_current_idx,
                /*new*/ lf_u_new_idx);

            d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        d_rho_p_integrator->integrate(dt);
    }

    const int rho_cc_new_idx = d_rho_p_integrator->getUpdatedCellCenteredDensityPatchDataIndex();
    d_div_U_F_idx = d_rho_p_integrator->getRHSOfContinuityEquation();
    d_hier_cc_data_ops->copyData(rho_new_idx,
                                 rho_cc_new_idx,
                                 /*interior_only*/ true);

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
MassConservationHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                   const double new_time,
                                                                   const bool skip_synchronize_new_state_data,
                                                                   const int num_cycles)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    AdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

int
MassConservationHierarchyIntegrator::getRHSOfContinuityEquation()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_div_U_F_idx >= 0);
#endif
    return d_div_U_F_idx;
} // getRHSOfContinuityEquation

void
MassConservationHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    AdvDiffSemiImplicitHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
        base_hierarchy, coarsest_level, finest_level);

    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();

    // Reset the solution and rhs vectors.
    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    d_rho_p_integrator->setHierarchyMathOps(d_hier_math_ops);
    d_rho_p_integrator->initializeTimeIntegrator(base_hierarchy);

    return;
}

void
MassConservationHierarchyIntegrator::getFromInput(Pointer<Database> input_db, bool is_from_restart)
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
