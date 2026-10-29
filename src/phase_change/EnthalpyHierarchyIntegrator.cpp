// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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
#include <ibamr/AdvDiffCUIConvectiveOperator.h>
#include <ibamr/AdvDiffConservativeMassScalarTransportRKIntegrator.h>
#include <ibamr/AdvDiffConvectiveOperatorManager.h>
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/EnthalpyHierarchyIntegrator.h>
#include <ibamr/ibamr_enums.h>
#include <ibamr/ibamr_utilities.h>

#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CartGridFunction.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LaplaceOperator.h>
#include <ibtk/PoissonSolver.h>

#include <tbox/Database.h>
#include <tbox/MathUtilities.h>
#include <tbox/MemoryDatabase.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>
#include <tbox/Utilities.h>

#include <BasePatchHierarchy.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellDataFactory.h>
#include <CellVariable.h>
#include <FaceData.h>
#include <FaceVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchyDataOpsManager.h>
#include <HierarchyFaceDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchFaceDataOpsReal.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <string>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

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

// FORTRAN ROUTINES
#if (NDIM == 2)
#define SC_NORMAL_FC IBAMR_FC_FUNC_(sc_normal_2d, SC_NORMAL_2D)
#define NAVIER_STOKES_SIDE_TO_FACE_FC IBAMR_FC_FUNC_(navier_stokes_side_to_face2d, NAVIER_STOKES_SIDE_TO_FACE2D)
#endif

#if (NDIM == 3)
#define SC_NORMAL_FC IBAMR_FC_FUNC_(sc_normal_3d, SC_NORMAL_3D)
#define NAVIER_STOKES_SIDE_TO_FACE_FC IBAMR_FC_FUNC_(navier_stokes_side_to_face3d, NAVIER_STOKES_SIDE_TO_FACE3D)
#endif

extern "C"
{
    void SC_NORMAL_FC(double* N00,
                      double* N01,
#if (NDIM == 3)
                      double* N02,
#endif
                      double* N10,
                      double* N11,
#if (NDIM == 3)
                      double* N12,
                      double* N20,
                      double* N21,
                      double* N22,
#endif
                      const int& N_gcw,
                      const double* U,
                      const int& U_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1,
#if (NDIM == 3)
                      const int& ilower2,
                      const int& iupper2,
#endif
                      const double* dx);

    void NAVIER_STOKES_SIDE_TO_FACE_FC(
#if (NDIM == 2)
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        const double*,
        const int&,
        double*,
        double*,
        const int&
#endif
#if (NDIM == 3)
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        const double*,
        const double*,
        const int&,
        double*,
        double*,
        double*,
        const int&
#endif
    );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
// Version of EnthalpyHierarchyIntegrator restart file data.
static const int ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION = 6;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
static const int FACEG = 1;
static const int NOGHOSTS = 0;

static const double H_LIM = 0.5;

const auto clamp = [](const double x, const double x_min, const double x_max)
{ return std::max(x_min, std::min(x, x_max)); };

// Copy data from a side-centered variable to a face-centered variable.
void
copy_side_to_face(const int U_fc_idx, const int U_sc_idx, Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const hier::Index<NDIM>& ilower = patch->getBox().lower();
            const hier::Index<NDIM>& iupper = patch->getBox().upper();
            Pointer<SideData<NDIM, double>> U_sc_data = patch->getPatchData(U_sc_idx);
            Pointer<FaceData<NDIM, double>> U_fc_data = patch->getPatchData(U_fc_idx);
#if !defined(NDEBUG)
            TBOX_ASSERT(U_sc_data->getGhostCellWidth().min() == U_sc_data->getGhostCellWidth().max());
            TBOX_ASSERT(U_fc_data->getGhostCellWidth().min() == U_fc_data->getGhostCellWidth().max());
#endif
            const int U_sc_gcw = U_sc_data->getGhostCellWidth().max();
            const int U_fc_gcw = U_fc_data->getGhostCellWidth().max();
            NAVIER_STOKES_SIDE_TO_FACE_FC(ilower(0),
                                          iupper(0),
                                          ilower(1),
                                          iupper(1),
#if (NDIM == 3)
                                          ilower(2),
                                          iupper(2),
#endif
                                          U_sc_data->getPointer(0),
                                          U_sc_data->getPointer(1),
#if (NDIM == 3)
                                          U_sc_data->getPointer(2),
#endif
                                          U_sc_gcw,
                                          U_fc_data->getPointer(0),
                                          U_fc_data->getPointer(1),
#if (NDIM == 3)
                                          U_fc_data->getPointer(2),
#endif
                                          U_fc_gcw);
        }
    }
    return;
} // copy_side_to_face
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

EnthalpyHierarchyIntegrator::EnthalpyHierarchyIntegrator(const std::string& object_name,
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

    return;
} // EnthalpyHierarchyIntegrator

void
EnthalpyHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                           Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_integrator_is_initialized) return;

    PhaseChangeHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> side_ghosts = SIDEG;
    const IntVector<NDIM> no_ghosts = NOGHOSTS;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    registerVariable(d_h_current_idx,
                     d_h_new_idx,
                     d_h_scratch_idx,
                     d_h_var,
                     cell_ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");

    if (d_visit_writer) d_visit_writer->registerPlotQuantity("enthalpy", "SCALAR", d_h_current_idx);

    if (d_lf_extrap_var)
    {
        registerVariable(d_lf_extrap_current_idx,
                         d_lf_extrap_new_idx,
                         d_lf_extrap_scratch_idx,
                         d_lf_extrap_var,
                         cell_ghosts,
                         "CONSERVATIVE_COARSEN",
                         "CONSERVATIVE_LINEAR_REFINE");
        if (d_visit_writer) d_visit_writer->registerPlotQuantity("lf_extrap", "SCALAR", d_lf_extrap_current_idx);

        const IntVector<NDIM> face_ghosts = FACEG;

        d_lf_extrap_rhs_var = new CellVariable<NDIM, double>(d_object_name + "::lf_extrap_rhs_var");
        registerVariable(d_lf_extrap_rhs_scratch_idx, d_lf_extrap_rhs_var, cell_ghosts, getScratchContext());

        d_u_adv_fc_lf_extrap_var = new FaceVariable<NDIM, double>(d_object_name + "::u_adv_fc_lf_extrap_var");
        d_u_adv_fc_lf_extrap_current_idx =
            var_db->registerVariableAndContext(d_u_adv_fc_lf_extrap_var, getCurrentContext(), face_ghosts);

        d_u_adv_sc_lf_extrap_var = new SideVariable<NDIM, double>(d_object_name + "::u_adv_sc_lf_extrap_var");
        d_u_adv_sc_lf_extrap_current_idx =
            var_db->registerVariableAndContext(d_u_adv_sc_lf_extrap_var, getCurrentContext(), side_ghosts);

        d_normal_lf_extrap_var =
            new SideVariable<NDIM, double>(d_object_name + "::normal_lf_extrap_var", /*depth*/ NDIM);
        d_normal_lf_extrap_current_idx =
            var_db->registerVariableAndContext(d_normal_lf_extrap_var, getCurrentContext(), /*gcw*/ IntVector<NDIM>(2));

        d_lf_extrap_interp_var = new FaceVariable<NDIM, double>(d_lf_extrap_var->getName() + "::interp");
        d_lf_extrap_interp_idx =
            var_db->registerVariableAndContext(d_lf_extrap_interp_var, getCurrentContext(), no_ghosts);
    }

    d_dT_dh_var = new CellVariable<NDIM, double>(d_object_name + "::dT_dh_var");
    d_dT_dh_scratch_idx = var_db->registerVariableAndContext(d_dT_dh_var, getScratchContext(), no_ghosts);

    d_delta_h_var = new CellVariable<NDIM, double>(d_object_name + "::delta_h_var");
    d_delta_h_idx = var_db->registerVariableAndContext(d_delta_h_var, getScratchContext(), cell_ghosts);

    d_delta_h_rhs_var = new CellVariable<NDIM, double>(d_object_name + "::delta_h_rhs_var");
    d_delta_h_rhs_idx = var_db->registerVariableAndContext(d_delta_h_rhs_var, getScratchContext(), cell_ghosts);

    d_h_residual_var = new CellVariable<NDIM, double>(d_object_name + "::h_residual_var");
    d_h_residual_idx = var_db->registerVariableAndContext(d_h_residual_var, getScratchContext(), no_ghosts);

    d_h_C_var = new CellVariable<NDIM, double>(d_object_name + "::h_C_var");
    d_h_C_idx = var_db->registerVariableAndContext(d_h_C_var, getScratchContext(), no_ghosts);

    d_h_diffusion_coef_var = new SideVariable<NDIM, double>(d_object_name + "::h_diffusion_coef_var");
    d_h_diffusion_coef_idx =
        var_db->registerVariableAndContext(d_h_diffusion_coef_var, getScratchContext(), side_ghosts);

    if (d_solve_mass_conservation && d_u_adv_var)
    {
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_integrator = d_rho_p_integrator;
        rho_integrator->setCellCenteredTransportQuantityBoundaryConditions(d_h_bc_coef);
    }

    if (d_lf_extrap_var) d_lf_extrap_convective_op = getLiquidFractionExtrapConvectiveOperator(d_lf_extrap_var);

    // Create the boundary condition object for the enthalpy correction equation.
    d_h_correction_bc_coef = new EnthalpyCorrectionBcCoefs(d_object_name + "::enthalpy_correction_bc");

    // Create the level set to VOF converter object
    d_ls_to_vof_converter =
        new VCINSVOFUtilities::LevelSetToVOFConverter(d_object_name + "::vof_converter", this, d_ls_var, d_H_var);

    d_integrator_is_initialized = true;

} // initializeHierarchyIntegrator

void
EnthalpyHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                          const double new_time,
                                                          const int num_cycles)
{
    PhaseChangeHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    if (initial_time)
    {
        TBOX_ASSERT(d_start_time == current_time);
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);

        if (!level->checkAllocated(d_dT_dh_scratch_idx)) level->allocatePatchData(d_dT_dh_scratch_idx, current_time);
        if (!level->checkAllocated(d_delta_h_idx)) level->allocatePatchData(d_delta_h_idx, current_time);
        if (!level->checkAllocated(d_delta_h_rhs_idx)) level->allocatePatchData(d_delta_h_rhs_idx, current_time);
        if (!level->checkAllocated(d_h_residual_idx)) level->allocatePatchData(d_h_residual_idx, current_time);
        if (!level->checkAllocated(d_h_C_idx)) level->allocatePatchData(d_h_C_idx, current_time);
        if (!level->checkAllocated(d_h_diffusion_coef_idx))
            level->allocatePatchData(d_h_diffusion_coef_idx, current_time);

        if (d_lf_extrap_var)
        {
            if (!level->checkAllocated(d_u_adv_fc_lf_extrap_current_idx))
                level->allocatePatchData(d_u_adv_fc_lf_extrap_current_idx, current_time);
            if (!level->checkAllocated(d_u_adv_sc_lf_extrap_current_idx))
                level->allocatePatchData(d_u_adv_sc_lf_extrap_current_idx, current_time);
            if (!level->checkAllocated(d_normal_lf_extrap_current_idx))
                level->allocatePatchData(d_normal_lf_extrap_current_idx, current_time);
            if (!level->checkAllocated(d_lf_extrap_interp_idx))
                level->allocatePatchData(d_lf_extrap_interp_idx, current_time);
        }
    }

    if (initial_time)
    {
        const double apply_time = current_time;

        // Convert level set to VOF data
        d_ls_to_vof_converter->computeVOFFromLevelSet(apply_time, current_time, new_time);

        // Compute initial enthalpy and liquid fraction
        const int H_current_idx = var_db->mapVariableAndContextToIndex(d_H_var, getCurrentContext());
        computeEnthalpyBasedOnTemperature(d_h_current_idx, d_T_current_idx, H_current_idx);
        computeLiquidFractionBasedOnEnthalpy(d_lf_current_idx, d_h_current_idx, H_current_idx);

        // Update material properties based on initial enthalpy and liquid fraction.
        for (unsigned k = 0; k < d_reset_rho_fcns.size(); ++k)
        {
            d_reset_rho_fcns[k](d_rho_current_idx,
                                d_rho_var,
                                d_hier_math_ops,
                                -1 /*cycle_num*/,
                                apply_time,
                                current_time,
                                new_time,
                                d_reset_rho_fcns_ctx[k]);
        }

        for (unsigned k = 0; k < d_reset_specific_heat_fcns.size(); ++k)
        {
            d_reset_specific_heat_fcns[k](d_specific_heat_current_idx,
                                          d_specific_heat_var,
                                          d_hier_math_ops,
                                          -1 /*cycle_num*/,
                                          apply_time,
                                          current_time,
                                          new_time,
                                          d_reset_specific_heat_fcns_ctx[k]);
        }

        for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
        {
            d_reset_kappa_fcns[k](d_T_diffusion_coef_cc_current_idx,
                                  d_T_diffusion_coef_cc_var,
                                  d_hier_math_ops,
                                  -1,
                                  apply_time,
                                  current_time,
                                  new_time,
                                  d_reset_kappa_fcns_ctx[k]);
        }

        d_hier_cc_data_ops->copyData(d_rho_new_idx, d_rho_current_idx);
        d_hier_cc_data_ops->copyData(d_specific_heat_new_idx, d_specific_heat_current_idx);
        d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_new_idx, d_T_diffusion_coef_cc_current_idx);
    }

    if (d_solve_mass_conservation && d_u_adv_var)
    {
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_integrator = d_rho_p_integrator;
        rho_integrator->setTransportQuantityPatchDataIndices(d_h_current_idx, -1);

        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            rho_integrator->setFluidVelocityPatchDataIndices(-1, d_u_adv_current_idx, -1);
        }
        else
        {
            rho_integrator->setFluidVelocityPatchDataIndices(d_U_old_current_idx, d_u_adv_current_idx, -1);
            rho_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }
        rho_integrator->integrate(new_time - current_time);
    }

    d_hier_cc_data_ops->copyData(d_h_new_idx, d_h_current_idx);
    d_hier_cc_data_ops->copyData(d_T_new_idx, d_T_current_idx);
    d_hier_cc_data_ops->copyData(d_lf_new_idx, d_lf_current_idx);

    if (d_lf_extrap_var)
    {
        d_lf_extrap_convective_op->initializeOperatorState(*d_lf_extrap_sol, *d_lf_extrap_rhs);
    }

    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
} // preprocessIntegrateHierarchy

void
EnthalpyHierarchyIntegrator::integrateHierarchySpecialized(const double current_time,
                                                           const double new_time,
                                                           const int cycle_num)
{
    // First advance all quantities managed by the base advection-diffusion
    // integrator. In particular, this produces ls^{n+1} in the new context
    AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchySpecialized(current_time, new_time, cycle_num);

    // Convert level set to VOF data at the new time
    d_ls_to_vof_converter->computeVOFFromLevelSet(/*apply_time*/ new_time, current_time, new_time);

    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE({
            pout << "EnthalpyHierarchyIntegrator::integrateHierarchy():\n"
                 << "  WARNING: num_cycles = " << d_current_num_cycles
                 << " but expected num_cycles = " << expected_num_cycles << ".\n";
        });
    }

    // Update rho^n+1 and find the convective derivative term N^n+1.
    Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_integrator = d_rho_p_integrator;
    if (cycle_num > 0 && d_solve_mass_conservation && d_u_adv_var)
    {
        rho_p_integrator->setSolutionTime(half_time);
        rho_p_integrator->setCycleNumber(cycle_num);
        rho_p_integrator->setConvectiveDerivativePatchDataIndex(d_T_N_scratch_idx);
        rho_p_integrator->setDensityPatchDataIndex(d_rho_current_idx);
        rho_p_integrator->setTransportQuantityPatchDataIndices(d_h_current_idx, d_h_new_idx);

        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            rho_p_integrator->setFluidVelocityPatchDataIndices(-1, d_u_adv_current_idx, d_u_adv_new_idx);
        }
        else
        {
            rho_p_integrator->setFluidVelocityPatchDataIndices(
                d_U_old_current_idx, d_u_adv_current_idx, d_u_adv_new_idx);
            rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        rho_p_integrator->integrate(dt);
    }
    d_updated_rho_idx = d_rho_p_integrator ? d_rho_p_integrator->getUpdatedDensityPatchDataIndex() : d_rho_new_idx;
    d_hier_cc_data_ops->copyData(d_rho_new_idx, d_updated_rho_idx, /*interior_only*/ true);

    if (d_lf_extrap_var)
    {
        computeAdvectionVelocityForExtrapolation(d_u_adv_fc_lf_extrap_current_idx, current_time, new_time);
    }

    // Inner iterations for the Newton-Ralphson scheme.
    double h_correction_error = 1.0;
    int inner_iterations = 0;

    while (h_correction_error > d_newton_error_tolerance && inner_iterations < d_max_newton_iterations)
    {
        // beta = dT/dh.
        computeTemperatureDerivative(d_dT_dh_scratch_idx, d_h_new_idx, H_new_idx);

        // Build the nonlinear residual R_h(h^m).
        computeEnthalpyResidual(d_h_residual_idx, current_time, new_time);

        // Construct the temperature-to-enthalpy Robin correction condition.
        updateEnthalpyCorrectionBoundaryCondition(new_time);

        // Right-hand side is -R_h.
        d_hier_cc_data_ops->scale(d_delta_h_rhs_idx, -1.0, d_h_residual_idx);

        // Temporal Jacobian coefficient: rho / dt.
        d_hier_cc_data_ops->copyData(d_h_C_idx, d_rho_new_idx);
        d_hier_cc_data_ops->scale(d_h_C_idx, 1.0 / dt, d_h_C_idx);

        // Form the cell-centered enthalpy diffusion coefficient:
        // kappa_h = kappa_T * dT/dh.
        d_hier_cc_data_ops->multiply(
            d_T_diffusion_coef_cc_scratch_idx, d_T_diffusion_coef_cc_new_idx, d_dT_dh_scratch_idx);
        d_k_bdry_bc_fill_op->fillData(new_time);

        // Interpolate the cell-centered coefficient to side centers.
        if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
        {
            interpolateCCToSCSimpleAveraging(d_h_diffusion_coef_idx, d_T_diffusion_coef_cc_scratch_idx);
        }
        else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
        {
            interpolateCCToSCHarmonicAveraging(d_h_diffusion_coef_idx, d_T_diffusion_coef_cc_scratch_idx);
        }
        else
        {
            TBOX_ERROR("Unsupported diffusion coefficient interpolation type.");
        }
        d_hier_sc_data_ops->scale(d_h_diffusion_coef_idx, -1.0, d_h_diffusion_coef_idx);

        // Solve for the enthalpy correction delta_h.
        PoissonSpecifications h_correction_spec(d_object_name + "::enthalpy_correction_spec");
        h_correction_spec.setCPatchDataId(d_h_C_idx);
        h_correction_spec.setDPatchDataId(d_h_diffusion_coef_idx);

        Pointer<PoissonSolver> h_correction_solver = getEnergyEquationHelmholtzSolver(d_h_var);
        h_correction_solver->setPoissonSpecifications(h_correction_spec);
        h_correction_solver->setPhysicalBcCoef(d_h_correction_bc_coef.getPointer());
        h_correction_solver->setHomogeneousBc(false);
        h_correction_solver->setSolutionTime(new_time);
        h_correction_solver->setTimeInterval(current_time, new_time);
        h_correction_solver->initializeSolverState(*d_delta_h_sol,
                                                   *d_delta_h_rhs); // Coefs are changing; update each time
        d_hier_cc_data_ops->setToScalar(d_delta_h_idx, 0.0);
        h_correction_solver->solveSystem(*d_delta_h_sol, *d_delta_h_rhs);

        // Update the enthalpy iterate: h^{m+1} = h^m + damping * delta_h.
        const double h_old_norm = d_hier_cc_data_ops->L2Norm(d_h_new_idx, wgt_cc_idx);
        const double damping = 1.0;
        d_hier_cc_data_ops->axpy(d_h_new_idx, damping, d_delta_h_idx, d_h_new_idx);

        // Update the temperature iterate: T^{m+1} = T(h^{m+1}).
        computeTemperatureBasedOnEnthalpy(d_T_new_idx, d_h_new_idx, H_new_idx);

        // Copy the old liquid fraction, and compute its norm
        d_hier_cc_data_ops->copyData(d_lf_pre_idx, d_lf_new_idx);
        const double lf_old_norm = d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx);

        // Update the liquid fraction iterate: lf^{m+1} = lf(h^{m+1}).
        computeLiquidFractionBasedOnEnthalpy(d_lf_new_idx, d_h_new_idx, H_new_idx);

        // Update specific heat based on new enthalpy, temperature, and liquid fraction.
        const double apply_time = new_time;
        for (unsigned k = 0; k < d_reset_specific_heat_fcns.size(); ++k)
        {
            d_reset_specific_heat_fcns[k](d_specific_heat_new_idx,
                                          d_specific_heat_var,
                                          d_hier_math_ops,
                                          -1 /*cycle_num*/,
                                          apply_time,
                                          current_time,
                                          new_time,
                                          d_reset_specific_heat_fcns_ctx[k]);
        }

        // Update conductivity based on new enthalpy, temperature, and liquid fraction.
        for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
        {
            d_reset_kappa_fcns[k](d_T_diffusion_coef_cc_new_idx,
                                  d_T_diffusion_coef_cc_var,
                                  d_hier_math_ops,
                                  -1 /*cycle_num*/,
                                  apply_time,
                                  current_time,
                                  new_time,
                                  d_reset_kappa_fcns_ctx[k]);
        }

        // Compute norms of the liquid fraction correction and enthalpy correction
        d_hier_cc_data_ops->subtract(d_lf_pre_idx, d_lf_new_idx, d_lf_pre_idx);
        const double lf_correction_error = d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_old_norm);

        h_correction_error = d_hier_cc_data_ops->L2Norm(d_delta_h_idx, wgt_cc_idx) / (1.0 + h_old_norm);

        plog << "At Newton iteration: " << inner_iterations << "\n"
             << "Relative liquid fraction correction error: " << lf_correction_error << "\n"
             << "Relative enthalpy correction error:        " << h_correction_error << "\n";

        if (h_correction_error < d_newton_error_tolerance)
        {
            break;
        }

        ++inner_iterations;
    }

    computeDivergenceVelocitySourceTerm(d_Div_U_F_idx, new_time);

    if (d_lf_extrap_var)
    {
        extrapolateLiquidFractionToGasRegion(d_lf_new_idx);
    }
} // integrateHierarchySpecialized

void
EnthalpyHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                           const double new_time,
                                                           const bool skip_synchronize_new_state_data,
                                                           const int num_cycles)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Deallocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (d_lf_extrap_var)
        {
            level->deallocatePatchData(d_u_adv_fc_lf_extrap_current_idx);
            level->deallocatePatchData(d_u_adv_sc_lf_extrap_current_idx);
            level->deallocatePatchData(d_normal_lf_extrap_current_idx);
            level->deallocatePatchData(d_lf_extrap_interp_idx);
        }
    }

    PhaseChangeHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
EnthalpyHierarchyIntegrator::registerSpecificEnthalpyVariable(Pointer<CellVariable<NDIM, double>> h_var,
                                                              const bool output_h_var)
{
    d_h_var = h_var;
    d_output_h = output_h_var;

    return;
} // registerLiquidFractionVariable

void
EnthalpyHierarchyIntegrator::setEnthalpyBcCoef(RobinBcCoefStrategy<NDIM>* h_bc_coef)
{
    d_h_bc_coef = h_bc_coef;
    return;
} // setEnthalpyBcCoef

void
EnthalpyHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    db->putInteger("ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION", ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("specific_heat_liquid", d_specific_heat_liquid);
    db->putDouble("specific_heat_solid", d_specific_heat_solid);
    db->putDouble("specific_heat_gas", d_specific_heat_gas);
    db->putDouble("liquidus_temperature", d_liquidus_temperature);
    db->putDouble("solidus_temperature", d_solidus_temperature);
    db->putDouble("reference_temperature", d_reference_temperature);
    db->putDouble("specific_heat_mushy", d_specific_heat_mushy);
    db->putDouble("gas_liquid_fraction", d_gas_liquid_fraction);
    db->putInteger("max_newton_iterations", d_max_newton_iterations);
    db->putDouble("newton_error_tolerance", d_newton_error_tolerance);
    db->putBool("require_lf_extrapolation", d_require_lf_extrapolation);
    if (d_require_lf_extrapolation)
    {
        db->putInteger("lf_extrap_max_num_time_steps", d_lf_extrap_max_num_time_steps);
        db->putDouble("lf_extrap_cell_size", d_lf_extrap_cell_size);
    }

    PhaseChangeHierarchyIntegrator::putToDatabaseSpecialized(db);
    return;
} // putToDatabaseSpecialized

void
EnthalpyHierarchyIntegrator::registerLiquidFractionVariableForExtrapolation(
    Pointer<CellVariable<NDIM, double>> lf_extrap_var)
{
    d_lf_extrap_var = lf_extrap_var;
    return;
} // registerLiquidFractionVariableForExtrapolation

/////////////////////////////// PROTECTED //////////////////////////////////////

void
EnthalpyHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM>> base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    PhaseChangeHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
        base_hierarchy, coarsest_level, finest_level);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Construct the Newton correction solution vector.
    // The unknown is delta_h, stored in d_delta_h_idx.
    d_delta_h_sol = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::sol_vec::" + d_delta_h_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_delta_h_sol->addComponent(d_delta_h_var, d_delta_h_idx, wgt_idx, d_hier_cc_data_ops);

    // Construct the Newton correction right-hand-side vector.
    // The right-hand side is -R_h and is stored in d_delta_h_rhs_idx.
    d_delta_h_rhs = new SAMRAIVectorReal<NDIM, double>(
        d_object_name + "::rhs_vec::" + d_delta_h_rhs_var->getName(), d_hierarchy, 0, finest_hier_level);
    d_delta_h_rhs->addComponent(d_delta_h_rhs_var, d_delta_h_rhs_idx, wgt_idx, d_hier_cc_data_ops);

    // Reinitialize any additional hierarchy-dependent objects.
    d_delta_h_solver_needs_init = true;

    if (d_lf_extrap_var)
    {
        const int lf_extrap_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_extrap_var, getScratchContext());
        d_lf_extrap_sol = new SAMRAIVectorReal<NDIM, double>(
            d_object_name + "::sol_vec::" + d_lf_var->getName(), d_hierarchy, 0, finest_hier_level);
        d_lf_extrap_sol->addComponent(d_lf_extrap_var, lf_extrap_scratch_idx, wgt_idx, d_hier_cc_data_ops);

        const int lf_extrap_rhs_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_lf_extrap_rhs_var, getScratchContext());
        d_lf_extrap_rhs = new SAMRAIVectorReal<NDIM, double>(
            d_object_name + "::rhs_vec::" + d_lf_extrap_var->getName(), d_hierarchy, 0, finest_hier_level);
        d_lf_extrap_rhs->addComponent(d_lf_extrap_rhs_var, lf_extrap_rhs_scratch_idx, wgt_idx, d_hier_cc_data_ops);

        d_lf_extrap_convective_op_needs_init = true;
    }
    return;
} // resetHierarchyConfigurationSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
EnthalpyHierarchyIntegrator::computeDivergenceVelocitySourceTerm(int Div_U_F_idx, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double current_time = d_integrator_time;
    const double half_time = current_time + 0.5 * (new_time - current_time);

    const double delta_T = d_liquidus_temperature - d_solidus_temperature;
    const double delta_rho = d_rho_liquid - d_rho_solid;
    const double hp_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double B = d_latent_heat + d_specific_heat_mushy * delta_T;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());

    // Build the physical side-centered conductivity.
    // This is kappa_T, not the correction coefficient kappa_T*dT/dh.
    d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_scratch_idx, d_T_diffusion_coef_cc_new_idx);
    d_k_bdry_bc_fill_op->fillData(new_time);

    if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
    {
        interpolateCCToSCSimpleAveraging(d_T_diffusion_coef_rhs_scratch_idx, d_T_diffusion_coef_cc_scratch_idx);
    }
    else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
    {
        interpolateCCToSCHarmonicAveraging(d_T_diffusion_coef_rhs_scratch_idx, d_T_diffusion_coef_cc_scratch_idx);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::computeDivergenceVelocitySourceTerm():\n"
                                 << "Unsupported conductivity interpolation type.\n");
    }

    // C = 0 and D = +kappa give div(kappa*grad(T)).
    PoissonSpecifications T_rhs_op_spec(d_object_name + "::divergence_source_spec");
    T_rhs_op_spec.setCZero();
    T_rhs_op_spec.setDPatchDataId(d_T_diffusion_coef_rhs_scratch_idx);

    Pointer<LaplaceOperator> T_rhs_op = d_T_rhs_op;
    T_rhs_op->setPoissonSpecifications(T_rhs_op_spec);
    T_rhs_op->setPhysicalBcCoef(d_T_bc_coef);
    T_rhs_op->setHomogeneousBc(false);
    T_rhs_op->setSolutionTime(new_time);
    T_rhs_op->setTimeInterval(current_time, new_time);

    if (d_T_rhs_op_needs_init)
    {
        T_rhs_op->initializeOperatorState(*d_T_sol, *d_T_rhs);
        d_T_rhs_op_needs_init = false;
    }

    d_hier_cc_data_ops->copyData(d_T_scratch_idx, d_T_new_idx, false);
    T_rhs_op->apply(*d_T_sol, *d_T_rhs);

    // Temporarily store Q = div(kappa*grad(T)) + S in Div_U_F_idx.
    d_hier_cc_data_ops->copyData(Div_U_F_idx, d_T_rhs_scratch_idx);

    // Use the same forcing evaluation time as computeEnthalpyResidual().
    // For now we are not adding the source term. We need to test that effect.
    // d_hier_cc_data_ops->axpy(Div_U_F_idx, 1.0, d_T_F_scratch_idx, Div_U_F_idx);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(d_h_new_idx);
            Pointer<CellData<NDIM, double>> T_data = patch->getPatchData(d_T_new_idx);
            Pointer<CellData<NDIM, double>> lf_data = patch->getPatchData(d_lf_new_idx);
            Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double>> Div_U_F_data = patch->getPatchData(Div_U_F_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double alpha_p = clamp((*H_data)(ci), 0.0, 1.0);
                const double alpha_g = 1.0 - alpha_p;
                const double h_cell = (*h_data)(ci);
                const double Q_cell = (*Div_U_F_data)(ci);

                (*Div_U_F_data)(ci) = 0.0;

                if (alpha_p <= 0.0) continue;

                const double R0 = alpha_g * d_rho_gas + alpha_p * d_rho_solid;
                const double R1 = alpha_p * delta_rho;

                const double H0 =
                    alpha_g * d_rho_gas * d_specific_heat_gas * (d_solidus_temperature - d_reference_temperature) +
                    alpha_p * d_rho_solid * hp_s;

                const double H1 = alpha_g * d_rho_gas * d_specific_heat_gas * delta_T +
                                  alpha_p * (delta_rho * hp_s + d_rho_liquid * B);

                const double hs = H0 / R0;
                const double hl = (H0 + H1) / (R0 + R1);

                // No phase-change divergence outside the active interval.
                if (h_cell < hs || h_cell > hl) continue;

                const double T_cell = (*T_data)(ci);
                const double phi = clamp((T_cell - d_solidus_temperature) / delta_T, 0.0, 1.0);
                const double rho_p = (1.0 - phi) * d_rho_solid + phi * d_rho_liquid;

                const double denominator = alpha_g * d_rho_gas * d_specific_heat_gas * rho_p * delta_T +
                                           alpha_p * d_rho_solid * d_rho_liquid * B;

                (*Div_U_F_data)(ci) = alpha_p * (d_rho_solid - d_rho_liquid) * Q_cell / denominator;
            }
        }
    }

    return;
} // computeDivergenceVelocitySourceTerm

void
EnthalpyHierarchyIntegrator::computeEnthalpyBasedOnTemperature(int h_idx, const int T_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double delta_T = d_liquidus_temperature - d_solidus_temperature;
    const double delta_rho = d_rho_liquid - d_rho_solid;

    const double hp_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double B = d_latent_heat + d_specific_heat_mushy * delta_T;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double>> T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
            Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double alpha_p = clamp((*H_data)(ci), 0.0, 1.0);
                const double alpha_g = 1.0 - alpha_p;

                const double R0 = alpha_g * d_rho_gas + alpha_p * d_rho_solid;
                const double R1 = alpha_p * delta_rho;

                const double H0 =
                    alpha_g * d_rho_gas * d_specific_heat_gas * (d_solidus_temperature - d_reference_temperature) +
                    alpha_p * d_rho_solid * hp_s;
                const double H1 = alpha_g * d_rho_gas * d_specific_heat_gas * delta_T +
                                  alpha_p * (delta_rho * hp_s + d_rho_liquid * B);

                const double c_gs =
                    (alpha_g * d_rho_gas * d_specific_heat_gas + alpha_p * d_rho_solid * d_specific_heat_solid) / R0;

                const double c_gl =
                    (alpha_g * d_rho_gas * d_specific_heat_gas + alpha_p * d_rho_liquid * d_specific_heat_liquid) /
                    (R0 + R1);

                const double T = (*T_data)(ci);
                if (T < d_solidus_temperature)
                {
                    (*h_data)(ci) = (T - d_reference_temperature) * c_gs;
                }
                else if (T > d_liquidus_temperature)
                {
                    const double hl = (H0 + H1) / (R0 + R1);

                    (*h_data)(ci) = hl + c_gl * (T - d_liquidus_temperature);
                }
                else
                {
                    const double phi_p = clamp(
                        (T - d_solidus_temperature) / delta_T, 0.0, 1.0); // bounded liquid fraction within the PCM
                    (*h_data)(ci) = (H0 + H1 * phi_p) / (R0 + R1 * phi_p);
                }
            }
        }
    }
    return;
} // computeEnthalpyBasedOnTemperature

void
EnthalpyHierarchyIntegrator::computeTemperatureBasedOnEnthalpy(const int T_idx, const int h_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double delta_rho = d_rho_liquid - d_rho_solid;
    const double delta_T = d_liquidus_temperature - d_solidus_temperature;
    const double hp_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double B = d_latent_heat + d_specific_heat_mushy * delta_T;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
            Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double alpha_p = clamp((*H_data)(ci), 0.0, 1.0);
                const double alpha_g = 1.0 - alpha_p;

                const double R0 = alpha_g * d_rho_gas + alpha_p * d_rho_solid;
                const double R1 = alpha_p * delta_rho;

                const double H0 =
                    alpha_g * d_rho_gas * d_specific_heat_gas * (d_solidus_temperature - d_reference_temperature) +
                    alpha_p * d_rho_solid * hp_s;

                const double H1 = alpha_g * d_rho_gas * d_specific_heat_gas * delta_T +
                                  alpha_p * (delta_rho * hp_s + d_rho_liquid * B);

                const double c_gs =
                    (alpha_g * d_rho_gas * d_specific_heat_gas + alpha_p * d_rho_solid * d_specific_heat_solid) / R0;

                const double c_gl =
                    (alpha_g * d_rho_gas * d_specific_heat_gas + alpha_p * d_rho_liquid * d_specific_heat_liquid) /
                    (R0 + R1);

                const double hs = H0 / R0;
                const double hl = (H0 + H1) / (R0 + R1);

                const double h = (*h_data)(ci);
                double& T = (*T_data)(ci);
                if (h < hs)
                {
                    T = d_reference_temperature + h / c_gs;
                }
                else if (h > hl)
                {
                    T = d_liquidus_temperature + (h - hl) / c_gl;
                }
                else
                {
                    const double phi = clamp((h * R0 - H0) / (H1 - h * R1), 0.0, 1.0);
                    T = d_solidus_temperature + delta_T * phi;
                }
            }
        }
    }
    return;
} // computeTemperatureBasedOnEnthalpy

void
EnthalpyHierarchyIntegrator::computeTemperatureDerivative(const int dT_dh_idx, const int h_idx, const int H_idx)
{
    const double delta_T = d_liquidus_temperature - d_solidus_temperature;
    const double delta_rho = d_rho_liquid - d_rho_solid;
    const double hp_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double B = d_latent_heat + d_specific_heat_mushy * delta_T;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> dT_dh_data = patch->getPatchData(dT_dh_idx);
            Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
            Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double alpha_p = clamp((*H_data)(ci), 0.0, 1.0);
                const double alpha_g = 1.0 - alpha_p;

                const double R0 = alpha_g * d_rho_gas + alpha_p * d_rho_solid;
                const double R1 = alpha_p * delta_rho;

                const double H0 =
                    alpha_g * d_rho_gas * d_specific_heat_gas * (d_solidus_temperature - d_reference_temperature) +
                    alpha_p * d_rho_solid * hp_s;

                const double H1 = alpha_g * d_rho_gas * d_specific_heat_gas * delta_T +
                                  alpha_p * (delta_rho * hp_s + d_rho_liquid * B);

                const double c_gs =
                    (alpha_g * d_rho_gas * d_specific_heat_gas + alpha_p * d_rho_solid * d_specific_heat_solid) / R0;

                const double c_gl =
                    (alpha_g * d_rho_gas * d_specific_heat_gas + alpha_p * d_rho_liquid * d_specific_heat_liquid) /
                    (R0 + R1);

                const double hs = H0 / R0;
                const double hl = (H0 + H1) / (R0 + R1);

                const double h = (*h_data)(ci);
                double& beta = (*dT_dh_data)(ci);
                if (h < hs)
                {
                    beta = 1.0 / c_gs;
                }
                else if (h > hl)
                {
                    beta = 1.0 / c_gl;
                }
                else
                {
                    const double denominator = H1 - h * R1;
                    const double D = H1 * R0 - H0 * R1;
                    beta = delta_T * D / (denominator * denominator);
                }
            }
        }
    }
} // computeTemperatureDerivative

void
EnthalpyHierarchyIntegrator::computeEnthalpyResidual(const int residual_idx,
                                                     const double current_time,
                                                     const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    const bool advection_is_computed = (d_solve_mass_conservation && d_u_adv_var);

    // Account for forcing terms.
    if (d_T_F_fcn)
    {
        d_T_F_fcn->setDataOnPatchHierarchy(d_T_F_scratch_idx, d_T_F_var, d_hierarchy, half_time);
    }
    else
    {
        d_hier_cc_data_ops->setToScalar(d_T_F_scratch_idx, 0.0);
    }

    /* Compute the temporal, advective, and source contributions:
     *
     *      R_h = (rho_new h_new - rho_current h_current) / dt + N_h - F.
     *
     */
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> residual_data = patch->getPatchData(residual_idx);
            Pointer<CellData<NDIM, double>> rho_new_data = patch->getPatchData(d_rho_new_idx);
            Pointer<CellData<NDIM, double>> rho_current_data = patch->getPatchData(d_rho_current_idx);
            Pointer<CellData<NDIM, double>> h_new_data = patch->getPatchData(d_h_new_idx);
            Pointer<CellData<NDIM, double>> h_current_data = patch->getPatchData(d_h_current_idx);
            Pointer<CellData<NDIM, double>> F_data = patch->getPatchData(d_T_F_scratch_idx);
            Pointer<CellData<NDIM, double>> N_data = patch->getPatchData(d_T_N_scratch_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double temporal_term =
                    ((*rho_new_data)(ci) * (*h_new_data)(ci) - (*rho_current_data)(ci) * (*h_current_data)(ci)) / dt;

                const double advection_term = advection_is_computed ? (*N_data)(ci) : 0.0;
                (*residual_data)(ci) = temporal_term + advection_term - (*F_data)(ci);
            }
        }
    }

    // Construct kappa(T) for the residual diffusion operator.
    // The cell-centered coefficient is interpolated to side centers.
    d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_scratch_idx, d_T_diffusion_coef_cc_new_idx);
    d_k_bdry_bc_fill_op->fillData(new_time);

    if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
    {
        interpolateCCToSCSimpleAveraging(d_T_diffusion_coef_rhs_scratch_idx, d_T_diffusion_coef_cc_scratch_idx);
    }
    else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
    {
        interpolateCCToSCHarmonicAveraging(d_T_diffusion_coef_rhs_scratch_idx, d_T_diffusion_coef_cc_scratch_idx);
    }
    else
    {
        TBOX_ERROR("Unsupported diffusion coefficient interpolation type.");
    }

    PoissonSpecifications diffusion_spec(d_object_name + "::enthalpy_residual_diffusion_spec");
    diffusion_spec.setCZero();
    diffusion_spec.setDPatchDataId(d_T_diffusion_coef_rhs_scratch_idx);

    Pointer<LaplaceOperator> diffusion_op = d_T_rhs_op;
    diffusion_op->setPoissonSpecifications(diffusion_spec);
    diffusion_op->setPhysicalBcCoef(d_T_bc_coef);
    diffusion_op->setHomogeneousBc(false);
    diffusion_op->setSolutionTime(new_time);
    diffusion_op->setTimeInterval(current_time, new_time);

    // d_T_sol uses d_T_scratch_idx as its solution component.
    // Copy the current nonlinear temperature iterate into it.
    d_hier_cc_data_ops->copyData(d_T_scratch_idx, d_T_new_idx, false);

    /*
     * The operator writes
     *
     *     div(kappa grad(T_new))
     *
     * into d_T_rhs_scratch_idx through d_T_rhs.
     */
    if (d_T_rhs_op_needs_init)
    {
        diffusion_op->initializeOperatorState(*d_T_sol, *d_T_rhs);
        d_T_rhs_op_needs_init = false;
    }
    diffusion_op->apply(*d_T_sol, *d_T_rhs);

    /*
     * Add the diffusion contribution to the residual.
     *
     *  R_h =  (rho_new h_new - rho_current h_current) / dt + N_h -div(kappa grad(T_new)) - F.
     */
    d_hier_cc_data_ops->axpy(residual_idx, -1.0, d_T_rhs_scratch_idx, residual_idx);

    return;
} // computeEnthalpyResidual

void
EnthalpyHierarchyIntegrator::updateEnthalpyCorrectionBoundaryCondition(const double new_time)
{
    // Preserve the current temperature iterate for boundary conversion.
    d_hier_cc_data_ops->copyData(d_T_scratch_idx, d_T_new_idx);

    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    std::vector<InterpolationTransactionComponent> transaction_comps(1);
    transaction_comps[0] = InterpolationTransactionComponent(
        d_T_scratch_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, d_T_bc_coef);

    HierarchyGhostCellInterpolation T_bdry_fill_op;
    T_bdry_fill_op.initializeOperatorState(transaction_comps, d_hierarchy);
    T_bdry_fill_op.setHomogeneousBc(false);
    T_bdry_fill_op.fillData(new_time);

    // beta = dT/dh is sampled from the adjacent interior cell.
    // Its physical ghost values are not needed by the correction BC.
    d_h_correction_bc_coef->setTemperatureBcCoef(d_T_bc_coef);
    d_h_correction_bc_coef->setTemperatureAndDerivativePatchDataIndices(d_T_scratch_idx, d_dT_dh_scratch_idx);

    return;
} // updateEnthalpyCorrectionBoundaryCondition

void
EnthalpyHierarchyIntegrator::computeLiquidFractionBasedOnEnthalpy(const int lf_idx, const int h_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double delta_T = d_liquidus_temperature - d_solidus_temperature;
    const double delta_rho = d_rho_liquid - d_rho_solid;

    const double hp_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double B = d_latent_heat + d_specific_heat_mushy * delta_T;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
            Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double alpha_p = clamp((*H_data)(ci), 0.0, 1.0);
                const double alpha_g = 1.0 - alpha_p;

                if (alpha_p <= 0.0)
                {
                    (*lf_data)(ci) = 0.0;
                    continue;
                }

                const double R0 = alpha_g * d_rho_gas + alpha_p * d_rho_solid;
                const double R1 = alpha_p * delta_rho;

                const double H0 =
                    alpha_g * d_rho_gas * d_specific_heat_gas * (d_solidus_temperature - d_reference_temperature) +
                    alpha_p * d_rho_solid * hp_s;

                const double H1 = alpha_g * d_rho_gas * d_specific_heat_gas * delta_T +
                                  alpha_p * (delta_rho * hp_s + d_rho_liquid * B);

                const double hs = H0 / R0;
                const double hl = (H0 + H1) / (R0 + R1);
                const double h = (*h_data)(ci);

                if (h <= hs)
                {
                    (*lf_data)(ci) = 0.0;
                }
                else if (h >= hl)
                {
                    (*lf_data)(ci) = alpha_p;
                }
                else
                {
                    const double phi = (h * R0 - H0) / (H1 - h * R1);
                    (*lf_data)(ci) = alpha_p * clamp(phi, 0.0, 1.0);
                }
            }
        }
    }
    return;
} // computeLiquidFractionBasedOnEnthalpy

void
EnthalpyHierarchyIntegrator::extrapolateLiquidFractionToGasRegion(int lf_new_idx)
{
    const double CFL_extrap = 0.3; // CFL is fixed to be 0.3.
    const double dt_adv = CFL_extrap * d_lf_extrap_cell_size;
    const double H_eps = 1.0e-6; // PCM volume fraction threshold in a cell.

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double>> lf_data = patch->getPatchData(lf_new_idx);
            Pointer<CellData<NDIM, double>> lf_extrap_data = patch->getPatchData(d_lf_extrap_current_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double H = (*H_data)(ci);    // PCM volume fraction in the cell
                const double phi = (*lf_data)(ci); // Liquid fraction in the cell

                double phi_p = 0.0;
                if ((*H_data)(ci) > H_eps)
                {
                    // Find liquid fraction in the PCM volume from liquid fraction in the cell.
                    // Note: The extrapolated liquid fraction is used to turn on or off Brinkman and surface tension
                    // forces. Therefore, liquid fraction in the PCM volume is the correct quantity to use for that
                    // purpose.
                    phi_p = phi / H;
                    phi_p = clamp(phi_p, 0.0, 1.0);
                }

                (*lf_extrap_data)(ci) = phi_p;
            }
        }
    }

    // Initialize variables for extrapolation.
    d_hier_cc_data_ops->copyData(d_lf_extrap_scratch_idx, d_lf_extrap_current_idx);
    d_hier_cc_data_ops->copyData(d_lf_extrap_new_idx, d_lf_extrap_current_idx);
    d_hier_cc_data_ops->setToScalar(d_lf_extrap_rhs_scratch_idx, 0.0);

    int current_time_step = 0;
    while (current_time_step < d_lf_extrap_max_num_time_steps)
    {
        // Using CUI limiter, limit \varphi.
        d_lf_extrap_convective_op->interpolateToFaceOnHierarchy(
            d_lf_extrap_interp_idx, d_lf_extrap_scratch_idx, d_u_adv_fc_lf_extrap_current_idx, /*synch_cf_bdry*/ false);

        // Compute N = u . \grad{\varphi}.
        d_lf_extrap_convective_op->computeAdvectiveDerivativeOnHierarchy(d_lf_extrap_rhs_scratch_idx,
                                                                         d_lf_extrap_interp_idx,
                                                                         d_u_adv_fc_lf_extrap_current_idx,
                                                                         /*synch_cf_bdry*/ false);

        d_hier_cc_data_ops->scale(d_lf_extrap_rhs_scratch_idx, -dt_adv, d_lf_extrap_rhs_scratch_idx);
        d_hier_cc_data_ops->add(d_lf_extrap_rhs_scratch_idx, d_lf_extrap_rhs_scratch_idx, d_lf_extrap_current_idx);
        d_hier_cc_data_ops->copyData(d_lf_extrap_new_idx, d_lf_extrap_rhs_scratch_idx);
        current_time_step += 1;

        // For visit writing and for next iteration.
        d_hier_cc_data_ops->copyData(d_lf_extrap_scratch_idx, d_lf_extrap_new_idx);
        d_hier_cc_data_ops->copyData(d_lf_extrap_current_idx, d_lf_extrap_new_idx);
    }
    return;
} // extrapolateLiquidFractionToGasRegion

void
EnthalpyHierarchyIntegrator::computeAdvectionVelocityForExtrapolation(int u_adv_fc_lf_extrap_current_idx,
                                                                      const double current_time,
                                                                      const double new_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_new_idx = var_db->mapVariableAndContextToIndex(d_ls_var, getNewContext());
    const int ls_scratch_idx = var_db->mapVariableAndContextToIndex(d_ls_var, getScratchContext());

    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    const int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());

    // Fill the ghost cell values of level set at new_time.
    std::vector<RobinBcCoefStrategy<NDIM>*> ls_bc_coef = getPhysicalBcCoefs(d_ls_var);
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> interface_transaction(2);
    interface_transaction[0] = InterpolationTransactionComponent(ls_scratch_idx,
                                                                 ls_new_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 true,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 ls_bc_coef);

    std::vector<RobinBcCoefStrategy<NDIM>*> H_bc_coef = getPhysicalBcCoefs(d_H_var);
    interface_transaction[1] = InterpolationTransactionComponent(H_scratch_idx,
                                                                 H_new_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 true,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 H_bc_coef);

    Pointer<HierarchyGhostCellInterpolation> interface_fill_op = new HierarchyGhostCellInterpolation();
    interface_fill_op->initializeOperatorState(interface_transaction, d_hierarchy);
    interface_fill_op->fillData(new_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            Pointer<SideData<NDIM, double>> u_sc_data = patch->getPatchData(d_u_adv_sc_lf_extrap_current_idx);
            Pointer<SideData<NDIM, double>> normal_data = patch->getPatchData(d_normal_lf_extrap_current_idx);
            Pointer<CellData<NDIM, double>> ls_data = patch->getPatchData(ls_scratch_idx);
            Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_scratch_idx);

            // Compute normal = grad(phi_data)
            // phi is -ve in the gas region and +ve in the pcm region.
            // grad(phi) points from gas to pcm region.
            SC_NORMAL_FC(normal_data->getPointer(0, 0),
                         normal_data->getPointer(0, 1),
#if (NDIM == 3)
                         normal_data->getPointer(0, 2),
#endif
                         normal_data->getPointer(1, 0),
                         normal_data->getPointer(1, 1),
#if (NDIM == 3)
                         normal_data->getPointer(1, 2),
                         normal_data->getPointer(2, 0),
                         normal_data->getPointer(2, 1),
                         normal_data->getPointer(2, 2),
#endif
                         normal_data->getGhostCellWidth().max(),
                         ls_data->getPointer(),
                         ls_data->getGhostCellWidth().max(),
                         patch_box.lower(0),
                         patch_box.upper(0),
                         patch_box.lower(1),
                         patch_box.upper(1),
#if (NDIM == 3)
                         patch_box.lower(2),
                         patch_box.upper(2),
#endif
                         dx);

            for (unsigned int axis = 0; axis < NDIM; axis++)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);
                    const double H_lower = (*H_data)(si.toCell(0));
                    const double H_upper = (*H_data)(si.toCell(1));

                    const double H = 0.5 * (H_lower + H_upper);

                    (*u_sc_data)(si) =
                        -(1.0 - H) * (*normal_data)(si, axis); // -ve sign to point the normal towards gas region;
                }
            }
        }
    }

    // copy side-centered H*n to face-centered.
    copy_side_to_face(u_adv_fc_lf_extrap_current_idx, d_u_adv_sc_lf_extrap_current_idx, d_hierarchy);

    return;
} // computeAdvectionVelocityForExtrapolation

Pointer<CellConvectiveOperator>
EnthalpyHierarchyIntegrator::getLiquidFractionExtrapConvectiveOperator(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> lf_extrap_var)
{
    // Allocate convective operator. // using H_bc for lf_var.
    std::vector<RobinBcCoefStrategy<NDIM>*> lf_bc_coef = getPhysicalBcCoefs(d_H_var);
    const ConvectiveDifferencingType lf_difference_form =
        IBAMR::string_to_enum<ConvectiveDifferencingType>("ADVECTIVE");

    // Set the convective operator.
    if (!d_lf_extrap_convective_op)
    {
        d_lf_extrap_convective_op = new AdvDiffCUIConvectiveOperator(d_object_name + "::lfExtrapConvectiveOperator",
                                                                     lf_extrap_var,
                                                                     d_lf_extrap_convective_op_input_db,
                                                                     lf_difference_form,
                                                                     lf_bc_coef);
    }

    return d_lf_extrap_convective_op;
} // getLiquidFractionExtrapConvectiveOperator

void
EnthalpyHierarchyIntegrator::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db);
#endif
    // Read in data members from input database.
    if (!is_from_restart)
    {
        d_specific_heat_liquid = input_db->getDouble("specific_heat_liquid");
        d_specific_heat_solid = input_db->getDouble("specific_heat_solid");
        d_specific_heat_gas = input_db->getDouble("specific_heat_gas");
        d_liquidus_temperature = input_db->getDouble("liquidus_temperature");
        d_solidus_temperature = input_db->getDouble("solidus_temperature");
        d_reference_temperature = input_db->getDouble("reference_temperature");

        if (input_db->keyExists("specific_heat_mushy"))
            d_specific_heat_mushy = input_db->getDouble("specific_heat_mushy");
        else
            d_specific_heat_mushy = 0.5 * (d_specific_heat_liquid + d_specific_heat_solid);
        if (input_db->keyExists("gas_liquid_fraction"))
            d_gas_liquid_fraction = input_db->getDouble("gas_liquid_fraction");
        if (input_db->keyExists("max_newton_iterations"))
            d_max_newton_iterations = input_db->getInteger("max_newton_iterations");
        if (input_db->keyExists("newton_error_tolerance"))
            d_newton_error_tolerance = input_db->getDouble("newton_error_tolerance");
        if (input_db->keyExists("require_lf_extrapolation"))
            d_require_lf_extrapolation = input_db->getBool("require_lf_extrapolation");
        if (d_require_lf_extrapolation)
        {
            d_lf_extrap_max_num_time_steps = input_db->getInteger("lf_extrap_max_num_time_steps");
            d_lf_extrap_cell_size = input_db->getDouble("lf_extrap_cell_size");
        }
    }
    return;
} // getFromInput

void
EnthalpyHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != ENTHALPY_PC_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }

    d_specific_heat_liquid = db->getDouble("specific_heat_liquid");
    d_specific_heat_solid = db->getDouble("specific_heat_solid");
    d_specific_heat_gas = db->getDouble("specific_heat_gas");
    d_liquidus_temperature = db->getDouble("liquidus_temperature");
    d_solidus_temperature = db->getDouble("solidus_temperature");
    d_reference_temperature = db->getDouble("reference_temperature");
    d_specific_heat_mushy = db->getDouble("specific_heat_mushy");
    d_gas_liquid_fraction = db->getDouble("gas_liquid_fraction");
    d_max_newton_iterations = db->getInteger("max_newton_iterations");
    d_newton_error_tolerance = db->getDouble("newton_error_tolerance");
    d_require_lf_extrapolation = db->getBool("require_lf_extrapolation");
    if (d_require_lf_extrapolation)
    {
        d_lf_extrap_max_num_time_steps = db->getInteger("lf_extrap_max_num_time_steps");
        d_lf_extrap_cell_size = db->getDouble("lf_extrap_cell_size");
    }

    return;
} // getFromRestart

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
