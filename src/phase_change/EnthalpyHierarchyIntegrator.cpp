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
#include "ibamr/AdvDiffCUIConvectiveOperator.h"
#include "ibamr/AdvDiffConservativeMassScalarTransportRKIntegrator.h"
#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/EnthalpyHierarchyIntegrator.h"
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

#include <string>
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

auto clamp = [](double x, double lower, double upper) {return (x < lower) ? lower : (x > upper ? upper : x);};

// Copy data from a side-centered variable to a face-centered variable.
void
copy_side_to_face(const int U_fc_idx, const int U_sc_idx, Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const hier::Index<NDIM>& ilower = patch->getBox().lower();
            const hier::Index<NDIM>& iupper = patch->getBox().upper();
            Pointer<SideData<NDIM, double> > U_sc_data = patch->getPatchData(U_sc_idx);
            Pointer<FaceData<NDIM, double> > U_fc_data = patch->getPatchData(U_fc_idx);
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
EnthalpyHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                           Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Perform hierarchy initialization operations common to all implementations
    // of PhaseChangeHierarchyIntegrator.
    PhaseChangeHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Register additional variables required for present time stepping algorithm.
    const IntVector<NDIM> cell_ghosts = CELLG;
    const IntVector<NDIM> no_ghosts = NOGHOSTS;

    // Register specific enthalpy.
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

        const IntVector<NDIM> side_ghosts = SIDEG;
        const IntVector<NDIM> face_ghosts = FACEG;

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
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

        d_u_adv_cc_lf_extrap_var = new CellVariable<NDIM, double>(d_object_name + "::rho_interp_cc", NDIM);
        registerVariable(d_u_adv_cc_lf_extrap_current_idx, d_u_adv_cc_lf_extrap_var, no_ghosts, getCurrentContext());

        if (d_visit_writer)
        {
            d_visit_writer->registerPlotQuantity("u_extrap_cc", "VECTOR", d_u_adv_cc_lf_extrap_current_idx, 0, 1.0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == 0)
                    d_visit_writer->registerPlotQuantity(
                        "u_extrap_cc_x", "SCALAR", d_u_adv_cc_lf_extrap_current_idx, d, 1.0);
                if (d == 1)
                    d_visit_writer->registerPlotQuantity(
                        "u_extrap_cc_y", "SCALAR", d_u_adv_cc_lf_extrap_current_idx, d, 1.0);
                if (d == 2)
                    d_visit_writer->registerPlotQuantity(
                        "u_extrap_cc_z", "SCALAR", d_u_adv_cc_lf_extrap_current_idx, d, 1.0);
            }
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_T_pre_var = new CellVariable<NDIM, double>(d_object_name + "::T_pre_var");
    d_T_pre_idx = var_db->registerVariableAndContext(d_T_pre_var, getCurrentContext());

    d_drhoh_dT_var = new CellVariable<NDIM, double>(d_object_name + "::drhoh_dT_var");
    d_drhoh_dT_scratch_idx = var_db->registerVariableAndContext(d_drhoh_dT_var, getCurrentContext(), no_ghosts);


    d_dh_dT_var = new CellVariable<NDIM, double>(d_object_name + "::dh_dT_var");
    d_dh_dT_scratch_idx = var_db->registerVariableAndContext(d_dh_dT_var, getCurrentContext(), no_ghosts);

    d_d2h_dT_var = new CellVariable<NDIM, double>(d_object_name + "::d2h_dT_var");
    d_d2h_dT_scratch_idx = var_db->registerVariableAndContext(d_d2h_dT_var, getCurrentContext(), no_ghosts);

    d_grad_T_var = new SideVariable<NDIM, double>(d_object_name + "::grad_T");
    d_grad_T_idx =
        var_db->registerVariableAndContext(d_grad_T_var, var_db->getContext(d_object_name + "grad_T::SCRATCH"));

    if (d_solve_mass_conservation)
    {
        // Set various objects with conservative time integrator.
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;
        rho_p_cc_integrator->setCellCenteredTransportQuantityBoundaryConditions(d_h_bc_coef);
    }

    // For extrapolating the lf_var from PCM into the gas phase.
    if (d_lf_extrap_var) d_lf_extrap_convective_op = getLiquidFractionExtrapConvectiveOperator(d_lf_extrap_var);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
EnthalpyHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                          const double new_time,
                                                          const int num_cycles)
{
    PhaseChangeHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_drhoh_dT_scratch_idx)) level->allocatePatchData(d_drhoh_dT_scratch_idx, current_time);
        if (!level->checkAllocated(d_dh_dT_scratch_idx)) level->allocatePatchData(d_dh_dT_scratch_idx, current_time);
        if (!level->checkAllocated(d_d2h_dT_scratch_idx)) level->allocatePatchData(d_d2h_dT_scratch_idx, current_time);
        if (!level->checkAllocated(d_T_pre_idx)) level->allocatePatchData(d_T_pre_idx, current_time);
        if (!level->checkAllocated(d_grad_T_idx)) level->allocatePatchData(d_grad_T_idx, current_time);
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

    const int H_current_idx = var_db->mapVariableAndContextToIndex(d_H_var, getCurrentContext());


    if (initial_time)
    computeEnthalpyBasedOnTemperature(
        d_h_current_idx, d_T_current_idx, d_rho_current_idx, d_lf_current_idx, H_current_idx);

    if (d_solve_mass_conservation)
    {
        Pointer<AdvDiffConservativeMassScalarTransportRKIntegrator> rho_p_cc_integrator = d_rho_p_integrator;

        // Set the velocities used to update the density and the previous time step
        // size
        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ -1, /*current*/ d_u_adv_current_idx, /*new*/ -1);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*current*/ d_h_current_idx, /*new*/ -1);
        }
        else
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ d_U_old_current_idx, /*current*/ d_u_adv_current_idx, /*new*/ -1);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*current*/ d_h_current_idx, /*new*/ -1);
            d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        // Integrate density and convective term of energy equation.
        d_rho_p_integrator->integrate(dt);
    }
    // Setup the operators and solvers and compute the right-hand-side terms.
    // Setup the problem coefficients for the linear solver.
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

    const double apply_time = current_time;
    for (unsigned k = 0; k < d_reset_kappa_fcns.size(); ++k)
    {
        d_reset_kappa_fcns[k](d_T_diffusion_coef_cc_current_idx,
                              d_T_diffusion_coef_cc_var,
                              d_hier_math_ops,
                              -1 /*cycle_num*/,
                              apply_time,
                              current_time,
                              new_time,
                              d_reset_kappa_fcns_ctx[k]);
    }

    // Interpolate the cell-centered diffusion coef to side-centered.
    d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_scratch_idx, d_T_diffusion_coef_cc_current_idx);
    d_k_bdry_bc_fill_op->fillData(current_time);

    if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
    {
        interpolateCCToSCSimpleAveraging(d_T_diffusion_coef_current_idx, d_T_diffusion_coef_cc_scratch_idx);
    }
    else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
    {
        interpolateCCToSCHarmonicAveraging(d_T_diffusion_coef_current_idx, d_T_diffusion_coef_cc_scratch_idx);
    }
    else
    {
        TBOX_ERROR("this statement should not be reached");
    }

    // For plotting purpose.
    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(d_D_cc_new_idx,
                            d_D_cc_var,
                            d_T_diffusion_coef_current_idx,
                            d_T_diffusion_coef_var,
                            d_no_fill_op,
                            d_integrator_time,
                            synch_cf_interface);

    d_hier_sc_data_ops->scale(d_T_diffusion_coef_rhs_scratch_idx, (1.0 - alpha), d_T_diffusion_coef_current_idx);
    T_rhs_op_spec.setDPatchDataId(d_T_diffusion_coef_rhs_scratch_idx);

    // Initialize the RHS operator and compute the RHS vector for the temperature equation.
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
    d_hier_cc_data_ops->copyData(d_T_scratch_idx, d_T_current_idx, false);
    T_rhs_op->apply(*d_T_sol, *d_T_rhs);

    d_hier_cc_data_ops->copyData(d_h_new_idx, d_h_current_idx);
    d_hier_cc_data_ops->copyData(d_T_new_idx, d_T_current_idx);
    d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_new_idx, d_T_diffusion_coef_cc_current_idx);
    d_hier_cc_data_ops->copyData(d_lf_new_idx, d_lf_current_idx);

    // Initializing operator state.
    if (d_lf_extrap_var) d_lf_extrap_convective_op->initializeOperatorState(*d_lf_extrap_sol, *d_lf_extrap_rhs);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
EnthalpyHierarchyIntegrator::integrateHierarchySpecialized(const double current_time,
                                                           const double new_time,
                                                           const int cycle_num)
{
    // First advance all quantities managed by the base advection-diffusion
    // integrator. In particular, this produces ls^{n+1} in the new context.
    AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchySpecialized(current_time, new_time, cycle_num);

    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

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

    // Perform a single step of fixed point iteration.

    
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());

    // computeEnthalpyBasedOnTemperature(
        // d_h_new_idx, d_T_new_idx, d_rho_new_idx, d_lf_new_idx, H_new_idx);

        // computeTemperatureBasedOnEnthalpy(d_T_new_idx, d_h_new_idx, H_new_idx);

        // Find lf^n+1, m+1 based on h^n+1, m+1.
        // computeLiquidFraction(d_lf_new_idx, d_h_new_idx, H_new_idx);

    // In the special case of conservative discretization, the updated
    // density is calculated by the mass and convective integrator.
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
        d_rho_p_integrator->setConvectiveDerivativePatchDataIndex(d_T_N_scratch_idx);

        // Always set to current because we want to update rho^{n} to rho^{n+1}
        d_rho_p_integrator->setDensityPatchDataIndex(d_rho_current_idx);

        // Set the velocities used to update the density
        if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ -1, /*current*/ d_u_adv_current_idx, /*new*/ d_u_adv_new_idx);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*current*/ d_h_current_idx, /*new*/ d_h_new_idx);
        }
        else
        {
            d_rho_p_integrator->setFluidVelocityPatchDataIndices(
                /*old*/ d_U_old_current_idx,
                /*current*/ d_u_adv_current_idx,
                /*new*/ d_u_adv_new_idx);
            rho_p_cc_integrator->setTransportQuantityPatchDataIndices(
                /*current*/ d_h_current_idx,
                /*new*/ d_h_new_idx);

            d_rho_p_integrator->setPreviousTimeStepSize(d_dt_previous[0]);
        }

        d_rho_p_integrator->integrate(dt);
    }

    d_updated_rho_idx = d_rho_p_integrator ? d_rho_p_integrator->getUpdatedDensityPatchDataIndex() : d_rho_new_idx; // this is the sts mass flux integrator density; 
    d_hier_cc_data_ops->copyData(d_rho_new_idx,
                                 d_updated_rho_idx,
                                 /*interior_only*/ true);

    // Heaviside is H^n+1, calculate temperautre and liquid fraction for this heaviside. 
    // computeTemperatureBasedOnEnthalpy(d_T_new_idx, d_h_new_idx, H_new_idx); // T(h^n, H^n+1); 
    // computeLiquidFraction(d_lf_new_idx, d_h_new_idx, H_new_idx); // lf(h^n, H^n+1); 

    // similarly need to update thermal conductivity 
const double StateApplyTime = new_time; 
    for (unsigned k =0; k< d_reset_kappa_fcns.size(); ++k)
    {
        d_reset_kappa_fcns[k](
            d_T_diffusion_coef_cc_new_idx, 
            d_T_diffusion_coef_cc_var, 
            d_hier_math_ops, 
            -1, 
            StateApplyTime, 
            current_time, 
            new_time, 
            d_reset_kappa_fcns_ctx[k]
        );
    }

    // Account for the convective acceleration term N_full.
    if (d_solve_mass_conservation && d_u_adv_var)
        d_hier_cc_data_ops->axpy(d_T_rhs_scratch_idx, -1.0, d_T_N_scratch_idx, d_T_rhs_scratch_idx);

    PoissonSpecifications T_solver_spec(d_object_name + "::solver_spec::" + d_T_var->getName());

    // compute u_adv = H*n based on phi^{n+1, k+1}.
    if (d_lf_extrap_var)
        computeAdvectionVelocityForExtrapolation(d_u_adv_fc_lf_extrap_current_idx, current_time, new_time);

    double lf_relative_iteration_error = 1.0;
    double inner_iterations = 1.0;

    // I can create rho^n h^n here. 
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex(); 
    
    Pointer<SAMRAIVectorReal<NDIM, double>> h_current_vec = 
        new SAMRAIVectorReal<NDIM, double > (
            d_object_name + "::h_current_wrapper",
            d_hierarchy, 
            coarsest_ln, 
            finest_ln
        );

    h_current_vec->addComponent(
        d_h_var, 
        d_h_current_idx, 
        wgt_idx, 
        d_hier_cc_data_ops
    ); 

    Pointer<SAMRAIVectorReal<NDIM, double>> rho_h_n_vec = h_current_vec->cloneVector(d_object_name + "::rho_h_n");

    rho_h_n_vec->allocateVectorData(current_time); 
    const int rho_h_n_idx = rho_h_n_vec->getComponentDescriptorIndex(0);

    d_hier_cc_data_ops->multiply(rho_h_n_idx,d_rho_current_idx,d_h_current_idx); // now rho h is present before the newton iteration loop starts
    double T_relative_iteration_error = 1.0; 
    // Inner iterations for the Newton-Ralphson scheme.
    while ((lf_relative_iteration_error >= d_lf_iteration_error_tolerance || T_relative_iteration_error >= d_lf_iteration_error_tolerance )&& inner_iterations <= d_max_inner_iterations)
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

        // Interpolate the cell-centered diffusion coef to side-centered.
        d_hier_cc_data_ops->copyData(d_T_diffusion_coef_cc_scratch_idx, d_T_diffusion_coef_cc_new_idx);
        d_k_bdry_bc_fill_op->fillData(new_time);

        if (d_k_vc_interp_type == VC_AVERAGE_INTERP)
        {
            interpolateCCToSCSimpleAveraging(d_T_diffusion_coef_new_idx, d_T_diffusion_coef_cc_scratch_idx);
        }
        else if (d_k_vc_interp_type == VC_HARMONIC_INTERP)
        {
            interpolateCCToSCHarmonicAveraging(d_T_diffusion_coef_new_idx, d_T_diffusion_coef_cc_scratch_idx);
        }
        else
        {
            TBOX_ERROR("this statement should not be reached");
        }

        // For plotting purpose.
        static const bool synch_cf_interface = true;
        d_hier_math_ops->interp(d_D_cc_new_idx,
                                d_D_cc_var,
                                d_T_diffusion_coef_new_idx,
                                d_T_diffusion_coef_var,
                                d_no_fill_op,
                                d_integrator_time,
                                synch_cf_interface);

        d_hier_sc_data_ops->scale(d_T_diffusion_coef_scratch_idx, -alpha, d_T_diffusion_coef_new_idx);
        T_solver_spec.setDPatchDataId(d_T_diffusion_coef_scratch_idx);

        computeRhoEnthalpyDerivative(d_drhoh_dT_scratch_idx, d_T_new_idx, H_new_idx);

        computeEnthalpyDerivative(d_dh_dT_scratch_idx, d_T_new_idx, H_new_idx);

        computeEnthalpySecondDerivative(d_d2h_dT_scratch_idx, d_T_new_idx, H_new_idx);

        // Set rho*Cp/dt.
        // d_hier_cc_data_ops->multiply(d_C_new_idx, d_rho_new_idx, d_dh_dT_scratch_idx); // THis does not need to happen now
        // d_hier_cc_data_ops->copyData(d_T_C_idx, d_drhoh_dT_scratch_idx);
        d_hier_cc_data_ops->multiply(d_C_new_idx, d_rho_new_idx, d_dh_dT_scratch_idx);
        d_hier_cc_data_ops->copyData(d_T_C_idx, d_C_new_idx);

        d_hier_cc_data_ops->scale(d_T_C_idx, 1.0 / dt, d_T_C_idx);
        T_solver_spec.setCPatchDataId(d_T_C_idx);

        // Initialize the linear solver for temperature equation.
        Pointer<PoissonSolver> T_solver = d_T_solver;
        T_solver->setPoissonSpecifications(T_solver_spec);
        T_solver->setPhysicalBcCoef(d_T_bc_coef);
        T_solver->setHomogeneousBc(false);
        T_solver->setSolutionTime(new_time);
        T_solver->setTimeInterval(current_time, new_time);
        // Initialize solver each time because the coefficients are changing.
        if (d_enable_logging)
        {
            plog << d_object_name << ": "
                 << "Initializing the solvers for" << d_T_var->getName() << "\n";
        }
        T_solver->initializeSolverState(*d_T_sol, *d_T_rhs);

        // Account for forcing terms.
        if (d_T_F_fcn)
        {
            d_T_F_fcn->setDataOnPatchHierarchy(d_T_F_scratch_idx, d_T_F_var, d_hierarchy, half_time);
        }
        else
            d_hier_cc_data_ops->setToScalar(d_T_F_scratch_idx, 0.0);

        // Compute and add temporal and linearized terms to the RHS of the energy equation.
        addTemporalAndLinearTermstoRHSOfEnergyEquation(d_T_F_scratch_idx, rho_h_n_idx, H_new_idx, dt);
        d_hier_cc_data_ops->axpy(d_T_rhs_scratch_idx, +1.0, d_T_F_scratch_idx, d_T_rhs_scratch_idx);

        // Storing T^n+1,m.
        d_hier_cc_data_ops->copyData(d_T_pre_idx, d_T_new_idx);
        // d_hier_cc_data_ops->copyData(d_lf_pre_idx, d_lf_new_idx); 

        // Solve for T(n+1, m+1).
        T_solver->solveSystem(*d_T_sol, *d_T_rhs);
        d_hier_cc_data_ops->copyData(d_T_new_idx, d_T_scratch_idx);

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
        updateEnthalpy(d_h_new_idx, d_T_new_idx, d_T_pre_idx);

        // computeEnthalpyBasedOnTemperatureOnly(d_T_new_idx,d_h_new_idx , H_new_idx); 

        // Find T^n+1,m+1 based on h^n+1, m+1.
        computeTemperatureBasedOnEnthalpy(d_T_new_idx, d_h_new_idx, H_new_idx);

        // Find lf^n+1, m+1 based on h^n+1, m+1.
        d_hier_cc_data_ops->copyData(d_lf_pre_idx, d_lf_new_idx);
        computeLiquidFraction(d_lf_new_idx, d_h_new_idx, H_new_idx);

        // Update specific heat
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

        // Update conductivity
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

        // Finding L2 norm of lf^m-1 iteration.
        const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
        lf_relative_iteration_error = d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx);

        const double T_old_norm = d_hier_cc_data_ops->L2Norm(d_T_pre_idx, wgt_cc_idx); 
        d_hier_cc_data_ops->subtract(d_T_pre_idx, d_T_new_idx, d_T_pre_idx); 
        T_relative_iteration_error = d_hier_cc_data_ops->L2Norm(d_T_pre_idx, wgt_cc_idx)/(1.0+T_old_norm); 



        // Finding lf^m - lf^m-1.
        d_hier_cc_data_ops->subtract(d_lf_pre_idx, d_lf_new_idx, d_lf_pre_idx);

        plog << "liquid fraction relative error norms at Newton iteration: " << inner_iterations << "\n"
             << "L1 : " << d_hier_cc_data_ops->L1Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error)
             << " || "
             << "L2 : " << d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error)
             << " || "
             << "L_oo : " << d_hier_cc_data_ops->maxNorm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error)
             << " || "
             << "L_2T : " << T_relative_iteration_error
             << "\n";

        lf_relative_iteration_error =
            d_hier_cc_data_ops->L2Norm(d_lf_pre_idx, wgt_cc_idx) / (1.0 + lf_relative_iteration_error);
        inner_iterations++;

        d_hier_cc_data_ops->axpy(d_T_rhs_scratch_idx, -1.0, d_T_F_scratch_idx, d_T_rhs_scratch_idx);
        d_hier_cc_data_ops->copyData(d_T_F_new_idx, d_T_F_scratch_idx);
    }

    rho_h_n_vec->freeVectorComponents(); 
    rho_h_n_vec.setNull(); 
    h_current_vec.setNull(); 


    // Reset the right-hand side vector.
    if (d_solve_mass_conservation && d_u_adv_var)
        d_hier_cc_data_ops->axpy(d_T_rhs_scratch_idx, +1.0, d_T_N_scratch_idx, d_T_rhs_scratch_idx);

    // Compute the source term for Div U equation.
    computeDivergenceVelocitySourceTerm(d_Div_U_F_idx, new_time);

    // extrapolate the liquid fraction to gas region.
    if (d_lf_extrap_var) extrapolateLiquidFractionToGasRegion(d_lf_new_idx);

    return;
} // integrateHierarchy

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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_T_pre_idx);
        level->deallocatePatchData(d_drhoh_dT_scratch_idx);
        level->deallocatePatchData(d_dh_dT_scratch_idx);
        level->deallocatePatchData(d_d2h_dT_scratch_idx);
        level->deallocatePatchData(d_grad_T_idx);
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
EnthalpyHierarchyIntegrator::registerSpecificEnthalpyVariable(Pointer<CellVariable<NDIM, double> > h_var,
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
EnthalpyHierarchyIntegrator::registerUpdateVOFFromLevelSetFcn(UpdateVOFFromLevelSetFcnPtr callback, void* ctx)
{
    d_update_vof_from_ls_fcn = callback;
    d_update_vof_from_ls_ctx = ctx;
    return; 
}

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
    db->putInteger("max_inner_iterations", d_max_inner_iterations);
    db->putDouble("lf_iteration_error_tolerance", d_lf_iteration_error_tolerance);
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
EnthalpyHierarchyIntegrator::registerLevelSetVariable(Pointer<CellVariable<NDIM, double> > phi_var)
{
    d_phi_var = phi_var;
    return;
} // registerLevelSetVariable

void
EnthalpyHierarchyIntegrator::registerLiquidFractionVariableForExtrapolation(
    Pointer<CellVariable<NDIM, double> > lf_extrap_var)
{
    d_lf_extrap_var = lf_extrap_var;
    return;
} // registerLiquidFractionVariableForExtrapolation

/////////////////////////////// PROTECTED //////////////////////////////////////

void
EnthalpyHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    PhaseChangeHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
        base_hierarchy, coarsest_level, finest_level);

    if (d_lf_extrap_var)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        // Setup the patch boundary filling objects.
        const int finest_hier_level = d_hierarchy->getFinestLevelNumber();

        // Reset the solution and rhs vectors.
        const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

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
EnthalpyHierarchyIntegrator::addTemporalAndLinearTermstoRHSOfEnergyEquation(int F_scratch_idx, int rho_h_n_idx, int H_new_idx, const double dt)
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
            Pointer<CellData<NDIM, double> > T_new_data = patch->getPatchData(d_T_new_idx);
            Pointer<CellData<NDIM, double> > h_new_data = patch->getPatchData(d_h_new_idx);
            Pointer<CellData<NDIM, double> > H_new_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double> > h_current_data = patch->getPatchData(d_h_current_idx);
            Pointer<CellData<NDIM, double> > rho_new_data = patch->getPatchData(d_rho_new_idx);
            Pointer<CellData<NDIM, double> > rho_current_data = patch->getPatchData(d_rho_current_idx);
            Pointer<CellData<NDIM, double> > dh_dT_data = patch->getPatchData(d_dh_dT_scratch_idx);
            Pointer<CellData<NDIM, double> > d2h_dT_data = patch->getPatchData(d_d2h_dT_scratch_idx);
            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_scratch_idx);
            Pointer<CellData<NDIM, double> > rho_h_n_data = patch->getPatchData(rho_h_n_idx);
            Pointer<CellData<NDIM, double> > drhoh_dT_data = patch->getPatchData(d_drhoh_dT_scratch_idx);
            



            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double T = (*T_new_data)(ci); 
                const double H = (*H_new_data)(ci); 

                const double rho = (*rho_new_data)(ci);
const double dhdT = (*dh_dT_data)(ci);
const double d2hdT2 = (*d2h_dT_data)(ci);

                double rho_thermo = 0.0; 
                double rho_h_m = 0.0;
                double rhoT; 
                
                
                const double delta_T =
                d_liquidus_temperature -
                d_solidus_temperature;
            
            
            
            if (T < d_solidus_temperature)
            {
                rho_thermo =(1.0 - H) *d_rho_gas +H *d_rho_solid; 
                rhoT = 0.0; 
                   
            }
            else if (T > d_liquidus_temperature)
            {
                rho_thermo =(1.0 - H) *d_rho_gas +H *d_rho_liquid; 
                rhoT = 0.0; 
            
            }
            else
            {
                const double theta =(T - d_solidus_temperature) /delta_T;
                const double R0 = (1.0 - H)*d_rho_gas + H*d_rho_solid; 
                const double R1 = H* ( d_rho_liquid - d_rho_solid); 
                rho_thermo = R0 + R1*theta; 
                rhoT = H*(d_rho_liquid - d_rho_solid)/delta_T; 

            }   

            // Actually I shoudl only construct rho^n+1 not rho h here. h new data is already present. 
                // (*F_data)(ci) +=  -1.0 / dt *
                //                     (rho_thermo*(*h_new_data)(ci) -  (*drhoh_dT_data)(ci)*(*T_new_data)(ci) - (*rho_h_n_data)(ci));
                
                
                // (*F_data)(ci) += -1.0 / dt *
                //                  ((*rho_new_data)(ci) * (*h_new_data)(ci) - (*drhoh_dT_data)(ci) * (*T_new_data)(ci) -
                //                   (*rho_h_n_data)(ci));//using rho breve in new method

                (*F_data)(ci) += -1.0 / dt *
                                 ((*rho_new_data)(ci) * (*h_new_data)(ci) - (*rho_new_data)(ci) *(*dh_dT_data)(ci) * (*T_new_data)(ci) -
                                  (*rho_h_n_data)(ci)); //original method using dh/dT


                // (*F_data)(ci) += -1.0 / dt *
                //                  ((*rho_new_data)(ci) * (*h_new_data)(ci) - (*rho_new_data)(ci) *(*dh_dT_data)(ci) * (*T_new_data)(ci) 
                //                         - (*h_new_data)(ci) * rhoT * (*T_new_data)(ci) - 
                //                   (*rho_h_n_data)(ci)); //original method using dh/dT multiplying with rho breve with splitting derivative. 


            }
        }
    }
    return;
} // addTemporalAndLinearTermstoRHSOfEnergyEquation


void
EnthalpyHierarchyIntegrator::computeDivergenceVelocitySourceTerm(int Div_U_F_idx, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());

    // Filling ghost cells for temperature.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> T_transaction_comps(1);
    T_transaction_comps[0] = InterpolationTransactionComponent(d_T_scratch_idx,
                                                               d_T_new_idx,
                                                               "CONSERVATIVE_LINEAR_REFINE",
                                                               false,
                                                               "CONSERVATIVE_COARSEN",
                                                               "LINEAR",
                                                               false,
                                                               d_T_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(T_transaction_comps, d_hierarchy);
    hier_bdry_fill->fillData(new_time);

    // Compute gradient of temperature.
    d_hier_math_ops->grad(d_grad_T_idx, d_grad_T_var, true, 1.0, d_T_scratch_idx, d_T_var, nullptr, new_time);

    // Compute k*grad_T.
    d_hier_sc_data_ops->multiply(d_grad_T_idx, d_grad_T_idx, d_T_diffusion_coef_new_idx);

    // Compute div(k*grad_T).
    d_hier_math_ops->div(Div_U_F_idx, d_Div_U_F_var, 1.0, d_grad_T_idx, d_grad_T_var, nullptr, new_time, false);

    const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;
    const double delta_T_m = d_liquidus_temperature - d_solidus_temperature;
    const double delta_rho = d_rho_liquid - d_rho_solid;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(d_h_new_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(d_rho_new_idx);
            Pointer<CellData<NDIM, double> > Div_U_F_data = patch->getPatchData(Div_U_F_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(d_T_new_idx);



            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                double hs_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_solid * h_s)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);
                 double hl_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_liquidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_liquid * h_l)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                 double R0 = (1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid; 
                 double R1 = (*H_data)(ci)*(d_rho_liquid - d_rho_solid); 

                 double H0 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) + (*H_data)(ci)*d_rho_solid * h_s; 
                 double H1 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas*(d_liquidus_temperature - d_solidus_temperature) + 
                         (*H_data)(ci)*((d_rho_liquid - d_rho_solid)*h_s + d_rho_liquid*( d_latent_heat + d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature)));

                // double Rho_T = (((*H_data)(ci) - 1.0)*d_solidus_temperature*d_rho_gas - (*H_data)(ci)* d_solidus_temperature*d_rho_liquid + (*H_data)(ci)*(*T_data)(ci)*(d_rho_liquid - d_rho_solid) + d_liquidus_temperature*(d_rho_gas - (*H_data)(ci)*d_rho_gas + (*H_data)(ci)*d_rho_solid))/ 
                // (d_liquidus_temperature - d_solidus_temperature); 

                double C_LG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_liquid*d_specific_heat_liquid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                double C_SG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_solid*d_specific_heat_solid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);




                double material_derivative = 0.0;
        //         if ((*h_data)(ci) >= hs_cell && (*h_data)(ci) <= hl_cell )
        //         {

        //             // material_derivative
        //             material_derivative = (*Div_U_F_data)(ci)*(-d_rho_solid + d_rho_liquid)/(H1*R0 - H0*R1);
        //         //  material_derivative = (*Div_U_F_data)(ci)*(-d_rho_solid + d_rho_liquid)/(d_rho_solid*d_rho_liquid*(h_l-h_s));
		// }

            // (*Div_U_F_data)(ci) =-material_derivative*(*H_data)(ci);


        const double Qmix = (*Div_U_F_data)(ci);

double H = clamp((*H_data)(ci), 0.0, 1.0);
const double alpha_p = H;
const double alpha_g = 1.0 - H;

const double T = (*T_data)(ci);
const double DelT = d_liquidus_temperature - d_solidus_temperature;

if (T >= d_solidus_temperature && T <= d_liquidus_temperature)
{
    double phi = (T - d_solidus_temperature) / DelT;
    phi = clamp(phi, 0.0, 1.0);

    const double rho_p =
        d_rho_solid*(1.0 - phi) + d_rho_liquid*phi;

    const double Cv_pm =
        d_rho_liquid*d_rho_solid*
        (d_specific_heat_mushy + d_latent_heat/DelT) / rho_p;

    const double Cv_mix =
        alpha_g*d_rho_gas*d_specific_heat_gas
      + alpha_p*Cv_pm;

    (*Div_U_F_data)(ci) =
        -alpha_p*(d_rho_liquid - d_rho_solid) /
        (rho_p*DelT*Cv_mix) * Qmix;
}
else
{
    (*Div_U_F_data)(ci) = 0.0;
}


	    }
        }
    }

    return;
} // computeDivergenceVelocitySourceTerm

void
EnthalpyHierarchyIntegrator::computeEnthalpyBasedOnTemperature(int h_idx,
                                                               const int T_idx,
                                                               const int rho_idx,
                                                               const int lf_idx,
                                                               const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;
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

                double hs_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_solid * h_s)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);
                 double hl_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_liquidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_liquid * h_l)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                 double R0 = (1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid; 
                 double R1 = (*H_data)(ci)*(d_rho_liquid - d_rho_solid); 

                 double H0 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) + (*H_data)(ci)*d_rho_solid * h_s; 
                 double H1 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas*(d_liquidus_temperature - d_solidus_temperature) + 
                         (*H_data)(ci)*((d_rho_liquid - d_rho_solid)*h_s + d_rho_liquid*( d_latent_heat + d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature)));

                double Rho_T = (((*H_data)(ci) - 1.0)*d_solidus_temperature*d_rho_gas - (*H_data)(ci)* d_solidus_temperature*d_rho_liquid + (*H_data)(ci)*(*T_data)(ci)*(d_rho_liquid - d_rho_solid) + d_liquidus_temperature*(d_rho_gas - (*H_data)(ci)*d_rho_gas + (*H_data)(ci)*d_rho_solid))/ 
                (d_liquidus_temperature - d_solidus_temperature); 

                double C_LG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_liquid*d_specific_heat_liquid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                double C_SG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_solid*d_specific_heat_solid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);

                // if ((*H_data)(ci) >= H_LIM)
                // {
                    if ((*T_data)(ci) < d_solidus_temperature)
                    {
                        (*h_data)(ci) =  ((*T_data)(ci)-d_reference_temperature)*C_SG;
                    }
                    else if ((*T_data)(ci) > d_liquidus_temperature)
                    {
                    (*h_data)(ci) = hl_cell + C_LG*((*T_data)(ci)-d_liquidus_temperature);
                    
                    // (d_rho_liquid/(*rho_data)(ci))*d_specific_heat_mushy * ((*T_data)(ci)-d_solidus_temperature) + h_s +
			        //                                     (*lf_data)(ci)*d_rho_liquid * d_latent_heat / (*rho_data)(ci);
		            }
                    else
                    {
                        (*h_data)(ci) = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas*((*T_data)(ci)-d_reference_temperature) + 
                        ((*H_data)(ci) - (*lf_data)(ci))*d_rho_solid*h_s + 
                        (*lf_data)(ci)*d_rho_liquid*(h_s + d_latent_heat + d_specific_heat_mushy*(d_liquidus_temperature - d_solidus_temperature)))/( (*rho_data)(ci));
                    }
                // }
                // else
                // {
                //     (*h_data)(ci) = d_specific_heat_gas * ((*T_data)(ci)-d_reference_temperature);
                // }
            }
        }
    }
    return;
} // computeEnthalpyBasedOnTemperature



void
EnthalpyHierarchyIntegrator::computeEnthalpyBasedOnTemperatureOnly(const int T_idx,
                                                                int h_idx,
                                                               const int H_idx)
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
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(h_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double H = (*H_data)(ci); 
                const double T = (*T_data)(ci); 

                const double theta =(T - d_solidus_temperature) / (d_liquidus_temperature -d_solidus_temperature);
            
                const double R0 =(1.0 - H) * d_rho_gas + H * d_rho_solid;
            
                const double R1 =H * (d_rho_liquid - d_rho_solid);
            
                // const double rho =R0 + R1 * theta;

                const double delta_T = d_liquidus_temperature - d_solidus_temperature; 

                const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
                const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;

                 
                //  double H0 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) + (*H_data)(ci)*d_rho_solid * h_s; 
                //  double H1 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas*(d_liquidus_temperature - d_solidus_temperature) + 
                        //  (*H_data)(ci)*((d_rho_liquid - d_rho_solid)*h_s + d_rho_liquid*( d_latent_heat + d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature)));

                double Rho_T = (((*H_data)(ci) - 1.0)*d_solidus_temperature*d_rho_gas - (*H_data)(ci)* d_solidus_temperature*d_rho_liquid + (*H_data)(ci)*(*T_data)(ci)*(d_rho_liquid - d_rho_solid) + d_liquidus_temperature*(d_rho_gas - (*H_data)(ci)*d_rho_gas + (*H_data)(ci)*d_rho_solid))/ 
                (d_liquidus_temperature - d_solidus_temperature); 

                double C_LG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_liquid*d_specific_heat_liquid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                double C_SG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_solid*d_specific_heat_solid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);


                double rho; 
                double rho_h; 
                if (T < d_solidus_temperature)
                {
                    rho =
                        (1.0 - H) * d_rho_gas +
                        H * d_rho_solid;
                
                    rho_h =
                        (1.0 - H) *
                            d_rho_gas *
                            d_specific_heat_gas *
                            (T - d_reference_temperature)
                        +
                        H *
                            d_rho_solid *
                            d_specific_heat_solid *
                            (T - d_reference_temperature);
                }
                else if (T > d_liquidus_temperature)
                {
                    rho =
                        (1.0 - H) * d_rho_gas +
                        H * d_rho_liquid;
                
                    // rho*h evaluated exactly at T_liquidus.
                    const double rho_h_liquidus =
                        (1.0 - H) *
                            d_rho_gas *
                            d_specific_heat_gas *
                            (d_liquidus_temperature -
                             d_reference_temperature)
                        +
                        H *
                            d_rho_liquid *
                            h_l;
                
                    // volumetric heat capacity of the liquid/gas mixture.
                    const double rho_cp_LG =
                        (1.0 - H) *
                            d_rho_gas *
                            d_specific_heat_gas
                        +
                        H *
                            d_rho_liquid *
                            d_specific_heat_liquid;
                
                    rho_h =
                        rho_h_liquidus +
                        rho_cp_LG *
                            (T - d_liquidus_temperature);
                }
                else
                {
                    const double theta =
                        (T - d_solidus_temperature) /
                        delta_T;
                
                    // rho = R0 + R1*theta
                    const double R0 =
                        (1.0 - H) * d_rho_gas +
                        H * d_rho_solid;
                
                    const double R1 =
                        H * (d_rho_liquid - d_rho_solid);
                
                    rho =
                        R0 + R1 * theta;
                
                
                    // rho*h = H0 + H1*theta
                    const double H0 =
                        (1.0 - H) *
                            d_rho_gas *
                            d_specific_heat_gas *
                            (d_solidus_temperature -
                             d_reference_temperature)
                        +
                        H *
                            d_rho_solid *
                            h_s;
                
                    const double H1 =
                        (1.0 - H) *
                            d_rho_gas *
                            d_specific_heat_gas *
                            delta_T
                        +
                        H *
                            (
                                (d_rho_liquid - d_rho_solid) * h_s
                                +
                                d_rho_liquid *
                                    (
                                        d_latent_heat +
                                        d_specific_heat_mushy * delta_T
                                    )
                            );
                
                    rho_h =
                        H0 + H1 * theta;
                }

            
                (*h_data)(ci) = rho_h / rho;
                // // double rho_h = H0 + H1*theta;

                // // if ((*H_data)(ci) >= H_LIM)
                // // {
                //     if ((*T_data)(ci) < d_solidus_temperature)
                //     {
                //         (*h_data)(ci) =  ((*T_data)(ci)-d_reference_temperature)*C_SG;
                //     }
                //     else if ((*T_data)(ci) > d_liquidus_temperature)
                //     {
                //     (*h_data)(ci) = hl_cell + C_LG*((*T_data)(ci)-d_liquidus_temperature);
                    
                //     // (d_rho_liquid/(*rho_data)(ci))*d_specific_heat_mushy * ((*T_data)(ci)-d_solidus_temperature) + h_s +
			    //     //                                     (*lf_data)(ci)*d_rho_liquid * d_latent_heat / (*rho_data)(ci);
		        //     }
                //     else
                //     {
                //         (*h_data)(ci) = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas*((*T_data)(ci)-d_reference_temperature) + 
                //         ((*H_data)(ci) - (*lf_data)(ci))*d_rho_solid*h_s + 
                //         (*lf_data)(ci)*d_rho_liquid*(h_s + d_latent_heat + d_specific_heat_mushy*(d_liquidus_temperature - d_solidus_temperature)))/( (*rho_data)(ci));
                //     }
                // // }
                // // else
                // // {
                // //     (*h_data)(ci) = d_specific_heat_gas * ((*T_data)(ci)-d_reference_temperature);
                // // }
            }
        }
    }
    return;
} // computeEnthalpyBasedOnTemperature


void
EnthalpyHierarchyIntegrator::computeTemperatureBasedOnEnthalpy(int T_idx, const int h_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;
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

                 double hs_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_solid * h_s)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);
                 double hl_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_liquidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_liquid * h_l)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                 double R0 = (1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid; 
                 double R1 = (*H_data)(ci)*(d_rho_liquid - d_rho_solid); 

                 double H0 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) + (*H_data)(ci)*d_rho_solid * h_s; 
                 double H1 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas*(d_liquidus_temperature - d_solidus_temperature) + 
                         (*H_data)(ci)*((d_rho_liquid - d_rho_solid)*h_s + d_rho_liquid*( d_latent_heat + d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature)));

                double Rho_T = (((*H_data)(ci) - 1.0)*d_solidus_temperature*d_rho_gas - (*H_data)(ci)* d_solidus_temperature*d_rho_liquid + (*H_data)(ci)*(*T_data)(ci)*(d_rho_liquid - d_rho_solid) + d_liquidus_temperature*(d_rho_gas - (*H_data)(ci)*d_rho_gas + (*H_data)(ci)*d_rho_solid))/ 
                                (d_liquidus_temperature - d_solidus_temperature); 

                double C_LG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_liquid*d_specific_heat_liquid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                double C_SG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_solid*d_specific_heat_solid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);


                // if ((*H_data)(ci) > 0.0)
                // {
                    if ((*h_data)(ci) < hs_cell)
                    {
                        (*T_data)(ci) = d_reference_temperature + ((*h_data)(ci))/C_SG;
                    }
                    else if ((*h_data)(ci) > hl_cell)
                    {
                        (*T_data)(ci) = d_liquidus_temperature + ((*h_data)(ci)-hl_cell) / C_LG;
		            }
                    else
                    {
                        (*T_data)(ci) = d_solidus_temperature + (d_liquidus_temperature - d_solidus_temperature)*(((*h_data)(ci)*R0 - H0)/(H1 - (*h_data)(ci)*R1));
                    }
                // }
                // // else
                // // {
                // //     (*T_data)(ci) = (*h_data)(ci) / d_specific_heat_gas + d_reference_temperature;
                // // }
            }
        }
    }
    return;
} // computeTemperatureBasedOnEnthalpy

void
EnthalpyHierarchyIntegrator::updateEnthalpy(int h_new_idx, const int T_new_idx, const int T_pre_idx)
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
            Pointer<CellData<NDIM, double> > h_new_data = patch->getPatchData(h_new_idx);
            Pointer<CellData<NDIM, double> > T_new_data = patch->getPatchData(T_new_idx);
            Pointer<CellData<NDIM, double> > T_pre_data = patch->getPatchData(T_pre_idx);
            Pointer<CellData<NDIM, double> > dh_dT_data = patch->getPatchData(d_dh_dT_scratch_idx);
            Pointer<CellData<NDIM, double> > drhoh_dT_data = patch->getPatchData(d_drhoh_dT_scratch_idx);
            Pointer<CellData<NDIM, double> > d2h_dT_data = patch->getPatchData(d_d2h_dT_scratch_idx);
            // Pointer<CellData<NDIM, double> > rho_h_n_data = patch->getPatchData(rho_h_n_idx);
            // Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);


            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                
                 (*h_new_data)(ci) +=  (*dh_dT_data)(ci) * ((*T_new_data)(ci) - (*T_pre_data)(ci)) ; //+  0.5*(*d2h_dT_data)(ci) * ((*T_new_data)(ci) - (*T_pre_data)(ci)) * ((*T_new_data)(ci) - (*T_pre_data)(ci)) ;
                // (*rho_h_n_data)(ci) +=  (*drhoh_dT_data)(ci) * ((*T_new_data)(ci) - (*T_pre_data)(ci)) ; 
                // (*h_new_data)(ci) = (*rho_h_n_data)(ci)/(*d_rho_new_data);

            }

        }
    }
    return;
} // updateEnthalpy

void
EnthalpyHierarchyIntegrator::computeRhoEnthalpyDerivative(int drhoh_dT_idx, const int T_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

     double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
         const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            Pointer<CellData<NDIM, double> > drhoh_dT_data = patch->getPatchData(drhoh_dT_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double T = (*T_data)(ci); 
                const double H = (*H_data)(ci); 

                const double delta_T = d_liquidus_temperature - d_solidus_temperature; 
                

                double hs_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_solid * h_s)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);
                 double hl_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_liquidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_liquid * h_l)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                 double R0 = (1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid; 
                 double R1 = (*H_data)(ci)*(d_rho_liquid - d_rho_solid); 

                 double H0 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) + (*H_data)(ci)*d_rho_solid * h_s; 
                 double H1 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas*(d_liquidus_temperature - d_solidus_temperature) + 
                         (*H_data)(ci)*((d_rho_liquid - d_rho_solid)*h_s + d_rho_liquid*( d_latent_heat + d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature)));

                double Rho_T = (((*H_data)(ci) - 1.0)*d_solidus_temperature*d_rho_gas - (*H_data)(ci)* d_solidus_temperature*d_rho_liquid + (*H_data)(ci)*(*T_data)(ci)*(d_rho_liquid - d_rho_solid) + d_liquidus_temperature*(d_rho_gas - (*H_data)(ci)*d_rho_gas + (*H_data)(ci)*d_rho_solid))/ 
                (d_liquidus_temperature - d_solidus_temperature); 

                double C_LG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_liquid*d_specific_heat_liquid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                double C_SG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_solid*d_specific_heat_solid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);


            //     // if ((*H_data)(ci) >= H_LIM)
            //     // {
            //         if ((*T_data)(ci) < d_solidus_temperature)
            //         {
            //             (*dh_dT_data)(ci) = //C_SG;
            //         }
            //         else if ((*T_data)(ci) > d_liquidus_temperature)
            //         {
            //             (*dh_dT_data)(ci) =  C_LG;
            //         }
                    
            // //         else if ((*T_data)(ci) >= d_solidus_temperature && (*T_data)(ci) <= d_liquidus_temperature)
            // //         {
            // //          (*dh_dT_data)(ci) =d_rho_liquid*d_rho_solid*(h_l - h_s)*(d_liquidus_temperature - d_solidus_temperature)/
			// //                                              (std::pow(((*T_data)(ci)*(d_rho_liquid - d_rho_solid) +d_liquidus_temperature*d_rho_solid - d_solidus_temperature*d_rho_liquid ),2.0)); 
		    // // }
            //         else
            //         {
            //             (*dh_dT_data)(ci) =  (H1*R0 - H0*R1)/((d_liquidus_temperature - d_solidus_temperature)*std::pow(Rho_T ,2.0)); 
            //         }
            //     // }
            //     // else
            //     // {
            //     //     (*dh_dT_data)(ci) = d_specific_heat_gas;
            //     // }

            
            if (T < d_solidus_temperature)
{
    (*drhoh_dT_data)(ci) =
        (1.0 - H) * d_rho_gas * d_specific_heat_gas +
        H * d_rho_solid * d_specific_heat_solid;
}
else if (T > d_liquidus_temperature)
{
    (*drhoh_dT_data)(ci) =
        (1.0 - H) * d_rho_gas * d_specific_heat_gas +
        H * d_rho_liquid * d_specific_heat_liquid;
}
else
{
    const double H1 =
        (1.0 - H) *
            d_rho_gas *
            d_specific_heat_gas *
            delta_T
        +
        H *
            ((d_rho_liquid - d_rho_solid) * h_s +
             d_rho_liquid *
                 (d_latent_heat +
                  d_specific_heat_mushy * delta_T));

    (*drhoh_dT_data)(ci) = H1 / delta_T;
}

            }
        }
    }
    return;
} // computeRhoEnthalpyDerivative



void
EnthalpyHierarchyIntegrator::computeEnthalpyDerivative(int dh_dT_idx, const int T_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

     double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
         const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;

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

                const double T = (*T_data)(ci); 
                const double H = (*H_data)(ci); 

                const double delta_T = d_liquidus_temperature - d_solidus_temperature; 
                

                double hs_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_solid * h_s)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);
                 double hl_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_liquidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_liquid * h_l)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                 double R0 = (1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid; 
                 double R1 = (*H_data)(ci)*(d_rho_liquid - d_rho_solid); 

                 double H0 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) + (*H_data)(ci)*d_rho_solid * h_s; 
                 double H1 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas*(d_liquidus_temperature - d_solidus_temperature) + 
                         (*H_data)(ci)*((d_rho_liquid - d_rho_solid)*h_s + d_rho_liquid*( d_latent_heat + d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature)));

                double Rho_T = (((*H_data)(ci) - 1.0)*d_solidus_temperature*d_rho_gas - (*H_data)(ci)* d_solidus_temperature*d_rho_liquid + (*H_data)(ci)*(*T_data)(ci)*(d_rho_liquid - d_rho_solid) + d_liquidus_temperature*(d_rho_gas - (*H_data)(ci)*d_rho_gas + (*H_data)(ci)*d_rho_solid))/ 
                (d_liquidus_temperature - d_solidus_temperature); 

                double C_LG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_liquid*d_specific_heat_liquid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                double C_SG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_solid*d_specific_heat_solid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);


                // if ((*H_data)(ci) >= H_LIM)
                // {
                    if ((*T_data)(ci) < d_solidus_temperature)
                    {
                        (*dh_dT_data)(ci) = C_SG;
                    }
                    else if ((*T_data)(ci) > d_liquidus_temperature)
                    {
                        (*dh_dT_data)(ci) =  C_LG;
                    }
                    
            //         else if ((*T_data)(ci) >= d_solidus_temperature && (*T_data)(ci) <= d_liquidus_temperature)
            //         {
            //          (*dh_dT_data)(ci) =d_rho_liquid*d_rho_solid*(h_l - h_s)*(d_liquidus_temperature - d_solidus_temperature)/
			//                                              (std::pow(((*T_data)(ci)*(d_rho_liquid - d_rho_solid) +d_liquidus_temperature*d_rho_solid - d_solidus_temperature*d_rho_liquid ),2.0)); 
		    // }
                    else
                    {
                        (*dh_dT_data)(ci) =  (H1*R0 - H0*R1)/((d_liquidus_temperature - d_solidus_temperature)*std::pow(Rho_T ,2.0)); 
                    }
                // }
                // else
                // {
                //     (*dh_dT_data)(ci) = d_specific_heat_gas;
                // }

            
        

            }
        }
    }
    return;
} // computeRhoEnthalpyDerivative



void
EnthalpyHierarchyIntegrator::computeEnthalpySecondDerivative(int d2h_dT_idx, const int T_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

     double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
         const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            Pointer<CellData<NDIM, double> > d2h_dT_data = patch->getPatchData(d2h_dT_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double T = (*T_data)(ci); 
                const double H = (*H_data)(ci); 

                const double delta_T = d_liquidus_temperature - d_solidus_temperature; 
                

                double hs_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_solid * h_s)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);
                 double hl_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_liquidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_liquid * h_l)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                 double R0 = (1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid; 
                 double R1 = (*H_data)(ci)*(d_rho_liquid - d_rho_solid); 

                 double H0 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) + (*H_data)(ci)*d_rho_solid * h_s; 
                 double H1 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas*(d_liquidus_temperature - d_solidus_temperature) + 
                         (*H_data)(ci)*((d_rho_liquid - d_rho_solid)*h_s + d_rho_liquid*( d_latent_heat + d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature)));

                double Rho_T = (((*H_data)(ci) - 1.0)*d_solidus_temperature*d_rho_gas - (*H_data)(ci)* d_solidus_temperature*d_rho_liquid + (*H_data)(ci)*(*T_data)(ci)*(d_rho_liquid - d_rho_solid) + d_liquidus_temperature*(d_rho_gas - (*H_data)(ci)*d_rho_gas + (*H_data)(ci)*d_rho_solid))/ 
                (d_liquidus_temperature - d_solidus_temperature); 

                double C_LG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_liquid*d_specific_heat_liquid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                double C_SG = ((1.0 - (*H_data)(ci))*d_rho_gas*d_specific_heat_gas + (*H_data)(ci)*d_rho_solid*d_specific_heat_solid)/
                ((1.0-(*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);


                // if ((*H_data)(ci) >= H_LIM)
                // {
                    if ((*T_data)(ci) < d_solidus_temperature)
                    {
                        (*d2h_dT_data)(ci) = 0.0;
                    }
                    else if ((*T_data)(ci) > d_liquidus_temperature)
                    {
                        (*d2h_dT_data)(ci) =  0.0;
                    }
                    
            //         else if ((*T_data)(ci) >= d_solidus_temperature && (*T_data)(ci) <= d_liquidus_temperature)
            //         {
            //          (*dh_dT_data)(ci) =d_rho_liquid*d_rho_solid*(h_l - h_s)*(d_liquidus_temperature - d_solidus_temperature)/
			//                                              (std::pow(((*T_data)(ci)*(d_rho_liquid - d_rho_solid) +d_liquidus_temperature*d_rho_solid - d_solidus_temperature*d_rho_liquid ),2.0)); 
		    // }
                    else
                    {
                        (*d2h_dT_data)(ci) =
                        -2.0 * R1 * (H1 * R0 - H0 * R1) /
                        (delta_T * delta_T * std::pow(Rho_T, 3.0));
                    }
                // }
                // else
                // {
                //     (*dh_dT_data)(ci) = d_specific_heat_gas;
                // }

            
        

            }
        }
    }
    return;
} // computeRhoEnthalpyDerivative



void
EnthalpyHierarchyIntegrator::computeLiquidFraction(int lf_idx, const int h_idx, const int H_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double h_s = d_specific_heat_solid * (d_solidus_temperature - d_reference_temperature);
    const double h_l = d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature) + h_s + d_latent_heat;

    
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

                 double hs_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_solid * h_s)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid);
                 double hl_cell = ((1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_liquidus_temperature - d_reference_temperature) +  (*H_data)(ci)*d_rho_liquid * h_l)/
                                        ((1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_liquid);

                 double R0 = (1.0 - (*H_data)(ci))*d_rho_gas + (*H_data)(ci)*d_rho_solid; 
                 double R1 = (*H_data)(ci)*(d_rho_liquid - d_rho_solid); 

                 double H0 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas* (d_solidus_temperature - d_reference_temperature) + (*H_data)(ci)*d_rho_solid * h_s; 
                 double H1 = (1.0 - (*H_data)(ci))*d_specific_heat_gas*d_rho_gas*(d_liquidus_temperature - d_solidus_temperature) + 
                         (*H_data)(ci)*((d_rho_liquid - d_rho_solid)*h_s + d_rho_liquid*( d_latent_heat + d_specific_heat_mushy * (d_liquidus_temperature - d_solidus_temperature)));
                // if ((*H_data)(ci) >= H_LIM)
                // {
                    if ((*h_data)(ci) < hs_cell)
                    {
                        (*lf_data)(ci) = 0.0;
                    }
                    else if ((*h_data)(ci) > hl_cell)
                    {
                        (*lf_data)(ci) = (*H_data)(ci);
                    }
                    else
                    {
                        (*lf_data)(ci) = (*H_data)(ci)*(((*h_data)(ci)*R0 - H0)/(H1 - (*h_data)(ci)*R1));
                            // d_rho_solid * (h_s - (*h_data)(ci)) /
                            // ((d_rho_liquid - d_rho_solid) * (*h_data)(ci)-d_rho_liquid * h_l + d_rho_solid * h_s);
                    }
                // }
                // else
                // {
                //     (*lf_data)(ci) = d_gas_liquid_fraction;
                // }
		 (*lf_data)(ci) = clamp((*lf_data)(ci),0.0,(*H_data)(ci));
            }
        }
    }
    return;
} // computeLiquidFraction

void
EnthalpyHierarchyIntegrator::extrapolateLiquidFractionToGasRegion(int lf_new_idx)
{
    // CFL is fixed to be 0.3.
    const double dt = 0.3 * d_lf_extrap_cell_size;
    int current_time_step = 0;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase(); 
    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    const double H_eps = 1.0e-6; 
    
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_new_idx);
            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_new_idx);

            Pointer<CellData<NDIM, double> > lf_extrap_data = patch->getPatchData(d_lf_extrap_current_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double H = (*H_data)(ci);
                const double L = (*lf_data)(ci);

                double varphi = 0.0; 

                if ((*H_data)(ci) > H_eps)
                {
                    varphi = L/H; 
                    varphi = clamp(varphi,0.0,1.0);

                    }
                    (*lf_extrap_data)(ci) = varphi;
            }
        }
    }

    // Initially, copy lf from pcm for extrapolation.
    // d_hier_cc_data_ops->copyData(d_lf_extrap_current_idx, lf_new_idx);
    d_hier_cc_data_ops->copyData(d_lf_extrap_scratch_idx, d_lf_extrap_current_idx);

    d_hier_cc_data_ops->copyData(d_lf_extrap_new_idx, d_lf_extrap_current_idx);




    // Initializing with zero.
    d_hier_cc_data_ops->setToScalar(d_lf_extrap_rhs_scratch_idx, 0.0);
    while (current_time_step < d_lf_extrap_max_num_time_steps)
    {
        // Using CUI limiter, limit \varphi.
        d_lf_extrap_convective_op->interpolateToFaceOnHierarchy(
            d_lf_extrap_interp_idx, d_lf_extrap_scratch_idx, d_u_adv_fc_lf_extrap_current_idx, /*synch_cf_bdry*/ false);

        // compute N = u . \grad{\varphi}.
        d_lf_extrap_convective_op->computeAdvectiveDerivativeOnHierarchy(d_lf_extrap_rhs_scratch_idx,
                                                                         d_lf_extrap_interp_idx,
                                                                         d_u_adv_fc_lf_extrap_current_idx,
                                                                         /*synch_cf_bdry*/ false);

        d_hier_cc_data_ops->scale(d_lf_extrap_rhs_scratch_idx, -dt, d_lf_extrap_rhs_scratch_idx);
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
    const int phi_new_idx = var_db->mapVariableAndContextToIndex(d_phi_var, getNewContext());
    const int phi_scratch_idx = var_db->mapVariableAndContextToIndex(d_phi_var, getScratchContext());

    const int H_new_idx = var_db->mapVariableAndContextToIndex(d_H_var, getNewContext());
    const int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, getScratchContext());

    // Fill the ghost cell values of level set at new_time.
    std::vector<RobinBcCoefStrategy<NDIM>*> phi_bc_coef = getPhysicalBcCoefs(d_phi_var);
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> interface_transaction(2);
    interface_transaction[0] = InterpolationTransactionComponent(phi_scratch_idx,
                                                                 phi_new_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 true,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 phi_bc_coef);

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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            Pointer<SideData<NDIM, double> > u_sc_data = patch->getPatchData(d_u_adv_sc_lf_extrap_current_idx);
            Pointer<SideData<NDIM, double> > normal_data = patch->getPatchData(d_normal_lf_extrap_current_idx);
            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_scratch_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_scratch_idx);

            // computes normal_data = grad(phi_data)
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
                         phi_data->getPointer(),
                         phi_data->getGhostCellWidth().max(),
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
                        -(1.0 - H) * (*normal_data)(si, axis); // -ive sign to point the normal towards gas region;
                }
            }
        }
    }

    // copy side-centered H*n to face-centered.
    copy_side_to_face(u_adv_fc_lf_extrap_current_idx, d_u_adv_sc_lf_extrap_current_idx, d_hierarchy);

    static const bool synch_cf_interface = true;
    d_hier_math_ops->interp(d_u_adv_cc_lf_extrap_current_idx,
                            d_u_adv_cc_lf_extrap_var,
                            d_u_adv_sc_lf_extrap_current_idx,
                            d_u_adv_sc_lf_extrap_var,
                            d_no_fill_op,
                            current_time,
                            synch_cf_interface);

    return;
} // computeAdvectionVelocityForExtrapolation

Pointer<CellConvectiveOperator>
EnthalpyHierarchyIntegrator::getLiquidFractionExtrapConvectiveOperator(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_extrap_var)
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
        if (input_db->keyExists("max_inner_iterations"))
            d_max_inner_iterations = input_db->getInteger("max_inner_iterations");
        if (input_db->keyExists("lf_iteration_error_tolerance"))
            d_lf_iteration_error_tolerance = input_db->getDouble("lf_iteration_error_tolerance");
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
    d_max_inner_iterations = db->getInteger("max_inner_iterations");
    d_lf_iteration_error_tolerance = db->getDouble("lf_iteration_error_tolerance");
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