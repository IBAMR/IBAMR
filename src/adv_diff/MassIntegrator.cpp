// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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

#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/MassIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"

#include "BasePatchHierarchy.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CoarseFineBoundary.h"
#include "FaceData.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchySideDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include <array>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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
#define VC_UPDATE_DENSITY_FC IBAMR_FC_FUNC_(vc_update_density2d, VC_UPDATE_DENSITY2D)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate2d, CUI_EXTRAPOLATE2D)
#define C_TO_F_CWISE_INTERP_2ND_FC IBAMR_FC_FUNC_(ctofcwiseinterp2nd2d, CTOFINTERP2ND2D)
#define GODUNOV_EXTRAPOLATE_FC IBAMR_FC_FUNC_(godunov_extrapolate2d, GODUNOV_EXTRAPOLATE2D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux2d, ADVECT_FLUX2D)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#endif

#if (NDIM == 3)
#define VC_UPDATE_DENSITY_FC IBAMR_FC_FUNC_(vc_update_density3d, VC_UPDATE_DENSITY3D)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate3d, CUI_EXTRAPOLATE3D)
#define C_TO_F_CWISE_INTERP_2ND_FC IBAMR_FC_FUNC_(ctofcwiseinterp2nd3d, CTOFINTERP2ND3D)
#define GODUNOV_EXTRAPOLATE_FC IBAMR_FC_FUNC_(godunov_extrapolate3d, GODUNOV_EXTRAPOLATE3D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux3d, ADVECT_FLUX3D)
#define ADVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#endif

extern "C"
{
    void VC_UPDATE_DENSITY_FC(const double*,
                              const double&,
                              const double&,
                              const double&,
                              const double&,
#if (NDIM == 2)
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,

#endif
#if (NDIM == 3)
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const int&,
                              const int&,
                              const int&,

#endif
                              double*);

    void VC_SSP_RK2_UPDATE_DENSITY_FC(const double*,
                                      const double&,
#if (NDIM == 2)
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const double*,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const double*,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const double*,
                                      const int&,
                                      const int&,
#endif
#if (NDIM == 3)
                                      const int&,
                                      const int&,
                                      const int&,
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
                                      const int&,
                                      const int&,
                                      const double*,
                                      const double*,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const int&,
#endif
                                      double*);

    void CUI_EXTRAPOLATE_FC(
#if (NDIM == 2)
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        double*,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        const double*,
        double*,
        double*
#endif
#if (NDIM == 3)
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        double*,
        double*,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        const double*,
        const double*,
        double*,
        double*,
        double*
#endif
    );

    void C_TO_F_CWISE_INTERP_2ND_FC(
#if (NDIM == 2)
        double*,
        double*,
        const int&,
        const double*,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&
#endif
#if (NDIM == 3)
        double*,
        double*,
        double*,
        const int&,
        const double*,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&
#endif
    );

    void GODUNOV_EXTRAPOLATE_FC(
#if (NDIM == 2)
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        double*,
        double*,
        double*,
        double*,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        const double*,
        double*,
        double*
#endif
#if (NDIM == 3)
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        double*,
        double*,
        double*,
        double*,
        double*,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double*,
        const double*,
        const double*,
        double*,
        double*,
        double*
#endif
    );

    void ADVECT_FLUX_FC(const double&,
#if (NDIM == 2)
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const double*,
                        const double*,
                        const double*,
                        const double*,
                        double*,
                        double*
#endif
#if (NDIM == 3)
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const int&,
                        const double*,
                        const double*,
                        const double*,
                        const double*,
                        const double*,
                        const double*,
                        double*,
                        double*,
                        double*
#endif
    );

    void ADVECT_DERIVATIVE_FC(const double*,
#if (NDIM == 2)
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
#endif
#if (NDIM == 3)
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const int&,
                              const double*,
                              const double*,
                              const double*,
                              const double*,
                              const double*,
                              const double*,
                              const int&,
                              const int&,
                              const int&,
#endif
                              double*);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// NOTE: The number of ghost cells required by the convection scheme depends
// on the chosen convective limiter, which will be set via input file
static const int GPPMG = 4;
static const int GCUIG = 3;
static const int NOGHOSTS = 0;

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_integrate;
static Timer* t_initialize_integrator;
static Timer* t_deallocate_integrator;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

MassIntegrator::MassIntegrator(std::string object_name, Pointer<Database> input_db)
    : d_object_name(std::move(object_name)), d_u_sc_bc_coefs(NDIM), d_rho_sc_bc_coefs(NDIM)
{
    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type"))
        {
            d_density_bdry_extrap_type = input_db->getString("bdry_extrap_type");
        }
        if (input_db->keyExists("density_bdry_extrap_type"))
        {
            d_density_bdry_extrap_type = input_db->getString("density_bdry_extrap_type");
        }
        if (input_db->keyExists("convective_limiter"))
        {
            d_density_convective_limiter =
                IBAMR::string_to_enum<LimiterType>(input_db->getString("convective_limiter"));
        }
        if (input_db->keyExists("density_convective_limiter"))
        {
            d_density_convective_limiter =
                IBAMR::string_to_enum<LimiterType>(input_db->getString("density_convective_limiter"));
        }

        if (input_db->keyExists("density_time_stepping_type"))
        {
            d_density_time_stepping_type =
                IBAMR::string_to_enum<TimeSteppingType>(input_db->getString("density_time_stepping_type"));
        }
        if (input_db->keyExists("enable_logging"))
        {
            d_enable_logging = input_db->getBool("enable_logging");
        }
    }

    switch (d_density_convective_limiter)
    {
    case PPM:
        d_density_limiter_gcw = GPPMG;
        break;
    case CUI:
        d_density_limiter_gcw = GCUIG;
        break;
    default:
        TBOX_ERROR(
            "MassIntegrator::"
            "MassIntegrator():\n"
            << "  unsupported density convective limiter: "
            << IBAMR::enum_to_string<LimiterType>(d_density_convective_limiter) << " \n"
            << "  valid choices are: PPM, CUI\n");
    }

    switch (d_density_time_stepping_type)
    {
    case FORWARD_EULER:
        d_num_steps = 1;
        break;
    case SSPRK2:
        d_num_steps = 2;
        break;
    case SSPRK3:
        d_num_steps = 3;
        break;
    default:
        TBOX_ERROR(
            "MassIntegrator::"
            "MassIntegrator():\n"
            << "  unsupported density time stepping type: "
            << IBAMR::enum_to_string<TimeSteppingType>(d_density_time_stepping_type) << " \n"
            << "  valid choices are: FORWARD_EULER, SSPRK2, SSPRK3\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("MassIntegrator::CONTEXT");

    const std::string V_var_name = "MassIntegrator::V";
    d_V_var = var_db->getVariable(V_var_name);
    if (d_V_var)
    {
        d_V_scratch_idx = var_db->mapVariableAndContextToIndex(d_V_var, var_db->getContext(V_var_name + "::SCRATCH"));
        d_V_composite_idx =
            var_db->mapVariableAndContextToIndex(d_V_var, var_db->getContext(V_var_name + "::COMPOSITE"));
    }
    else
    {
        d_V_var = new FaceVariable<NDIM, double>(V_var_name);
        d_V_scratch_idx = var_db->registerVariableAndContext(
            d_V_var, var_db->getContext(V_var_name + "::SCRATCH"), IntVector<NDIM>(NOGHOSTS));
        d_V_composite_idx = var_db->registerVariableAndContext(
            d_V_var, var_db->getContext(V_var_name + "::COMPOSITE"), IntVector<NDIM>(NOGHOSTS));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_V_scratch_idx >= 0);
    TBOX_ASSERT(d_V_composite_idx >= 0);
#endif

    const std::string rho_cc_name = "MassIntegrator::RHO_CELL_CENTERED";
    d_rho_cc_var = var_db->getVariable(rho_cc_name);
    if (d_rho_cc_var)
    {
        d_rho_cc_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_rho_cc_var, var_db->getContext(rho_cc_name + "::SCRATCH"));
        d_rho_cc_new_idx =
            var_db->mapVariableAndContextToIndex(d_rho_cc_var, var_db->getContext(rho_cc_name + "::NEW"));
    }
    else
    {
        d_rho_cc_var = new CellVariable<NDIM, double>(rho_cc_name);
        d_rho_cc_scratch_idx = var_db->registerVariableAndContext(
            d_rho_cc_var, var_db->getContext(rho_cc_name + "::SCRATCH"), IntVector<NDIM>(d_density_limiter_gcw));
        d_rho_cc_new_idx = var_db->registerVariableAndContext(
            d_rho_cc_var, var_db->getContext(rho_cc_name + "::NEW"), IntVector<NDIM>(NOGHOSTS));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_cc_scratch_idx >= 0);
    TBOX_ASSERT(d_rho_cc_new_idx >= 0);
#endif

    const std::string S_var_name = "MassIntegrator::S";
    d_S_var = var_db->getVariable(S_var_name);
    if (d_S_var)
    {
        d_S_scratch_idx = var_db->mapVariableAndContextToIndex(d_S_var, context);
    }
    else
    {
        d_S_var = new CellVariable<NDIM, double>(S_var_name);
        d_S_scratch_idx = var_db->registerVariableAndContext(d_S_var, context, IntVector<NDIM>(NOGHOSTS));
    }

    const std::string rho_extrap_var_name = d_object_name + "::rho_extrap";
    d_rho_extrap_var = var_db->getVariable(rho_extrap_var_name);
    if (d_rho_extrap_var)
    {
        d_rho_extrap_idx = var_db->mapVariableAndContextToIndex(d_rho_extrap_var, context);
    }
    else
    {
        d_rho_extrap_var = new FaceVariable<NDIM, double>(rho_extrap_var_name, 1);
        d_rho_extrap_idx = var_db->registerVariableAndContext(d_rho_extrap_var, context, IntVector<NDIM>(0));
    }

    d_div_U_N_var = new CellVariable<NDIM, double>(d_object_name + "::div_U_N");
    d_div_U_N_idx = var_db->registerVariableAndContext(
        d_div_U_N_var, var_db->getContext(d_object_name + "div_U_N::SCRATCH"), IntVector<NDIM>(0));

#if !defined(NDEBUG)
    TBOX_ASSERT(d_S_scratch_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator = TimerManager::getManager()->getTimer("IBAMR::MassIntegrator::"
                                                                                     "applyConvectiveOperator()");
                  t_integrate = TimerManager::getManager()->getTimer("IBAMR::MassIntegrator::integrate("
                                                                     ")");
                  t_initialize_integrator = TimerManager::getManager()->getTimer("IBAMR::MassIntegrator::"
                                                                                 "initializeTimeIntegrator()");
                  t_deallocate_integrator = TimerManager::getManager()->getTimer("IBAMR::MassIntegrator::"
                                                                                 "deallocateTimeIntegrator()"););
    return;
} // MassIntegrator

MassIntegrator::~MassIntegrator()
{
    deallocateTimeIntegrator();
    return;
} // ~MassIntegrator

void
MassIntegrator::integrate(double dt)
{
    // Get hierarchy operation object
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_cc_data_ops =
        hier_ops_manager->getOperationsDouble(new CellVariable<NDIM, double>("cc_var"), d_hierarchy, true);
    d_hier_fc_data_ops =
        hier_ops_manager->getOperationsDouble(new FaceVariable<NDIM, double>("fc_var"), d_hierarchy, true);

    IBAMR_TIMER_START(t_integrate)
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("MassIntegrator::integrate():\n"
                   << "  time integrator must be initialized prior to call to "
                      "integrate()\n");
    }

    TBOX_ASSERT(d_rho_cc_current_idx >= 0);
    TBOX_ASSERT(d_V_old_idx >= 0);
    TBOX_ASSERT(d_V_current_idx >= 0);
    TBOX_ASSERT(d_V_new_idx >= 0);
#endif

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(dt, getDt()));
#endif

    if (d_V_old_idx == d_V_current_idx)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_dt_prev <= 0.0);
#endif
        // Ensure that previous time step is set for initial times
        d_dt_prev = dt;
    }
#if !defined(NDEBUG)
    if (!(dt > 0.0))
    {
        TBOX_ERROR("MassIntegrator::integrate():\n"
                   << " invalid time step size dt = " << dt << "\n");
    }
#endif

// Assertions for velocity interpolation and extrapolation
#if !defined(NDEBUG)
    if (d_cycle_num < 0)
    {
        TBOX_ERROR("MassIntegrator::integrate():\n"
                   << "  invalid cycle number = " << d_cycle_num << "\n");
    }
    if (d_dt_prev <= 0.0 && d_density_time_stepping_type != FORWARD_EULER)
    {
        TBOX_ERROR("MassIntegrator::integrate():\n"
                   << "  invalid previous time step size = " << d_dt_prev << "\n");
    }
#endif

    // Fill ghost cell values
    static const bool homogeneous_bc = false;
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    // Fill ghost cells for current density
    std::vector<InterpolationTransactionComponent> rho_transaction_comps(1);

    std::vector<RobinBcCoefStrategy<NDIM>*> rho_cc_bc_coefs(1, d_rho_cc_bc_coefs);
    rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_cc_scratch_idx,
                                                                 d_rho_cc_current_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 d_density_bdry_extrap_type,
                                                                 false,
                                                                 rho_cc_bc_coefs);

    d_hier_rho_bdry_fill->fillData(d_current_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    // Fill ghost cells for the velocity used to compute the density update
    // Note, enforce divergence free condition on all physical boundaries to
    // ensure boundedness of density update
    d_hier_fc_data_ops->copyData(d_V_composite_idx,
                                 d_V_current_idx,
                                 /*interior_only*/ true);

    // Compute the old mass
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double old_mass = d_hier_cc_data_ops->integral(d_rho_cc_current_idx, wgt_cc_idx);

    if (d_enable_logging)
    {
        plog << "MassIntegrator::integrate(): "
                "old mass in the domain = "
             << old_mass << "\n";
    }

    // Compute the convective derivative.
    for (int step = 0; step < d_num_steps; ++step)
    {
        double eval_time = std::numeric_limits<double>::quiet_NaN();
        double w0 = std::numeric_limits<double>::quiet_NaN();
        double w1 = std::numeric_limits<double>::quiet_NaN();
        double w2 = std::numeric_limits<double>::quiet_NaN();
        const double omega = dt / d_dt_prev;
        const double sum_dt = dt + d_dt_prev;

        switch (step)
        {
        case 0:
            eval_time = d_current_time;
            break;
        case 1:
            eval_time = d_current_time + dt;
            if (d_cycle_num > 0)
            {
                w0 = 0.0, w1 = 0.0, w2 = 1.0;
            }
            else
            {
                w0 = -1.0 * omega, w1 = 1.0 + omega, w2 = 0.0;
            }
            break;
        case 2:
            eval_time = d_current_time + dt / 2.0;
            if (d_cycle_num > 0)
            {
                w0 = -0.25 * dt * dt / (d_dt_prev * sum_dt);
                w1 = 0.25 * (2.0 + omega);
                w2 = 0.25 * (dt + 2.0 * d_dt_prev) / sum_dt;
            }
            else
            {
                w0 = -0.5 * omega, w1 = 1.0 + 0.5 * omega, w2 = 0.0;
            }
            break;
        default:
            TBOX_ERROR("This statement should not be reached");
        }
        // Fill ghost cells for new density, if needed
        if (step > 0)
        {
            std::vector<InterpolationTransactionComponent> update_rho_transaction_comps(1);
            update_rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_cc_scratch_idx,
                                                                                d_rho_cc_new_idx,
                                                                                "CONSERVATIVE_LINEAR_REFINE",
                                                                                false,
                                                                                "CONSERVATIVE_COARSEN",
                                                                                d_density_bdry_extrap_type,
                                                                                false,
                                                                                rho_cc_bc_coefs);

            //
            d_hier_rho_bdry_fill->resetTransactionComponents(rho_transaction_comps);
            d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
            d_hier_rho_bdry_fill->fillData(eval_time);
            d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

            // Compute an approximation to velocity at eval_time.
            d_hier_fc_data_ops->linearSum(
                d_V_composite_idx, w0, d_V_old_idx, w1, d_V_current_idx, /*interior_only*/ true);
            d_hier_fc_data_ops->axpy(d_V_composite_idx, w2, d_V_new_idx, d_V_composite_idx, /*interior_only*/ true);
        }

        // Compute the source term
        if (d_S_fcn)
        {
            d_S_fcn->setDataOnPatchHierarchy(d_S_scratch_idx, d_S_var, d_hierarchy, eval_time);
        }
        else
        {
            d_hier_cc_data_ops->setToScalar(d_S_scratch_idx, 0.0);
        }

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const dx = patch_geom->getDx();

                const Box<NDIM>& patch_box = patch->getBox();
                const IntVector<NDIM>& patch_lower = patch_box.lower();
                const IntVector<NDIM>& patch_upper = patch_box.upper();

                Pointer<FaceData<NDIM, double> > V_adv_data = patch->getPatchData(d_V_composite_idx);
                Pointer<CellData<NDIM, double> > R_cur_data = patch->getPatchData(d_rho_cc_current_idx);
                Pointer<CellData<NDIM, double> > R_pre_data = patch->getPatchData(d_rho_cc_scratch_idx);
                Pointer<CellData<NDIM, double> > R_new_data = patch->getPatchData(d_rho_cc_new_idx);

                Pointer<CellData<NDIM, double> > R_src_data = patch->getPatchData(d_S_scratch_idx);

                // Define variables that live on the "faces" of control
                // volumes centered about side-centered staggered velocity
                // components
                const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
                Pointer<FaceData<NDIM, double> > C_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
                Pointer<FaceData<NDIM, double> > T_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
                Pointer<FaceData<NDIM, double> > R_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
                Pointer<FaceData<NDIM, double> > P_half_data =
                    new FaceData<NDIM, double>(patch_box, 1, ghosts); // to store (rho*C*T)^n+half

                std::vector<RobinBcCoefStrategy<NDIM>*> rho_cc_bc_coefs(1, d_rho_cc_bc_coefs);
                // Enforce physical boundary conditions at inflow boundaries.
                AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                    R_pre_data,
                    V_adv_data,
                    patch,
                    rho_cc_bc_coefs,
                    d_solution_time,
                    /*inflow_boundary_only*/ d_density_bdry_extrap_type != "NONE",
                    homogeneous_bc);

                // Upwind cell-centered densities onto faces.
                interpolateCellQuantity(R_half_data,
                                        V_adv_data,
                                        R_pre_data,
                                        patch_lower,
                                        patch_upper,
                                        patch_box,
                                        d_density_convective_limiter);
                //                std::ofstream rho_lim;
                //                rho_lim.open("rho_lim.txt");
                //                R_half_data->print(patch_box, rho_lim);
                //                rho_lim.close();

                // Compute the updated density
                double a0, a1, a2;
                switch (step)
                {
                case 0:
                    a0 = 0.5;
                    a1 = 0.5;
                    a2 = 1.0;
                    break;
                case 1:
                    if (d_density_time_stepping_type == SSPRK2)
                    {
                        a0 = 0.5;
                        a1 = 0.5;
                        a2 = 0.5;
                        break;
                    }
                    else if (d_density_time_stepping_type == SSPRK3)
                    {
                        a0 = 0.75;
                        a1 = 0.25;
                        a2 = 0.25;
                        break;
                    }
                    else
                    {
                        TBOX_ERROR("This statement should not be reached");
                        break;
                    }
                case 2:
                    a0 = 1.0 / 3.0;
                    a1 = 2.0 / 3.0;
                    a2 = 2.0 / 3.0;
                    break;
                default:
                    TBOX_ERROR("This statement should not be reached");
                }
                computeDensityUpdate(R_new_data,
                                     a0,
                                     R_cur_data,
                                     a1,
                                     R_pre_data,
                                     a2,
                                     V_adv_data,
                                     R_half_data,
                                     R_src_data,
                                     patch_box,
                                     dt,
                                     dx);
            }
        }
    }

    // Refill boundary values of newest density
    const double new_time = d_current_time + dt;
    std::vector<InterpolationTransactionComponent> new_rho_transaction_comps(1);
    new_rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_cc_scratch_idx,
                                                                     d_rho_cc_new_idx,
                                                                     "CONSERVATIVE_LINEAR_REFINE",
                                                                     false,
                                                                     "CONSERVATIVE_COARSEN",
                                                                     d_density_bdry_extrap_type,
                                                                     false,
                                                                     d_rho_cc_bc_coefs);
    d_hier_rho_bdry_fill->resetTransactionComponents(new_rho_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(new_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    d_hier_cc_data_ops->copyData(d_rho_cc_new_idx,
                                 d_rho_cc_scratch_idx,
                                 /*interior_only*/ true);

    computeRHSOfContinuityEquation(dt);

    //    std::cout << "L2 norm of update density after" << d_cycle_num << "\t" <<
    //    d_hier_cc_data_ops->L2Norm(d_rho_cc_new_idx) << std::endl;

    // Compute the new mass
    const double new_mass = d_hier_cc_data_ops->integral(d_rho_cc_new_idx, wgt_cc_idx);
    if (d_enable_logging)
    {
        plog << "MassIntegrator::integrate(): "
                "new mass in the domain = "
             << new_mass << "\n";
        plog << "MassIntegrator::integrate(): "
                "change in mass = "
             << new_mass - old_mass << "\n";
    }

    // Reset select options
    d_N_idx = -1;
    d_rho_cc_current_idx = -1;
    d_V_old_idx = -1;
    d_V_current_idx = -1;
    d_V_new_idx = -1;
    d_cycle_num = -1;
    d_dt_prev = -1.0;

    IBAMR_TIMER_STOP(t_integrate);
    return;
} // integrate

void
MassIntegrator::initializeTimeIntegrator(Pointer<BasePatchHierarchy<NDIM> > base_hierarchy)
{
    IBAMR_TIMER_START(t_initialize_integrator);

    if (d_is_initialized) deallocateTimeIntegrator();

    // Get the hierarchy configuration.
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    d_hierarchy = hierarchy;
    d_coarsest_ln = 0;
    d_finest_ln = d_hierarchy->getFinestLevelNumber();

    // Setup the interpolation transaction information.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_rho_transaction_comps.resize(1);
    d_rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_cc_scratch_idx,
                                                                   d_rho_cc_new_idx,
                                                                   "CONSERVATIVE_LINEAR_REFINE",
                                                                   false,
                                                                   "CONSERVATIVE_COARSEN",
                                                                   d_density_bdry_extrap_type,
                                                                   false,
                                                                   d_rho_cc_bc_coefs);

    //    d_v_transaction_comps.resize(1);
    //    d_v_transaction_comps[0] = InterpolationTransactionComponent(d_V_scratch_idx,
    //                                                                 d_V_composite_idx,
    //                                                                 "CONSERVATIVE_LINEAR_REFINE",
    //                                                                 false,
    //                                                                 "CONSERVATIVE_COARSEN",
    //                                                                 d_velocity_bdry_extrap_type,
    //                                                                 false,
    //                                                                 d_u_sc_bc_coefs);

    // Initialize the interpolation operators.
    d_hier_rho_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_rho_bdry_fill->initializeOperatorState(d_rho_transaction_comps, d_hierarchy);

    // Do I need this here?
    //    // Initialize the BC helper.
    //    d_bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    //    d_bc_helper->cacheBcCoefData(d_u_sc_bc_coefs, d_solution_time, d_hierarchy);
    //
    //    // Create the coarse-fine boundary boxes.
    //    d_cf_boundary.resize(d_finest_ln + 1);
    //    const IntVector<NDIM>& max_ghost_width = CF_GHOST_WIDTH;
    //    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    //    {
    //        d_cf_boundary[ln] = CoarseFineBoundary<NDIM>(*d_hierarchy, ln, max_ghost_width);
    //    }

    // Allocate data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_V_scratch_idx)) level->allocatePatchData(d_V_scratch_idx);
        if (!level->checkAllocated(d_V_composite_idx)) level->allocatePatchData(d_V_composite_idx);
        if (!level->checkAllocated(d_rho_cc_scratch_idx)) level->allocatePatchData(d_rho_cc_scratch_idx);
        if (!level->checkAllocated(d_rho_cc_new_idx)) level->allocatePatchData(d_rho_cc_new_idx);
        if (!level->checkAllocated(d_S_scratch_idx)) level->allocatePatchData(d_S_scratch_idx);
        if (!level->checkAllocated(d_rho_extrap_idx)) level->allocatePatchData(d_rho_extrap_idx);
        if (!level->checkAllocated(d_div_U_N_idx)) level->allocatePatchData(d_div_U_N_idx);
    }

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops =
            new HierarchyMathOps("MassIntegrator::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_hier_math_ops);
#endif
    }

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_integrator);
    return;
} // initializeTimeIntegrator

void
MassIntegrator::deallocateTimeIntegrator()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_integrator);

    // Deallocate the communications operators and BC helpers.
    d_hier_rho_bdry_fill.setNull();
    d_bc_helper.setNull();

    // Deallocate data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_V_scratch_idx)) level->deallocatePatchData(d_V_scratch_idx);
        if (level->checkAllocated(d_V_composite_idx)) level->deallocatePatchData(d_V_composite_idx);
        if (level->checkAllocated(d_rho_cc_scratch_idx)) level->deallocatePatchData(d_rho_cc_scratch_idx);
        if (level->checkAllocated(d_rho_cc_new_idx)) level->deallocatePatchData(d_rho_cc_new_idx);
        if (level->checkAllocated(d_S_scratch_idx)) level->deallocatePatchData(d_S_scratch_idx);
        if (level->checkAllocated(d_rho_extrap_idx)) level->deallocatePatchData(d_rho_extrap_idx);
        if (level->checkAllocated(d_div_U_N_idx)) level->deallocatePatchData(d_div_U_N_idx);
    }

    // Deallocate coarse-fine boundary object.
    d_cf_boundary.clear();

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_integrator);
    return;
} // deallocateOperatorState

void
MassIntegrator::setCellCenteredDensityPatchDataIndex(int rho_cc_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(rho_cc_idx >= 0);
#endif
    d_rho_cc_current_idx = rho_cc_idx;
} // setCellCenteredDensityPatchDataIndex

void
MassIntegrator::setCellCenteredDensityBoundaryConditions(RobinBcCoefStrategy<NDIM>*& rho_cc_bc_coefs)
{
    d_rho_cc_bc_coefs = rho_cc_bc_coefs;
    return;
} // setCellCenteredDensityBoundaryConditions

int
MassIntegrator::getUpdatedCellCenteredDensityPatchDataIndex()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_cc_new_idx >= 0);
#endif
    return d_rho_cc_new_idx;
} // getUpdatedCellCenteredDensityPatchDataIndex

void
MassIntegrator::setMassDensitySourceTerm(const Pointer<CartGridFunction> S_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(S_fcn);
#endif
    d_S_fcn = S_fcn;
    return;
} // setMassDensitySourceTerm

void
MassIntegrator::setFluidVelocityPatchDataIndices(int V_old_idx, int V_current_idx, int V_new_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(V_current_idx >= 0);
#endif

    // Set the old velocity if it has been set, otherwise set to current.
    if (V_old_idx >= 0)
    {
        d_V_old_idx = V_old_idx;
    }
    else
    {
        d_V_old_idx = V_current_idx;
    }

    // Set the current velocity
    d_V_current_idx = V_current_idx;

    // Set the new velocity if it has been set, otherwise set to current.
    if (V_new_idx >= 0)
    {
        d_V_new_idx = V_new_idx;
    }
    else
    {
        d_V_new_idx = V_current_idx;
    }
    return;
} // setFluidVelocityPatchDataIndices

void
MassIntegrator::setCycleNumber(int cycle_num)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(cycle_num >= 0);
#endif
    d_cycle_num = cycle_num;
    return;
} // setCycleNumber

void
MassIntegrator::setSolutionTime(double solution_time)
{
    d_solution_time = solution_time;
} // setSolutionTime

void
MassIntegrator::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
} // setTimeInterval

std::pair<double, double>
MassIntegrator::getTimeInterval() const
{
    return std::make_pair(d_current_time, d_new_time);
} // getTimeInterval

double
MassIntegrator::getDt() const
{
    return d_new_time - d_current_time;
} // getDt

void
MassIntegrator::setHierarchyMathOps(Pointer<HierarchyMathOps> hier_math_ops)
{
    d_hier_math_ops = hier_math_ops;
    d_hier_math_ops_external = d_hier_math_ops;
    return;
} // setHierarchyMathOps

Pointer<HierarchyMathOps>
MassIntegrator::getHierarchyMathOps() const
{
    return d_hier_math_ops;
} // getHierarchyMathOps

void
MassIntegrator::setPreviousTimeStepSize(double dt_prev)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(dt_prev > 0.0);
#endif
    d_dt_prev = dt_prev;
    return;
} // setPreviousTimeStepSize

int
MassIntegrator::getRHSOfContinuityEquation()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_div_U_N_idx >= 0);
#endif
    return d_div_U_N_idx;
} // getRHSOfContinuityEquatio
/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////
void
MassIntegrator::interpolateCellQuantity(Pointer<FaceData<NDIM, double> > Q_half_data,
                                        Pointer<FaceData<NDIM, double> > U_adv_data,
                                        const Pointer<CellData<NDIM, double> > Q_data,
                                        const IntVector<NDIM>& patch_lower,
                                        const IntVector<NDIM>& patch_upper,
                                        const Box<NDIM>& patch_box,
                                        const LimiterType& convective_limiter)
{
    const IntVector<NDIM>& Q_data_gcw = Q_data->getGhostCellWidth();
    const IntVector<NDIM>& U_adv_data_gcw = U_adv_data->getGhostCellWidth();
    const IntVector<NDIM>& Q_half_data_gcw = Q_half_data->getGhostCellWidth();

    CellData<NDIM, double>& Q0_data = *Q_data;
    CellData<NDIM, double> Q1_data(patch_box, 1, Q_data_gcw);
#if (NDIM == 3)
    CellData<NDIM, double> Q2_data(patch_box, 1, Q_data_gcw);
#endif
    CellData<NDIM, double> dQ_data(patch_box, 1, Q_data_gcw);
    CellData<NDIM, double> Q_L_data(patch_box, 1, Q_data_gcw);
    CellData<NDIM, double> Q_R_data(patch_box, 1, Q_data_gcw);

    switch (convective_limiter)
    {
    case PPM:
        // Upwind cell-centered densities onto faces.
        GODUNOV_EXTRAPOLATE_FC(
#if (NDIM == 2)
            patch_lower(0),
            patch_upper(0),
            patch_lower(1),
            patch_upper(1),
            Q_data_gcw(0),
            Q_data_gcw(1),
            Q0_data.getPointer(0),
            Q1_data.getPointer(),
            dQ_data.getPointer(),
            Q_L_data.getPointer(),
            Q_R_data.getPointer(),
            U_adv_data_gcw(0),
            U_adv_data_gcw(1),
            Q_half_data_gcw(0),
            Q_half_data_gcw(1),
            U_adv_data->getPointer(0),
            U_adv_data->getPointer(1),
            Q_half_data->getPointer(0),
            Q_half_data->getPointer(1)
#endif
#if (NDIM == 3)
                patch_lower(0),
            patch_upper(0),
            patch_lower(1),
            patch_upper(1),
            patch_lower(2),
            patch_upper(2),
            Q_data_gcw(0),
            Q_data_gcw(1),
            Q_data_gcw(2),
            Q0_data.getPointer(0),
            Q1_data.getPointer(),
            Q2_data.getPointer(),
            dQ_data.getPointer(),
            Q_L_data.getPointer(),
            Q_R_data.getPointer(),
            U_adv_data_gcw(0),
            U_adv_data_gcw(1),
            U_adv_data_gcw(2),
            Q_half_data_gcw(0),
            Q_half_data_gcw(1),
            Q_half_data_gcw(2),
            U_adv_data->getPointer(0),
            U_adv_data->getPointer(1),
            U_adv_data->getPointer(2),
            Q_half_data->getPointer(0),
            Q_half_data->getPointer(1),
            Q_half_data->getPointer(2)
#endif
        );
        break;

    case CUI:
        // Upwind cell-centered densities onto faces.
        CUI_EXTRAPOLATE_FC(
#if (NDIM == 2)
            patch_lower(0),
            patch_upper(0),
            patch_lower(1),
            patch_upper(1),
            Q_data_gcw(0),
            Q_data_gcw(1),
            Q0_data.getPointer(0),
            Q1_data.getPointer(),
            U_adv_data_gcw(0),
            U_adv_data_gcw(1),
            Q_half_data_gcw(0),
            Q_half_data_gcw(1),
            U_adv_data->getPointer(0),
            U_adv_data->getPointer(1),
            Q_half_data->getPointer(0),
            Q_half_data->getPointer(1)
#endif
#if (NDIM == 3)
                patch_lower(0),
            patch_upper(0),
            patch_lower(1),
            patch_upper(1),
            patch_lower(2),
            patch_upper(2),
            Q_data_gcw(0),
            Q_data_gcw(1),
            Q_data_gcw(2),
            Q0_data.getPointer(0),
            Q1_data.getPointer(),
            Q2_data.getPointer(),
            U_adv_data_gcw(0),
            U_adv_data_gcw(1),
            U_adv_data_gcw(2),
            Q_half_data_gcw(0),
            Q_half_data_gcw(1),
            Q_half_data_gcw(2),
            U_adv_data->getPointer(0),
            U_adv_data->getPointer(1),
            U_adv_data->getPointer(2),
            Q_half_data->getPointer(0),
            Q_half_data->getPointer(1),
            Q_half_data->getPointer(2)
#endif
        );

        break;

    default:
        TBOX_ERROR(
            "MassIntegrator::"
            "interpolateCellQuantity():\n"
            << "  unsupported convective limiter: " << IBAMR::enum_to_string<LimiterType>(convective_limiter) << " \n"
            << "  valid choices are: PPM, CUI\n");
    }
} // interpolateCellQuantity

void
MassIntegrator::computeDensityUpdate(Pointer<CellData<NDIM, double> > R_data,
                                     const double& a0,
                                     const Pointer<CellData<NDIM, double> > R0_data,
                                     const double& a1,
                                     const Pointer<CellData<NDIM, double> > R1_data,
                                     const double& a2,
                                     const Pointer<FaceData<NDIM, double> > U_adv_data,
                                     const Pointer<FaceData<NDIM, double> > R_half_data,
                                     const Pointer<CellData<NDIM, double> > S_data,
                                     const Box<NDIM>& patch_box,
                                     const double& dt,
                                     const double* const dx)
{
    const IntVector<NDIM>& R_data_gcw = R_data->getGhostCellWidth();
    const IntVector<NDIM>& R0_data_gcw = R0_data->getGhostCellWidth();
    const IntVector<NDIM>& R1_data_gcw = R1_data->getGhostCellWidth();
    const IntVector<NDIM>& U_adv_data_gcw = U_adv_data->getGhostCellWidth();
    const IntVector<NDIM>& R_half_data_gcw = R_half_data->getGhostCellWidth();
    const IntVector<NDIM>& S_data_gcw = S_data->getGhostCellWidth();

#if (NDIM == 2)
    VC_UPDATE_DENSITY_FC(dx,
                         dt,
                         a0,
                         a1,
                         a2,
                         patch_box.lower(0),
                         patch_box.upper(0),
                         patch_box.lower(1),
                         patch_box.upper(1),
                         R0_data_gcw(0),
                         R0_data_gcw(1),
                         R0_data->getPointer(0),
                         R1_data_gcw(0),
                         R1_data_gcw(1),
                         R1_data->getPointer(0),
                         U_adv_data_gcw(0),
                         U_adv_data_gcw(1),
                         U_adv_data->getPointer(0),
                         U_adv_data->getPointer(1),
                         R_half_data_gcw(0),
                         R_half_data_gcw(1),
                         R_half_data->getPointer(0),
                         R_half_data->getPointer(1),
                         S_data_gcw(0),
                         S_data_gcw(1),
                         S_data->getPointer(0),
                         R_data_gcw(0),
                         R_data_gcw(1),
                         R_data->getPointer(0));
#endif
#if (NDIM == 3)
    VC_UPDATE_DENSITY_FC(dx,
                         dt,
                         a0,
                         a1,
                         a2,
                         patch_box.lower(0),
                         patch_box.upper(0),
                         patch_box.lower(1),
                         patch_box.upper(1),
                         patch_box.lower(2),
                         patch_box.upper(2),
                         R0_data_gcw(0),
                         R0_data_gcw(1),
                         R0_data_gcw(2),
                         R0_data->getPointer(0),
                         R1_data_gcw(0),
                         R1_data_gcw(1),
                         R1_data_gcw(2),
                         R1_data->getPointer(0),
                         U_adv_data_gcw(0),
                         U_adv_data_gcw(1),
                         U_adv_data_gcw(2),
                         U_adv_data->getPointer(0),
                         U_adv_data->getPointer(1),
                         U_adv_data->getPointer(2),
                         R_half_data_gcw(0),
                         R_half_data_gcw(1),
                         R_half_data_gcw(2),
                         R_half_data->getPointer(0),
                         R_half_data->getPointer(1),
                         R_half_data->getPointer(2),
                         S_data_gcw(0),
                         S_data_gcw(1),
                         S_data_gcw(2),
                         S_data->getPointer(0),
                         R_data_gcw(0),
                         R_data_gcw(1),
                         R_data_gcw(2),
                         R_data->getPointer(0));
#endif
} // computeDensityUpdate

void
MassIntegrator::computeRHSOfContinuityEquation(const double dt)
{
    // Get hierarchy operation object
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_cc_data_ops =
        hier_ops_manager->getOperationsDouble(new CellVariable<NDIM, double>("cc_var"), d_hierarchy, true);
    d_hier_fc_data_ops =
        hier_ops_manager->getOperationsDouble(new FaceVariable<NDIM, double>("fc_var"), d_hierarchy, true);

    static const bool homogeneous_bc = false;
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            Pointer<FaceData<NDIM, double> > V_adv_data = patch->getPatchData(d_V_new_idx);
            Pointer<CellData<NDIM, double> > R_cur_data = patch->getPatchData(d_rho_cc_current_idx);
            Pointer<CellData<NDIM, double> > R_pre_data = patch->getPatchData(d_rho_cc_scratch_idx);
            Pointer<CellData<NDIM, double> > R_new_data = patch->getPatchData(d_rho_cc_new_idx);
            Pointer<CellData<NDIM, double> > div_U_N_data = patch->getPatchData(d_div_U_N_idx);
            const IntVector<NDIM>& N_data_gcw = div_U_N_data->getGhostCellWidth();
            const IntVector<NDIM>& V_adv_data_gcw = V_adv_data->getGhostCellWidth();

            // Define variables that live on the "faces" of control
            // volumes centered about side-centered staggered velocity
            // components
            const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
            Pointer<FaceData<NDIM, double> > C_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
            Pointer<FaceData<NDIM, double> > T_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
            Pointer<FaceData<NDIM, double> > R_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
            Pointer<FaceData<NDIM, double> > P_half_data =
                new FaceData<NDIM, double>(patch_box, 1, ghosts); // to store (rho*C*T)^n+half

            std::vector<RobinBcCoefStrategy<NDIM>*> rho_cc_bc_coefs(1, d_rho_cc_bc_coefs);
            // Enforce physical boundary conditions at inflow boundaries.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                R_pre_data,
                V_adv_data,
                patch,
                rho_cc_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_density_bdry_extrap_type != "NONE",
                homogeneous_bc);

            // Upwind cell-centered densities onto faces.
            interpolateCellQuantity(
                R_half_data, V_adv_data, R_pre_data, patch_lower, patch_upper, patch_box, d_density_convective_limiter);

            const IntVector<NDIM>& R_half_data_gcw = R_half_data->getGhostCellWidth();

            ADVECT_DERIVATIVE_FC(dx,
#if (NDIM == 2)
                                 patch_lower(0),
                                 patch_upper(0),
                                 patch_lower(1),
                                 patch_upper(1),
                                 V_adv_data_gcw(0),
                                 V_adv_data_gcw(1),
                                 R_half_data_gcw(0),
                                 R_half_data_gcw(1),
                                 V_adv_data->getPointer(0),
                                 V_adv_data->getPointer(1),
                                 R_half_data->getPointer(0),
                                 R_half_data->getPointer(1),
                                 N_data_gcw(0),
                                 N_data_gcw(1),
#endif
#if (NDIM == 3)
                                 patch_lower(0),
                                 patch_upper(0),
                                 patch_lower(1),
                                 patch_upper(1),
                                 patch_lower(2),
                                 patch_upper(2),
                                 V_adv_data_gcw(0),
                                 V_adv_data_gcw(1),
                                 V_adv_data_gcw(2),
                                 R_half_data_gcw(0),
                                 R_half_data_gcw(1),
                                 R_half_data_gcw(2),
                                 V_adv_data->getPointer(0),
                                 V_adv_data->getPointer(1),
                                 V_adv_data->getPointer(2),
                                 R_half_data->getPointer(0),
                                 R_half_data->getPointer(1),
                                 R_half_data->getPointer(2),
                                 N_data_gcw(0),
                                 N_data_gcw(1),
                                 N_data_gcw(2),
#endif
                                 div_U_N_data->getPointer());

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                (*div_U_N_data)(ci) = (*div_U_N_data)(ci) + ((*R_new_data)(ci) - (*R_cur_data)(ci)) / dt;
                (*div_U_N_data)(ci) /= (*R_new_data)(ci);
            }
        }
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
