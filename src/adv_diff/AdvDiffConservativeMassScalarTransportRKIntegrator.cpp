// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2023 by the IBAMR developers
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
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
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
#define CONVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(convect_derivative2d, CONVECT_DERIVATIVE2D)
#define VC_UPDATE_DENSITY_FC IBAMR_FC_FUNC_(vc_update_density2d, VC_UPDATE_DENSITY2D)
#define VC_MASS_CONSERVATION_ERROR_FC IBAMR_FC_FUNC_(vc_mass_conservation_error2d, VC_MASS_CONSERVATION_ERROR2D)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate2d, CUI_EXTRAPOLATE2D)
#define C_TO_F_CWISE_INTERP_2ND_FC IBAMR_FC_FUNC_(ctofcwiseinterp2nd2d, CTOFINTERP2ND2D)
#define GODUNOV_EXTRAPOLATE_FC IBAMR_FC_FUNC_(godunov_extrapolate2d, GODUNOV_EXTRAPOLATE2D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux2d, ADVECT_FLUX2D)
#endif

#if (NDIM == 3)
#define CONVECT_DERIVATIVE_FC IBAMR_FC_FUNC_(convect_derivative3d, CONVECT_DERIVATIVE3D)
#define VC_UPDATE_DENSITY_FC IBAMR_FC_FUNC_(vc_update_density3d, VC_UPDATE_DENSITY3D)
#define VC_MASS_CONSERVATION_ERROR_FC IBAMR_FC_FUNC_(vc_mass_conservation_error3d, VC_MASS_CONSERVATION_ERROR3D)
#define CUI_EXTRAPOLATE_FC IBAMR_FC_FUNC_(cui_extrapolate3d, CUI_EXTRAPOLATE3D)
#define C_TO_F_CWISE_INTERP_2ND_FC IBAMR_FC_FUNC_(ctofcwiseinterp2nd3d, CTOFINTERP2ND3D)
#define GODUNOV_EXTRAPOLATE_FC IBAMR_FC_FUNC_(godunov_extrapolate3d, GODUNOV_EXTRAPOLATE3D)
#define ADVECT_FLUX_FC IBAMR_FC_FUNC_(advect_flux3d, ADVECT_FLUX3D)
#endif

extern "C"
{
    void CONVECT_DERIVATIVE_FC(const double*,
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

    void VC_MASS_CONSERVATION_ERROR_FC(const double*,
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

AdvDiffConservativeMassScalarTransportRKIntegrator::AdvDiffConservativeMassScalarTransportRKIntegrator(
    std::string object_name,
    Pointer<Database> input_db)
    : STSMassFluxIntegrator(object_name, input_db)
{
    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type"))
        {
            d_transport_quantity_bdry_extrap_type = input_db->getString("bdry_extrap_type");
            d_material_property_bdry_extrap_type = input_db->getString("bdry_extrap_type");
        }
        if (input_db->keyExists("transport_quantity_bdry_extrap_type"))
        {
            d_transport_quantity_bdry_extrap_type = input_db->getString("transport_quantity_bdry_extrap_type");
        }
        if (input_db->keyExists("material_property_bdry_extrap_type"))
        {
            d_material_property_bdry_extrap_type = input_db->getString("material_property_bdry_extrap_type");
        }
        if (input_db->keyExists("convective_limiter"))
        {
            d_transport_quantity_convective_limiter =
                IBAMR::string_to_enum<LimiterType>(input_db->getString("convective_limiter"));
            d_material_property_convective_limiter =
                IBAMR::string_to_enum<LimiterType>(input_db->getString("convective_limiter"));
        }
        if (input_db->keyExists("transport_quantity_convective_limiter"))
        {
            d_transport_quantity_convective_limiter =
                IBAMR::string_to_enum<LimiterType>(input_db->getString("transport_quantity_convective_limiter"));
        }
        if (input_db->keyExists("material_property_convective_limiter"))
        {
            d_material_property_convective_limiter =
                IBAMR::string_to_enum<LimiterType>(input_db->getString("material_property_convective_limiter"));
        }
    }

    switch (d_transport_quantity_convective_limiter)
    {
    case PPM:
        d_transport_quantity_limiter_gcw = GPPMG;
        break;
    case CUI:
        d_transport_quantity_limiter_gcw = GCUIG;
        break;
    default:
        TBOX_ERROR(
            "AdvDiffConservativeMassScalarTransportRKIntegrator::"
            "AdvDiffConservativeMassScalarTransportRKIntegrator():\n"
            << "  unsupported transport_quantity convective limiter: "
            << IBAMR::enum_to_string<LimiterType>(d_transport_quantity_convective_limiter) << " \n"
            << "  valid choices are: PPM, CUI\n");
    }

    switch (d_material_property_convective_limiter)
    {
    case PPM:
        d_material_property_limiter_gcw = GPPMG;
        break;
    case CUI:
        d_material_property_limiter_gcw = GCUIG;
        break;
    default:
        TBOX_ERROR(
            "AdvDiffConservativeMassScalarTransportRKIntegrator::"
            "AdvDiffConservativeMassScalarTransportRKIntegrator():\n"
            << "  unsupported specific heat convective limiter: "
            << IBAMR::enum_to_string<LimiterType>(d_material_property_convective_limiter) << " \n"
            << "  valid choices are: PPM, CUI\n");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context =
        var_db->getContext("AdvDiffConservativeMassScalarTransportRKIntegrator::CONTEXT");

    const std::string V_var_name = "AdvDiffConservativeMassScalarTransportRKIntegrator::V";
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

    const std::string rho_cc_name = "AdvDiffConservativeMassScalarTransportRKIntegrator::RHO_CELL_CENTERED";
    d_rho_var = var_db->getVariable(rho_cc_name);
    if (d_rho_var)
    {
        d_rho_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_rho_var, var_db->getContext(rho_cc_name + "::SCRATCH"));
        d_rho_new_idx = var_db->mapVariableAndContextToIndex(d_rho_var, var_db->getContext(rho_cc_name + "::NEW"));
    }
    else
    {
        d_rho_var = new CellVariable<NDIM, double>(rho_cc_name);
        d_rho_scratch_idx = var_db->registerVariableAndContext(
            d_rho_var, var_db->getContext(rho_cc_name + "::SCRATCH"), IntVector<NDIM>(d_density_limiter_gcw));
        d_rho_new_idx = var_db->registerVariableAndContext(
            d_rho_var, var_db->getContext(rho_cc_name + "::NEW"), IntVector<NDIM>(NOGHOSTS));
        d_rho_composite_idx = var_db->registerVariableAndContext(
            d_rho_var, var_db->getContext(rho_cc_name + "::COMPOSITE"), IntVector<NDIM>(NOGHOSTS));
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_scratch_idx >= 0);
    TBOX_ASSERT(d_rho_new_idx >= 0);
#endif

    const std::string Q_cc_name = "AdvDiffConservativeMassScalarTransportRKIntegrator::Q_CELL_CENTERED";
    d_Q_cc_var = var_db->getVariable(Q_cc_name);
    if (d_Q_cc_var)
    {
        d_Q_cc_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_Q_cc_var, var_db->getContext(Q_cc_name + "::SCRATCH"));
        d_Q_cc_composite_idx =
            var_db->mapVariableAndContextToIndex(d_Q_cc_var, var_db->getContext(Q_cc_name + "::COMPOSITE"));
    }
    else
    {
        d_Q_cc_var = new CellVariable<NDIM, double>(Q_cc_name);
        d_Q_cc_scratch_idx = var_db->registerVariableAndContext(
            d_Q_cc_var, var_db->getContext(Q_cc_name + "::SCRATCH"), IntVector<NDIM>(d_transport_quantity_limiter_gcw));
        d_Q_cc_composite_idx = var_db->registerVariableAndContext(
            d_Q_cc_var, var_db->getContext(Q_cc_name + "::COMPOSITE"), IntVector<NDIM>(NOGHOSTS));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_cc_scratch_idx >= 0);
    TBOX_ASSERT(d_Q_cc_composite_idx >= 0);
#endif

    const std::string S_var_name = "AdvDiffConservativeMassScalarTransportRKIntegrator::S";
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

#if !defined(NDEBUG)
    TBOX_ASSERT(d_S_scratch_idx >= 0);
#endif

    const std::string E_var_name = "AdvDiffConservativeMassScalarTransportRKIntegrator::E";
    d_E_var = var_db->getVariable(E_var_name);
    if (d_E_var)
    {
        d_E_scratch_idx = var_db->mapVariableAndContextToIndex(d_E_var, context);
    }
    else
    {
        d_E_var = new CellVariable<NDIM, double>(E_var_name);
        d_E_scratch_idx = var_db->registerVariableAndContext(d_E_var, context, IntVector<NDIM>(NOGHOSTS));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(d_E_scratch_idx >= 0);
#endif

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator =
                      TimerManager::getManager()->getTimer("IBAMR::AdvDiffConservativeMassScalarTransportRKIntegrator::"
                                                           "applyConvectiveOperator()");
                  t_integrate = TimerManager::getManager()->getTimer(
                      "IBAMR::AdvDiffConservativeMassScalarTransportRKIntegrator::integrate("
                      ")");
                  t_initialize_integrator =
                      TimerManager::getManager()->getTimer("IBAMR::AdvDiffConservativeMassScalarTransportRKIntegrator::"
                                                           "initializeSTSIntegrator()");
                  t_deallocate_integrator =
                      TimerManager::getManager()->getTimer("IBAMR::AdvDiffConservativeMassScalarTransportRKIntegrator::"
                                                           "deallocateSTSIntegrator()"););
    return;
} // AdvDiffConservativeMassScalarTransportRKIntegrator

AdvDiffConservativeMassScalarTransportRKIntegrator::~AdvDiffConservativeMassScalarTransportRKIntegrator()
{
    deallocateSTSIntegrator();
    return;
} // ~AdvDiffConservativeMassScalarTransportRKIntegrator

void
AdvDiffConservativeMassScalarTransportRKIntegrator::integrate(double dt)
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
        TBOX_ERROR("AdvDiffConservativeMassScalarTransportRKIntegrator::integrate():\n"
                   << "  time integrator must be initialized prior to call to "
                      "integrate()\n");
    }

    TBOX_ASSERT(d_rho_current_idx >= 0);
    TBOX_ASSERT(d_Q_cc_current_idx >= 0);
    TBOX_ASSERT(d_Q_cc_new_idx >= 0);
    TBOX_ASSERT(d_V_old_idx >= 0);
    TBOX_ASSERT(d_V_current_idx >= 0);
    TBOX_ASSERT(d_V_new_idx >= 0);
#endif

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(dt, getTimeStepSize()));
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
        TBOX_ERROR("AdvDiffConservativeMassScalarTransportRKIntegrator::integrate():\n"
                   << " invalid time step size dt = " << dt << "\n");
    }
#endif

// Assertions for velocity interpolation and extrapolation
#if !defined(NDEBUG)
    if (d_cycle_num < 0)
    {
        TBOX_ERROR("AdvDiffConservativeMassScalarTransportRKIntegrator::integrate():\n"
                   << "  invalid cycle number = " << d_cycle_num << "\n");
    }
    if (d_dt_prev <= 0.0 && d_density_time_stepping_type != FORWARD_EULER)
    {
        TBOX_ERROR("AdvDiffConservativeMassScalarTransportRKIntegrator::integrate():\n"
                   << "  invalid previous time step size = " << d_dt_prev << "\n");
    }
#endif

    // Fill ghost cell values
    static const bool homogeneous_bc = false;
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    // Fill ghost cells for current density
    std::vector<InterpolationTransactionComponent> rho_transaction_comps(1);

    d_hier_cc_data_ops->copyData(d_rho_composite_idx,
                                 d_rho_current_idx,
                                 /*interior_only*/ true);
    std::vector<RobinBcCoefStrategy<NDIM>*> rho_cc_bc_coefs(1, d_rho_cc_bc_coefs);
    rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_scratch_idx,
                                                                 d_rho_composite_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 d_density_bdry_extrap_type,
                                                                 false,
                                                                 d_rho_bc_coefs);

    d_hier_cc_data_ops->copyData(d_Q_cc_composite_idx,
                                 d_Q_cc_current_idx,
                                 /*interior_only*/ true);
    std::vector<InterpolationTransactionComponent> Q_transaction_comps(1);
    std::vector<RobinBcCoefStrategy<NDIM>*> Q_cc_bc_coefs(1, d_Q_cc_bc_coefs);
    Q_transaction_comps[0] = InterpolationTransactionComponent(d_Q_cc_scratch_idx,
                                                               d_Q_cc_composite_idx,
                                                               "CONSERVATIVE_LINEAR_REFINE",
                                                               false,
                                                               "CONSERVATIVE_COARSEN",
                                                               d_density_bdry_extrap_type,
                                                               false,
                                                               Q_cc_bc_coefs);

    double eval_time = d_current_time;
    d_hier_rho_bdry_fill->resetTransactionComponents(rho_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(eval_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    if (d_gamma_cc_var)
    {
        d_hier_cc_data_ops->copyData(d_gamma_cc_composite_idx,
                                     d_gamma_cc_current_idx,
                                     /*interior_only*/ true);

        std::vector<InterpolationTransactionComponent> gamma_transaction_comps(1);
        std::vector<RobinBcCoefStrategy<NDIM>*> gamma_cc_bc_coefs(1, d_gamma_cc_bc_coefs);
        gamma_transaction_comps[0] = InterpolationTransactionComponent(d_gamma_cc_scratch_idx,
                                                                       d_gamma_cc_composite_idx,
                                                                       "CONSERVATIVE_LINEAR_REFINE",
                                                                       false,
                                                                       "CONSERVATIVE_COARSEN",
                                                                       d_density_bdry_extrap_type,
                                                                       false,
                                                                       gamma_cc_bc_coefs);

        d_hier_gamma_bdry_fill->resetTransactionComponents(gamma_transaction_comps);
        d_hier_gamma_bdry_fill->setHomogeneousBc(homogeneous_bc);
        d_hier_gamma_bdry_fill->fillData(eval_time);
        d_hier_gamma_bdry_fill->resetTransactionComponents(d_gamma_transaction_comps);
    }

    d_hier_Q_bdry_fill->resetTransactionComponents(Q_transaction_comps);
    d_hier_Q_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_Q_bdry_fill->fillData(eval_time);
    d_hier_Q_bdry_fill->resetTransactionComponents(d_Q_transaction_comps);

    // Fill ghost cells for the velocity used to compute the density update
    // Note, enforce divergence free condition on all physical boundaries to
    // ensure boundedness of density update
    d_hier_fc_data_ops->copyData(d_V_composite_idx,
                                 d_V_current_idx,
                                 /*interior_only*/ true);

    // Compute the old mass
    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const double old_mass = d_hier_cc_data_ops->integral(d_rho_current_idx, wgt_cc_idx);

    if (d_enable_logging)
    {
        plog << "AdvDiffConservativeMassScalarTransportRKIntegrator::integrate(): "
                "old mass in the domain = "
             << old_mass << "\n";
    }

    // Compute the convective derivative.
    // Fill ghost cells for new density specific heat and transport_quantity, if needed
    if (d_cycle_num > 0)
    {
        eval_time = d_current_time + dt / 2.0;
        // compute rho^n+1/2
        d_hier_cc_data_ops->linearSum(
            d_rho_composite_idx, 0.5, d_rho_new_idx, 0.5, d_rho_current_idx, /*interior_only*/ true);
        std::vector<InterpolationTransactionComponent> update_rho_transaction_comps(1);
        update_rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_scratch_idx,
                                                                            d_rho_composite_idx,
                                                                            "CONSERVATIVE_LINEAR_REFINE",
                                                                            false,
                                                                            "CONSERVATIVE_COARSEN",
                                                                            d_density_bdry_extrap_type,
                                                                            false,
                                                                            d_rho_bc_coefs);

        d_hier_cc_data_ops->linearSum(
            d_Q_cc_composite_idx, 0.5, d_Q_cc_new_idx, 0.5, d_Q_cc_current_idx, /*interior_only*/ true);
        std::vector<InterpolationTransactionComponent> update_Q_transaction_comps(1);
        update_Q_transaction_comps[0] = InterpolationTransactionComponent(d_Q_cc_scratch_idx,
                                                                          d_Q_cc_composite_idx,
                                                                          "CONSERVATIVE_LINEAR_REFINE",
                                                                          false,
                                                                          "CONSERVATIVE_COARSEN",
                                                                          d_density_bdry_extrap_type,
                                                                          false,
                                                                          Q_cc_bc_coefs);
        //
        d_hier_rho_bdry_fill->resetTransactionComponents(update_rho_transaction_comps);
        d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
        d_hier_rho_bdry_fill->fillData(eval_time);
        d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

        if (d_gamma_cc_var)
        {
            d_hier_cc_data_ops->linearSum(
                d_gamma_cc_composite_idx, 0.5, d_gamma_cc_new_idx, 0.5, d_gamma_cc_current_idx, /*interior_only*/ true);
            std::vector<InterpolationTransactionComponent> update_gamma_transaction_comps(1);
            std::vector<RobinBcCoefStrategy<NDIM>*> gamma_cc_bc_coefs(1, d_gamma_cc_bc_coefs);
            update_gamma_transaction_comps[0] = InterpolationTransactionComponent(d_gamma_cc_scratch_idx,
                                                                                  d_gamma_cc_composite_idx,
                                                                                  "CONSERVATIVE_LINEAR_REFINE",
                                                                                  false,
                                                                                  "CONSERVATIVE_COARSEN",
                                                                                  d_density_bdry_extrap_type,
                                                                                  false,
                                                                                  gamma_cc_bc_coefs);
            d_hier_gamma_bdry_fill->resetTransactionComponents(update_gamma_transaction_comps);
            d_hier_gamma_bdry_fill->setHomogeneousBc(homogeneous_bc);
            d_hier_gamma_bdry_fill->fillData(eval_time);
            d_hier_gamma_bdry_fill->resetTransactionComponents(d_gamma_transaction_comps);
        }

        d_hier_Q_bdry_fill->resetTransactionComponents(update_Q_transaction_comps);
        d_hier_Q_bdry_fill->setHomogeneousBc(homogeneous_bc);
        d_hier_Q_bdry_fill->fillData(eval_time);
        d_hier_Q_bdry_fill->resetTransactionComponents(d_Q_transaction_comps);

        // Compute an approximation to velocity at eval_time.
        d_hier_fc_data_ops->linearSum(
            d_V_composite_idx, 0.5, d_V_new_idx, 0.5, d_V_current_idx, /*interior_only*/ true);
    }

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

            Pointer<CellData<NDIM, double> > N_data = patch->getPatchData(d_N_idx);
            Pointer<FaceData<NDIM, double> > V_adv_data = patch->getPatchData(d_V_composite_idx);
            Pointer<CellData<NDIM, double> > R_cur_data = patch->getPatchData(d_rho_current_idx);
            Pointer<CellData<NDIM, double> > R_pre_data = patch->getPatchData(d_rho_scratch_idx);
            Pointer<CellData<NDIM, double> > R_new_data = patch->getPatchData(d_rho_new_idx);
            Pointer<CellData<NDIM, double> > Q_pre_data = patch->getPatchData(d_Q_cc_scratch_idx);

            Pointer<CellData<NDIM, double> > R_src_data = patch->getPatchData(d_S_scratch_idx);

            Pointer<CellData<NDIM, double> > E_data = patch->getPatchData(d_E_scratch_idx);

            // Define variables that live on the "faces" of control
            // volumes centered about side-centered staggered velocity
            // components
            const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
            Pointer<FaceData<NDIM, double> > C_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
            Pointer<FaceData<NDIM, double> > Q_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
            Pointer<FaceData<NDIM, double> > R_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);
            Pointer<FaceData<NDIM, double> > P_half_data = new FaceData<NDIM, double>(patch_box, 1, ghosts);

            std::vector<RobinBcCoefStrategy<NDIM>*> rho_cc_bc_coefs(1, d_rho_cc_bc_coefs);

            // Enforce physical boundary conditions at inflow boundaries.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                R_pre_data,
                V_adv_data,
                patch,
                d_rho_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_density_bdry_extrap_type != "NONE",
                homogeneous_bc);

            // Upwind cell-centered densities onto faces.
            interpolateCellQuantity(
                R_half_data, V_adv_data, R_pre_data, patch_lower, patch_upper, patch_box, d_density_convective_limiter);

            if (d_gamma_cc_var)
            {
                Pointer<CellData<NDIM, double> > C_pre_data = patch->getPatchData(d_gamma_cc_scratch_idx);
                Pointer<CellData<NDIM, double> > C_new_data = patch->getPatchData(d_gamma_cc_new_idx);

                std::vector<RobinBcCoefStrategy<NDIM>*> gamma_cc_bc_coefs(1, d_gamma_cc_bc_coefs);

                // Enforce physical boundary conditions at inflow boundaries.
                AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                    C_pre_data,
                    V_adv_data,
                    patch,
                    gamma_cc_bc_coefs,
                    d_solution_time,
                    /*inflow_boundary_only*/ d_material_property_bdry_extrap_type != "NONE",
                    homogeneous_bc);

                interpolateCellQuantity(C_half_data,
                                        V_adv_data,
                                        C_pre_data,
                                        patch_lower,
                                        patch_upper,
                                        patch_box,
                                        d_material_property_convective_limiter);
            }

            std::vector<RobinBcCoefStrategy<NDIM>*> Q_cc_bc_coefs(1, d_Q_cc_bc_coefs);

            // Enforce physical boundary conditions at inflow boundaries.
            AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
                Q_pre_data,
                V_adv_data,
                patch,
                Q_cc_bc_coefs,
                d_solution_time,
                /*inflow_boundary_only*/ d_transport_quantity_bdry_extrap_type != "NONE",
                homogeneous_bc);

            interpolateCellQuantity(Q_half_data,
                                    V_adv_data,
                                    Q_pre_data,
                                    patch_lower,
                                    patch_upper,
                                    patch_box,
                                    d_transport_quantity_convective_limiter);

            IBAMR_TIMER_START(t_apply_convective_operator);

            // Compute the convective derivative with the penultimate density and
            // velocity, if necessary
            computeConvectiveDerivative(
                N_data, P_half_data, V_adv_data, R_half_data, Q_half_data, C_half_data, patch_box, dx);

            IBAMR_TIMER_STOP(t_apply_convective_operator);

            // Compute the updated density rho_new = a0*rho_cur + a1*rho_pre - a2*dt*R(rho_lim, u_adv)
            const double a0 = 1.0;
            const double a1 = 0.0;
            const double a2 = 1.0;
            computeDensityUpdate(
                R_new_data, a0, R_cur_data, a1, R_pre_data, a2, V_adv_data, R_half_data, R_src_data, patch_box, dt, dx);

            computeErrorOfMassConservationEquation(
                E_data, R_new_data, R_cur_data, V_adv_data, R_half_data, patch_box, dt, dx);

            // subtract Error*Q*gamma from the convective operator.
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if (d_gamma_cc_var)
                {
                    Pointer<CellData<NDIM, double> > C_pre_data = patch->getPatchData(d_gamma_cc_scratch_idx);
                    (*N_data)(ci) -= (*C_pre_data)(ci) * (*Q_pre_data)(ci) * (*E_data)(ci);
                }
                else
                    (*N_data)(ci) -= (*Q_pre_data)(ci) * (*E_data)(ci);
            }
        }
    }

    // Refill boundary values of newest density
    const double new_time = d_current_time + dt;
    std::vector<InterpolationTransactionComponent> new_rho_transaction_comps(1);
    new_rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_scratch_idx,
                                                                     d_rho_new_idx,
                                                                     "CONSERVATIVE_LINEAR_REFINE",
                                                                     false,
                                                                     "CONSERVATIVE_COARSEN",
                                                                     d_density_bdry_extrap_type,
                                                                     false,
                                                                     d_rho_bc_coefs);
    d_hier_rho_bdry_fill->resetTransactionComponents(new_rho_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(new_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    d_hier_cc_data_ops->copyData(d_rho_new_idx,
                                 d_rho_scratch_idx,
                                 /*interior_only*/ true);

    // Compute the new mass
    const double new_mass = d_hier_cc_data_ops->integral(d_rho_new_idx, wgt_cc_idx);
    if (d_enable_logging)
    {
        plog << "AdvDiffConservativeMassScalarTransportRKIntegrator::integrate(): "
                "new mass in the domain = "
             << new_mass << "\n";
        plog << "AdvDiffConservativeMassScalarTransportRKIntegrator::integrate(): "
                "change in mass = "
             << new_mass - old_mass << "\n";
    }

    // Reset select options
    d_N_idx = -1;
    d_rho_current_idx = -1;
    d_V_old_idx = -1;
    d_V_current_idx = -1;
    d_V_new_idx = -1;
    d_Q_cc_current_idx = -1;
    d_Q_cc_new_idx = -1;
    d_gamma_cc_current_idx = -1;
    d_gamma_cc_new_idx = -1;
    d_cycle_num = -1;
    d_dt_prev = -1.0;

    IBAMR_TIMER_STOP(t_integrate);
    return;
} // integrate

void
AdvDiffConservativeMassScalarTransportRKIntegrator::initializeSTSIntegrator(
    Pointer<BasePatchHierarchy<NDIM> > base_hierarchy)
{
    IBAMR_TIMER_START(t_initialize_integrator);

    // Get the hierarchy configuration.
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    d_hierarchy = hierarchy;

    if (d_is_initialized) deallocateSTSIntegrator();

    // Update the level numbers.
    d_coarsest_ln = 0;
    d_finest_ln = d_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    if (d_gamma_cc_var)
    {
        d_gamma_cc_scratch_idx =
            var_db->registerVariableAndContext(d_gamma_cc_var,
                                               var_db->getContext(d_gamma_cc_var->getName() + "::SCRATCH"),
                                               IntVector<NDIM>(d_material_property_limiter_gcw));
        d_gamma_cc_composite_idx = var_db->registerVariableAndContext(
            d_gamma_cc_var, var_db->getContext(d_gamma_cc_var->getName() + "::COMPOSITE"), IntVector<NDIM>(NOGHOSTS));
    }

    // Setup the interpolation transaction information.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_rho_transaction_comps.resize(1);
    d_rho_transaction_comps[0] = InterpolationTransactionComponent(d_rho_scratch_idx,
                                                                   d_rho_composite_idx,
                                                                   "CONSERVATIVE_LINEAR_REFINE",
                                                                   false,
                                                                   "CONSERVATIVE_COARSEN",
                                                                   d_density_bdry_extrap_type,
                                                                   false,
                                                                   d_rho_bc_coefs);

    if (d_gamma_cc_var)
    {
        d_gamma_transaction_comps.resize(1);
        d_gamma_transaction_comps[0] = InterpolationTransactionComponent(d_gamma_cc_scratch_idx,
                                                                         d_gamma_cc_composite_idx,
                                                                         "CONSERVATIVE_LINEAR_REFINE",
                                                                         false,
                                                                         "CONSERVATIVE_COARSEN",
                                                                         d_density_bdry_extrap_type,
                                                                         false,
                                                                         d_gamma_cc_bc_coefs);
    }

    d_Q_transaction_comps.resize(1);
    d_Q_transaction_comps[0] = InterpolationTransactionComponent(d_Q_cc_scratch_idx,
                                                                 d_Q_cc_composite_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 d_density_bdry_extrap_type,
                                                                 false,
                                                                 d_Q_cc_bc_coefs);

    // Initialize the interpolation operators.
    d_hier_rho_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_rho_bdry_fill->initializeOperatorState(d_rho_transaction_comps, d_hierarchy);
    d_hier_gamma_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_gamma_bdry_fill->initializeOperatorState(d_gamma_transaction_comps, d_hierarchy);
    d_hier_Q_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_Q_bdry_fill->initializeOperatorState(d_Q_transaction_comps, d_hierarchy);

    // Allocate data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_V_scratch_idx)) level->allocatePatchData(d_V_scratch_idx);
        if (!level->checkAllocated(d_V_composite_idx)) level->allocatePatchData(d_V_composite_idx);
        if (!level->checkAllocated(d_rho_scratch_idx)) level->allocatePatchData(d_rho_scratch_idx);
        if (!level->checkAllocated(d_rho_new_idx)) level->allocatePatchData(d_rho_new_idx);
        if (!level->checkAllocated(d_rho_composite_idx)) level->allocatePatchData(d_rho_composite_idx);
        if (d_gamma_cc_var)
        {
            if (!level->checkAllocated(d_gamma_cc_scratch_idx)) level->allocatePatchData(d_gamma_cc_scratch_idx);
            if (!level->checkAllocated(d_gamma_cc_composite_idx)) level->allocatePatchData(d_gamma_cc_composite_idx);
        }
        if (!level->checkAllocated(d_Q_cc_scratch_idx)) level->allocatePatchData(d_Q_cc_scratch_idx);
        if (!level->checkAllocated(d_Q_cc_composite_idx)) level->allocatePatchData(d_Q_cc_composite_idx);
        if (!level->checkAllocated(d_S_scratch_idx)) level->allocatePatchData(d_S_scratch_idx);
        if (!level->checkAllocated(d_E_scratch_idx)) level->allocatePatchData(d_E_scratch_idx);
    }

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops = new HierarchyMathOps("AdvDiffConservativeMassScalarTransportRKIntegrator::HierarchyMathOps",
                                               d_hierarchy,
                                               d_coarsest_ln,
                                               d_finest_ln);
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
} // initializeSTSIntegrator

void
AdvDiffConservativeMassScalarTransportRKIntegrator::deallocateSTSIntegrator()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_integrator);

    // Deallocate the communications operators and BC helpers.
    d_hier_rho_bdry_fill.setNull();
    d_hier_gamma_bdry_fill.setNull();
    d_hier_Q_bdry_fill.setNull();
    d_bc_helper.setNull();

    // Deallocate data.
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_V_scratch_idx)) level->deallocatePatchData(d_V_scratch_idx);
        if (level->checkAllocated(d_V_composite_idx)) level->deallocatePatchData(d_V_composite_idx);
        if (level->checkAllocated(d_rho_scratch_idx)) level->deallocatePatchData(d_rho_scratch_idx);
        if (level->checkAllocated(d_rho_new_idx)) level->deallocatePatchData(d_rho_new_idx);
        if (level->checkAllocated(d_rho_composite_idx)) level->deallocatePatchData(d_rho_composite_idx);
        if (d_gamma_cc_var)
        {
            if (level->checkAllocated(d_gamma_cc_scratch_idx)) level->deallocatePatchData(d_gamma_cc_scratch_idx);
            if (level->checkAllocated(d_gamma_cc_composite_idx)) level->deallocatePatchData(d_gamma_cc_composite_idx);
        }
        if (level->checkAllocated(d_Q_cc_scratch_idx)) level->deallocatePatchData(d_Q_cc_scratch_idx);
        if (level->checkAllocated(d_Q_cc_composite_idx)) level->deallocatePatchData(d_Q_cc_composite_idx);
        if (level->checkAllocated(d_S_scratch_idx)) level->deallocatePatchData(d_S_scratch_idx);
        if (level->checkAllocated(d_E_scratch_idx)) level->deallocatePatchData(d_E_scratch_idx);
    }

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_integrator);
    return;
} // deallocateSTSIntegrator

void
AdvDiffConservativeMassScalarTransportRKIntegrator::setCellCenteredMaterialPropertyPatchDataIndex(int gamma_cc_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(gamma_cc_idx >= 0);
#endif
    d_gamma_cc_current_idx = gamma_cc_idx;
} // setCellCenteredMaterialPropertyPatchDataIndex

void
AdvDiffConservativeMassScalarTransportRKIntegrator::setCellCenteredTransportQuantityPatchDataIndex(int Q_cc_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_cc_idx >= 0);
#endif
    d_Q_cc_current_idx = Q_cc_idx;
} // setCellCenteredTransportQuantityPatchDataIndex

void
AdvDiffConservativeMassScalarTransportRKIntegrator::setCellCenteredDensityBoundaryConditions(
    RobinBcCoefStrategy<NDIM>*& rho_cc_bc_coefs)
{
    d_rho_cc_bc_coefs = rho_cc_bc_coefs;
    return;
} // setCellCenteredDensityBoundaryConditions

void
AdvDiffConservativeMassScalarTransportRKIntegrator::setCellCenteredMaterialPropertyBoundaryConditions(
    RobinBcCoefStrategy<NDIM>*& gamma_cc_bc_coefs)
{
    d_gamma_cc_bc_coefs = gamma_cc_bc_coefs;
    return;
} // setCellCenteredMaterialPropertyBoundaryConditions

void
AdvDiffConservativeMassScalarTransportRKIntegrator::setCellCenteredTransportQuantityBoundaryConditions(
    RobinBcCoefStrategy<NDIM>*& Q_cc_bc_coefs)
{
    d_Q_cc_bc_coefs = Q_cc_bc_coefs;
    return;
} // setCellCenteredTransportQuantityBoundaryConditions

int
AdvDiffConservativeMassScalarTransportRKIntegrator::getUpdatedCellCenteredDensityPatchDataIndex()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_new_idx >= 0);
#endif
    return d_rho_new_idx;
} // getUpdatedCellCenteredDensityPatchDataIndex

void
AdvDiffConservativeMassScalarTransportRKIntegrator::setMaterialPropertyPatchDataIndices(int gamma_current_idx,
                                                                                        int gamma_new_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(gamma_current_idx >= 0);
#endif

    // Set the current specific heat
    d_gamma_cc_current_idx = gamma_current_idx;

    // Set the new specific heat if it has been set, otherwise set to current.
    if (gamma_new_idx >= 0)
    {
        d_gamma_cc_new_idx = gamma_new_idx;
    }
    else
    {
        d_gamma_cc_new_idx = gamma_current_idx;
    }
    return;
} // setMaterialPropertyPatchDataIndices

void
AdvDiffConservativeMassScalarTransportRKIntegrator::setTransportQuantityPatchDataIndices(int Q_current_idx,
                                                                                         int Q_new_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_current_idx >= 0);
#endif

    // Set the current transport_quantity
    d_Q_cc_current_idx = Q_current_idx;

    // Set the new transport_quantity if it has been set, otherwise set to current.
    if (Q_new_idx >= 0)
    {
        d_Q_cc_new_idx = Q_new_idx;
    }
    else
    {
        d_Q_cc_new_idx = Q_current_idx;
    }
    return;
} // settransport_quantityPatchDataIndices

void
AdvDiffConservativeMassScalarTransportRKIntegrator::setMaterialPropertyVariable(
    Pointer<CellVariable<NDIM, double> > gamma_var)
{
    d_gamma_cc_var = gamma_var;
    return;
}
/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////
void
AdvDiffConservativeMassScalarTransportRKIntegrator::interpolateCellQuantity(
    Pointer<FaceData<NDIM, double> > Q_half_data,
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
            "AdvDiffConservativeMassScalarTransportRKIntegrator::"
            "interpolateCellQuantity():\n"
            << "  unsupported convective limiter: " << IBAMR::enum_to_string<LimiterType>(convective_limiter) << " \n"
            << "  valid choices are: PPM, CUI\n");
    }
} // interpolateCellQuantity

void
AdvDiffConservativeMassScalarTransportRKIntegrator::computeConvectiveDerivative(
    Pointer<CellData<NDIM, double> > N_data,
    Pointer<FaceData<NDIM, double> > P_half_data,
    const Pointer<FaceData<NDIM, double> > U_adv_data,
    const Pointer<FaceData<NDIM, double> > R_half_data,
    const Pointer<FaceData<NDIM, double> > Q_half_data,
    const Pointer<FaceData<NDIM, double> > G_half_data,
    const Box<NDIM>& patch_box,
    const double* const dx)
{
    static const double dt = 1.0;
    Pointer<FaceData<NDIM, double> > GQ_half_data =
        new FaceData<NDIM, double>(patch_box, 1, 1); // to store (G*Q)^n+half
    const IntVector<NDIM>& U_adv_data_gcw = U_adv_data->getGhostCellWidth();
    const IntVector<NDIM>& P_half_data_gcw = P_half_data->getGhostCellWidth();
    const IntVector<NDIM>& R_half_data_gcw = R_half_data->getGhostCellWidth();
    const IntVector<NDIM>& Q_half_data_gcw = Q_half_data->getGhostCellWidth();
    const IntVector<NDIM>& G_half_data_gcw = G_half_data->getGhostCellWidth();
    const IntVector<NDIM>& GQ_half_data_gcw = GQ_half_data->getGhostCellWidth();
    const IntVector<NDIM>& N_data_gcw = N_data->getGhostCellWidth();

    if (d_gamma_cc_var)
    {
        // compute GQ^(n+half) = G^(n+half) * Q^(n+half).
        ADVECT_FLUX_FC(dt,
#if (NDIM == 2)
                       patch_box.lower(0),
                       patch_box.upper(0),
                       patch_box.lower(1),
                       patch_box.upper(1),
                       G_half_data_gcw(0),
                       G_half_data_gcw(1),
                       Q_half_data_gcw(0),
                       Q_half_data_gcw(1),
                       GQ_half_data_gcw(0),
                       GQ_half_data_gcw(1),
                       G_half_data->getPointer(0),
                       G_half_data->getPointer(1),
                       Q_half_data->getPointer(0),
                       Q_half_data->getPointer(1),
                       GQ_half_data->getPointer(0),
                       GQ_half_data
                           ->getPointer(1)
#endif
#if (NDIM == 3)
                               patch_box.lower(0),
                       patch_box.upper(0),
                       patch_box.lower(1),
                       patch_box.upper(1),
                       patch_box.lower(2),
                       patch_box.upper(2),
                       G_half_data_gcw(0),
                       G_half_data_gcw(1),
                       G_half_data_gcw(2),
                       Q_half_data_gcw(0),
                       Q_half_data_gcw(1),
                       Q_half_data_gcw(2),
                       GQ_half_data_gcw(0),
                       GQ_half_data_gcw(1),
                       GQ_half_data_gcw(2),
                       G_half_data->getPointer(0),
                       G_half_data->getPointer(1),
                       G_half_data->getPointer(2),
                       Q_half_data->getPointer(0),
                       Q_half_data->getPointer(1),
                       Q_half_data->getPointer(2),
                       GQ_half_data->getPointer(0),
                       GQ_half_data->getPointer(1),
                       GQ_half_data->getPointer(2)
#endif
        );
    }
    else
        GQ_half_data->copy(*Q_half_data); // When gamma is not registered, copy Q_half_data into GQ_half_data.

    // compute RGQ^(n+half) = R^(n+half) * GQ^(n+half).
    ADVECT_FLUX_FC(dt,
#if (NDIM == 2)
                   patch_box.lower(0),
                   patch_box.upper(0),
                   patch_box.lower(1),
                   patch_box.upper(1),
                   GQ_half_data_gcw(0),
                   GQ_half_data_gcw(1),
                   R_half_data_gcw(0),
                   R_half_data_gcw(1),
                   P_half_data_gcw(0),
                   P_half_data_gcw(1),
                   GQ_half_data->getPointer(0),
                   GQ_half_data->getPointer(1),
                   R_half_data->getPointer(0),
                   R_half_data->getPointer(1),
                   P_half_data->getPointer(0),
                   P_half_data
                       ->getPointer(1)
#endif
#if (NDIM == 3)
                           patch_box.lower(0),
                   patch_box.upper(0),
                   patch_box.lower(1),
                   patch_box.upper(1),
                   patch_box.lower(2),
                   patch_box.upper(2),
                   GQ_half_data_gcw(0),
                   GQ_half_data_gcw(1),
                   GQ_half_data_gcw(2),
                   R_half_data_gcw(0),
                   R_half_data_gcw(1),
                   R_half_data_gcw(2),
                   P_half_data_gcw(0),
                   P_half_data_gcw(1),
                   P_half_data_gcw(2),
                   GQ_half_data->getPointer(0),
                   GQ_half_data->getPointer(1),
                   GQ_half_data->getPointer(2),
                   R_half_data->getPointer(0),
                   R_half_data->getPointer(1),
                   R_half_data->getPointer(2),
                   P_half_data->getPointer(0),
                   P_half_data->getPointer(1),
                   P_half_data->getPointer(2)
#endif
    );

    // compute div (rho Q gamma u)
#if (NDIM == 2)
    CONVECT_DERIVATIVE_FC(dx,
                          patch_box.lower(0),
                          patch_box.upper(0),
                          patch_box.lower(1),
                          patch_box.upper(1),
                          U_adv_data_gcw(0),
                          U_adv_data_gcw(1),
                          P_half_data_gcw(0),
                          P_half_data_gcw(1),
                          U_adv_data->getPointer(0),
                          U_adv_data->getPointer(1),
                          P_half_data->getPointer(0),
                          P_half_data->getPointer(1),
                          N_data_gcw(0),
                          N_data_gcw(1),
                          N_data->getPointer(0));
#endif
#if (NDIM == 3)
    CONVECT_DERIVATIVE_FC(dx,
                          patch_box.lower(0),
                          patch_box.upper(0),
                          patch_box.lower(1),
                          patch_box.upper(1),
                          patch_box.lower(2),
                          patch_box.upper(2),
                          U_adv_data_gcw(0),
                          U_adv_data_gcw(1),
                          U_adv_data_gcw(2),
                          P_half_data_gcw(0),
                          P_half_data_gcw(1),
                          P_half_data_gcw(2),
                          U_adv_data->getPointer(0),
                          U_adv_data->getPointer(1),
                          U_adv_data->getPointer(2),
                          P_half_data->getPointer(0),
                          P_half_data->getPointer(1),
                          P_half_data->getPointer(2),
                          N_data_gcw(0),
                          N_data_gcw(1),
                          N_data_gcw(2),
                          N_data->getPointer(0));
#endif
} // computeConvectiveDerivative

void
AdvDiffConservativeMassScalarTransportRKIntegrator::computeDensityUpdate(
    Pointer<CellData<NDIM, double> > R_data,
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
AdvDiffConservativeMassScalarTransportRKIntegrator::computeErrorOfMassConservationEquation(
    Pointer<CellData<NDIM, double> > E_data,
    const Pointer<CellData<NDIM, double> > Rnew_data,
    const Pointer<CellData<NDIM, double> > Rold_data,
    const Pointer<FaceData<NDIM, double> > U_adv_data,
    const Pointer<FaceData<NDIM, double> > R_half_data,
    const Box<NDIM>& patch_box,
    const double& dt,
    const double* const dx)
{
    const IntVector<NDIM>& E_data_gcw = E_data->getGhostCellWidth();
    const IntVector<NDIM>& Rnew_data_gcw = Rnew_data->getGhostCellWidth();
    const IntVector<NDIM>& Rold_data_gcw = Rold_data->getGhostCellWidth();
    const IntVector<NDIM>& U_adv_data_gcw = U_adv_data->getGhostCellWidth();
    const IntVector<NDIM>& R_half_data_gcw = R_half_data->getGhostCellWidth();
#if (NDIM == 2)
    VC_MASS_CONSERVATION_ERROR_FC(dx,
                                  dt,
                                  patch_box.lower(0),
                                  patch_box.upper(0),
                                  patch_box.lower(1),
                                  patch_box.upper(1),
                                  Rnew_data_gcw(0),
                                  Rnew_data_gcw(1),
                                  Rnew_data->getPointer(0),
                                  Rold_data_gcw(0),
                                  Rold_data_gcw(1),
                                  Rold_data->getPointer(0),
                                  U_adv_data_gcw(0),
                                  U_adv_data_gcw(1),
                                  U_adv_data->getPointer(0),
                                  U_adv_data->getPointer(1),
                                  R_half_data_gcw(0),
                                  R_half_data_gcw(1),
                                  R_half_data->getPointer(0),
                                  R_half_data->getPointer(1),
                                  E_data_gcw(0),
                                  E_data_gcw(1),
                                  E_data->getPointer(0));
#endif
#if (NDIM == 3)
    VC_MASS_CONSERVATION_ERROR_FC(dx,
                                  dt,
                                  patch_box.lower(0),
                                  patch_box.upper(0),
                                  patch_box.lower(1),
                                  patch_box.upper(1),
                                  patch_box.lower(2),
                                  patch_box.upper(2),
                                  Rnew_data_gcw(0),
                                  Rnew_data_gcw(1),
                                  Rnew_data_gcw(2),
                                  Rnew_data->getPointer(0),
                                  Rold_data_gcw(0),
                                  Rold_data_gcw(1),
                                  Rold_data_gcw(2),
                                  Rold_data->getPointer(0),
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
                                  E_data_gcw(0),
                                  E_data_gcw(1),
                                  E_data_gcw(2),
                                  E_data->getPointer(0));
#endif
} // computeErrorOfMassConservationEquation
//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
