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
#include "ibamr/INSVCStaggeredConservativeMassMomentumSSPRKIntegrator.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IndexUtilities.h"

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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_integrate;
static Timer* t_initialize_integrator;
static Timer* t_deallocate_integrator;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator(
    std::string object_name,
    Pointer<Database> input_db)
    : INSVCStaggeredConservativeMassMomentumRKIntegrator(object_name, input_db)
{
    switch (d_density_time_stepping_type)
    {
    case SSPRK2:
        d_num_steps = 2;
        break;
    case SSPRK3:
        d_num_steps = 3;
        break;
    default:
        TBOX_ERROR(
            "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
            "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator():\n"
            << "  unsupported density time stepping type: "
            << IBAMR::enum_to_string<TimeSteppingType>(d_density_time_stepping_type) << " \n"
            << "  valid choices are: SSPRK2 and SSPRK3\n");
    }

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
                      "applyConvectiveOperator()");
                  t_integrate = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate("
                      ")");
                  t_initialize_integrator = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
                      "initializeSTSIntegrator()");
                  t_deallocate_integrator = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
                      "deallocateSTSIntegrator()"););
    return;
} // INSVCStaggeredConservativeMassMomentumSSPRKIntegrator

void
INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(double dt)
{
    // Get hierarchy operation object
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(new SideVariable<NDIM, double>("sc_var"), d_hierarchy, true);

    IBAMR_TIMER_START(t_integrate)
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate():\n"
                   << "  time integrator must be initialized prior to call to "
                      "integrate()\n");
    }

    TBOX_ASSERT(d_rho_current_idx >= 0);
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
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate():\n"
                   << " invalid time step size dt = " << dt << "\n");
    }
#endif

// Assertions for velocity interpolation and extrapolation
#if !defined(NDEBUG)
    if (d_cycle_num < 0)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate():\n"
                   << "  invalid cycle number = " << d_cycle_num << "\n");
    }
    if (d_dt_prev <= 0.0)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate():\n"
                   << "  invalid previous time step size = " << d_dt_prev << "\n");
    }
#endif

    // Fill ghost cell values
    static const bool homogeneous_bc = false;
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    // Fill ghost cells for current density
    std::vector<InterpolationTransactionComponent> rho_transaction_comps({ { d_rho_scratch_idx,
                                                                             d_rho_current_idx,
                                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                                             false,
                                                                             "CONSERVATIVE_COARSEN",
                                                                             d_density_bdry_extrap_type,
                                                                             false,
                                                                             d_rho_bc_coefs } });
    d_hier_rho_bdry_fill->resetTransactionComponents(rho_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(d_current_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    // Fill ghost cells for the velocity used to compute the density update
    // Note, enforce divergence free condition on all physical boundaries to
    // ensure boundedness of density update
    d_hier_sc_data_ops->copyData(d_V_composite_idx,
                                 d_V_current_idx,
                                 /*interior_only*/ true);
    std::vector<InterpolationTransactionComponent> v_transaction_comps({ { d_V_scratch_idx,
                                                                           d_V_composite_idx,
                                                                           "CONSERVATIVE_LINEAR_REFINE",
                                                                           false,
                                                                           "CONSERVATIVE_COARSEN",
                                                                           d_velocity_bdry_extrap_type,
                                                                           false,
                                                                           d_u_bc_coefs } });
    d_hier_v_bdry_fill->resetTransactionComponents(v_transaction_comps);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_u_bc_coefs, nullptr, d_V_scratch_idx, -1, homogeneous_bc);
    d_hier_v_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_v_bdry_fill->fillData(d_current_time);
    d_bc_helper->enforceDivergenceFreeConditionAtBoundary(
        d_V_scratch_idx, d_coarsest_ln, d_finest_ln, StaggeredStokesPhysicalBoundaryHelper::ALL_BDRY);
    enforceDivergenceFreeConditionAtCoarseFineInterface(d_V_scratch_idx);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);
    d_hier_v_bdry_fill->resetTransactionComponents(d_v_transaction_comps);

    // Compute the old mass
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    const double old_mass = d_hier_sc_data_ops->integral(d_rho_current_idx, wgt_sc_idx);
    if (d_enable_logging)
    {
        plog << "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(): "
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
        // Fill ghost cells for new density and velocity, if needed
        if (step > 0)
        {
            std::vector<InterpolationTransactionComponent> update_transaction_comps({ { d_rho_scratch_idx,
                                                                                        d_rho_new_idx,
                                                                                        "CONSERVATIVE_LINEAR_REFINE",
                                                                                        false,
                                                                                        "CONSERVATIVE_COARSEN",
                                                                                        d_density_bdry_extrap_type,
                                                                                        false,
                                                                                        d_rho_bc_coefs } });
            d_hier_rho_bdry_fill->resetTransactionComponents(update_transaction_comps);
            d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
            d_hier_rho_bdry_fill->fillData(eval_time);
            d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

            // Compute an approximation to velocity at eval_time Note, enforce
            // divergence free condition on all physical boundaries to ensure
            // boundedness of density update
            d_hier_sc_data_ops->linearSum(
                d_V_composite_idx, w0, d_V_old_idx, w1, d_V_current_idx, /*interior_only*/ true);
            d_hier_sc_data_ops->axpy(d_V_composite_idx, w2, d_V_new_idx, d_V_composite_idx, /*interior_only*/ true);
            std::vector<InterpolationTransactionComponent> v_update_transaction_comps({ { d_V_scratch_idx,
                                                                                          d_V_composite_idx,
                                                                                          "CONSERVATIVE_LINEAR_REFINE",
                                                                                          false,
                                                                                          "CONSERVATIVE_COARSEN",
                                                                                          d_velocity_bdry_extrap_type,
                                                                                          false,
                                                                                          d_u_bc_coefs } });
            d_hier_v_bdry_fill->resetTransactionComponents(v_update_transaction_comps);
            StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
                d_u_bc_coefs, nullptr, d_V_scratch_idx, -1, homogeneous_bc);
            d_hier_v_bdry_fill->setHomogeneousBc(homogeneous_bc);
            d_hier_v_bdry_fill->fillData(eval_time);
            d_bc_helper->enforceDivergenceFreeConditionAtBoundary(
                d_V_scratch_idx, d_coarsest_ln, d_finest_ln, StaggeredStokesPhysicalBoundaryHelper::ALL_BDRY);
            enforceDivergenceFreeConditionAtCoarseFineInterface(d_V_scratch_idx);
            StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);
            d_hier_v_bdry_fill->resetTransactionComponents(d_v_transaction_comps);
        }

        // Compute the source term
        if (d_S_fcn)
        {
            d_S_fcn->setDataOnPatchHierarchy(d_S_scratch_idx, d_S_var, d_hierarchy, eval_time);
        }
        else
        {
            d_hier_sc_data_ops->setToScalar(d_S_scratch_idx, 0.0);
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

                Pointer<SideData<NDIM, double> > N_data = patch->getPatchData(d_N_idx);
                Pointer<SideData<NDIM, double> > V_data = patch->getPatchData(d_V_scratch_idx);
                Pointer<SideData<NDIM, double> > R_cur_data = patch->getPatchData(d_rho_current_idx);
                Pointer<SideData<NDIM, double> > R_pre_data = patch->getPatchData(d_rho_scratch_idx);
                Pointer<SideData<NDIM, double> > R_new_data = patch->getPatchData(d_rho_new_idx);
                Pointer<SideData<NDIM, double> > R_src_data = patch->getPatchData(d_S_scratch_idx);
                Pointer<SideData<NDIM, double> > E_data = patch->getPatchData(d_E_scratch_idx);

                // Define variables that live on the "faces" of control
                // volumes centered about side-centered staggered velocity
                // components
                const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
                std::array<Box<NDIM>, NDIM> side_boxes;
                std::array<Pointer<FaceData<NDIM, double> >, NDIM> V_adv_data;
                std::array<Pointer<FaceData<NDIM, double> >, NDIM> V_half_data;
                std::array<Pointer<FaceData<NDIM, double> >, NDIM> R_half_data;
                std::array<Pointer<FaceData<NDIM, double> >, NDIM> P_half_data;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                    V_adv_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                    V_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                    R_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                    P_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                }
                // Interpolate velocity components onto "faces" using simple averages.
                computeAdvectionVelocity(V_adv_data, V_data, patch_lower, patch_upper, side_boxes);

                // Upwind side-centered densities onto faces.
                interpolateSideQuantity(R_half_data,
                                        V_adv_data,
                                        R_pre_data,
                                        patch_lower,
                                        patch_upper,
                                        side_boxes,
                                        d_density_convective_limiter);

                // Compute the convective derivative with the penultimate density and
                // velocity, if necessary
                if ((d_density_time_stepping_type == SSPRK2 && step == 1) ||
                    (d_density_time_stepping_type == SSPRK3 && step == 2))
                {
                    interpolateSideQuantity(V_half_data,
                                            V_adv_data,
                                            V_data,
                                            patch_lower,
                                            patch_upper,
                                            side_boxes,
                                            d_velocity_convective_limiter);

                    IBAMR_TIMER_START(t_apply_convective_operator);

                    computeConvectiveDerivative(
                        N_data, P_half_data, V_adv_data, R_half_data, V_half_data, side_boxes, dx, patch);

                    IBAMR_TIMER_STOP(t_apply_convective_operator);
                }

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
                                     side_boxes,
                                     dt,
                                     dx);

                Pointer<SideData<NDIM, double> > V_cur_data = patch->getPatchData(d_V_current_idx);
                if ((d_density_time_stepping_type == SSPRK2 && step == 1) ||
                    (d_density_time_stepping_type == SSPRK3 && step == 2))
                {
                    computeErrorOfMassConservationEquation(
                        E_data, R_new_data, R_cur_data, V_adv_data, R_half_data, side_boxes, dt, dx);

                    // subtract Error*U from the convective operator.
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                        {
                            SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

                            (*N_data)(si) -= (*V_cur_data)(si) * (*E_data)(si);
                        }
                    }
                }
            }
        }
    }

    // Refill boundary values of newest density
    const double new_time = d_current_time + dt;
    std::vector<InterpolationTransactionComponent> new_transaction_comps({ { d_rho_scratch_idx,
                                                                             d_rho_new_idx,
                                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                                             false,
                                                                             "CONSERVATIVE_COARSEN",
                                                                             d_density_bdry_extrap_type,
                                                                             false,
                                                                             d_rho_bc_coefs } });
    d_hier_rho_bdry_fill->resetTransactionComponents(new_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(new_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    d_hier_sc_data_ops->copyData(d_rho_new_idx,
                                 d_rho_scratch_idx,
                                 /*interior_only*/ true);

    // Compute the new mass
    const double new_mass = d_hier_sc_data_ops->integral(d_rho_new_idx, wgt_sc_idx);
    if (d_enable_logging)
    {
        plog << "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(): "
                "new mass in the domain = "
             << new_mass << "\n";
        plog << "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(): "
                "change in mass = "
             << new_mass - old_mass << "\n";
    }

    // Reset select options
    d_N_idx = -1;
    d_rho_current_idx = -1;
    d_V_old_idx = -1;
    d_V_current_idx = -1;
    d_V_new_idx = -1;
    d_cycle_num = -1;
    d_dt_prev = -1.0;

    IBAMR_TIMER_STOP(t_integrate);
    return;
} // integrate

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
