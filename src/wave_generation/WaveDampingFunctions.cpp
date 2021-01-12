// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/WaveDampingFunctions.h"
#include "ibamr/WaveUtilities.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <boost/math/tools/roots.hpp>

#include <cmath>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
/*!
 * A struct holding the required information used by the variable alpha damping
 * method
 */
struct MassConservationFunctor
{
    std::pair<double, double> operator()(const double& alpha);
    double d_dt;
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > d_patch_hierarchy;
    int d_u_idx;
    int d_phi_new_idx;
    int d_phi_current_idx;
    int d_I_idx;
    int d_dI_idx;
    WaveDampingData* d_ptr_wave_damper;
    static double s_newton_guess, s_newton_min, s_newton_max;
};

double MassConservationFunctor::s_newton_guess = 0.0;
double MassConservationFunctor::s_newton_min = 0.0;
double MassConservationFunctor::s_newton_max = 3.0;

// Functor returning the value and 1st derviative.
std::pair<double, double>
MassConservationFunctor::operator()(const double& alpha)
{
    // Zone and wave parameters
    const double x_zone_start = d_ptr_wave_damper->d_x_zone_start;
    const double x_zone_end = d_ptr_wave_damper->d_x_zone_end;
    const double depth = d_ptr_wave_damper->d_depth;

    // Compute the integrad
    const int coarsest_ln = 0;
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    static const int dir = NDIM == 2 ? 1 : 2;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
        const double* const grid_X_lower = grid_geom->getXLower();
        const double* const grid_X_upper = grid_geom->getXUpper();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(d_u_idx);
            Pointer<CellData<NDIM, double> > phi_new_data = patch->getPatchData(d_phi_new_idx);
            Pointer<CellData<NDIM, double> > phi_current_data = patch->getPatchData(d_phi_current_idx);
            Pointer<CellData<NDIM, double> > I_data = patch->getPatchData(d_I_idx);
            Pointer<CellData<NDIM, double> > dI_data = patch->getPatchData(d_dI_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                SAMRAI::hier::Index<NDIM> i = it();
                CellIndex<NDIM> ci(it());
                CellIndex<NDIM> cl = ci;
                CellIndex<NDIM> cu = ci;
                cl(0) -= 1;
                cu(0) += 1;
                const SideIndex<NDIM> sl(ci, /*axis*/ 0, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> su(ci, /*axis*/ 0, SideIndex<NDIM>::Upper);

                // Note: assumes (0,0) in bottom left corner, and water phase denoted by
                // negative level set
                const double phi_new = (*phi_new_data)(ci);
                const double phi_current = (*phi_current_data)(ci);
                const double phi_new_sc_lower = 0.5 * (phi_new + (*phi_new_data)(cl));
                const double phi_new_sc_upper = 0.5 * (phi_new + (*phi_new_data)(cu));
                const double u_sc_lower = (*u_data)(sl);
                const double u_sc_upper = (*u_data)(su);

                const auto x_posn = patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5);

                const auto z_posn =
                    patch_x_lower[dir] + patch_dx[dir] * (static_cast<double>(i(dir) - patch_lower(dir)) + 0.5);

                // Indicator for absorption region
                const bool within_damping = (x_zone_start <= x_posn && x_posn <= x_zone_end);

                if (!within_damping)
                {
                    (*I_data)(ci) = 0.0;
                    (*dI_data)(ci) = 0.0;
                }
                else
                {
                    // 1 - Heaviside
                    double theta_new_sc_lower, theta_new_sc_upper;
                    double vol_cell = 1.0;
                    for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                    auto beta = 1.0 * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));
                    if (phi_new_sc_lower < -beta)
                    {
                        theta_new_sc_lower = 1.0;
                    }
                    else if (std::abs(phi_new_sc_lower) <= beta)
                    {
                        theta_new_sc_lower = 1.0 - (0.5 + 0.5 * phi_new_sc_lower / beta +
                                                    1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_new_sc_lower / beta));
                    }
                    else
                    {
                        theta_new_sc_lower = 0.0;
                    }
                    if (phi_new_sc_upper < -beta)
                    {
                        theta_new_sc_upper = 1.0;
                    }
                    else if (std::abs(phi_new_sc_upper) <= beta)
                    {
                        theta_new_sc_upper = 1.0 - (0.5 + 0.5 * phi_new_sc_upper / beta +
                                                    1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_new_sc_upper / beta));
                    }
                    else
                    {
                        theta_new_sc_upper = 0.0;
                    }

                    // NDIM - 1 area
                    double area = 1.0;
                    for (int d = 1; d < NDIM; ++d) area *= (grid_X_upper[d] - grid_X_lower[d]);

                    const double fac = 5;
                    const double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                    const double xalph = std::pow(1.0 - xtilde, alpha);
                    const double gamma = 1.0 - exp(-fac * xalph);
                    const double eta_new = -phi_new + (z_posn - depth);
                    const double eta_current = -phi_current + (z_posn - depth);
                    const double dtu_dx =
                        (theta_new_sc_upper * u_sc_upper - theta_new_sc_lower * u_sc_lower) / patch_dx[0];
                    const double dgamma_dalpha = fac * exp(-fac * xalph) * xalph * log(1.0 - xtilde);
                    (*I_data)(ci) = d_dt * gamma * dtu_dx + (gamma * eta_new - eta_current) / area;
                    (*dI_data)(ci) = dgamma_dalpha * (d_dt * dtu_dx + eta_new / area);
                }
            }
        }
    }
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_patch_hierarchy, coarsest_ln, finest_ln);
    HierarchyMathOps hier_math_ops("HierarchyMathOps", d_patch_hierarchy);
    hier_math_ops.setPatchHierarchy(d_patch_hierarchy);
    hier_math_ops.resetLevels(coarsest_ln, finest_ln);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
    const double val = hier_cc_data_ops.integral(d_I_idx, wgt_cc_idx);
    const double dval = hier_cc_data_ops.integral(d_dI_idx, wgt_cc_idx);
    TBOX_ASSERT(!std::isnan(val));
    TBOX_ASSERT(!std::isnan(dval));
    return std::make_pair(val, dval);
} // MassConservationFunctor::operator()

} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

namespace WaveDampingFunctions
{
void
callRelaxationZoneCallbackFunction(double /*current_time*/,
                                   double /*new_time*/,
                                   bool /*skip_synchronize_new_state_data*/,
                                   int /*num_cycles*/,
                                   void* ctx)
{
    auto ptr_wave_damper = static_cast<WaveDampingData*>(ctx);
    const double x_zone_start = ptr_wave_damper->d_x_zone_start;
    const double x_zone_end = ptr_wave_damper->d_x_zone_end;
    const double depth = ptr_wave_damper->d_depth;
    const double alpha = ptr_wave_damper->d_alpha;
    const double sign_gas = ptr_wave_damper->d_sign_gas_phase;

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = ptr_wave_damper->d_ins_hier_integrator->getPatchHierarchy();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int u_new_idx = var_db->mapVariableAndContextToIndex(ptr_wave_damper->d_ins_hier_integrator->getVelocityVariable(),
                                                         ptr_wave_damper->d_ins_hier_integrator->getNewContext());

    Pointer<CellVariable<NDIM, double> > phi_cc_var = ptr_wave_damper->d_phi_var;
    if (!phi_cc_var)
        TBOX_ERROR(
            "WaveDampingStrategy::callRelaxationZoneCallbackFunction: Level "
            "set variable must be cell centered!");
    int phi_new_idx =
        var_db->mapVariableAndContextToIndex(phi_cc_var, ptr_wave_damper->d_adv_diff_hier_integrator->getNewContext());
    static const int dir = NDIM == 2 ? 1 : 2;

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_new_idx);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SAMRAI::hier::Index<NDIM> i = it();
                    SideIndex<NDIM> i_side(i, axis, SideIndex<NDIM>::Lower);
                    double x_posn = patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)));
                    const double shift_x = (axis == 0 ? 0.0 : 0.5);
                    x_posn += patch_dx[0] * shift_x;

                    if (x_posn >= x_zone_start && x_posn <= x_zone_end)
                    {
                        const double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                        const double gamma = 1.0 - std::expm1(std::pow(xtilde, alpha)) / std::expm1(1.0);
                        const double target = 0.0;

                        if (axis == 0 || axis == dir)
                        {
                            (*u_data)(i_side, 0) = gamma * (*u_data)(i_side, 0) + (1.0 - gamma) * target;
                        }
                    }
                }
            }
        }
    }

    // Modify the level set in the generation zone.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_new_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                SAMRAI::hier::Index<NDIM> i = it();

                const auto x_posn = patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5);

                const auto dir_posn =
                    patch_x_lower[dir] + patch_dx[dir] * (static_cast<double>(i(dir) - patch_lower(dir)) + 0.5);

                const double dir_surface = sign_gas * (dir_posn - depth);

                if (x_posn >= x_zone_start && x_posn <= x_zone_end)
                {
                    const double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                    const double gamma = 1.0 - std::expm1(std::pow(xtilde, alpha)) / std::expm1(1.0);
                    const double target = dir_surface;

                    (*phi_data)(i, 0) = gamma * (*phi_data)(i, 0) + (1.0 - gamma) * target;
                }
            }
        }
    }

    return;
} // callRelaxationZoneCallbackFunction

void
callConservedWaveAbsorbingCallbackFunction(double current_time,
                                           double new_time,
                                           bool /*skip_synchronize_new_state_data*/,
                                           int /*num_cycles*/,
                                           void* ctx)
{
    // Reference: Hu Zhe, et al., "Numerical Wave Tank Based on A Conserved
    // Wave-Absorbing Method"
    auto ptr_wave_damper = static_cast<WaveDampingData*>(ctx);
    const double x_zone_start = ptr_wave_damper->d_x_zone_start;
    const double x_zone_end = ptr_wave_damper->d_x_zone_end;
    const double depth = ptr_wave_damper->d_depth;
    const double sign_gas = ptr_wave_damper->d_sign_gas_phase;

    // Initialize the required variables
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const std::string var_name = "callConservedWaveAbsorbingCallbackFunction::I_var";
    IntVector<NDIM> no_ghosts = 0;
    Pointer<CellVariable<NDIM, double> > I_var = var_db->getVariable(var_name);
    if (!I_var) I_var = new CellVariable<NDIM, double>(var_name);
    const int I_idx =
        var_db->registerVariableAndContext(I_var, var_db->getContext(var_name + "::INTEGRAND"), no_ghosts);
    const int dI_idx =
        var_db->registerVariableAndContext(I_var, var_db->getContext(var_name + "::DERIVATIVE"), no_ghosts);
#if !defined(NDEBUG)
    TBOX_ASSERT(I_idx >= 0);
    TBOX_ASSERT(dI_idx >= 0);
#endif

    // Allocate patch data
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = ptr_wave_damper->d_ins_hier_integrator->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(I_idx, new_time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(dI_idx, new_time);
    }

    // Get the required patch data
    int u_new_idx = var_db->mapVariableAndContextToIndex(ptr_wave_damper->d_ins_hier_integrator->getVelocityVariable(),
                                                         ptr_wave_damper->d_ins_hier_integrator->getNewContext());
    Pointer<CellVariable<NDIM, double> > phi_cc_var = ptr_wave_damper->d_phi_var;
    if (!phi_cc_var)
        TBOX_ERROR(
            "WaveDampingStrategy::"
            "callConservedWaveAbsorbingCallbackFunction: Level set variable "
            "must be cell "
            "centered!");
    int phi_current_idx = var_db->mapVariableAndContextToIndex(
        phi_cc_var, ptr_wave_damper->d_adv_diff_hier_integrator->getCurrentContext());
    int phi_new_idx =
        var_db->mapVariableAndContextToIndex(phi_cc_var, ptr_wave_damper->d_adv_diff_hier_integrator->getNewContext());

    // Allocate scratch data
    IntVector<NDIM> cell_ghosts = 1;
    const int phi_scratch_idx = var_db->registerVariableAndContext(
        phi_cc_var, var_db->getContext("callConservedWaveAbsorbingCallbackFunction::phi_cc_var::SCRATCH"), cell_ghosts);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(phi_scratch_idx, new_time);
    }

    // Fill ghost cells for level set
    RobinBcCoefStrategy<NDIM>* phi_bc_coef =
        ptr_wave_damper->d_adv_diff_hier_integrator->getPhysicalBcCoefs(phi_cc_var).front();
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent phi_transaction(phi_scratch_idx,
                                                      phi_new_idx,
                                                      "CONSERVATIVE_LINEAR_REFINE",
                                                      false,
                                                      "CONSERVATIVE_COARSEN",
                                                      "LINEAR",
                                                      false,
                                                      phi_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(phi_transaction, patch_hierarchy);
    hier_bdry_fill->fillData(new_time);

    // Create the functor and set its options
    MassConservationFunctor mass_eval_functor;
    mass_eval_functor.d_dt = new_time - current_time;
    mass_eval_functor.d_patch_hierarchy = patch_hierarchy;
    mass_eval_functor.d_u_idx = u_new_idx;
    mass_eval_functor.d_phi_new_idx = phi_scratch_idx;
    mass_eval_functor.d_phi_current_idx = phi_current_idx;
    mass_eval_functor.d_I_idx = I_idx;
    mass_eval_functor.d_dI_idx = dI_idx;
    mass_eval_functor.d_ptr_wave_damper = ptr_wave_damper;

    // Set up Newton itertion
    const double guess = MassConservationFunctor::s_newton_guess;
    const double min = MassConservationFunctor::s_newton_min;
    const double max = MassConservationFunctor::s_newton_max;

    // Accuracy doubles at each step, so stop when just over half of the
    // digits are correct, and rely on that step to polish off the remainder:
    auto get_digits = static_cast<int>(std::numeric_limits<double>::digits * 0.6);
    const boost::uintmax_t maxit = 50;
    boost::uintmax_t it = maxit;

    double alpha = boost::math::tools::newton_raphson_iterate(mass_eval_functor, guess, min, max, get_digits, it);
    if (alpha > max - 0.1 || alpha < 0.0) alpha = 0.0;
    plog << "alpha relaxation dynamic = " << alpha << std::endl;

    // Carry out the damping
    static const int dir = NDIM == 2 ? 1 : 2;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_new_idx);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SAMRAI::hier::Index<NDIM> i = it();
                    SideIndex<NDIM> i_side(i, axis, SideIndex<NDIM>::Lower);
                    double x_posn = patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)));
                    const double shift_x = (axis == 0 ? 0.0 : 0.5);
                    x_posn += patch_dx[0] * shift_x;
                    if (x_posn >= x_zone_start && x_posn <= x_zone_end)
                    {
                        double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                        double gamma = 1.0 - exp(-5.0 * std::pow(1.0 - xtilde, alpha));
                        const double target = 0.0;

                        if (axis == 0 || axis == dir)
                        {
                            (*u_data)(i_side, 0) = gamma * (*u_data)(i_side, 0) + (1.0 - gamma) * target;
                        }
                    }
                }
            }
        }
    }

    // Modify the level set in the generation zone.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_new_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                SAMRAI::hier::Index<NDIM> i = it();

                const auto x_posn = patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5);

                const auto dir_posn =
                    patch_x_lower[dir] + patch_dx[dir] * (static_cast<double>(i(dir) - patch_lower(dir)) + 0.5);
                const double dir_surface = sign_gas * (dir_posn - depth);

                if (x_posn >= x_zone_start && x_posn <= x_zone_end)
                {
                    double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                    double gamma = 1.0 - exp(-5.0 * std::pow(1.0 - xtilde, alpha));
                    double target = dir_surface;

                    (*phi_data)(i, 0) = gamma * (*phi_data)(i, 0) + (1.0 - gamma) * target;
                }
            }
        }
    }

    // Deallocate patch data and remove indices
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(I_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(dI_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(phi_scratch_idx);
    }
    var_db->removePatchDataIndex(I_idx);
    var_db->removePatchDataIndex(dI_idx);
    var_db->removePatchDataIndex(phi_scratch_idx);

    return;

} // callConservedWaveAbsorbingCallbackFunction

} // namespace WaveDampingFunctions

} // namespace IBAMR
