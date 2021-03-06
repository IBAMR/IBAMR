// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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
#include "ibamr/BrinkmanAdvDiffBcHelper.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include <cmath>
#include <string>
#include <utility>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
static const int SIDEG = 1;
static const int NOGHOSTS = 0;
} // namespace
/////////////////////////////// PUBLIC //////////////////////////////////////
BrinkmanAdvDiffBcHelper::BrinkmanAdvDiffBcHelper(std::string object_name,
                                                 Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver)
    : d_object_name(object_name), d_adv_diff_solver(adv_diff_solver)
{
    // Initialize variables
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> no_ghosts = NOGHOSTS;
    const IntVector<NDIM> side_ghosts = SIDEG;
    d_B_var = new SideVariable<NDIM, double>(d_object_name + "::");
    d_B_scratch_idx =
        var_db->registerVariableAndContext(d_B_var, var_db->getContext(d_object_name + "B::SCRATCH"), no_ghosts);
    d_B_chi_scratch_idx =
        var_db->registerVariableAndContext(d_B_var, var_db->getContext(d_object_name + "BCHI::SCRATCH"), no_ghosts);

    d_div_var = new CellVariable<NDIM, double>(d_object_name + "::DIV");
    d_div_B_scratch_idx =
        var_db->registerVariableAndContext(d_div_var, var_db->getContext(d_object_name + "B::SCRATCH"), no_ghosts);
    d_div_B_chi_scratch_idx =
        var_db->registerVariableAndContext(d_div_var, var_db->getContext(d_object_name + "BCHI::SCRATCH"), no_ghosts);
    d_chi_scratch_idx =
        var_db->registerVariableAndContext(d_div_var, var_db->getContext(d_object_name + "CHI::SCRATCH"), no_ghosts);

    d_n_var = new SideVariable<NDIM, double>(d_object_name + "::N");
    d_n_scratch_idx = var_db->registerVariableAndContext(d_n_var, var_db->getContext(d_object_name + "N::SCRATCH"));
    d_n_chi_scratch_idx =
        var_db->registerVariableAndContext(d_n_var, var_db->getContext(d_object_name + "NCHI::SCRATCH"));

    d_div_n_scratch_idx =
        var_db->registerVariableAndContext(d_div_var, var_db->getContext(d_object_name + "N::SCRATCH"));
    d_div_n_chi_scratch_idx =
        var_db->registerVariableAndContext(d_div_var, var_db->getContext(d_object_name + "NCHI::SCRATCH"));

    d_variable_g_scratch_idx =
        var_db->registerVariableAndContext(d_div_var, var_db->getContext(d_object_name + "G::SCRATCH"));

    return;
} // BrinkmanAdvDiffBcHelper

void
BrinkmanAdvDiffBcHelper::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
} // setTimeInterval

void
BrinkmanAdvDiffBcHelper::preprocessBrinkmanAdvDiffBcHelper(double current_time, double /*new_time*/, int /*num_cycles*/)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Allocate required patch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_B_scratch_idx, current_time);
        level->allocatePatchData(d_B_chi_scratch_idx, current_time);
        level->allocatePatchData(d_div_B_scratch_idx, current_time);
        level->allocatePatchData(d_div_B_chi_scratch_idx, current_time);
        level->allocatePatchData(d_chi_scratch_idx, current_time);
        level->allocatePatchData(d_n_scratch_idx, current_time);
        level->allocatePatchData(d_n_chi_scratch_idx, current_time);
        level->allocatePatchData(d_div_n_scratch_idx, current_time);
        level->allocatePatchData(d_div_n_chi_scratch_idx, current_time);
        level->allocatePatchData(d_variable_g_scratch_idx, current_time);
    }

    return;
} // preprocessBrinkmanAdvDiffBcHelper

void
BrinkmanAdvDiffBcHelper::postprocessBrinkmanAdvDiffBcHelper(double /*current_time*/,
                                                            double /*new_time*/,
                                                            int /*num_cycles*/)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Deallocate patch data
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_B_scratch_idx);
        level->deallocatePatchData(d_B_chi_scratch_idx);
        level->deallocatePatchData(d_div_B_scratch_idx);
        level->deallocatePatchData(d_div_B_chi_scratch_idx);
        level->deallocatePatchData(d_chi_scratch_idx);
        level->deallocatePatchData(d_n_scratch_idx);
        level->deallocatePatchData(d_n_chi_scratch_idx);
        level->deallocatePatchData(d_div_n_scratch_idx);
        level->deallocatePatchData(d_div_n_chi_scratch_idx);
        level->deallocatePatchData(d_variable_g_scratch_idx);
    }
    return;
} // postprocessBrinkmanAdvDiffBcHelper

void
BrinkmanAdvDiffBcHelper::setPenaltyCoefficient(double eta_penalty_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(eta_penalty_coef > 0);
#endif
    for (auto& Q_bc_pair : d_Q_bc)
    {
        for (auto& bc_prop : Q_bc_pair.second)
        {
            bc_prop.eta = eta_penalty_coef;
        }
    }
    return;
} // setPenaltyCoefficient

void
BrinkmanAdvDiffBcHelper::setNumInterfaceCells(double num_interface_cells)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(num_interface_cells > 0);
#endif
    for (auto& Q_bc_pair : d_Q_bc)
    {
        for (auto& bc_prop : Q_bc_pair.second)
        {
            bc_prop.num_interface_cells = num_interface_cells;
        }
    }
    return;
} // setNumInterfaceCells

void
BrinkmanAdvDiffBcHelper::registerHomogeneousBC(Pointer<CellVariable<NDIM, double> > Q_var,
                                               Pointer<CellVariable<NDIM, double> > ls_solid_var,
                                               std::string bc_type,
                                               std::string indicator_func_type,
                                               double num_interface_cells,
                                               double eta_penalty_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_var);
    TBOX_ASSERT(ls_solid_var);
    TBOX_ASSERT(num_interface_cells > 0);
    TBOX_ASSERT(eta_penalty_coef > 0);
#endif
    auto bc_type_enum = string_to_enum<AdvDiffBrinkmanPenalizationBcType>(bc_type);
    if (bc_type_enum == UNKNOWN_BRINKMAN_BC_TYPE)
    {
        TBOX_ERROR("BrinkmanAdvDiffBcHelper::registerHomogeneousBC\n"
                   << "  unsupported Brinkman penalization BC type: "
                   << enum_to_string<AdvDiffBrinkmanPenalizationBcType>(bc_type_enum) << " \n"
                   << "  valid choices are: DIRICHLET, NEUMANN and ROBIN\n");
    }

    // Store BC options within the map.
    for (const auto& Q_bc_pair : d_Q_bc)
    {
        // Ensuring that a solid level set is only registered with each transported quantity once.
        if (Q_bc_pair.first == Q_var)
        {
            for (const auto& bc_prop : Q_bc_pair.second)
            {
                if (bc_prop.ls_solid_var == ls_solid_var)
                {
                    TBOX_ERROR("BrinkmanAdvDiffBcHelper::registerHomogeneousBC\n"
                               << "  two separate boundary conditions on level set variable " << ls_solid_var->getName()
                               << " \n"
                               << "  are being registered for transported quantity " << Q_var->getName() << "\n");
                }
            }
        }
    }

    // No user supplied callbacks are needed.
    BrinkmanInhomogeneousBCsFcnPtr callback = nullptr;
    void* ctx = nullptr;

    // Get the user specified indicator function.
    auto indicator_func_type_enum = string_to_enum<IndicatorFunctionType>(indicator_func_type);
    if (indicator_func_type_enum == UNKNOWN_INDICATOR_FUNC_TYPE)
    {
        TBOX_ERROR("BrinkmanAdvDiffBcHelper::registerHomogeneousBC\n"
                   << "  unsupported Indicator Function Type: "
                   << enum_to_string<IndicatorFunctionType>(indicator_func_type_enum) << " \n"
                   << "  valid choices are: SMOOTH and DISCONTINUOUS\n");
    }

    // Store BC properties in a struct.
    BCProperties bc_prop;
    bc_prop.ls_solid_var = ls_solid_var;
    bc_prop.bc_type = bc_type_enum;
    bc_prop.indicator_func_type = indicator_func_type_enum;
    bc_prop.bc_val = 0.0;
    bc_prop.num_interface_cells = num_interface_cells;
    bc_prop.eta = eta_penalty_coef;
    bc_prop.callback = callback;
    bc_prop.ctx = ctx;
    d_Q_bc[Q_var].push_back(bc_prop);

    return;
} // registerHomogeneousBC

void
BrinkmanAdvDiffBcHelper::registerInhomogeneousBC(Pointer<CellVariable<NDIM, double> > Q_var,
                                                 Pointer<CellVariable<NDIM, double> > ls_solid_var,
                                                 std::string bc_type,
                                                 BrinkmanInhomogeneousBCsFcnPtr callback,
                                                 void* ctx,
                                                 std::string indicator_func_type,
                                                 double num_interface_cells,
                                                 double eta_penalty_coef,
                                                 double bc_val)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_var);
    TBOX_ASSERT(ls_solid_var);
    TBOX_ASSERT(callback);
    TBOX_ASSERT(num_interface_cells > 0);
    TBOX_ASSERT(eta_penalty_coef > 0);
#endif
    auto bc_type_enum = string_to_enum<AdvDiffBrinkmanPenalizationBcType>(bc_type);
    if (bc_type_enum == UNKNOWN_BRINKMAN_BC_TYPE)
    {
        TBOX_ERROR("BrinkmanAdvDiffBcHelper::registerInhomogeneousBC\n"
                   << "  unsupported Brinkman penalization BC type: "
                   << enum_to_string<AdvDiffBrinkmanPenalizationBcType>(bc_type_enum) << " \n"
                   << "  valid choices are: DIRICHLET, NEUMANN and ROBIN\n");
    }

    // Store BC options within the map.
    for (const auto& Q_bc_pair : d_Q_bc)
    {
        // Ensuring that a solid level set is only registered with each transported quantity once.
        if (Q_bc_pair.first == Q_var)
        {
            for (const auto& bc_prop : Q_bc_pair.second)
            {
                if (bc_prop.ls_solid_var == ls_solid_var)
                {
                    TBOX_ERROR("BrinkmanAdvDiffBcHelper::registerInhomogeneousBC\n"
                               << "  two separate boundary conditions on level set variable " << ls_solid_var->getName()
                               << " \n"
                               << "  are being registered for transported quantity " << Q_var->getName() << "\n");
                }
            }
        }
    }

    // Specifying a dummy nonzero value for bc_val when inhomogeneous bcs are imposed except
    // for Dirichlet spatially constant g value.
    if (callback != nullptr) bc_val = std::numeric_limits<double>::max();

    // Get the user specified indicator function.
    auto indicator_func_type_enum = string_to_enum<IndicatorFunctionType>(indicator_func_type);
    if (indicator_func_type_enum == UNKNOWN_INDICATOR_FUNC_TYPE)
    {
        TBOX_ERROR("BrinkmanAdvDiffBcHelper::registerInhomogeneousBC\n"
                   << "  unsupported Indicator Function Type: "
                   << enum_to_string<IndicatorFunctionType>(indicator_func_type_enum) << " \n"
                   << "  valid choices are: SMOOTH and DISCONTINUOUS\n");
    }

    // Store BC properties in a struct.
    BCProperties bc_prop;
    bc_prop.ls_solid_var = ls_solid_var;
    bc_prop.bc_type = bc_type_enum;
    bc_prop.indicator_func_type = indicator_func_type_enum;
    bc_prop.bc_val = bc_val;
    bc_prop.num_interface_cells = num_interface_cells;
    bc_prop.eta = eta_penalty_coef;
    bc_prop.callback = callback;
    bc_prop.ctx = ctx;
    d_Q_bc[Q_var].push_back(bc_prop);

    return;
} // registerInhomogeneousBC

void
BrinkmanAdvDiffBcHelper::computeDampingCoefficient(int C_idx, Pointer<CellVariable<NDIM, double> > Q_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_bc.find(Q_var) != d_Q_bc.end());
#endif
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    std::vector<BCProperties> brinkman_zones = d_Q_bc[Q_var];

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // First, deal with Dirichlet BCs, which can be computed pointwise.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const double* patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];

            Pointer<CellData<NDIM, double> > C_data = patch->getPatchData(C_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                // Loop over all the level sets and add up their contributions.
                // TODO: Pointwise check that level sets don't overlap
                double brinkman_coef = 0.0;
                for (const auto& bc_prop : brinkman_zones)
                {
                    // Get the BC specifications for each zone.
                    Pointer<CellVariable<NDIM, double> > ls_solid_var = bc_prop.ls_solid_var;
                    AdvDiffBrinkmanPenalizationBcType bc_type = bc_prop.bc_type;
                    if (bc_type != DIRICHLET) continue;

                    double eta = bc_prop.eta;
                    double num_interface_cells = bc_prop.num_interface_cells;
                    const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));
                    const int phi_idx =
                        var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getNewContext());
                    Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(phi_idx);
                    double phi = (*ls_solid_data)(ci);

                    IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                    double Hphi = 0.0;
                    if (indicator_func_type == SMOOTH)
                        Hphi = IBTK::smooth_heaviside(phi, alpha);
                    else if (indicator_func_type == DISCONTINUOUS)
                        Hphi = IBTK::discontinuous_heaviside(phi);

                    // Note: assumes that chi is positive when phi is negative.
                    const double chi = 1.0 - Hphi;

                    brinkman_coef += (chi / eta);
                }

                // Set the Brinkman damping coefficient.
                (*C_data)(ci) = brinkman_coef;
            }
        }
    }

    // Next, deal with Robin bcs which can be computed hierarchy-wise.
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
    Pointer<HierarchyMathOps> hier_math_ops = d_adv_diff_solver->getHierarchyMathOps();
    for (const auto& bc_prop : brinkman_zones)
    {
        AdvDiffBrinkmanPenalizationBcType bc_type = bc_prop.bc_type;
        if (bc_type != ROBIN) continue;

        // Get the BC specifications for each zone.
        Pointer<CellVariable<NDIM, double> > ls_solid_var = bc_prop.ls_solid_var;
        double num_interface_cells = bc_prop.num_interface_cells;
        const int phi_new_idx = var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getNewContext());
        const int phi_scratch_idx =
            var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getScratchContext());

        // Fill the ghost cells of phi_scratch.
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> phi_transaction_comps(1);
        phi_transaction_comps[0] =
            InterpolationTransactionComponent(phi_scratch_idx,
                                              phi_new_idx,
                                              "CONSERVATIVE_LINEAR_REFINE",
                                              false,
                                              "CONSERVATIVE_COARSEN",
                                              "LINEAR",
                                              false,
                                              d_adv_diff_solver->getPhysicalBcCoefs(ls_solid_var));
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(phi_transaction_comps, patch_hierarchy);
        hier_bdry_fill->fillData(d_new_time);

        // Computing the normal at the edges (resp. faces) for 2D (resp. 3D) cases.
        hier_math_ops->grad(d_n_scratch_idx, d_n_var, true, -1.0, phi_scratch_idx, ls_solid_var, nullptr, d_new_time);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const Box<NDIM>& patch_box = patch->getBox();
                const double* patch_dx = patch_geom->getDx();
                double vol_cell = 1.0;
                for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(phi_scratch_idx);
                Pointer<SideData<NDIM, double> > n_data = patch->getPatchData(d_n_scratch_idx);
                Pointer<SideData<NDIM, double> > n_chi_data = patch->getPatchData(d_n_chi_scratch_idx);
                Pointer<CellData<NDIM, double> > chi_data = patch->getPatchData(d_chi_scratch_idx);

                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);
                        const double phi_lower = (*ls_solid_data)(s_i.toCell(0));
                        const double phi_upper = (*ls_solid_data)(s_i.toCell(1));
                        const double phi = 0.5 * (phi_lower + phi_upper);

                        IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                        double Hphi = 0.0;
                        if (indicator_func_type == SMOOTH)
                            Hphi = IBTK::smooth_heaviside(phi, alpha);
                        else if (indicator_func_type == DISCONTINUOUS)
                            Hphi = IBTK::discontinuous_heaviside(phi);

                        // Note: assumes that chi is positive when phi is negative.
                        const double chi = 1.0 - Hphi;

                        // Set n_chi = n*chi.
                        (*n_chi_data)(s_i) = chi * (*n_data)(s_i);
                    }
                }
                // Compute persistent cell patch data for chi.
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    double phi = (*ls_solid_data)(ci);

                    IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                    double Hphi = 0.0;
                    if (indicator_func_type == SMOOTH)
                        Hphi = IBTK::smooth_heaviside(phi, alpha);
                    else if (indicator_func_type == DISCONTINUOUS)
                        Hphi = IBTK::discontinuous_heaviside(phi);

                    (*chi_data)(ci) = (1.0 - Hphi);
                }
            }
        }
        // Compute C += div(chi * n) - chi * div(n) throughout the hierarchy.
        hier_math_ops->div(
            d_div_n_chi_scratch_idx, d_div_var, 1.0, d_n_chi_scratch_idx, d_n_var, NULL, d_current_time, false);
        hier_math_ops->div(d_div_n_scratch_idx, d_div_var, 1.0, d_n_scratch_idx, d_n_var, NULL, d_current_time, false);
        hier_cc_data_ops.multiply(d_div_n_scratch_idx, d_chi_scratch_idx, d_div_n_scratch_idx);
        // Ensure that we don't overwrite contributions from other BCs.
        hier_cc_data_ops.add(C_idx, C_idx, d_div_n_chi_scratch_idx);
        hier_cc_data_ops.subtract(C_idx, C_idx, d_div_n_scratch_idx);
#if !defined(NDEBUG)
        // Nullify patch data just to be safe.
        hier_sc_data_ops.setToScalar(d_n_scratch_idx, std::numeric_limits<double>::quiet_NaN());
        hier_sc_data_ops.setToScalar(d_n_chi_scratch_idx, std::numeric_limits<double>::quiet_NaN());
        hier_cc_data_ops.setToScalar(d_div_n_chi_scratch_idx, std::numeric_limits<double>::quiet_NaN());
        hier_cc_data_ops.setToScalar(d_div_n_scratch_idx, std::numeric_limits<double>::quiet_NaN());
        hier_cc_data_ops.setToScalar(d_chi_scratch_idx, std::numeric_limits<double>::quiet_NaN());
#endif
    }
    return;
} // computeDampingCoefficient

void
BrinkmanAdvDiffBcHelper::computeDiffusionCoefficient(int D_idx,
                                                     Pointer<CellVariable<NDIM, double> > Q_var,
                                                     int kappa_idx,
                                                     double kappa)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_bc.find(Q_var) != d_Q_bc.end());
#endif
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    std::vector<BCProperties> brinkman_zones = d_Q_bc[Q_var];
    const bool variable_kappa = (kappa_idx != IBTK::invalid_index);
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();

    // Fill the ghost cells for each solid level set, which is required to obtain appropriate side-centered data.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> phi_transaction_comps(brinkman_zones.size());
    for (std::size_t i = 0; i < brinkman_zones.size(); ++i)
    {
        const auto& bc_prop = brinkman_zones[i];
        Pointer<CellVariable<NDIM, double> > ls_solid_var = bc_prop.ls_solid_var;
        const int phi_idx = var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getNewContext());
        const int phi_scratch_idx =
            var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getScratchContext());
        phi_transaction_comps[i] =
            InterpolationTransactionComponent(phi_scratch_idx,
                                              phi_idx,
                                              "CONSERVATIVE_LINEAR_REFINE",
                                              false,
                                              "CONSERVATIVE_COARSEN",
                                              "LINEAR",
                                              false,
                                              d_adv_diff_solver->getPhysicalBcCoefs(ls_solid_var));
    }
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(phi_transaction_comps, patch_hierarchy);
    hier_bdry_fill->fillData(d_new_time);

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const double* patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];

            Pointer<SideData<NDIM, double> > D_data = patch->getPatchData(D_idx);
            Pointer<SideData<NDIM, double> > kappa_data = variable_kappa ? patch->getPatchData(kappa_idx) : nullptr;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                    // Loop over all the level sets and add up their contributions.
                    // TODO: Pointwise check that level sets don't overlap
                    const double kp = variable_kappa ? (*kappa_data)(s_i) : kappa;
                    double brinkman_coef = kp;
                    for (const auto& bc_prop : brinkman_zones)
                    {
                        // Get the BC specifications for each zone.
                        Pointer<CellVariable<NDIM, double> > ls_solid_var = bc_prop.ls_solid_var;
                        AdvDiffBrinkmanPenalizationBcType bc_type = bc_prop.bc_type;
                        double eta = bc_prop.eta;
                        double num_interface_cells = bc_prop.num_interface_cells;
                        const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                        // Ghost cells for scratch data filled above.
                        const int phi_idx =
                            var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getScratchContext());
                        Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(phi_idx);
                        const double phi_lower = (*ls_solid_data)(s_i.toCell(0));
                        const double phi_upper = (*ls_solid_data)(s_i.toCell(1));
                        const double phi = 0.5 * (phi_lower + phi_upper);

                        IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                        double Hphi = 0.0;
                        if (indicator_func_type == SMOOTH)
                            Hphi = IBTK::smooth_heaviside(phi, alpha);
                        else if (indicator_func_type == DISCONTINUOUS)
                            Hphi = IBTK::discontinuous_heaviside(phi);

                        // Note: assumes that chi is positive when phi is negative.
                        const double chi = 1.0 - Hphi;

                        // Dirichlet BCs have no contribution to the diffusion coefficient.
                        if (bc_type == NEUMANN || bc_type == ROBIN)
                        {
                            brinkman_coef += (-kp * chi) + (eta * chi);
                        }
                    }

                    // Set the Brinkman diffusion coefficient.
                    (*D_data)(s_i) = brinkman_coef;
                }
            }
        }
    }
    return;
} // computeDiffusionCoefficient

void
BrinkmanAdvDiffBcHelper::computeForcing(int F_idx, Pointer<CellVariable<NDIM, double> > Q_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_bc.find(Q_var) != d_Q_bc.end());
#endif
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    std::vector<BCProperties> brinkman_zones = d_Q_bc[Q_var];

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
    hier_cc_data_ops.setToScalar(F_idx, 0.0);

    // First, deal with Dirichlet, spatially constant g, which can be computed pointwise.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const double* patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];

            Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                // Loop over all the level sets and add up their contributions.
                // TODO: Pointwise check that level sets don't overlap
                double brinkman_forcing = 0.0;
                for (const auto& bc_prop : brinkman_zones)
                {
                    // Get the BC specifications for each zone.
                    Pointer<CellVariable<NDIM, double> > ls_solid_var = bc_prop.ls_solid_var;
                    AdvDiffBrinkmanPenalizationBcType bc_type = bc_prop.bc_type;
                    double bc_val = bc_prop.bc_val;
                    double eta = bc_prop.eta;

                    bool requires_callback = bc_prop.callback != nullptr;
                    double num_interface_cells = bc_prop.num_interface_cells;
                    const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                    if (requires_callback) continue;
#if !defined(NDEBUG)
                    TBOX_ASSERT(bc_type == DIRICHLET);
                    TBOX_ASSERT(!requires_callback);
#endif
                    const int phi_idx =
                        var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getNewContext());
                    Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(phi_idx);
                    double phi = (*ls_solid_data)(ci);

                    IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                    double Hphi = 0.0;
                    if (indicator_func_type == SMOOTH)
                        Hphi = IBTK::smooth_heaviside(phi, alpha);
                    else if (indicator_func_type == DISCONTINUOUS)
                        Hphi = IBTK::discontinuous_heaviside(phi);

                    // Note: assumes that chi is positive when phi is negative.
                    const double chi = (1.0 - Hphi);

                    // Homogeneous Neumann BCs have no contribution to the forcing term.
                    if (bc_type == DIRICHLET)
                    {
                        brinkman_forcing += (chi / eta) * bc_val;
                    }
                }
                // Set the Brinkman body force.
                (*F_data)(ci) += brinkman_forcing;
            }
        }
    }

    // Next, deal with inhomogeneous Neumann/Robin BCs and Dirichlet with spatially varying g, which need
    // to be computed 'hierarchy-wise'.
    for (const auto& bc_prop : brinkman_zones)
    {
        Pointer<CellVariable<NDIM, double> > ls_solid_var = bc_prop.ls_solid_var;
        AdvDiffBrinkmanPenalizationBcType bc_type = bc_prop.bc_type;
        bool requires_callback = (bc_prop.callback != nullptr);
        if (!requires_callback) continue;

        // Fill in patch data from user-prescribed callback.
        BrinkmanInhomogeneousBCsFcnPtr bp_fcn = bc_prop.callback;
        void* bp_ctx = bc_prop.ctx;

#if !defined(NDEBUG)
        TBOX_ASSERT(bp_fcn);
#endif
        Pointer<HierarchyMathOps> hier_math_ops = d_adv_diff_solver->getHierarchyMathOps();
        if (bc_type == DIRICHLET)
        {
            bp_fcn(d_variable_g_scratch_idx, ls_solid_var, hier_math_ops, d_current_time, bp_ctx);

            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                    const Box<NDIM>& patch_box = patch->getBox();
                    const double* patch_dx = patch_geom->getDx();
                    double vol_cell = 1.0;
                    for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                    double num_interface_cells = bc_prop.num_interface_cells;
                    double eta = bc_prop.eta;
                    const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                    Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_idx);
                    Pointer<CellData<NDIM, double> > g_data = patch->getPatchData(d_variable_g_scratch_idx);

                    for (Box<NDIM>::Iterator it(patch_box); it; it++)
                    {
                        CellIndex<NDIM> ci(it());

                        const int phi_idx =
                            var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getNewContext());
                        Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(phi_idx);
                        double phi = (*ls_solid_data)(ci);

                        IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                        double Hphi = 0.0;
                        if (indicator_func_type == SMOOTH)
                            Hphi = IBTK::smooth_heaviside(phi, alpha);
                        else if (indicator_func_type == DISCONTINUOUS)
                            Hphi = IBTK::discontinuous_heaviside(phi);

                        // Note: assumes that chi is positive when phi is negative
                        const double chi = (1.0 - Hphi);

                        (*F_data)(ci) += (chi / eta) * (*g_data)(ci);
                    }
                }
            }
        }
        else
        {
            bp_fcn(d_B_scratch_idx, ls_solid_var, hier_math_ops, d_current_time, bp_ctx);

            // Compute chi*B throughout the hierarchy
            const int phi_idx = var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getNewContext());
            const int phi_scratch_idx =
                var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getScratchContext());

            // Fill the ghost cells of phi_scratch
            using InterpolationTransactionComponent =
                HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<InterpolationTransactionComponent> phi_transaction_comps(1);
            phi_transaction_comps[0] =
                InterpolationTransactionComponent(phi_scratch_idx,
                                                  phi_idx,
                                                  "CONSERVATIVE_LINEAR_REFINE",
                                                  false,
                                                  "CONSERVATIVE_COARSEN",
                                                  "LINEAR",
                                                  false,
                                                  d_adv_diff_solver->getPhysicalBcCoefs(ls_solid_var));
            Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
            hier_bdry_fill->initializeOperatorState(phi_transaction_comps, patch_hierarchy);
            hier_bdry_fill->fillData(d_new_time);

            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                    const Box<NDIM>& patch_box = patch->getBox();
                    const double* patch_dx = patch_geom->getDx();
                    double vol_cell = 1.0;
                    for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                    double num_interface_cells = bc_prop.num_interface_cells;
                    const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                    Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(phi_scratch_idx);
                    Pointer<SideData<NDIM, double> > B_data = patch->getPatchData(d_B_scratch_idx);
                    Pointer<SideData<NDIM, double> > B_chi_data = patch->getPatchData(d_B_chi_scratch_idx);
                    Pointer<CellData<NDIM, double> > chi_data = patch->getPatchData(d_chi_scratch_idx);

                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                        {
                            SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);
                            const double phi_lower = (*ls_solid_data)(s_i.toCell(0));
                            const double phi_upper = (*ls_solid_data)(s_i.toCell(1));
                            const double phi = 0.5 * (phi_lower + phi_upper);

                            IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                            double Hphi = 0.0;
                            if (indicator_func_type == SMOOTH)
                                Hphi = IBTK::smooth_heaviside(phi, alpha);
                            else if (indicator_func_type == DISCONTINUOUS)
                                Hphi = IBTK::discontinuous_heaviside(phi);

                            // Note: assumes that chi is positive when phi is negative
                            const double chi = 1.0 - Hphi;

                            // Set B_chi = chi*B
                            (*B_chi_data)(s_i) = chi * (*B_data)(s_i);
                        }
                    }
                    // Compute persistent cell patch data for chi
                    for (Box<NDIM>::Iterator it(patch_box); it; it++)
                    {
                        CellIndex<NDIM> ci(it());
                        double phi = (*ls_solid_data)(ci);

                        IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                        double Hphi = 0.0;
                        if (indicator_func_type == SMOOTH)
                            Hphi = IBTK::smooth_heaviside(phi, alpha);
                        else if (indicator_func_type == DISCONTINUOUS)
                            Hphi = IBTK::discontinuous_heaviside(phi);

                        (*chi_data)(ci) = (1.0 - Hphi);
                    }
                }
            }
            // Compute F += div(chi * B) - chi * div(B) throughout the hierarchy.
            HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            hier_math_ops->div(
                d_div_B_chi_scratch_idx, d_div_var, 1.0, d_B_chi_scratch_idx, d_B_var, NULL, d_current_time, false);
            hier_math_ops->div(
                d_div_B_scratch_idx, d_div_var, 1.0, d_B_scratch_idx, d_B_var, NULL, d_current_time, false);
            hier_cc_data_ops.multiply(d_div_B_scratch_idx, d_chi_scratch_idx, d_div_B_scratch_idx);
            // Ensure that we don't overwrite contributions from other BCs.
            hier_cc_data_ops.add(F_idx, F_idx, d_div_B_chi_scratch_idx);
            hier_cc_data_ops.subtract(F_idx, F_idx, d_div_B_scratch_idx);

#if !defined(NDEBUG)
            // Nullify patch data just to be safe
            hier_sc_data_ops.setToScalar(d_B_scratch_idx, std::numeric_limits<double>::quiet_NaN());
            hier_sc_data_ops.setToScalar(d_B_chi_scratch_idx, std::numeric_limits<double>::quiet_NaN());
            hier_cc_data_ops.setToScalar(d_div_B_chi_scratch_idx, std::numeric_limits<double>::quiet_NaN());
            hier_cc_data_ops.setToScalar(d_div_B_scratch_idx, std::numeric_limits<double>::quiet_NaN());
            hier_cc_data_ops.setToScalar(d_chi_scratch_idx, std::numeric_limits<double>::quiet_NaN());
#endif
        }
    }
    return;
} // computeForcing

void
BrinkmanAdvDiffBcHelper::maskForcingTerm(int N_idx,
                                         Pointer<CellVariable<NDIM, double> > Q_var,
                                         const bool mask_smeared_region)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_Q_bc.find(Q_var) != d_Q_bc.end());
#endif
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    std::vector<BCProperties> brinkman_zones = d_Q_bc[Q_var];

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const double* patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];

            Pointer<CellData<NDIM, double> > N_data = patch->getPatchData(N_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                // Loop over all the level sets and add up their contributions to the masking term
                // TODO: Pointwise check that level sets don't overlap
                double brinkman_mask = 1.0;
                for (const auto& bc_prop : brinkman_zones)
                {
                    // Get the BC specifications for each zone
                    Pointer<CellVariable<NDIM, double> > ls_solid_var = bc_prop.ls_solid_var;
                    AdvDiffBrinkmanPenalizationBcType bc_type = bc_prop.bc_type;
                    double num_interface_cells = bc_prop.num_interface_cells;
                    const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));
                    const int phi_idx =
                        var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_solver->getNewContext());
                    Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(phi_idx);
                    double phi = (*ls_solid_data)(ci);

                    IndicatorFunctionType indicator_func_type = bc_prop.indicator_func_type;
                    double Hphi = 0.0;
                    if (indicator_func_type == SMOOTH)
                        Hphi = IBTK::smooth_heaviside(phi, alpha);
                    else if (indicator_func_type == DISCONTINUOUS)
                        Hphi = IBTK::discontinuous_heaviside(phi);

                    // Note: assumes that chi is positive when phi is negative
                    double chi = 1.0 - Hphi;
                    if (mask_smeared_region && chi > 0.0) chi = 1.0;

                    // Dirichlet BCs have no contribution to the masking term
                    if (bc_type == NEUMANN || bc_type == ROBIN)
                    {
                        brinkman_mask -= chi;
                    }
                }

                // Mask the forcing term
                (*N_data)(ci) = brinkman_mask * (*N_data)(ci);
            }
        }
    }
    return;
} // maskForcingTerm
/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
