// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/AdvDiffStochasticForcing.h"
#include "ibamr/RNG.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/samrai_compatibility_names.h"

#include "MultiblockDataTranslator.h"
#include "SAMRAIArray.h"
#include "SAMRAIArrayData.h"
#include "SAMRAIBoundaryBox.h"
#include "SAMRAIBox.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellDataFactory.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIHierarchyDataOpsManager.h"
#include "SAMRAIHierarchyDataOpsReal.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAISideData.h"
#include "SAMRAISideGeometry.h"
#include "SAMRAISideIndex.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableContext.h"
#include "SAMRAIVariableDatabase.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
void
genrandn(SAMRAIArrayData<double>& data, const SAMRAIBox& box)
{
    for (int depth = 0; depth < data.getDepth(); ++depth)
    {
        for (SAMRAIBox::Iterator i(box); i; i++)
        {
            RNG::genrandn(&data(i(), depth));
        }
    }
    return;
} // genrandn
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffStochasticForcing::AdvDiffStochasticForcing(std::string object_name,
                                                   SAMRAIPointer<SAMRAIDatabase> input_db,
                                                   SAMRAIPointer<SAMRAICellVariable<double>> C_var,
                                                   const AdvDiffSemiImplicitHierarchyIntegrator* const adv_diff_solver)
    : d_object_name(std::move(object_name)), d_C_var(C_var), d_adv_diff_solver(adv_diff_solver)
{
    std::string f_expression = "1.0";
    if (input_db)
    {
        if (input_db->keyExists("std")) d_std = input_db->getDouble("std");
        if (input_db->keyExists("num_rand_vals")) d_num_rand_vals = input_db->getInteger("num_rand_vals");
        int k = 0;
        std::string key_name = "weights_0";
        while (input_db->keyExists(key_name))
        {
            d_weights.push_back(input_db->getDoubleArray(key_name));
#if !defined(NDEBUG)
            TBOX_ASSERT(d_weights.back().size() == d_num_rand_vals);
#endif
            ++k;
            key_name = "weights_" + std::to_string(k);
        }
        if (input_db->keyExists("dirichlet_bc_scaling"))
            d_dirichlet_bc_scaling = input_db->getDouble("dirichlet_bc_scaling");
        if (input_db->keyExists("neumann_bc_scaling")) d_neumann_bc_scaling = input_db->getDouble("neumann_bc_scaling");
        if (input_db->keyExists("f_expression")) f_expression = input_db->getString("f_expression");
    }
    d_f_parser.SetExpr(f_expression);

    // Determine the number of components that need to be allocated.
    SAMRAIPointer<SAMRAICellDataFactory<double>> C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();

    // Setup variables and variable context objects.
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    d_C_cc_var = new SAMRAICellVariable<double>(d_object_name + "::C_cc", C_depth);
    static const SAMRAIIntVector ghosts_cc = 1;
    d_C_current_cc_idx = var_db->registerVariableAndContext(d_C_cc_var, d_context, ghosts_cc);
    d_C_half_cc_idx = var_db->registerClonedPatchDataIndex(d_C_cc_var, d_C_current_cc_idx);
    d_C_new_cc_idx = var_db->registerClonedPatchDataIndex(d_C_cc_var, d_C_current_cc_idx);
    d_F_sc_var = new SAMRAISideVariable<double>(d_object_name + "::F_sc", C_depth);
    static const SAMRAIIntVector ghosts_sc = 0;
    d_F_sc_idx = var_db->registerVariableAndContext(d_F_sc_var, d_context, ghosts_sc);
    for (int k = 0; k < d_num_rand_vals; ++k)
        d_F_sc_idxs.push_back(var_db->registerClonedPatchDataIndex(d_F_sc_var, d_F_sc_idx));
    return;
} // AdvDiffStochasticForcing

bool
AdvDiffStochasticForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
AdvDiffStochasticForcing::setDataOnPatchHierarchy(const int data_idx,
                                                  SAMRAIPointer<SAMRAIVariable> var,
                                                  SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                                                  const double data_time,
                                                  const bool initial_time,
                                                  const int coarsest_ln_in,
                                                  const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == IBTK::invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    const int cycle_num = d_adv_diff_solver->getCurrentCycleNumber();
    if (!initial_time)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(cycle_num >= 0);
#endif
        // Allocate data to store components of the stochastic stress components.
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(level_num);
            if (!level->checkAllocated(d_C_current_cc_idx)) level->allocatePatchData(d_C_current_cc_idx);
            if (!level->checkAllocated(d_C_half_cc_idx)) level->allocatePatchData(d_C_half_cc_idx);
            if (!level->checkAllocated(d_C_new_cc_idx)) level->allocatePatchData(d_C_new_cc_idx);
            if (!level->checkAllocated(d_F_sc_idx)) level->allocatePatchData(d_F_sc_idx);
            for (int k = 0; k < d_num_rand_vals; ++k)
                if (!level->checkAllocated(d_F_sc_idxs[k])) level->allocatePatchData(d_F_sc_idxs[k]);
        }

        // Set concentration value used to compute concentration-dependent flux
        // scaling.
        const double dt = d_adv_diff_solver->getCurrentTimeStepSize();
        const double current_time = d_adv_diff_solver->getIntegratorTime();
        const double half_time = current_time + 0.5 * dt;
        const double new_time = current_time + dt;
        SAMRAIHierarchyDataOpsManager* hier_data_ops_manager = SAMRAIHierarchyDataOpsManager::getManager();
        SAMRAIPointer<SAMRAIHierarchyDataOpsReal<double>> hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_C_cc_var,
                                                       hierarchy,
                                                       /*get_unique*/ true);
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        const int C_current_idx = var_db->mapVariableAndContextToIndex(d_C_var, d_adv_diff_solver->getCurrentContext());
        const int C_new_idx = var_db->mapVariableAndContextToIndex(d_C_var, d_adv_diff_solver->getNewContext());
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_fill_components(1);
        HierarchyGhostCellInterpolation ghost_fill_op;
        const std::vector<SAMRAIRobinBcCoefStrategy*>& C_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_C_var);
        const TimeSteppingType convective_time_stepping_type =
            d_adv_diff_solver->getConvectiveTimeSteppingType(d_C_var);
        switch (convective_time_stepping_type)
        {
        case FORWARD_EULER:
            hier_cc_data_ops->copyData(d_C_current_cc_idx, C_current_idx);
            ghost_fill_components[0] =
                InterpolationTransactionComponent(d_C_current_cc_idx, "NONE", false, "NONE", "NONE", false, C_bc_coef);
            ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
            ghost_fill_op.fillData(current_time);
            break;
        case MIDPOINT_RULE:
            hier_cc_data_ops->linearSum(d_C_half_cc_idx, 0.5, C_current_idx, 0.5, C_new_idx);
            ghost_fill_components[0] =
                InterpolationTransactionComponent(d_C_half_cc_idx, "NONE", false, "NONE", "NONE", false, C_bc_coef);
            ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
            ghost_fill_op.fillData(half_time);
            break;
        case TRAPEZOIDAL_RULE:
            if (cycle_num == 0)
            {
                hier_cc_data_ops->copyData(d_C_current_cc_idx, C_current_idx);
                ghost_fill_components[0] = InterpolationTransactionComponent(
                    d_C_current_cc_idx, "NONE", false, "NONE", "NONE", false, C_bc_coef);
                ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
                ghost_fill_op.fillData(current_time);
            }
            else
            {
                hier_cc_data_ops->copyData(d_C_new_cc_idx, C_new_idx);
                ghost_fill_components[0] =
                    InterpolationTransactionComponent(d_C_new_cc_idx, "NONE", false, "NONE", "NONE", false, C_bc_coef);
                ghost_fill_op.initializeOperatorState(ghost_fill_components, hierarchy);
                ghost_fill_op.fillData(new_time);
            }
            break;
        default:
            TBOX_ERROR(d_object_name << "::setDataOnPatchHierarchy():\n"
                                     << "  unsupported default convective time stepping type: "
                                     << enum_to_string<TimeSteppingType>(convective_time_stepping_type) << " \n"
                                     << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE, "
                                        "TRAPEZOIDAL_RULE\n");
        }

        // Generate random components.
        if (cycle_num == 0)
        {
            for (int k = 0; k < d_num_rand_vals; ++k)
            {
                for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
                {
                    SAMRAIPointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(level_num);
                    for (SAMRAIPatchLevel::Iterator p(level); p; p++)
                    {
                        SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                        SAMRAIPointer<SAMRAISideData<double>> F_sc_data = patch->getPatchData(d_F_sc_idxs[k]);
                        for (int d = 0; d < NDIM; ++d)
                        {
                            genrandn(F_sc_data->getArrayData(d), SAMRAISideGeometry::toSideBox(F_sc_data->getBox(), d));
                        }
                    }
                }
            }
        }

// Set random values for the present cycle as weighted combinations of
// the generated random values.
#if !defined(NDEBUG)
        TBOX_ASSERT(cycle_num >= 0 && cycle_num < static_cast<int>(d_weights.size()));
#endif
        const SAMRAIArray<double>& weights = d_weights[cycle_num];
        SAMRAIPointer<SAMRAIHierarchyDataOpsReal<double>> hier_sc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_F_sc_var,
                                                       hierarchy,
                                                       /*get_unique*/ true);
        hier_sc_data_ops->setToScalar(d_F_sc_idx, 0.0);
        for (int k = 0; k < d_num_rand_vals; ++k)
            hier_sc_data_ops->axpy(d_F_sc_idx, weights[k], d_F_sc_idxs[k], d_F_sc_idx);

        // Modify the flux values (if necessary).
        SAMRAIPointer<SAMRAICellDataFactory<double>> C_factory = d_C_var->getPatchDataFactory();
        const int C_depth = C_factory->getDefaultDepth();
        const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs = d_adv_diff_solver->getPhysicalBcCoefs(d_C_var);
        for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(level_num);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAISideData<double>> F_sc_data = patch->getPatchData(d_F_sc_idx);

                const SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch->getPatchGeometry();
                if (!pgeom->getTouchesRegularBoundary()) continue;

                const SAMRAIBox& patch_box = patch->getBox();
                SAMRAIBox side_boxes[NDIM];
                for (int d = 0; d < NDIM; ++d)
                {
                    side_boxes[d] = SAMRAISideGeometry::toSideBox(patch_box, d);
                }
                const SAMRAIArray<SAMRAIBoundaryBox> physical_codim1_boxes =
                    PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const SAMRAIBoundaryBox& bdry_box = physical_codim1_boxes[n];
                    const SAMRAIIntVector gcw_to_fill = 1;
                    const SAMRAIBox bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const int location_index = bdry_box.getLocationIndex();
                    const int bdry_normal_axis = location_index / 2;
                    const SAMRAIBoundaryBox trimmed_bdry_box(
                        bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    const SAMRAIBox bc_coef_box =
                        PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                    SAMRAIPointer<SAMRAIArrayData<double>> acoef_data = new SAMRAIArrayData<double>(bc_coef_box, 1);
                    SAMRAIPointer<SAMRAIArrayData<double>> bcoef_data = new SAMRAIArrayData<double>(bc_coef_box, 1);
                    SAMRAIPointer<SAMRAIArrayData<double>> gcoef_data = new SAMRAIArrayData<double>(bc_coef_box, 1);

                    // Set the boundary condition coefficients and use them to
                    // rescale the stochastic fluxes.
                    for (int d = 0; d < C_depth; ++d)
                    {
                        SAMRAIRobinBcCoefStrategy* bc_coef = bc_coefs[d];
                        bc_coef->setBcCoefs(
                            acoef_data, bcoef_data, gcoef_data, var, *patch, trimmed_bdry_box, data_time);
                        for (SAMRAIBox::Iterator it(bc_coef_box * side_boxes[bdry_normal_axis]); it; it++)
                        {
                            const SAMRAIIndex& i = it();
                            const double& alpha = (*acoef_data)(i, 0);
                            const double& beta = (*bcoef_data)(i, 0);
                            const bool dirichlet_bc = (alpha != 0.0 && beta == 0.0);

                            SAMRAISideIndex s_i(i, bdry_normal_axis, 0);
                            if (dirichlet_bc)
                            {
                                (*F_sc_data)(s_i, d) *= d_dirichlet_bc_scaling;
                            }
                            else
                            {
                                (*F_sc_data)(s_i, d) = d_neumann_bc_scaling;
                            }
                        }
                    }
                }
            }
        }

        // Synchronize side-centered values.
        using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
        SynchronizationTransactionComponent synch_component(d_F_sc_idx);
        SideDataSynchronization synch_data_op;
        synch_data_op.initializeOperatorState(synch_component, hierarchy);
        synch_data_op.synchronizeData(data_time);
    }

    // Compute div F on each patch level.
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }
    return;
} // setDataOnPatchHierarchy

void
AdvDiffStochasticForcing::setDataOnPatch(const int data_idx,
                                         SAMRAIPointer<SAMRAIVariable> /*var*/,
                                         SAMRAIPointer<SAMRAIPatch> patch,
                                         const double /*data_time*/,
                                         const bool initial_time,
                                         SAMRAIPointer<SAMRAIPatchLevel> /*patch_level*/)
{
    SAMRAIPointer<SAMRAICellData<double>> divF_cc_data = patch->getPatchData(data_idx);
    divF_cc_data->fillAll(0.0);
    if (initial_time) return;
    const SAMRAIBox& patch_box = patch->getBox();
    const SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    double dV = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d) dV *= dx[d];
    const double kappa = d_adv_diff_solver->getDiffusionCoefficient(d_C_var);
    const double dt = d_adv_diff_solver->getCurrentTimeStepSize();
    const double scale = d_std * std::sqrt(2.0 * kappa / (dt * dV));
    double C;
    d_f_parser.DefineVar("c", &C);
    d_f_parser.DefineVar("C", &C);
    SAMRAIPointer<SAMRAICellData<double>> C_current_cc_data = patch->getPatchData(d_C_current_cc_idx);
    SAMRAIPointer<SAMRAICellData<double>> C_half_cc_data = patch->getPatchData(d_C_half_cc_idx);
    SAMRAIPointer<SAMRAICellData<double>> C_new_cc_data = patch->getPatchData(d_C_new_cc_idx);
    SAMRAIPointer<SAMRAISideData<double>> F_sc_data = patch->getPatchData(d_F_sc_idx);
    SAMRAIPointer<SAMRAICellDataFactory<double>> C_factory = d_C_var->getPatchDataFactory();
    const int C_depth = C_factory->getDefaultDepth();
    SAMRAISideData<double> f_scale_sc_data(patch_box, C_depth, SAMRAIIntVector(0));
    const TimeSteppingType convective_time_stepping_type = d_adv_diff_solver->getConvectiveTimeSteppingType(d_C_var);
    const int cycle_num = d_adv_diff_solver->getCurrentCycleNumber();
    for (int d = 0; d < C_depth; ++d)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (BoxIterator<NDIM> i(SAMRAISideGeometry::toSideBox(patch_box, axis)); i; i++)
            {
                const SAMRAIIndex& ic = i();
                SAMRAIIndex ic_lower(ic);
                ic_lower(axis) -= 1;
                SAMRAISideIndex is(ic, axis, SAMRAISideIndex::Lower);
                double f;
                switch (convective_time_stepping_type)
                {
                case FORWARD_EULER:
                {
                    C = 0.5 * ((*C_current_cc_data)(ic, d) + (*C_current_cc_data)(ic_lower, d));
                    f = d_f_parser.Eval();
                    break;
                }
                case MIDPOINT_RULE:
                {
                    C = 0.5 * ((*C_half_cc_data)(ic, d) + (*C_half_cc_data)(ic_lower, d));
                    f = d_f_parser.Eval();
                    break;
                }
                case TRAPEZOIDAL_RULE:
                {
                    if (cycle_num == 0)
                    {
                        C = 0.5 * ((*C_current_cc_data)(ic, d) + (*C_current_cc_data)(ic_lower, d));
                        f = d_f_parser.Eval();
                    }
                    else
                    {
                        C = 0.5 * ((*C_current_cc_data)(ic, d) + (*C_current_cc_data)(ic_lower, d));
                        f = 0.5 * d_f_parser.Eval();
                        C = 0.5 * ((*C_new_cc_data)(ic, d) + (*C_new_cc_data)(ic_lower, d));
                        f += 0.5 * d_f_parser.Eval();
                    }
                    break;
                }
                default:
                    f = std::numeric_limits<double>::quiet_NaN();
                    TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                                             << "  unsupported default convective time stepping type: "
                                             << enum_to_string<TimeSteppingType>(convective_time_stepping_type) << " \n"
                                             << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE, "
                                                "TRAPEZOIDAL_RULE\n");
                }
                f_scale_sc_data(is, d) = std::sqrt(f) * scale;
            }
            for (BoxIterator<NDIM> i(patch_box); i; i++)
            {
                const SAMRAIIndex& ic = i();
                SAMRAISideIndex is_lower(ic, axis, SAMRAISideIndex::Lower);
                SAMRAISideIndex is_upper(ic, axis, SAMRAISideIndex::Upper);
                const double scale_lower = f_scale_sc_data(is_lower, d);
                const double scale_upper = f_scale_sc_data(is_upper, d);
                (*divF_cc_data)(ic, d) +=
                    (scale_upper * (*F_sc_data)(is_upper, d) - scale_lower * (*F_sc_data)(is_lower, d)) / dx[axis];
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
