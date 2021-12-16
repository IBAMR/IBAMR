// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
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
#include "ibamr/IEPSemiImplicitHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/LaserBeamForceFunction.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"

#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "IBAMR_config.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <string>
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
////////////////////////////// PUBLIC ///////////////////////////////////////

LaserBeamForceFunction::LaserBeamForceFunction(const std::string& object_name,
                                               Pointer<Database> input_db,
                                               // AdvDiffHierarchyIntegrator* adv_diff_solver,
                                               // Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                               Pointer<IEPSemiImplicitHierarchyIntegrator> iep_solver,
                                               // IEPSemiImplicitHierarchyIntegrator* iep_solver,
                                               Pointer<Variable<NDIM> > H_var,
                                               const double rho_liquid,
                                               const double rho_gas,
                                               const double cp_liquid,
                                               const double cp_gas)
    : CartGridFunction(object_name),
      d_iep_solver(iep_solver),
      d_H_var(H_var),
      d_rho_liquid(rho_liquid),
      d_rho_gas(rho_gas),
      d_cp_liquid(cp_liquid),
      d_cp_gas(cp_gas)
{
    // Set some default values
    d_ts_type = MIDPOINT_RULE;

    if (input_db)
    {
        if (input_db->keyExists("time_stepping_type"))
        {
            d_ts_type = string_to_enum<TimeSteppingType>(input_db->getString("time_stepping_type"));
        }
    }
    return;
} // LaserBeamForceFunction

bool
LaserBeamForceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LaserBeamForceFunction::setDataOnPatchHierarchy(const int data_idx,
                                                Pointer<Variable<NDIM> > var,
                                                Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                const double data_time,
                                                const bool initial_time,
                                                const int coarsest_ln_in,
                                                const int finest_ln_in)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif

    // Get the newest patch data index for the level set variable
    Pointer<CellVariable<NDIM, double> > H_cc_var = d_H_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(!H_cc_var.isNull());
#endif
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int H_new_idx = var_db->mapVariableAndContextToIndex(H_cc_var, d_iep_solver->getNewContext());
    const int H_current_idx = var_db->mapVariableAndContextToIndex(H_cc_var, d_iep_solver->getCurrentContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(H_new_idx >= 0);
    TBOX_ASSERT(H_current_idx >= 0);
#endif

    IntVector<NDIM> cell_ghosts = 1;
    IntVector<NDIM> no_ghosts = 0;

    d_H_scratch_idx =
        var_db->registerVariableAndContext(H_cc_var, var_db->getContext(d_object_name + "::H"), cell_ghosts);
    d_grad_H_var = new CellVariable<NDIM, double>(d_object_name + "::grad_H_var", NDIM);
    d_grad_H_scratch_idx =
        var_db->registerVariableAndContext(d_grad_H_var, var_db->getContext(d_object_name + "::grad_H"), no_ghosts);

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_H_scratch_idx, data_time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_grad_H_scratch_idx, data_time);
    }

    // Copy level set into phi and C.
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_ts_type == MIDPOINT_RULE)
    {
        hier_cc_data_ops.linearSum(d_H_scratch_idx,
                                   0.5,
                                   H_new_idx,
                                   0.5,
                                   H_current_idx,
                                   /*interior_only*/ true);
    }
    else if (d_ts_type == BACKWARD_EULER)
    {
        hier_cc_data_ops.copyData(d_H_scratch_idx,
                                  H_new_idx,
                                  /*interior_only*/ true);
    }
    else if (d_ts_type == FORWARD_EULER)
    {
        hier_cc_data_ops.copyData(d_H_scratch_idx,
                                  H_current_idx,
                                  /*interior_only*/ true);
    }
    else
    {
        TBOX_ERROR("LaserBeamForceFunction::setDataOnPatchHierarchy : "
                   << "The class only supports BACKWARD_EULER, FORWARD_EULER, and "
                      "MIDPOINT_RULE"
                   << std::endl);
    }

    // Fill ghost cells
    RobinBcCoefStrategy<NDIM>* H_bc_coef = (d_iep_solver->getPhysicalBcCoefs(H_cc_var)).front();
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent H_transaction(
        d_H_scratch_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, H_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> H_fill_op = new HierarchyGhostCellInterpolation();
    H_fill_op->initializeOperatorState(H_transaction, hierarchy, coarsest_ln, finest_ln);

    H_fill_op->fillData(data_time);

    // Find |grad H|
    HierarchyMathOps* hier_math_ops = new HierarchyMathOps("HierarchyMathOps", hierarchy);
    hier_math_ops->grad(d_grad_H_scratch_idx, d_grad_H_var, 1.0, d_H_scratch_idx, H_cc_var, nullptr, data_time);

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    // compute the heat influx.
    const double apply_time = data_time;
    const double current_time = data_time;
    const double new_time = data_time;

    for (unsigned k = 0; k < d_heat_influx.size(); ++k)
    {
        d_heat_influx[k](
            data_idx, hier_math_ops, -1 /*cycle_num*/, apply_time, current_time, new_time, d_heat_influx_ctx[k]);
    }

    // Deallocate and remove scratch variables.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_H_scratch_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_grad_H_scratch_idx);
    }
    var_db->removePatchDataIndex(d_H_scratch_idx);
    var_db->removePatchDataIndex(d_grad_H_scratch_idx);

    return;
} // setDataOnPatchHierarchy

void
LaserBeamForceFunction::setDataOnPatch(const int data_idx,
                                       Pointer<Variable<NDIM> > /*var*/,
                                       Pointer<Patch<NDIM> > patch,
                                       const double data_time,
                                       const bool initial_time,
                                       Pointer<PatchLevel<NDIM> > level)
{
    Pointer<PatchData<NDIM> > f_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(f_data);
#endif

    const int rho_cc_idx = d_iep_solver->getUpdatedDensityIndex();
    const int cp_cc_idx = d_iep_solver->getUpdatedSpecificHeatIndex();

    const Box<NDIM>& patch_box = patch->getBox();

    Pointer<CellData<NDIM, double> > grad_H_data = patch->getPatchData(d_grad_H_scratch_idx);
    Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_cc_idx);
    Pointer<CellData<NDIM, double> > cp_data = patch->getPatchData(cp_cc_idx);
    Pointer<CellData<NDIM, double> > f_cc_data = f_data;

    for (Box<NDIM>::Iterator it(patch_box); it; it++)
    {
        CellIndex<NDIM> ci(it());

        const double grad_H_mag =
            std::sqrt(std::pow((*grad_H_data)(ci, 0), 2.0) + std::pow((*grad_H_data)(ci, 1), 2.0));
        const double multiplier =
            (*rho_data)(ci) * (*cp_data)(ci) / (d_rho_liquid * d_cp_liquid + d_rho_gas * d_cp_gas);
        IBTK::Vector coord = IBTK::Vector::Zero();

        (*f_cc_data)(ci) = grad_H_mag * multiplier;
    }

    return;
} // setDataOnPatch

void
LaserBeamForceFunction::registerHeatInflux(HeatInfluxPtr callback, void* ctx)
{
    d_heat_influx.push_back(callback);
    d_heat_influx_ctx.push_back(ctx);

    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
