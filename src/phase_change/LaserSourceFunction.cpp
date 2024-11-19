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

#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/LaserSourceFunction.h"
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

// FORTRAN ROUTINES
#if (NDIM == 2)
#define MOLLIFY_IB_4_FC IBAMR_FC_FUNC_(mollify_ib_4_2d, MOLLIFY_IB_4_2D)
#endif

#if (NDIM == 3)
#define MOLLIFY_IB_4_FC IBAMR_FC_FUNC_(mollify_ib_4_3d, MOLLIFY_IB_4_3D)
#endif

extern "C"
{
    void MOLLIFY_IB_4_FC(double* V,
                         const int& V_gcw,
                         const double* U,
                         const int& U_gcw,
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
////////////////////////////// PUBLIC ///////////////////////////////////////

LaserSourceFunction::LaserSourceFunction(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         Pointer<PhaseChangeHierarchyIntegrator> phase_change_solver,
                                         Pointer<Variable<NDIM> > phi_var)
    : CartGridFunction(object_name), d_phase_change_solver(phase_change_solver), d_phi_var(phi_var)
{
    // Set some default values
    d_ts_type = BACKWARD_EULER;
    d_kernel_fcn = "none";
    d_grad_H_var = new CellVariable<NDIM, double>(d_object_name + "::grad_H_var", NDIM);
    d_num_interface_cells = 2.0;

    if (input_db)
    {
        if (input_db->keyExists("time_stepping_type"))
        {
            d_ts_type = string_to_enum<TimeSteppingType>(input_db->getString("time_stepping_type"));
        }

        d_kernel_fcn = input_db->getStringWithDefault("kernel", d_kernel_fcn);
        d_kernel_fcn = input_db->getStringWithDefault("smoother", d_kernel_fcn);
        d_kernel_fcn = input_db->getStringWithDefault("kernel_fcn", d_kernel_fcn);
        d_kernel_fcn = input_db->getStringWithDefault("smoother_fcn", d_kernel_fcn);

        d_num_interface_cells = input_db->getDoubleWithDefault("num_interface_cells", d_num_interface_cells);
    }
    return;
} // LaserSourceFunction

bool
LaserSourceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
LaserSourceFunction::setDataOnPatchHierarchy(const int data_idx,
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

    // Get the newest patch data index for the Heaviside variable
    Pointer<CellVariable<NDIM, double> > phi_cc_var = d_phi_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(!phi_cc_var.isNull());
#endif
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int phi_new_idx = var_db->mapVariableAndContextToIndex(phi_cc_var, d_phase_change_solver->getNewContext());
    const int phi_scratch_idx =
        var_db->mapVariableAndContextToIndex(phi_cc_var, d_phase_change_solver->getScratchContext());
    const int phi_current_idx =
        var_db->mapVariableAndContextToIndex(phi_cc_var, d_phase_change_solver->getCurrentContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(phi_new_idx >= 0);
    TBOX_ASSERT(phi_current_idx >= 0);
#endif

    IntVector<NDIM> no_ghosts = 0;
    IntVector<NDIM> cell_ghosts = getMinimumGhostWidth(d_kernel_fcn);
    d_H_scratch_idx =
        var_db->registerVariableAndContext(phi_cc_var, var_db->getContext(d_object_name + "::H"), cell_ghosts);
    d_grad_H_scratch_idx =
        var_db->registerVariableAndContext(d_grad_H_var, var_db->getContext(d_object_name + "::grad_H"), no_ghosts);

    const int coarsest_ln = (coarsest_ln_in == IBTK::invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_H_scratch_idx, data_time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_grad_H_scratch_idx, data_time);
    }

    // Based on time stepping copy H index into H_scratch index.
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_ts_type == MIDPOINT_RULE)
    {
        hier_cc_data_ops.linearSum(phi_scratch_idx,
                                   0.5,
                                   phi_new_idx,
                                   0.5,
                                   phi_current_idx,
                                   /*interior_only*/ true);
    }
    else if (d_ts_type == BACKWARD_EULER)
    {
        hier_cc_data_ops.copyData(phi_scratch_idx,
                                  phi_new_idx,
                                  /*interior_only*/ true);
    }
    else if (d_ts_type == FORWARD_EULER)
    {
        hier_cc_data_ops.copyData(phi_scratch_idx,
                                  phi_current_idx,
                                  /*interior_only*/ true);
    }
    else
    {
        TBOX_ERROR("LaserSourceFunction::setDataOnPatchHierarchy : "
                   << "The class only supports BACKWARD_EULER, FORWARD_EULER, and "
                      "MIDPOINT_RULE"
                   << std::endl);
    }

    convertToHeaviside(d_H_scratch_idx, phi_scratch_idx, coarsest_ln, finest_ln, hierarchy);

    // Fill ghost cells
    RobinBcCoefStrategy<NDIM>* H_bc_coef = (d_phase_change_solver->getPhysicalBcCoefs(phi_cc_var)).front();
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent H_transaction(
        d_H_scratch_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, H_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> H_fill_op = new HierarchyGhostCellInterpolation();
    H_fill_op->initializeOperatorState(H_transaction, hierarchy, coarsest_ln, finest_ln);

    H_fill_op->fillData(data_time);

    // Mollify H.
    mollifyData(d_H_scratch_idx, coarsest_ln, finest_ln, data_time, hierarchy, H_fill_op);

    // Find grad H.
    HierarchyMathOps* hier_math_ops = new HierarchyMathOps("HierarchyMathOps", hierarchy, coarsest_ln, finest_ln);
    hier_math_ops->grad(d_grad_H_scratch_idx, d_grad_H_var, 1.0, d_H_scratch_idx, phi_cc_var, nullptr, data_time);

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    // compute the imposed heat flux.
    const double apply_time = data_time;
    const double current_time = data_time;
    const double new_time = data_time;

    if (d_heat_flux)
    {
        (*d_heat_flux)(data_idx, hier_math_ops, -1 /*cycle_num*/, apply_time, current_time, new_time, d_heat_flux_ctx);
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
LaserSourceFunction::setDataOnPatch(const int data_idx,
                                    Pointer<Variable<NDIM> > /*var*/,
                                    Pointer<Patch<NDIM> > patch,
                                    const double /*data_time*/,
                                    const bool initial_time,
                                    Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<CellData<NDIM, double> > f_cc_data = patch->getPatchData(data_idx);

#if !defined(NDEBUG)
    TBOX_ASSERT(f_cc_data);
#endif
    f_cc_data->fillAll(0.0);

    if (initial_time) return;

    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CellData<NDIM, double> > grad_H_data = patch->getPatchData(d_grad_H_scratch_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(grad_H_data);
#endif
    for (Box<NDIM>::Iterator it(patch_box); it; it++)
    {
        CellIndex<NDIM> ci(it());
        double grad_H_mag = 0.0;
#if (NDIM == 2)
        grad_H_mag = std::sqrt(std::pow((*grad_H_data)(ci, 0), 2.0) + std::pow((*grad_H_data)(ci, 1), 2.0));
#elif (NDIM == 3)
        grad_H_mag = std::sqrt(std::pow((*grad_H_data)(ci, 0), 2.0) + std::pow((*grad_H_data)(ci, 1), 2.0) +
                               std::pow((*grad_H_data)(ci, 2), 2.0));
#endif
        (*f_cc_data)(ci) = grad_H_mag;
    }

    return;
} // setDataOnPatch

void
LaserSourceFunction::registerHeatFlux(HeatFluxPtr callback, void* ctx)
{
    d_heat_flux = callback;
    d_heat_flux_ctx = ctx;

    return;
} // registerHeatflux

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////
void
LaserSourceFunction::convertToHeaviside(int H_idx,
                                        int phi_idx,
                                        int coarsest_ln,
                                        int finest_ln,
                                        Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double phi = (*phi_data)(ci);
                (*H_data)(ci) = IBTK::smooth_heaviside(phi, alpha);
            }
        }
    }
    return;
} // convertToHeaviside

void
LaserSourceFunction::mollifyData(int smooth_H_idx,
                                 int coarsest_ln,
                                 int finest_ln,
                                 double data_time,
                                 Pointer<PatchHierarchy<NDIM> > hierarchy,
                                 Pointer<HierarchyGhostCellInterpolation> fill_op)
{
    if (d_kernel_fcn == "none") return;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > smooth_H_data = patch->getPatchData(smooth_H_idx);
            CellData<NDIM, double> H_data(patch_box, /*depth*/ 1, smooth_H_data->getGhostCellWidth());

            H_data.copy(*smooth_H_data);

            double* const V = smooth_H_data->getPointer(0);
            const double* const U = H_data.getPointer(0);
            const int V_gcw = (smooth_H_data->getGhostCellWidth()).max();
            const int U_gcw = (H_data.getGhostCellWidth()).max();

            if (d_kernel_fcn == "IB_4")
            {
                MOLLIFY_IB_4_FC(V,
                                V_gcw,
                                U,
                                U_gcw,
                                patch_box.lower(0),
                                patch_box.upper(0),
                                patch_box.lower(1),
                                patch_box.upper(1)
#if (NDIM == 3)
                                    ,
                                patch_box.lower(2),
                                patch_box.upper(2)
#endif
                );
            }

            else
            {
                TBOX_ERROR("this statement should not be reached");
            }
        }
    }

    if (fill_op) fill_op->fillData(data_time);

    return;
} // mollifyData

int
LaserSourceFunction::getStencilSize(const std::string& kernel_fcn)
{
    if (kernel_fcn == "IB_4") return 4;
    TBOX_ERROR("LaserSourceFunction::getStencilSize()\n"
               << "  Unknown kernel function " << kernel_fcn << std::endl);
    return -1;

} // getStencilSize

int
LaserSourceFunction::getMinimumGhostWidth(const std::string& kernel_fcn)
{
    if (kernel_fcn == "none")
    {
        return 2;
    }
    else
    {
        return std::max(2, static_cast<int>(floor(0.5 * getStencilSize(kernel_fcn))) + 1);
    }
} // getMinimumGhostWidth

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
