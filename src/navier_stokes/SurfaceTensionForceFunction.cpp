// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
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
#include "ibamr/SurfaceTensionForceFunction.h"
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
#define SC_NORMAL_FC IBAMR_FC_FUNC_(sc_normal_2d, SC_NORMAL_2D)
#define CC_CURVATURE_FC IBAMR_FC_FUNC_(cc_curvature_2d, CC_CURVATURE_2D)
#define SC_SURFACE_TENSION_FORCE_FC IBAMR_FC_FUNC_(sc_surface_tension_force_2d, SC_SURFACE_TENSION_FORCE_2D)
#endif

#if (NDIM == 3)
#define MOLLIFY_IB_4_FC IBAMR_FC_FUNC_(mollify_ib_4_3d, MOLLIFY_IB_4_3D)
#define SC_NORMAL_FC IBAMR_FC_FUNC_(sc_normal_3d, SC_NORMAL_3D)
#define CC_CURVATURE_FC IBAMR_FC_FUNC_(cc_curvature_3d, CC_CURVATURE_3D)
#define SC_SURFACE_TENSION_FORCE_FC IBAMR_FC_FUNC_(sc_surface_tension_force_3d, SC_SURFACE_TENSION_FORCE_3D)
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

    void CC_CURVATURE_FC(double* K,
                         const int& K_gcw,
                         const double* N00,
                         const double* N01,
#if (NDIM == 3)
                         const double* N02,
#endif
                         const double* N10,
                         const double* N11,
#if (NDIM == 3)
                         const double* N12,
                         const double* N20,
                         const double* N21,
                         const double* N22,
#endif
                         const int& N_gcw,
                         const int& ilower0,
                         const int& iupper0,
                         const int& ilower1,
                         const int& iupper1,
#if (NDIM == 3)
                         const int& ilower2,
                         const int& iupper2,
#endif
                         const double* dx);

    void SC_SURFACE_TENSION_FORCE_FC(
#if (NDIM == 2)
        double* F0,
        double* F1,
        const int& F_gcw,
        const double* K,
        const int& K_gcw,
        const double* N00,
        const double* N11,
        const int& N_gcw,
        const int& ilower0,
        const int& iupper0,
        const int& ilower1,
        const int& iupper1
#endif
#if (NDIM == 3)
        double* F0,
        double* F1,
        double* F2,
        const int& F_gcw,
        const double* K,
        const int& K_gcw,
        const double* N00,
        const double* N11,
        const double* N22,
        const int& N_gcw,
        const int& ilower0,
        const int& iupper0,
        const int& ilower1,
        const int& iupper1,
        const int& ilower2,
        const int& iupper2
#endif
    );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
////////////////////////////// PUBLIC ///////////////////////////////////////

SurfaceTensionForceFunction::SurfaceTensionForceFunction(const std::string& object_name,
                                                         const Pointer<Database> input_db,
                                                         const AdvDiffHierarchyIntegrator* adv_diff_solver,
                                                         const Pointer<Variable<NDIM> > level_set_var)
    : CartGridFunction(object_name), d_adv_diff_solver(adv_diff_solver), d_ls_var(level_set_var)
{
    // Set some default values
    d_ts_type = MIDPOINT_RULE;
    d_kernel_fcn = "none";
    d_sigma = 1.0;
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

        d_sigma = input_db->getDoubleWithDefault("sigma", d_sigma);
        d_sigma = input_db->getDoubleWithDefault("surface_tension_coef", d_sigma);

        d_num_interface_cells = input_db->getDoubleWithDefault("num_interface_cells", d_num_interface_cells);
    }
    return;
} // SurfaceTensionForceFunction

void
SurfaceTensionForceFunction::setSmoother(const std::string& kernel_fcn)
{
    d_kernel_fcn = kernel_fcn;
    return;
} // setSmoother

void
SurfaceTensionForceFunction::setSurfaceTensionCoef(double sigma)
{
    d_sigma = sigma;
    return;
} // setSurfaceTensionCoef

void
SurfaceTensionForceFunction::setNumberOfInterfaceCells(double m)
{
    d_num_interface_cells = m;
    return;
} // setNumberOfInterfaceCells

bool
SurfaceTensionForceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
SurfaceTensionForceFunction::setDataOnPatchHierarchy(const int data_idx,
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

    const int coarsest_ln = (coarsest_ln_in == IBTK::invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln =
        (finest_ln_in == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    d_hier_math_ops = new HierarchyMathOps("HierarchyMathOps", hierarchy, coarsest_ln, finest_ln);
    d_hier_sc_data_ops = new HierarchySideDataOpsReal<NDIM, double>(hierarchy, coarsest_ln, finest_ln);

    // Get the newest patch data index for the level set variable
    Pointer<CellVariable<NDIM, double> > phi_cc_var = d_ls_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(!phi_cc_var.isNull());
#endif
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int phi_adv_diff_new_idx =
        var_db->mapVariableAndContextToIndex(phi_cc_var, d_adv_diff_solver->getNewContext());
    const int phi_adv_diff_current_idx =
        var_db->mapVariableAndContextToIndex(phi_cc_var, d_adv_diff_solver->getCurrentContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(phi_adv_diff_new_idx >= 0);
    TBOX_ASSERT(phi_adv_diff_current_idx >= 0);
#endif

    IntVector<NDIM> cell_ghosts = getMinimumGhostWidth(d_kernel_fcn);
    d_C_idx = var_db->registerVariableAndContext(phi_cc_var, var_db->getContext(d_object_name + "::C"), cell_ghosts);
    d_phi_idx =
        var_db->registerVariableAndContext(phi_cc_var, var_db->getContext(d_object_name + "::Phi"), cell_ghosts);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_C_idx, data_time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_phi_idx, data_time);
    }

    // Copy level set into phi and C.
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_ts_type == MIDPOINT_RULE)
    {
        hier_cc_data_ops.linearSum(d_phi_idx,
                                   0.5,
                                   phi_adv_diff_new_idx,
                                   0.5,
                                   phi_adv_diff_current_idx,
                                   /*interior_only*/ true);
    }
    else if (d_ts_type == BACKWARD_EULER)
    {
        hier_cc_data_ops.copyData(d_phi_idx,
                                  phi_adv_diff_new_idx,
                                  /*interior_only*/ true);
    }
    else if (d_ts_type == FORWARD_EULER)
    {
        hier_cc_data_ops.copyData(d_phi_idx,
                                  phi_adv_diff_current_idx,
                                  /*interior_only*/ true);
    }
    else
    {
        TBOX_ERROR("SurfaceTensionForceFunction::setDataOnPatchHierarchy : "
                   << "The class only supports BACKWARD_EULER, FORWARD_EULER, and "
                      "MIDPOINT_RULE"
                   << std::endl);
    }
    hier_cc_data_ops.copyData(d_C_idx, d_phi_idx, /*interior_only*/ true);

    // Convert C to a smoothed heaviside to ensure that the force is only
    // applied near the interface
    convertToHeaviside(d_C_idx, coarsest_ln, finest_ln, hierarchy);

    // Fill ghost cells
    RobinBcCoefStrategy<NDIM>* phi_bc_coef = (d_adv_diff_solver->getPhysicalBcCoefs(phi_cc_var)).front();
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent C_transaction(
        d_C_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, phi_bc_coef);
    InterpolationTransactionComponent phi_transaction(
        d_phi_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, phi_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> C_fill_op = new HierarchyGhostCellInterpolation();
    Pointer<HierarchyGhostCellInterpolation> phi_fill_op = new HierarchyGhostCellInterpolation();
    C_fill_op->initializeOperatorState(C_transaction, hierarchy, coarsest_ln, finest_ln);
    phi_fill_op->initializeOperatorState(phi_transaction, hierarchy, coarsest_ln, finest_ln);

    C_fill_op->fillData(data_time);
    phi_fill_op->fillData(data_time);

    // Mollify C.
    mollifyData(d_C_idx, coarsest_ln, finest_ln, data_time, hierarchy, C_fill_op);

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    // Limit the surface tension force if necessary. This is mainly used in phase change problems
    // to activate the surface tension force only at the liquid-gas interface.
    if (d_mask_surface_tension_force)
    {
        const double apply_time = data_time;
        const double current_time = data_time;
        const double new_time = data_time;
        d_mask_surface_tension_force(data_idx,
                                     d_hier_math_ops,
                                     -1 /*cycle_num*/,
                                     apply_time,
                                     current_time,
                                     new_time,
                                     d_mask_surface_tension_force_ctx);
    }

    // Deallocate and remove scratch/smooth phi.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_phi_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_C_idx);
    }
    var_db->removePatchDataIndex(d_phi_idx);
    var_db->removePatchDataIndex(d_C_idx);

    return;
} // setDataOnPatchHierarchy

void
SurfaceTensionForceFunction::setDataOnPatch(const int data_idx,
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
    Pointer<CellData<NDIM, double> > f_cc_data = f_data;
    Pointer<SideData<NDIM, double> > f_sc_data = f_data;
#if !defined(NDEBUG)
    TBOX_ASSERT(f_cc_data || f_sc_data);
#endif
    if (f_cc_data) f_cc_data->fillAll(0.0);
    if (f_sc_data) f_sc_data->fillAll(0.0);

    if (initial_time) return;

    if (f_cc_data) setDataOnPatchCell(f_cc_data, patch, data_time, initial_time, level);
    if (f_sc_data)
    {
        setDataOnPatchSide(f_sc_data, patch, data_time, initial_time, level);

        if (d_compute_surface_tension_coef)
        {
            // Compute variable surface tension coefficient sigma as F = sigma*F.
            const double apply_time = data_time;
            const double current_time = data_time;
            const double new_time = data_time;
            d_compute_surface_tension_coef(data_idx,
                                           d_hier_math_ops,
                                           -1 /*cycle_num*/,
                                           apply_time,
                                           current_time,
                                           new_time,
                                           d_compute_surface_tension_coef_ctx);
        }
        else
        {
            d_hier_sc_data_ops->scale(data_idx, d_sigma, data_idx, /*interior_only*/ true);
        }
    }
    return;
} // setDataOnPatch

void
SurfaceTensionForceFunction::registerSurfaceTensionForceMasking(MaskSurfaceTensionForcePtr callback, void* ctx)
{
    d_mask_surface_tension_force = callback;
    d_mask_surface_tension_force_ctx = ctx;

    return;
} // registerSurfaceTensionForceMasking

void
SurfaceTensionForceFunction::registerSurfaceTensionCoefficientFunction(ComputeSurfaceTensionCoefficientPtr callback,
                                                                       void* ctx)
{
    d_compute_surface_tension_coef = callback;
    d_compute_surface_tension_coef_ctx = ctx;

    return;
} // registerSurfaceTensionCoefficientFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SurfaceTensionForceFunction::convertToHeaviside(int phi_idx,
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
            double eps = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double phi = (*phi_data)(ci);
                (*phi_data)(ci) = IBTK::smooth_heaviside(phi, eps);
            }
        }
    }
    return;
} // convertToHeaviside

void
SurfaceTensionForceFunction::mollifyData(int smooth_C_idx,
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

            Pointer<CellData<NDIM, double> > smooth_C_data = patch->getPatchData(smooth_C_idx);
            CellData<NDIM, double> C_data(patch_box, /*depth*/ 1, smooth_C_data->getGhostCellWidth());

            C_data.copy(*smooth_C_data);

            double* const V = smooth_C_data->getPointer(0);
            const double* const U = C_data.getPointer(0);
            const int V_gcw = (smooth_C_data->getGhostCellWidth()).max();
            const int U_gcw = (C_data.getGhostCellWidth()).max();

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

void
SurfaceTensionForceFunction::setDataOnPatchCell(Pointer<CellData<NDIM, double> > /*F_data*/,
                                                Pointer<Patch<NDIM> > /*patch*/,
                                                const double /*data_time*/,
                                                const bool /*initial_time*/,
                                                Pointer<PatchLevel<NDIM> > /*level*/)
{
    TBOX_ERROR(
        "SurfaceTensionForceFunction::setDataOnPatchCell() Cell centered "
        "surface force tension is not implemented yet."
        << std::endl);

    return;
} // setDataOnPatchCell

void
SurfaceTensionForceFunction::setDataOnPatchSide(Pointer<SideData<NDIM, double> > F_data,
                                                Pointer<Patch<NDIM> > patch,
                                                const double /*data_time*/,
                                                const bool /*initial_time*/,
                                                Pointer<PatchLevel<NDIM> > /*level*/)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // First find normal in terms of gradient of phi.
    // N = grad(phi)
    SideData<NDIM, double> N(patch_box,
                             /*depth*/ NDIM,
                             /*gcw*/ IntVector<NDIM>(2));
    Pointer<CellData<NDIM, double> > Phi = patch->getPatchData(d_phi_idx);

    SC_NORMAL_FC(N.getPointer(0, 0),
                 N.getPointer(0, 1),
#if (NDIM == 3)
                 N.getPointer(0, 2),
#endif
                 N.getPointer(1, 0),
                 N.getPointer(1, 1),
#if (NDIM == 3)
                 N.getPointer(1, 2),
                 N.getPointer(2, 0),
                 N.getPointer(2, 1),
                 N.getPointer(2, 2),
#endif
                 N.getGhostCellWidth().max(),
                 Phi->getPointer(),
                 Phi->getGhostCellWidth().max(),
                 patch_box.lower(0),
                 patch_box.upper(0),
                 patch_box.lower(1),
                 patch_box.upper(1),
#if (NDIM == 3)
                 patch_box.lower(2),
                 patch_box.upper(2),
#endif
                 dx);

    // Next find the cell centered curvature.
    // K = -div (N/|N|)
    CellData<NDIM, double> K(patch_box, /*depth*/ 1, /*gcw*/ IntVector<NDIM>(1));
    CC_CURVATURE_FC(K.getPointer(),
                    K.getGhostCellWidth().max(),
                    N.getPointer(0, 0),
                    N.getPointer(0, 1),
#if (NDIM == 3)
                    N.getPointer(0, 2),
#endif
                    N.getPointer(1, 0),
                    N.getPointer(1, 1),
#if (NDIM == 3)
                    N.getPointer(1, 2),
                    N.getPointer(2, 0),
                    N.getPointer(2, 1),
                    N.getPointer(2, 2),
#endif
                    N.getGhostCellWidth().max(),
                    patch_box.lower(0),
                    patch_box.upper(0),
                    patch_box.lower(1),
                    patch_box.upper(1),
#if (NDIM == 3)
                    patch_box.lower(2),
                    patch_box.upper(2),
#endif
                    dx);

    // Compute N = grad(C)
    Pointer<CellData<NDIM, double> > C = patch->getPatchData(d_C_idx);

    SC_NORMAL_FC(N.getPointer(0, 0),
                 N.getPointer(0, 1),
#if (NDIM == 3)
                 N.getPointer(0, 2),
#endif
                 N.getPointer(1, 0),
                 N.getPointer(1, 1),
#if (NDIM == 3)
                 N.getPointer(1, 2),
                 N.getPointer(2, 0),
                 N.getPointer(2, 1),
                 N.getPointer(2, 2),
#endif
                 N.getGhostCellWidth().max(),
                 C->getPointer(),
                 C->getGhostCellWidth().max(),
                 patch_box.lower(0),
                 patch_box.upper(0),
                 patch_box.lower(1),
                 patch_box.upper(1),
#if (NDIM == 3)
                 patch_box.lower(2),
                 patch_box.upper(2),
#endif
                 dx);

    // Compute the surface tension force
    // F = K * grad(C)
    SC_SURFACE_TENSION_FORCE_FC(
#if (NDIM == 2)
        F_data->getPointer(0),
        F_data->getPointer(1),
        F_data->getGhostCellWidth().max(),
        K.getPointer(),
        K.getGhostCellWidth().max(),
        N.getPointer(0, 0),
        N.getPointer(1, 1),
        N.getGhostCellWidth().max(),
        patch_box.lower(0),
        patch_box.upper(0),
        patch_box.lower(1),
        patch_box
            .upper(1)
#endif
#if (NDIM == 3)
                F_data->getPointer(0),
        F_data->getPointer(1),
        F_data->getPointer(2),
        F_data->getGhostCellWidth().max(),
        K.getPointer(),
        K.getGhostCellWidth().max(),
        N.getPointer(0, 0),
        N.getPointer(1, 1),
        N.getPointer(2, 2),
        N.getGhostCellWidth().max(),
        patch_box.lower(0),
        patch_box.upper(0),
        patch_box.lower(1),
        patch_box.upper(1),
        patch_box.lower(2),
        patch_box.upper(2)
#endif
    );

    return;
} // setDataOnPatchSide

int
SurfaceTensionForceFunction::getStencilSize(const std::string& kernel_fcn)
{
    if (kernel_fcn == "IB_4") return 4;
    TBOX_ERROR("SurfaceTensionForceFunction::getStencilSize()\n"
               << "  Unknown kernel function " << kernel_fcn << std::endl);
    return -1;

} // getStencilSize

int
SurfaceTensionForceFunction::getMinimumGhostWidth(const std::string& kernel_fcn)
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
