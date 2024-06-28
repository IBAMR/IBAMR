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
#include "ibamr/MarangoniSurfaceTensionForceFunction.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"

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
#define SC_NORMAL_FC IBAMR_FC_FUNC(sc_normal_2d, SC_NORMAL_2D)
#define SC_MARANGONI_FORCE_FC IBAMR_FC_FUNC(sc_marangoni_force_2d, SC_MARANGONI_FORCE_2D)
#endif

#if (NDIM == 3)
#define SC_NORMAL_FC IBAMR_FC_FUNC(sc_normal_3d, SC_NORMAL_3D)
#define SC_MARANGONI_FORCE_FC IBAMR_FC_FUNC(sc_marangoni_force_3d, SC_MARANGONI_FORCE_3D)
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

    void SC_MARANGONI_FORCE_FC(
#if (NDIM == 2)
        double* F0,
        double* F1,
        const int& F_gcw,
        double* gradT00,
        double* gradT01,
        double* gradT10,
        double* gradT11,
        const int& gradT_gcw,
        double* N00,
        double* N01,
        double* N10,
        double* N11,
        const int& N_gcw,
        double* gradC00,
        double* gradC01,
        double* gradC10,
        double* gradC11,
        const int& gradC_gcw,
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
        double* gradT00,
        double* gradT01,
        double* gradT02,
        double* gradT10,
        double* gradT11,
        double* gradT12,
        double* gradT20,
        double* gradT21,
        double* gradT22,
        const int& gradT_gcw,
        double* N00,
        double* N01,
        double* N02,
        double* N10,
        double* N11,
        double* N12,
        double* N20,
        double* N21,
        double* N22,
        const int& N_gcw,
        double* gradC00,
        double* gradC01,
        double* gradC02,
        double* gradC10,
        double* gradC11,
        double* gradC12,
        double* gradC20,
        double* gradC21,
        double* gradC22,
        const int& gradC_gcw,
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

MarangoniSurfaceTensionForceFunction::MarangoniSurfaceTensionForceFunction(const std::string& object_name,
                                                                           const Pointer<Database> input_db,
                                                                           AdvDiffHierarchyIntegrator* adv_diff_solver,
                                                                           Pointer<Variable<NDIM> > level_set_var,
                                                                           Pointer<Variable<NDIM> > T_var,
                                                                           RobinBcCoefStrategy<NDIM>*& T_bc_coef)
    : SurfaceTensionForceFunction(object_name, input_db, adv_diff_solver, level_set_var)
{
    d_T_var = T_var;
    d_T_bc_coef = T_bc_coef;
    if (input_db->keyExists("marangoni_coefficient"))
        d_marangoni_coefficient = input_db->getDouble("marangoni_coefficient");

    return;
} // SurfaceTensionForceFunction

bool
MarangoniSurfaceTensionForceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
MarangoniSurfaceTensionForceFunction::setDataOnPatchHierarchy(const int data_idx,
                                                              Pointer<Variable<NDIM> > var,
                                                              Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                              const double data_time,
                                                              const bool initial_time,
                                                              const int coarsest_ln_in,
                                                              const int finest_ln_in)
{
    // Patch index to store F_marangoni.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_F_cloned_idx = var_db->registerClonedPatchDataIndex(var, data_idx);

    // Get the patch data index for the temperature variable
    Pointer<CellVariable<NDIM, double> > T_var = d_T_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(!T_var.isNull());
#endif
    const int T_new_idx = var_db->mapVariableAndContextToIndex(T_var, d_adv_diff_solver->getNewContext());
    const int T_current_idx = var_db->mapVariableAndContextToIndex(T_var, d_adv_diff_solver->getCurrentContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(T_new_idx >= 0);
    TBOX_ASSERT(T_current_idx >= 0);
#endif

    IntVector<NDIM> cell_ghosts = getMinimumGhostWidth(d_kernel_fcn);
    d_T_idx = var_db->registerVariableAndContext(T_var, var_db->getContext(d_object_name + "::T"), cell_ghosts);

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_F_cloned_idx, data_time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_T_idx, data_time);
    }

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);
    if (d_ts_type == MIDPOINT_RULE)
    {
        hier_cc_data_ops.linearSum(d_T_idx,
                                   0.5,
                                   T_new_idx,
                                   0.5,
                                   T_current_idx,
                                   /*interior_only*/ true);
    }
    else if (d_ts_type == BACKWARD_EULER)
    {
        hier_cc_data_ops.copyData(d_T_idx,
                                  T_new_idx,
                                  /*interior_only*/ true);
    }
    else if (d_ts_type == FORWARD_EULER)
    {
        hier_cc_data_ops.copyData(d_T_idx,
                                  T_current_idx,
                                  /*interior_only*/ true);
    }
    else
    {
        TBOX_ERROR("MarangoniSurfaceTensionForceFunction::setDataOnPatchHierarchy : "
                   << "The class only supports BACKWARD_EULER, FORWARD_EULER, and "
                      "MIDPOINT_RULE"
                   << std::endl);
    }

    // Fill ghost cells
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent T_transaction(
        d_T_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, d_T_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> T_fill_op = new HierarchyGhostCellInterpolation();
    T_fill_op->initializeOperatorState(T_transaction, hierarchy, coarsest_ln, finest_ln);

    T_fill_op->fillData(data_time);

    SurfaceTensionForceFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    // Deallocate and remove scratch/smooth phi.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_F_cloned_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_T_idx);
    }
    var_db->removePatchDataIndex(d_F_cloned_idx);
    var_db->removePatchDataIndex(d_T_idx);

    return;
} // setDataOnPatchHierarchy

void
MarangoniSurfaceTensionForceFunction::setDataOnPatch(const int data_idx,
                                                     Pointer<Variable<NDIM> > var,
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

    SurfaceTensionForceFunction::setDataOnPatch(data_idx, var, patch, data_time, initial_time, level);

    if (f_cc_data) setDataOnPatchCell(f_cc_data, patch, data_time, initial_time, level);
    if (f_sc_data)
    {
        Pointer<SideData<NDIM, double> > f_marangoni_sc_data = patch->getPatchData(d_F_cloned_idx);
        setDataOnPatchSide(f_marangoni_sc_data, patch, data_time, initial_time, level);

        if (d_compute_marangoni_coef)
        {
            // Use the callback function to compute F = marangoni_coef*F.
            const double apply_time = data_time;
            const double current_time = data_time;
            const double new_time = data_time;
            d_compute_marangoni_coef(d_F_cloned_idx,
                                     d_hier_math_ops,
                                     -1 /*cycle_num*/,
                                     apply_time,
                                     current_time,
                                     new_time,
                                     d_compute_marangoni_coef_ctx);
        }
        else
        {
            d_hier_sc_data_ops->scale(d_F_cloned_idx, d_marangoni_coefficient, d_F_cloned_idx, /*interior_only*/ true);
        }

        // Add Marangoni force with the surface tension force.
        d_hier_sc_data_ops->add(data_idx, data_idx, d_F_cloned_idx, /*interior_only*/ true);
    }
    return;
} // setDataOnPatch

void
MarangoniSurfaceTensionForceFunction::registerMarangoniCoefficientFunction(ComputeMarangoniCoefPtr callback, void* ctx)
{
    d_compute_marangoni_coef = callback;
    d_compute_marangoni_coef_ctx = ctx;
    return;
} // registerMarangoniCoefficientFunction

/////////////////////////////// PRIVATE ////////////////////////////////////

void
MarangoniSurfaceTensionForceFunction::setDataOnPatchCell(Pointer<CellData<NDIM, double> > /*F_data*/,
                                                         Pointer<Patch<NDIM> > /*patch*/,
                                                         const double /*data_time*/,
                                                         const bool /*initial_time*/,
                                                         Pointer<PatchLevel<NDIM> > /*level*/)
{
    TBOX_ERROR(
        "MarangoniSurfaceTensionForceFunction::setDataOnPatchCell() Cell centered "
        "surface force tension is not implemented yet."
        << std::endl);

    return;
} // setDataOnPatchCell

void
MarangoniSurfaceTensionForceFunction::setDataOnPatchSide(Pointer<SideData<NDIM, double> > F_data,
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

    // Find the gradient of T at the side-center.
    SideData<NDIM, double> grad_T(patch_box, /*depth*/ NDIM, /*gcw*/ IntVector<NDIM>(2));
    Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(d_T_idx);
    SC_NORMAL_FC(grad_T.getPointer(0, 0),
                 grad_T.getPointer(0, 1),
#if (NDIM == 3)
                 grad_T.getPointer(0, 2),
#endif
                 grad_T.getPointer(1, 0),
                 grad_T.getPointer(1, 1),
#if (NDIM == 3)
                 grad_T.getPointer(1, 2),
                 grad_T.getPointer(2, 0),
                 grad_T.getPointer(2, 1),
                 grad_T.getPointer(2, 2),
#endif
                 grad_T.getGhostCellWidth().max(),
                 T_data->getPointer(),
                 T_data->getGhostCellWidth().max(),
                 patch_box.lower(0),
                 patch_box.upper(0),
                 patch_box.lower(1),
                 patch_box.upper(1),
#if (NDIM == 3)
                 patch_box.lower(2),
                 patch_box.upper(2),
#endif
                 dx);

    // Find the gradient of C at the side-center.
    SideData<NDIM, double> grad_C(patch_box, /*depth*/ NDIM, /*gcw*/ IntVector<NDIM>(2));
    Pointer<CellData<NDIM, double> > C_data = patch->getPatchData(d_C_idx);
    SC_NORMAL_FC(grad_C.getPointer(0, 0),
                 grad_C.getPointer(0, 1),
#if (NDIM == 3)
                 grad_C.getPointer(0, 2),
#endif
                 grad_C.getPointer(1, 0),
                 grad_C.getPointer(1, 1),
#if (NDIM == 3)
                 grad_C.getPointer(1, 2),
                 grad_C.getPointer(2, 0),
                 grad_C.getPointer(2, 1),
                 grad_C.getPointer(2, 2),
#endif
                 grad_C.getGhostCellWidth().max(),
                 C_data->getPointer(),
                 C_data->getGhostCellWidth().max(),
                 patch_box.lower(0),
                 patch_box.upper(0),
                 patch_box.lower(1),
                 patch_box.upper(1),
#if (NDIM == 3)
                 patch_box.lower(2),
                 patch_box.upper(2),
#endif
                 dx);

    // Compute F_marangoni = (grad T |grad C| - (grad T dot grad \phi) grad C).
    SC_MARANGONI_FORCE_FC(
#if (NDIM == 2)
        F_data->getPointer(0),
        F_data->getPointer(1),
        F_data->getGhostCellWidth().max(),
        grad_T.getPointer(0, 0),
        grad_T.getPointer(0, 1),
        grad_T.getPointer(1, 0),
        grad_T.getPointer(1, 1),
        grad_T.getGhostCellWidth().max(),
        N.getPointer(0, 0),
        N.getPointer(0, 1),
        N.getPointer(1, 0),
        N.getPointer(1, 1),
        N.getGhostCellWidth().max(),
        grad_C.getPointer(0, 0),
        grad_C.getPointer(0, 1),
        grad_C.getPointer(1, 0),
        grad_C.getPointer(1, 1),
        grad_C.getGhostCellWidth().max(),
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
        grad_T.getPointer(0, 0),
        grad_T.getPointer(0, 1),
        grad_T.getPointer(0, 2),
        grad_T.getPointer(1, 0),
        grad_T.getPointer(1, 1),
        grad_T.getPointer(1, 2),
        grad_T.getPointer(2, 0),
        grad_T.getPointer(2, 1),
        grad_T.getPointer(2, 2),
        grad_T.getGhostCellWidth().max(),
        N.getPointer(0, 0),
        N.getPointer(0, 1),
        N.getPointer(0, 2),
        N.getPointer(1, 0),
        N.getPointer(1, 1),
        N.getPointer(1, 2),
        N.getPointer(2, 0),
        N.getPointer(2, 1),
        N.getPointer(2, 2),
        N.getGhostCellWidth().max(),
        grad_C.getPointer(0, 0),
        grad_C.getPointer(0, 1),
        grad_C.getPointer(0, 2),
        grad_C.getPointer(1, 0),
        grad_C.getPointer(1, 1),
        grad_C.getPointer(1, 2),
        grad_C.getPointer(2, 0),
        grad_C.getPointer(2, 1),
        grad_C.getPointer(2, 2),
        grad_C.getGhostCellWidth().max(),
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
