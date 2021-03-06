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

#include "ibamr/FastSweepingLSMethod.h"
#include "ibamr/LSInitStrategy.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"

#include "BasePatchLevel.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define FAST_SWEEP_1ST_ORDER_FC IBAMR_FC_FUNC(fastsweep1storder2d, FASTSWEEP1STORDER2D)
#endif

#if (NDIM == 3)
#define FAST_SWEEP_1ST_ORDER_FC IBAMR_FC_FUNC(fastsweep1storder3d, FASTSWEEP1STORDER3D)
#endif

extern "C"
{
    void FAST_SWEEP_1ST_ORDER_FC(double* U,
                                 const int& U_gcw,
                                 const int& ilower0,
                                 const int& iupper0,
                                 const int& ilower1,
                                 const int& iupper1,
#if (NDIM == 3)
                                 const int& ilower2,
                                 const int& iupper2,
#endif
                                 const int& dlower0,
                                 const int& dupper0,
                                 const int& dlower1,
                                 const int& dupper1,
#if (NDIM == 3)

                                 const int& dlower2,
                                 const int& dupper2,
#endif
                                 const double* dx,
                                 const int& patch_touches_bdry,
                                 const int* touches_wall_loc_idx);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FastSweepingLSMethod::FastSweepingLSMethod(std::string object_name, Pointer<Database> db, bool register_for_restart)
    : LSInitStrategy(std::move(object_name), register_for_restart)
{
    for (int& wall_idx : d_wall_location_idx) wall_idx = 0;

    if (d_registered_for_restart) getFromRestart();
    if (!db.isNull()) getFromInput(db);

    return;
} // FastSweepingLSMethod

void
FastSweepingLSMethod::initializeLSData(int D_idx,
                                       Pointer<HierarchyMathOps> hier_math_ops,
                                       int integrator_step,
                                       double time,
                                       bool initial_time)
{
    bool initialize_ls =
        d_reinitialize_ls || initial_time || (d_reinit_interval && integrator_step % d_reinit_interval == 0);
    if (!initialize_ls) return;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > data_var;
    var_db->mapIndexToVariable(D_idx, data_var);
    Pointer<CellVariable<NDIM, double> > D_var = data_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(!D_var.isNull());
#endif

    Pointer<PatchHierarchy<NDIM> > hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Create a temporary variable to hold previous iteration values with
    // appropriate ghost cell width since it is not guaranteed that D_idx will
    // have proper ghost cell width.
    IntVector<NDIM> cell_ghosts;
    if (d_ls_order == FIRST_ORDER_LS)
    {
        cell_ghosts = 1;
    }
    else
    {
        TBOX_ERROR("FastSweepLSMethod does not support " << enum_to_string(d_ls_order) << std::endl);
    }
    const int D_scratch_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);
    const int D_iter_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::ITER"), cell_ghosts);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(D_scratch_idx, time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(D_iter_idx, time);
    }

    // First, fill cells with some large positive/negative values
    // away from the interface and actual distance value near the interface.
    for (unsigned k = 0; k < d_locate_interface_fcns.size(); ++k)
    {
        (*d_locate_interface_fcns[k])(D_scratch_idx, hier_math_ops, time, initial_time, d_locate_interface_fcns_ctx[k]);
    }

    // Set hierarchy objects.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent D_transaction(
        D_scratch_idx, "LINEAR_REFINE", true, "NONE", "LINEAR", false, d_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> fill_op = new HierarchyGhostCellInterpolation();
    fill_op->initializeOperatorState(D_transaction, hierarchy);
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);

    // Carry out iterations
    double diff_L2_norm = 1.0e12;
    int outer_iter = 0;
    const int cc_wgt_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

    while (diff_L2_norm > d_abs_tol && outer_iter < d_max_its)
    {
        hier_cc_data_ops.copyData(D_iter_idx, D_scratch_idx);
        fill_op->fillData(time);

        fastSweep(hier_math_ops, D_scratch_idx);

        hier_cc_data_ops.axmy(D_iter_idx, 1.0, D_iter_idx, D_scratch_idx);
        diff_L2_norm = hier_cc_data_ops.L2Norm(D_iter_idx, cc_wgt_idx);

        outer_iter += 1;

        if (d_enable_logging)
        {
            plog << d_object_name << "::initializeLSData(): After iteration # " << outer_iter << std::endl;
            plog << d_object_name << "::initializeLSData(): L2-norm between successive iterations = " << diff_L2_norm
                 << std::endl;
        }

        if (diff_L2_norm <= d_abs_tol)
        {
            plog << d_object_name
                 << "::initializeLSData(): Fast sweeping algorithm "
                    "converged for entire domain"
                 << std::endl;
        }
    }

    if (outer_iter >= d_max_its)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << "::initializeLSData(): Reached maximum allowable outer iterations" << std::endl;
            plog << d_object_name << "::initializeLSData(): ||distance_new - distance_old||_2 = " << diff_L2_norm
                 << std::endl;
        }
    }

    // Copy signed distance into supplied patch data index
    hier_cc_data_ops.copyData(D_idx, D_scratch_idx);

    // Deallocate the temporary variable.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(D_scratch_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(D_iter_idx);
    }
    var_db->removePatchDataIndex(D_scratch_idx);
    var_db->removePatchDataIndex(D_iter_idx);

    // Indicate that the LS has been initialized.
    d_reinitialize_ls = false;

    return;
} // initializeLSData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FastSweepingLSMethod::fastSweep(Pointer<HierarchyMathOps> hier_math_ops, int dist_idx) const
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        const BoxArray<NDIM>& domain_boxes = level->getPhysicalDomain();
#if !defined(NDEBUG)
        TBOX_ASSERT(domain_boxes.size() == 1);
#endif

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > dist_data = patch->getPatchData(dist_idx);
            fastSweep(dist_data, patch, domain_boxes[0]);
        }
    }
    return;

} // fastSweep

void
FastSweepingLSMethod::fastSweep(Pointer<CellData<NDIM, double> > dist_data,
                                const Pointer<Patch<NDIM> > patch,
                                const Box<NDIM>& domain_box) const
{
    double* const D = dist_data->getPointer(0);
    const int D_ghosts = (dist_data->getGhostCellWidth()).max();

    // Check if the patch touches physical domain.
    int touches_wall_loc_idx[NDIM * 2] = { 0 };
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const bool patch_touches_bdry = pgeom->getTouchesRegularBoundary() || pgeom->getTouchesPeriodicBoundary();
    if (patch_touches_bdry)
    {
        int loc_idx = 0;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int upperlower = 0; upperlower < 2; ++upperlower, ++loc_idx)
            {
                touches_wall_loc_idx[loc_idx] = d_consider_phys_bdry_wall &&
                                                pgeom->getTouchesRegularBoundary(axis, upperlower) &&
                                                d_wall_location_idx[loc_idx];
            }
        }
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(dist_data->getDepth() == 1);
    if (d_ls_order == FIRST_ORDER_LS) TBOX_ASSERT(D_ghosts >= 1);
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    const double* const dx = pgeom->getDx();
    if (d_ls_order == FIRST_ORDER_LS)
    {
        FAST_SWEEP_1ST_ORDER_FC(D,
                                D_ghosts,
                                patch_box.lower(0),
                                patch_box.upper(0),
                                patch_box.lower(1),
                                patch_box.upper(1),
#if (NDIM == 3)
                                patch_box.lower(2),
                                patch_box.upper(2),
#endif
                                domain_box.lower(0),
                                domain_box.upper(0),
                                domain_box.lower(1),
                                domain_box.upper(1),
#if (NDIM == 3)
                                domain_box.lower(2),
                                domain_box.upper(2),
#endif
                                dx,
                                patch_touches_bdry,
                                touches_wall_loc_idx);
    }
    else
    {
        TBOX_ERROR("FastSweepingLSMethod does not support " << enum_to_string(d_ls_order) << std::endl);
    }

    return;
} // fastSweep

void
FastSweepingLSMethod::getFromInput(Pointer<Database> input_db)
{
    std::string ls_order = "FIRST_ORDER";
    ls_order = input_db->getStringWithDefault("order", ls_order);
    d_ls_order = string_to_enum<LevelSetOrder>(ls_order);

    d_max_its = input_db->getIntegerWithDefault("max_iterations", d_max_its);
    d_max_its = input_db->getIntegerWithDefault("max_its", d_max_its);

    d_abs_tol = input_db->getDoubleWithDefault("abs_tol", d_abs_tol);

    d_enable_logging = input_db->getBoolWithDefault("enable_logging", d_enable_logging);

    d_reinit_interval = input_db->getIntegerWithDefault("reinit_interval", d_reinit_interval);

    d_consider_phys_bdry_wall = input_db->getBoolWithDefault("physical_bdry_wall", d_consider_phys_bdry_wall);
    Array<int> wall_loc_idices;
    if (input_db->keyExists("physical_bdry_wall_loc_idx"))
    {
        input_db->getArray("physical_bdry_wall_loc_idx", wall_loc_idices);
    }
    for (int k = 0; k < wall_loc_idices.size(); ++k)
    {
        d_wall_location_idx[wall_loc_idices[k]] = 1;
    }

    return;
} // getFromInput

void
FastSweepingLSMethod::getFromRestart()
{
    // intentionally left-blank.
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
