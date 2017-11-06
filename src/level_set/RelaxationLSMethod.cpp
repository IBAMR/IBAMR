// Filename: RelaxationLSMethod.cpp
// Created on 10 Oct 2017 by Nishant Nangia and Amneet Bhalla
//
// Copyright (c) 2002-2017, Nishant Nangia and Amneet Bhalla
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/RelaxationLSMethod.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "IBAMR_config.h"
#include "VariableDatabase.h"
#include "ibamr/namespaces.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "tbox/RestartManager.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define RELAXATION_LS_1ST_ORDER_FC IBAMR_FC_FUNC(relaxationls1storder2d, RELAXATIONLS1STORDER2D)
#define RELAXATION_LS_3RD_ORDER_FC IBAMR_FC_FUNC(relaxationls3rdorder2d, RELAXATIONLS3RDORDER2D)
#endif

#if (NDIM == 3)
#define RELAXATION_LS_1ST_ORDER_FC IBAMR_FC_FUNC(relaxationls1storder3d, RELAXATIONLS1STORDER3D)
#define RELAXATION_LS_3RD_ORDER_FC IBAMR_FC_FUNC(relaxationls3rdorder3d, RELAXATIONLS3RDORDER3D)
#endif

extern "C" {
void RELAXATION_LS_1ST_ORDER_FC(double* U,
                                const int& U_gcw,
                                const double* V,
                                const int& V_gcw,
                                const int& ilower0,
                                const int& iupper0,
                                const int& ilower1,
                                const int& iupper1,
#if (NDIM == 3)
                                const int& ilower2,
                                const int& iupper2,
#endif
                                const double* dx,
                                const int& dir);

void RELAXATION_LS_3RD_ORDER_FC(double* U,
                                const int& U_gcw,
                                const double* V,
                                const int& V_gcw,
                                const int& ilower0,
                                const int& iupper0,
                                const int& ilower1,
                                const int& iupper1,
#if (NDIM == 3)
                                const int& ilower2,
                                const int& iupper2,
#endif
                                const double* dx,
                                const int& dir);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RelaxationLSMethod::RelaxationLSMethod(const std::string& object_name, Pointer<Database> db, bool register_for_restart)
    : LSInitStrategy(object_name, register_for_restart)
{
    // Some default values.
    d_ls_order = FIRST_ORDER_LS;
    d_max_its = 100;
    d_abs_tol = 1e-5;
    d_enable_logging = false;

    // Get any additional or overwrite base class options.
    if (d_registered_for_restart) getFromRestart();
    if (!db.isNull()) getFromInput(db);

    return;
} // RelaxationLSMethod

RelaxationLSMethod::~RelaxationLSMethod()
{
    // intentionally-left blank.
    return;
} // ~RelaxationLSMethod

void
RelaxationLSMethod::initializeLSData(int D_idx, Pointer<HierarchyMathOps> hier_math_ops, double time, bool initial_time)
{
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

    // Create a temporary variable to hold previous iteration values.
    const int D_iter_idx = var_db->registerClonedPatchDataIndex(D_var, D_idx);
    const int D_init_idx = var_db->registerClonedPatchDataIndex(D_var, D_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(D_iter_idx, time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(D_init_idx, time);
    }

    // First, fill cells with some large positive/negative values
    // away from the interface and actual distance value near the interface.
    for (unsigned k = 0; k < d_locate_interface_fcns.size(); ++k)
    {
        (*d_locate_interface_fcns[k])(D_idx, hier_math_ops, time, initial_time, d_locate_interface_fcns_ctx[k]);
    }

    // Set hierarchy objects.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent D_transaction(
        D_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "QUADRATIC", false, d_bc_coef);
    InterpolationTransactionComponent D_init_transaction(
        D_init_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "QUADRATIC", false, d_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> fill_op = new HierarchyGhostCellInterpolation();
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);

    // Carry out relaxation
    double diff_L2_norm = 1.0e12;
    int outer_iter = 0;
    const int cc_wgt_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Copy initial condition
    hier_cc_data_ops.copyData(D_init_idx, D_idx);
    fill_op->initializeOperatorState(D_init_transaction, hierarchy);
    fill_op->fillData(time);

    fill_op->resetTransactionComponent(D_transaction);
    while (diff_L2_norm > d_abs_tol && outer_iter < d_max_its)
    {
        hier_cc_data_ops.copyData(D_iter_idx, D_idx);
        fill_op->fillData(time);
        relax(hier_math_ops, D_idx, D_init_idx, outer_iter);

        // Compute error
        hier_cc_data_ops.axmy(D_iter_idx, 1.0, D_iter_idx, D_idx);
        diff_L2_norm = hier_cc_data_ops.L2Norm(D_iter_idx, cc_wgt_idx);

        outer_iter += 1;

        if (d_enable_logging)
        {
            plog << d_object_name << "::initializeLSData(): After iteration # " << outer_iter << std::endl;
            plog << d_object_name << "::initializeLSData(): L2-norm between successive iterations = " << diff_L2_norm
                 << std::endl;
        }

        if (diff_L2_norm <= d_abs_tol && d_enable_logging)
        {
            plog << d_object_name << "::initializeLSData(): Relaxation converged for entire domain" << std::endl;
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

    // Deallocate the temporary variable.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(D_iter_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(D_init_idx);
    }
    var_db->removePatchDataIndex(D_iter_idx);
    var_db->removePatchDataIndex(D_init_idx);

    return;
} // initializeLSData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
RelaxationLSMethod::relax(Pointer<HierarchyMathOps> hier_math_ops,
                          int dist_idx,
                          int dist_init_idx,
                          const int iter) const
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > dist_data = patch->getPatchData(dist_idx);
            const Pointer<CellData<NDIM, double> > dist_init_data = patch->getPatchData(dist_init_idx);
            relax(dist_data, dist_init_data, patch, iter);
        }
    }
    return;

} // relax

void
RelaxationLSMethod::relax(Pointer<CellData<NDIM, double> > dist_data,
                          const Pointer<CellData<NDIM, double> > dist_init_data,
                          const Pointer<Patch<NDIM> > patch,
                          const int iter) const
{
    double* const D = dist_data->getPointer(0);
    const double* const P = dist_init_data->getPointer(0);
    const int D_ghosts = (dist_data->getGhostCellWidth()).max();
    const int P_ghosts = (dist_init_data->getGhostCellWidth()).max();

#if !defined(NDEBUG)
    TBOX_ASSERT(dist_data->getDepth() == 1);
    TBOX_ASSERT(dist_init_data->getDepth() == 1);
    if (d_ls_order == FIRST_ORDER_LS)
    {
        TBOX_ASSERT(D_ghosts >= 1);
        TBOX_ASSERT(P_ghosts >= 1);
    }
    if (d_ls_order == THIRD_ORDER_LS)
    {
        TBOX_ASSERT(D_ghosts >= 2);
        TBOX_ASSERT(P_ghosts >= 2);
    }
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Get the direction of sweeping (alternates according to iteration number)
    const int num_dirs = (NDIM < 3) ? 4 : 8;
    const int dir = iter % num_dirs;

    if (d_ls_order == FIRST_ORDER_LS)
    {
        RELAXATION_LS_1ST_ORDER_FC(D,
                                   D_ghosts,
                                   P,
                                   P_ghosts,
                                   patch_box.lower(0),
                                   patch_box.upper(0),
                                   patch_box.lower(1),
                                   patch_box.upper(1),
#if (NDIM == 3)
                                   patch_box.lower(2),
                                   patch_box.upper(2),
#endif
                                   dx,
                                   dir);
    }
    else if (d_ls_order == THIRD_ORDER_LS)
    {
        RELAXATION_LS_3RD_ORDER_FC(D,
                                   D_ghosts,
                                   P,
                                   P_ghosts,
                                   patch_box.lower(0),
                                   patch_box.upper(0),
                                   patch_box.lower(1),
                                   patch_box.upper(1),
#if (NDIM == 3)
                                   patch_box.lower(2),
                                   patch_box.upper(2),
#endif
                                   dx,
                                   dir);
    }
    else
    {
        TBOX_ERROR("RelaxationLSMethod does not support " << enum_to_string(d_ls_order) << std::endl);
    }

    return;
} // relax

void
RelaxationLSMethod::getFromInput(Pointer<Database> input_db)
{
    std::string ls_order = "FIRST_ORDER";
    ls_order = input_db->getStringWithDefault("order", ls_order);
    d_ls_order = string_to_enum<LevelSetOrder>(ls_order);

    d_max_its = input_db->getIntegerWithDefault("max_iterations", d_max_its);
    d_max_its = input_db->getIntegerWithDefault("max_its", d_max_its);

    d_abs_tol = input_db->getDoubleWithDefault("abs_tol", d_abs_tol);

    d_enable_logging = input_db->getBoolWithDefault("enable_logging", d_enable_logging);

    return;
} // getFromInput

void
RelaxationLSMethod::getFromRestart()
{
    // intentionally left-blank.
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
