// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "BoussinesqForcing.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

BoussinesqForcing::BoussinesqForcing(Pointer<Variable<NDIM> > T_var,
                                     Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
                                     int gamma)
    : d_T_var(T_var), d_adv_diff_hier_integrator(adv_diff_hier_integrator), d_gamma(gamma)
{
    // intentionally blank
    return;
} // BoussinesqForcing

BoussinesqForcing::~BoussinesqForcing()
{
    // intentionally blank
    return;
} // ~BoussinesqForcing

bool
BoussinesqForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
BoussinesqForcing::setDataOnPatchHierarchy(const int data_idx,
                                           Pointer<Variable<NDIM> > var,
                                           Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           const double data_time,
                                           const bool initial_time,
                                           const int coarsest_ln_in,
                                           const int finest_ln_in)
{
    // Allocate scratch data when needed.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, d_adv_diff_hier_integrator->getScratchContext());
    const bool T_scratch_is_allocated = d_adv_diff_hier_integrator->isAllocatedPatchData(T_scratch_idx);
    if (!T_scratch_is_allocated)
    {
        d_adv_diff_hier_integrator->allocatePatchData(T_scratch_idx, data_time);
    }

    // Communicate ghost-cell data.
    if (!initial_time)
    {
        int T_current_idx =
            var_db->mapVariableAndContextToIndex(d_T_var, d_adv_diff_hier_integrator->getCurrentContext());
        int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, d_adv_diff_hier_integrator->getNewContext());
        const bool T_new_is_allocated = d_adv_diff_hier_integrator->isAllocatedPatchData(T_new_idx);
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(d_T_var, hierarchy, /*get_unique*/ true);
        if (d_adv_diff_hier_integrator->getCurrentCycleNumber() == 0 || !T_new_is_allocated)
        {
            hier_cc_data_ops->copyData(T_scratch_idx, T_current_idx);
        }
        else
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(d_adv_diff_hier_integrator->getCurrentCycleNumber() > 0);
#endif
            hier_cc_data_ops->linearSum(T_scratch_idx, 0.5, T_current_idx, 0.5, T_new_idx);
        }
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent ghost_fill_component(T_scratch_idx,
                                                               "CONSERVATIVE_LINEAR_REFINE",
                                                               false,
                                                               "CONSERVATIVE_COARSEN",
                                                               "LINEAR",
                                                               false,
                                                               d_adv_diff_hier_integrator->getPhysicalBcCoefs(d_T_var));
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_fill_component, hierarchy);
        ghost_fill_op.fillData(data_time);
    }

    // Compute a staggered-grid approximation to -gamma*T on each patch level.
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }

    // Deallocate scratch data when needed.
    if (!T_scratch_is_allocated)
    {
        d_adv_diff_hier_integrator->deallocatePatchData(T_scratch_idx);
    }
    return;
} // setDataOnPatchHierarchy

void
BoussinesqForcing::setDataOnPatch(const int data_idx,
                                  Pointer<Variable<NDIM> > /*var*/,
                                  Pointer<Patch<NDIM> > patch,
                                  const double /*data_time*/,
                                  const bool initial_time,
                                  Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);
    if (initial_time) return;
    Pointer<CellData<NDIM, double> > T_scratch_data =
        patch->getPatchData(d_T_var, d_adv_diff_hier_integrator->getScratchContext());
    const Box<NDIM>& patch_box = patch->getBox();
    const int axis = NDIM - 1;
    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
    {
        SideIndex<NDIM> s_i(it(), axis, 0);
        (*F_data)(s_i) = -d_gamma * 0.5 * ((*T_scratch_data)(s_i.toCell(1)) + (*T_scratch_data)(s_i.toCell(0)));
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
