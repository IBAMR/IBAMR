// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibamr/AdvDiffHierarchyIntegrator.h"

#include "BasePatchHierarchy.h"
#include "PatchHierarchy.h"
#include "TagInterfaceRefinementCells.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

TagInterfaceRefinementCells::TagInterfaceRefinementCells(Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                         Pointer<CellVariable<NDIM, double> > scalar_var,
                                                         double tag_min_val,
                                                         double tag_max_val)
    : d_adv_diff_solver(adv_diff_solver),
      d_scalar_var(scalar_var),
      d_tag_min_val(tag_min_val),
      d_tag_max_val(tag_max_val)
{
    return;
} // TagInterfaceRefinementCells

void
callTagInterfaceRefinementCellsCallbackFunction(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                const int level_number,
                                                const double /*error_data_time*/,
                                                const int tag_index,
                                                const bool /*initial_time*/,
                                                const bool /*uses_richardson_extrapolation_too*/,
                                                void* ctx)
{
    TagInterfaceRefinementCells* ptr_tagger = static_cast<TagInterfaceRefinementCells*>(ctx);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int scalar_current_idx = var_db->mapVariableAndContextToIndex(
        ptr_tagger->d_scalar_var, ptr_tagger->d_adv_diff_solver->getCurrentContext());

    // Tag cells based on the value of the level set variable
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
        Pointer<CellData<NDIM, double> > scalar_data = patch->getPatchData(scalar_current_idx);

        for (CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            const double scalar = (*scalar_data)(i);

            // if (scalar >= ptr_tagger->d_tag_min_val && scalar <= ptr_tagger->d_tag_max_val)
            if (scalar >= ptr_tagger->d_tag_min_val)
            {
                (*tags_data)(i) = 1;
            }
        }
    }

    return;

} // callTagInterfaceRefinementCellsCallBackFunction
