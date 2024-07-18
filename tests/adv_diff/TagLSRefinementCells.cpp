// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// APPLICATION INCLUDES
#include "TagLSRefinementCells.h"

#include <CartesianGridGeometry.h>

#include <fstream>

#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callTagSolidLSRefinementCellsCallbackFunction(const SAMRAIPointer<BasePatchHierarchyNd> hierarchy,
                                              const int level_number,
                                              const double /*error_data_time*/,
                                              const int tag_index,
                                              const bool /*initial_time*/,
                                              const bool /*uses_richardson_extrapolation_too*/,
                                              void* ctx)
{
    TagLSRefinementCells* ptr_ls_tagger = static_cast<TagLSRefinementCells*>(ctx);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    // Get the current level set information
    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
    const int ls_current_idx = var_db->mapVariableAndContextToIndex(
        ptr_ls_tagger->d_ls_var, ptr_ls_tagger->d_adv_diff_solver->getCurrentContext());

    // Tag cells based on the value of the level set variable
    SAMRAIPointer<PatchLevelNd> level = hierarchy->getPatchLevel(level_number);
    for (PatchLevelNd::Iterator p(level); p; p++)
    {
        SAMRAIPointer<PatchNd> patch = level->getPatch(p());
        const BoxNd& patch_box = patch->getBox();
        SAMRAIPointer<CellDataNd<int> > tags_data = patch->getPatchData(tag_index);
        SAMRAIPointer<CellDataNd<double> > ls_data = patch->getPatchData(ls_current_idx);

        for (CellIteratorNd ic(patch_box); ic; ic++)
        {
            const hier::IndexNd& i = ic();
            const double dist = (*ls_data)(i);

            if (std::abs(dist) <=
                (level_number > 0 ? ptr_ls_tagger->d_tag_abs_thresh : 2 * ptr_ls_tagger->d_tag_abs_thresh))
            {
                (*tags_data)(i) = 1;
            }
        }
    }

    return;

} // callTagSolidLSRefinementCellsCallBackFunction
