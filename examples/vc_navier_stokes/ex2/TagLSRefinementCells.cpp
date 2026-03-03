// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2021 by the IBAMR developers
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
#include <ibtk/samrai_compatibility_names.h>

#include "TagLSRefinementCells.h"

#include <SAMRAIBasePatchHierarchy.h>
#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIterator.h>
#include <SAMRAIIndex.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIVariableDatabase.h>

#include <fstream>

#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callTagLSRefinementCellsCallbackFunction(const Pointer<SAMRAIBasePatchHierarchy> hierarchy,
                                         const int level_number,
                                         const double /*error_data_time*/,
                                         const int tag_index,
                                         const bool initial_time,
                                         const bool /*uses_richardson_extrapolation_too*/,
                                         void* ctx)
{
    if (initial_time || level_number == hierarchy->getFinestLevelNumber()) return;

    TagLSRefinementCells* ptr_ls_tagger = static_cast<TagLSRefinementCells*>(ctx);

    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number < hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));

    // Get the current level set information
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    const int ls_current_idx = var_db->mapVariableAndContextToIndex(
        ptr_ls_tagger->d_ls_gas_var, ptr_ls_tagger->d_adv_diff_solver->getCurrentContext());

    // Tag cells based on the value of the level set variable
    Pointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(level_number);
    for (SAMRAIPatchLevel::Iterator p(level); p; p++)
    {
        Pointer<SAMRAIPatch> patch = level->getPatch(p());
        const SAMRAIBox& patch_box = patch->getBox();
        Pointer<SAMRAICellData<int>> tags_data = patch->getPatchData(tag_index);
        Pointer<SAMRAICellData<double>> ls_data = patch->getPatchData(ls_current_idx);

        for (SAMRAICellIterator ic(patch_box); ic; ic++)
        {
            const SAMRAIIndex& i = ic();
            const double dist = (*ls_data)(i);
            if (dist - ptr_ls_tagger->d_tag_value <= ptr_ls_tagger->d_tag_abs_thresh)
            {
                (*tags_data)(i) = 1;
            }
        }
    }

    return;

} // callTagLSRefinementCellsCallBackFunction
