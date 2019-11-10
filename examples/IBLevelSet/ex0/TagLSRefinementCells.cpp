// Filename: TagLSRefinementCells.cpp
// Created on Dec 28, 2017 by Nishant Nangia

// APPLICATION INCLUDES
#include "TagLSRefinementCells.h"
#include <fstream>

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callTagLSRefinementCellsCallbackFunction(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_current_idx = var_db->mapVariableAndContextToIndex(
        ptr_ls_tagger->d_ls_gas_var, ptr_ls_tagger->d_adv_diff_solver->getCurrentContext());

    // Tag cells based on the value of the level set variable
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
        Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_current_idx);

        for (CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            const double dist_norm = std::abs((*ls_data)(i)-ptr_ls_tagger->d_tag_value);

            if (dist_norm <= ptr_ls_tagger->d_tag_abs_thresh)
            {
                (*tags_data)(i) = 1;
            }
        }
    }

    return;

} // callTagLSRefinementCellsCallBackFunction
