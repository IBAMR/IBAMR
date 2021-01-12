// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <SAMRAI_config.h>

// Headers for basic SAMRAI objects
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <CartesianGridGeometry.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects

#ifndef included_IBAMR_multiphase_flow_TagLSRefinementCells
#define included_IBAMR_multiphase_flow_TagLSRefinementCells

struct TagLSRefinementCells
{
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_gas_var;
    double d_tag_value;
    double d_tag_abs_thresh;
};

inline void
callTagLSRefinementCellsCallbackFunction(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
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
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    const int ls_current_idx = var_db->mapVariableAndContextToIndex(
        ptr_ls_tagger->d_ls_gas_var, ptr_ls_tagger->d_adv_diff_solver->getCurrentContext());

    // Tag cells based on the value of the level set variable
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > ls_data = patch->getPatchData(ls_current_idx);

        for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const SAMRAI::hier::Index<NDIM>& i = ic();
            const double dist_norm = std::abs((*ls_data)(i)-ptr_ls_tagger->d_tag_value);

            if (dist_norm <= ptr_ls_tagger->d_tag_abs_thresh)
            {
                (*tags_data)(i) = 1;
            }
        }
    }

    return;

} // callTagLSRefinementCellsCallBackFunction

#endif
