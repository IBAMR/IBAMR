// Filename: LNodeInitStrategy.C
// Created on 11 Jul 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "LNodeInitStrategy.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeInitStrategy::LNodeInitStrategy()
{
    // intentionally blank
    return;
}// LNodeInitStrategy

LNodeInitStrategy::~LNodeInitStrategy()
{
    // intentionally blank
    return;
}// ~LNodeInitStrategy


void
LNodeInitStrategy::initializeStructureIndexingOnPatchLevel(
    std::map<int,std::string>& strct_id_to_strct_name_map,
    std::map<int,std::pair<int,int> >& strct_id_to_lag_idx_range_map,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    TBOX_WARNING("LNodeInitStrategy::initializeStructureIndexingOnPatchLevel()\n"
                 << "  default implementation employed, no indexing data provided.\n");
    return;
}// initializeStructureIndexingOnPatchLevel

int
LNodeInitStrategy::initializeMassDataOnPatchLevel(
    const int global_index_offset,
    const int local_index_offset,
    Pointer<LData>& M_data,
    Pointer<LData>& K_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    TBOX_WARNING("LNodeInitStrategy::initializeMassDataOnPatchLevel()\n"
                 << "  default implementation employed, no mass data initialized.\n");
    return 0;
}// initializeMassDataOnPatchLevel

int
LNodeInitStrategy::initializeDirectorDataOnPatchLevel(
    const int global_index_offset,
    const int local_index_offset,
    Pointer<LData>& D_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    TBOX_WARNING("LNodeInitStrategy::initializeDirectorDataOnPatchLevel()\n"
                 << "  default implementation employed, no director data initialized.\n");
    return 0;
}// initializeDirectorDataOnPatchLevel

void
LNodeInitStrategy::tagCellsForInitialRefinement(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    TBOX_WARNING("LNodeInitStrategy::tagCellsForInitialRefinement()\n"
                 << "  default implementation employed, no cells tagged for refinement.\n");
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::LNodeInitStrategy>;

//////////////////////////////////////////////////////////////////////////////
