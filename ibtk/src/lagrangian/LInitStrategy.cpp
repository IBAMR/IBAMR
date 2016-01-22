// Filename: LInitStrategy.cpp
// Created on 11 Jul 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "ibtk/LInitStrategy.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace IBTK
{
class LDataManager;
} // namespace IBTK

namespace IBTK
{
class LData;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LInitStrategy::LInitStrategy()
{
    // intentionally blank
    return;
} // LInitStrategy

LInitStrategy::~LInitStrategy()
{
    // intentionally blank
    return;
} // ~LInitStrategy

void
LInitStrategy::initializeStructureIndexingOnPatchLevel(
    std::map<int, std::string>& /*strct_id_to_strct_name_map*/,
    std::map<int, std::pair<int, int> >& /*strct_id_to_lag_idx_range_map*/,
    const int /*level_number*/,
    const double /*init_data_time*/,
    const bool /*can_be_refined*/,
    const bool /*initial_time*/,
    LDataManager* const /*l_data_manager*/)
{
    TBOX_WARNING("LInitStrategy::initializeStructureIndexingOnPatchLevel()\n"
                 << "  default implementation employed, no indexing data provided.\n");
    return;
} // initializeStructureIndexingOnPatchLevel

unsigned int
LInitStrategy::initializeMassDataOnPatchLevel(const unsigned int /*global_index_offset*/,
                                              const unsigned int /*local_index_offset*/,
                                              Pointer<LData> /*M_data*/,
                                              Pointer<LData> /*K_data*/,
                                              const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                              const int /*level_number*/,
                                              const double /*init_data_time*/,
                                              const bool /*can_be_refined*/,
                                              const bool /*initial_time*/,
                                              LDataManager* const /*l_data_manager*/)
{
    TBOX_WARNING("LInitStrategy::initializeMassDataOnPatchLevel()\n"
                 << "  default implementation employed, no mass data initialized.\n");
    return 0;
} // initializeMassDataOnPatchLevel

unsigned int
LInitStrategy::initializeDirectorDataOnPatchLevel(const unsigned int /*global_index_offset*/,
                                                  const unsigned int /*local_index_offset*/,
                                                  Pointer<LData> /*D_data*/,
                                                  const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                  const int /*level_number*/,
                                                  const double /*init_data_time*/,
                                                  const bool /*can_be_refined*/,
                                                  const bool /*initial_time*/,
                                                  LDataManager* const /*l_data_manager*/)
{
    TBOX_WARNING("LInitStrategy::initializeDirectorDataOnPatchLevel()\n"
                 << "  default implementation employed, no director data initialized.\n");
    return 0;
} // initializeDirectorDataOnPatchLevel

void
LInitStrategy::tagCellsForInitialRefinement(const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                            const int /*level_number*/,
                                            const double /*error_data_time*/,
                                            const int /*tag_index*/)
{
    TBOX_WARNING("LInitStrategy::tagCellsForInitialRefinement()\n"
                 << "  default implementation employed, no cells tagged for refinement.\n");
    return;
} // tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
