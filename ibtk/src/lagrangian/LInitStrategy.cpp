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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LInitStrategy.h"

#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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
