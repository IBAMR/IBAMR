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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/SAMRAIDataCache.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/samrai_compatibility_names.h>

#include <MultiblockDataTranslator.h>
#include <SAMRAICellDataFactory.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIEdgeDataFactory.h>
#include <SAMRAIEdgeVariable.h>
#include <SAMRAIFaceDataFactory.h>
#include <SAMRAIFaceVariable.h>
#include <SAMRAINodeDataFactory.h>
#include <SAMRAINodeVariable.h>
#include <SAMRAIOuteredgeDataFactory.h>
#include <SAMRAIOuteredgeVariable.h>
#include <SAMRAIOuterfaceDataFactory.h>
#include <SAMRAIOuterfaceVariable.h>
#include <SAMRAIOuternodeDataFactory.h>
#include <SAMRAIOuternodeVariable.h>
#include <SAMRAIOutersideDataFactory.h>
#include <SAMRAIOutersideVariable.h>
#include <SAMRAIPatchDescriptor.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAISideDataFactory.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIUtilities.h>
#include <SAMRAIVariable.h>
#include <SAMRAIVariableDatabase.h>

#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
template <typename U, typename T>
inline bool
can_convert_to(T* t)
{
    return dynamic_cast<U*>(t) != nullptr;
}

template <typename U, typename T>
inline bool
can_convert_to(const SAMRAIPointer<T>& t)
{
    return dynamic_cast<U*>(t.getPointer()) != nullptr;
}

template <typename FactoryType>
std::pair</*depth*/ int, /*ghost_width*/ int>
get_characteristics(SAMRAIPointer<SAMRAIPatchDescriptor> patch_descriptor, const int idx)
{
    SAMRAIPointer<FactoryType> pdat_fac = patch_descriptor->getPatchDataFactory(idx);
    const int depth = pdat_fac->getDefaultDepth();
    const int ghost_width = pdat_fac->getGhostCellWidth().max();
#if !defined NDEBUG
    TBOX_ASSERT(ghost_width == pdat_fac->getGhostCellWidth().min());
#endif
    return std::make_pair(depth, ghost_width);
}

template <typename T>
inline bool
get_data_characteristics(const int idx, int& depth, int& ghost_width)
{
    auto var_db = SAMRAIVariableDatabase::getDatabase();
    SAMRAIPointer<SAMRAIVariable> var;
    var_db->mapIndexToVariable(idx, var);
    auto patch_descriptor = var_db->getPatchDescriptor();

    // Determine the data depth and ghost width from the patch data factory.
    //
    // TODO: Make this more generic and extensible.
    bool convertable = false;
    std::pair<int, int> characteristics;
    if (can_convert_to<SAMRAICellVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAICellDataFactory<T>>(patch_descriptor, idx);
    }
    else if (can_convert_to<SAMRAIEdgeVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAIEdgeDataFactory<T>>(patch_descriptor, idx);
    }
    else if (can_convert_to<SAMRAIFaceVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAIFaceDataFactory<T>>(patch_descriptor, idx);
    }
    else if (can_convert_to<SAMRAINodeVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAINodeDataFactory<T>>(patch_descriptor, idx);
    }
    else if (can_convert_to<SAMRAIOuteredgeVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAIOuteredgeDataFactory<T>>(patch_descriptor, idx);
    }
    else if (can_convert_to<SAMRAIOuterfaceVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAIOuterfaceDataFactory<T>>(patch_descriptor, idx);
    }
    else if (can_convert_to<SAMRAIOuternodeVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAIOuternodeDataFactory<T>>(patch_descriptor, idx);
    }
    else if (can_convert_to<SAMRAIOutersideVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAIOutersideDataFactory<T>>(patch_descriptor, idx);
    }
    else if (can_convert_to<SAMRAISideVariable<T>>(var))
    {
        convertable = true;
        characteristics = get_characteristics<SAMRAISideDataFactory<T>>(patch_descriptor, idx);
    }
    depth = characteristics.first;
    ghost_width = characteristics.second;
    return convertable;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

SAMRAIDataCache::~SAMRAIDataCache()
{
    setPatchHierarchy(nullptr);
    auto var_db = SAMRAIVariableDatabase::getDatabase();
    for (auto cloned_idx : d_all_cloned_patch_data_idxs)
    {
        var_db->removePatchDataIndex(cloned_idx);
    }
}

void
SAMRAIDataCache::setPatchHierarchy(SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy)
{
    if (hierarchy != d_hierarchy && d_hierarchy && (d_coarsest_ln != IBTK::invalid_level_number) &&
        (d_finest_ln != IBTK::invalid_level_number))
    {
        // Clean up allocated patch data on the old hierarchy.
        for (auto cloned_idx : d_all_cloned_patch_data_idxs)
        {
            for (auto ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(cloned_idx)) level->deallocatePatchData(cloned_idx);
            }
        }
    }
    d_hierarchy = hierarchy;
    if (!hierarchy) resetLevels(IBTK::invalid_level_number, IBTK::invalid_level_number);
}

void
SAMRAIDataCache::resetLevels(const int coarsest_ln, const int finest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(
        (d_hierarchy &&
         ((0 <= coarsest_ln) && (coarsest_ln <= finest_ln) && (finest_ln <= d_hierarchy->getFinestLevelNumber()))) ||
        (!d_hierarchy && (coarsest_ln == IBTK::invalid_level_number) && (finest_ln == IBTK::invalid_level_number)));
#endif
    if (d_hierarchy)
    {
        // Clean up allocated patch data on the old range of levels.
        for (auto cloned_idx : d_all_cloned_patch_data_idxs)
        {
            for (auto ln = d_coarsest_ln; ln < coarsest_ln; ++ln)
            {
                SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(cloned_idx)) level->deallocatePatchData(cloned_idx);
            }
            for (auto ln = finest_ln + 1; ln <= d_finest_ln; ++ln)
            {
                SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(cloned_idx)) level->deallocatePatchData(cloned_idx);
            }
        }
    }
    d_coarsest_ln = coarsest_ln;
    d_finest_ln = finest_ln;
}

int
SAMRAIDataCache::getCoarsestLevelNumber() const
{
    return d_coarsest_ln;
} // getCoarsestLevelNumber

int
SAMRAIDataCache::getFinestLevelNumber() const
{
    return d_finest_ln;
} // getFinestLevelNumber

/////////////////////////////// PRIVATE //////////////////////////////////////

int
SAMRAIDataCache::lookupCachedPatchDataIndex(const int idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy && (0 <= d_coarsest_ln) && (d_coarsest_ln <= d_finest_ln) &&
                (d_finest_ln <= d_hierarchy->getFinestLevelNumber()));
#endif
    // Find the patch data index corresponding to the data descriptor.
    int cloned_idx;
    const SAMRAIDataCache::key_type data_descriptor = construct_data_descriptor(idx);
    auto it = d_available_data_idx_map.find(data_descriptor);
    if (it == d_available_data_idx_map.end())
    {
        auto var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<SAMRAIVariable> var;
        var_db->mapIndexToVariable(idx, var);
        cloned_idx = var_db->registerClonedPatchDataIndex(var, idx);
        d_all_cloned_patch_data_idxs.insert(cloned_idx);
        it = d_available_data_idx_map.emplace(data_descriptor, idx);
    }
    else
    {
        cloned_idx = it->second;
    }

    // Remove the index from the collection of available data indices and store it in the collection of checked-out
    // indices.
    d_available_data_idx_map.erase(it);
    d_unavailable_data_idx_map.emplace(data_descriptor, cloned_idx);

    // Allocate data if needed.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(cloned_idx)) level->allocatePatchData(cloned_idx);
    }
    return cloned_idx;
}

void
SAMRAIDataCache::restoreCachedPatchDataIndex(const int cached_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy && (0 <= d_coarsest_ln) && (d_coarsest_ln <= d_finest_ln) &&
                (d_finest_ln <= d_hierarchy->getFinestLevelNumber()));
#endif
    // Find the index in the collection of checked-out indices.
    const SAMRAIDataCache::key_type data_descriptor = construct_data_descriptor(cached_idx);
    const auto data_range = d_unavailable_data_idx_map.equal_range(data_descriptor);
    auto it = data_range.first;
    for (; it != data_range.second; ++it)
    {
        if (it->second == cached_idx) break;
    }
    if (it == data_range.second)
    {
        TBOX_ERROR("could not find cached index: " << cached_idx << "\n");
    }
    d_available_data_idx_map.emplace(data_descriptor, cached_idx);
    d_unavailable_data_idx_map.erase(it);
    return;
}

SAMRAIDataCache::key_type
SAMRAIDataCache::construct_data_descriptor(const int idx)
{
    // TODO: Make this more generic and extensible.
    SAMRAIPointer<SAMRAIVariable> var;
    auto var_db = SAMRAIVariableDatabase::getDatabase();
    var_db->mapIndexToVariable(idx, var);
    SAMRAIVariable& var_ref = *var;
    const std::type_index var_type_id = typeid(var_ref);
    int depth = 0, ghost_width = 0;
    bool convertable = get_data_characteristics<double>(idx, depth, ghost_width) ||
                       get_data_characteristics<float>(idx, depth, ghost_width) ||
                       get_data_characteristics<int>(idx, depth, ghost_width);
    if (!convertable)
    {
        TBOX_ERROR("unsupported data alignment for SAMRAIDataCache: " << var_type_id.name() << "\n");
    }
    return SAMRAIDataCache::key_type{ var_type_id, depth, ghost_width };
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
