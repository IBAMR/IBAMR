// Filename: SAMRAIDataCache.cpp
// Created on 23 Apr 2019 by Boyce Griffith
//
// Copyright (c) 2002-2019, Boyce Griffith
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

#include "ibtk/SAMRAIDataCache.h"
#include "CellDataFactory.h"
#include "CellVariable.h"
#include "EdgeDataFactory.h"
#include "EdgeVariable.h"
#include "FaceDataFactory.h"
#include "FaceVariable.h"
#include "NodeDataFactory.h"
#include "NodeVariable.h"
#include "OuteredgeDataFactory.h"
#include "OuteredgeVariable.h"
#include "OuterfaceDataFactory.h"
#include "OuterfaceVariable.h"
#include "OuternodeDataFactory.h"
#include "OuternodeVariable.h"
#include "OutersideDataFactory.h"
#include "OutersideVariable.h"
#include "SideDataFactory.h"
#include "SideVariable.h"
#include "VariableDatabase.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

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
can_convert_to(const SAMRAI::tbox::Pointer<T>& t)
{
    return dynamic_cast<U*>(t.getPointer()) != nullptr;
}

template <typename FactoryType>
std::pair</*depth*/ int, /*ghost_width*/ int>
get_characteristics(Pointer<PatchDescriptor<NDIM> > patch_descriptor, const int idx)
{
    Pointer<FactoryType> pdat_fac = patch_descriptor->getPatchDataFactory(idx);
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
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;
    var_db->mapIndexToVariable(idx, var);
    auto patch_descriptor = var_db->getPatchDescriptor();

    // Determine the data depth and ghost width from the patch data factory.
    //
    // TODO: Make this more generic and extensible.
    bool convertable = false;
    std::pair<int, int> characteristics;
    if (can_convert_to<CellVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<CellDataFactory<NDIM, T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<EdgeVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<EdgeDataFactory<NDIM, T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<FaceVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<FaceDataFactory<NDIM, T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<NodeVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<NodeDataFactory<NDIM, T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<OuteredgeVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<OuteredgeDataFactory<NDIM, T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<OuterfaceVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<OuterfaceDataFactory<NDIM, T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<OuternodeVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<OuternodeDataFactory<NDIM, T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<OutersideVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<OutersideDataFactory<NDIM, T> >(patch_descriptor, idx);
    }
    else if (can_convert_to<SideVariable<NDIM, T> >(var))
    {
        convertable = true;
        characteristics = get_characteristics<SideDataFactory<NDIM, T> >(patch_descriptor, idx);
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
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    for (auto cloned_idx : d_all_cloned_patch_data_idxs)
    {
        var_db->removePatchDataIndex(cloned_idx);
    }
}

void
SAMRAIDataCache::setPatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    if (hierarchy != d_hierarchy && d_hierarchy && (d_coarsest_ln != IBTK::invalid_level_number) &&
        (d_finest_ln != IBTK::invalid_level_number))
    {
        // Clean up allocated patch data on the old hierarchy.
        for (auto cloned_idx : d_all_cloned_patch_data_idxs)
        {
            for (auto ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
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
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(cloned_idx)) level->deallocatePatchData(cloned_idx);
            }
            for (auto ln = finest_ln + 1; ln <= d_finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(cloned_idx)) level->deallocatePatchData(cloned_idx);
            }
        }
    }
    d_coarsest_ln = coarsest_ln;
    d_finest_ln = finest_ln;
}

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
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<Variable<NDIM> > var;
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
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
    Pointer<Variable<NDIM> > var;
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->mapIndexToVariable(idx, var);
    Variable<NDIM>& var_ref = *var;
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
