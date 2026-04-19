// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_ibtk_utilities_inl
#define included_IBTK_ibtk_utilities_inl

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

namespace IBTK
{
inline void
allocate_patch_data(const int idx,
                    const double data_time,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(idx)) level->allocatePatchData(idx, data_time);
    }
}

template <class IntegerContainer>
inline void
allocate_patch_data(const IntegerContainer& idxs,
                    const double data_time,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (const auto idx : idxs) allocate_patch_data(idx, data_time, hierarchy);
}

inline void
allocate_patch_data(const std::initializer_list<int> idxs,
                    const double data_time,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (const int idx : idxs) allocate_patch_data(idx, data_time, hierarchy);
}

inline void
allocate_patch_data(const SAMRAI::hier::ComponentSelector& idxs,
                    const double data_time,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (int idx = 0; idx < idxs.getSize(); ++idx)
    {
        if (idxs.isSet(idx)) allocate_patch_data(idx, data_time, hierarchy);
    }
}

inline void
deallocate_patch_data(const int idx, const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(idx)) level->deallocatePatchData(idx);
    }
}

template <class IntegerContainer>
inline void
deallocate_patch_data(const IntegerContainer& idxs,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (const auto idx : idxs) deallocate_patch_data(idx, hierarchy);
}

inline void
deallocate_patch_data(const std::initializer_list<int> idxs,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (const int idx : idxs) deallocate_patch_data(idx, hierarchy);
}

inline void
deallocate_patch_data(const SAMRAI::hier::ComponentSelector& idxs,
                      const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (int idx = 0; idx < idxs.getSize(); ++idx)
    {
        if (idxs.isSet(idx)) deallocate_patch_data(idx, hierarchy);
    }
}
} // namespace IBTK

#endif // #ifndef included_IBTK_ibtk_utilities_inl
