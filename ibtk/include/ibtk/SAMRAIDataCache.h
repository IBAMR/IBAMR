// Filename: SAMRAIDataCache.h
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

#ifndef included_IBTK_SAMRAIDataCache
#define included_IBTK_SAMRAIDataCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <typeindex>

#include "PatchHierarchy.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class SAMRAIDataCache is a utility class for caching cloned SAMRAI patch data.  Patch data are allocated as
 * needed and should not be deallocated by the caller.
 */
class SAMRAIDataCache : public SAMRAI::tbox::DescribedClass
{
public:
    /// \brief Default constructor.
    SAMRAIDataCache() = default;

    /// \brief Destructor.
    ~SAMRAIDataCache();

    /// \name Methods to set the hierarchy and range of levels.
    //\{

    /**
     * @brief      Set the patch hierarchy to use in allocating patch data.
     *
     * @param[in]  hierarchy  The patch hierarchy
     */
    void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /**
     * @brief      Reset the range of patch levels over which operations occur.
     *
     * @param[in]  coarsest_ln  The coarsest level number
     * @param[in]  finest_ln    The finest level number
     */
    void resetLevels(int coarsest_ln, int finest_ln);

    //\}

    /**
     * @brief      Class for accessing cached patch data indices.
     *
     * @note       Cached indices will be released when all references to them are deleted.
     */
    class CachedPatchDataIndex
    {
    public:
        friend SAMRAIDataCache;

        CachedPatchDataIndex() = delete;

        ~CachedPatchDataIndex() = default;

        inline operator int() const
        {
            return getPatchDataIndex();
        }

        inline int getPatchDataIndex() const
        {
            return d_idx_data->d_idx;
        }

    private:
        inline CachedPatchDataIndex(const int idx, SAMRAIDataCache* const cache)
            : d_idx_data(new PatchDataIndexHandle(idx, cache))
        {
        }

        struct PatchDataIndexHandle
        {
            PatchDataIndexHandle() = delete;

            PatchDataIndexHandle(const PatchDataIndexHandle&) = delete;

            PatchDataIndexHandle& operator=(const PatchDataIndexHandle&) = delete;

            inline PatchDataIndexHandle(const int idx, SAMRAIDataCache* const cache) : d_idx(idx), d_cache(cache)
            {
            }

            inline ~PatchDataIndexHandle()
            {
                d_cache->restoreCachedPatchDataIndex(d_idx);
            }

            int d_idx;
            SAMRAIDataCache* d_cache;
        };

        std::shared_ptr<PatchDataIndexHandle> d_idx_data;
    };

    /**
     * @brief      Gets the cached patch data index.
     *
     * @param[in]  idx  The patch data index clone or lookup.
     *
     * @return     The cached patch data index.
     */
    inline CachedPatchDataIndex getCachedPatchDataIndex(int idx)
    {
        return CachedPatchDataIndex(lookupCachedPatchDataIndex(idx), this);
    }

private:
    /**
     * @brief      Lookup the cached patch data index.
     *
     * @param[in]  idx  The patch data index clone or lookup.
     *
     * @return     The cached patch data index.
     *
     * @note       Patch data indices should be restored to the caching object via restoreCachedPatchDataIndex().
     */
    int lookupCachedPatchDataIndex(int idx);

    /**
     * @brief      Restore the cached patch data index, to allow other to use the corresponding data.
     *
     * @param[in]  cached_idx  The cached index
     *
     * @note       cached_idx must have been obtained via getCachedPatchDataIndex().
     */
    void restoreCachedPatchDataIndex(int cached_idx);

    /// \brief Disable the copy constructor.
    SAMRAIDataCache(const SAMRAIDataCache& from) = delete;

    /// \brief Disable the assignment operator.
    SAMRAIDataCache& operator=(const SAMRAIDataCache& that) = delete;

    /// \brief The patch hierarchy and range of levels to use in allocating/deallocating patch data.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    /// \brief Key type for looking up cached data.
    using key_type = std::tuple<std::type_index, /*data_depth*/ int, /*ghost_cell_width*/ int>;

    /// \brief Mappings from data type to registered patch data indices.
    std::multimap<key_type, int> d_available_data_idx_map, d_unavailable_data_idx_map;

    /// \brief Set of all patch data indices cloned by this object.
    std::set<int> d_all_cloned_patch_data_idxs;

    /// \brief Construct the data descriptor for a given variable and patch data index.
    static key_type construct_data_descriptor(int idx);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_SAMRAIDataCache
