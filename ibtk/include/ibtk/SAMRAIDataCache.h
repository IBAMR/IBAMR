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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_SAMRAIDataCache
#define included_IBTK_SAMRAIDataCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <typeindex>

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

    /**
     * Return the coarsest patch level for which scratch data is allocated.
     */
    int getCoarsestLevelNumber() const;

    /**
     * Return the finest patch level for which scratch data is allocated.
     */
    int getFinestLevelNumber() const;

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

    /// \brief The patch hierarchy under consideration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    /// Coarsest level of allocated patch data.
    int d_coarsest_ln = IBTK::invalid_level_number;

    /// Coarsest level of allocated patch data.
    int d_finest_ln = IBTK::invalid_level_number;

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
