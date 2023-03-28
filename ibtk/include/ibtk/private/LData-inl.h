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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_LData_inl_h
#define included_IBTK_LData_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline const std::string&
LData::getName() const
{
    return d_name;
} // getName

inline unsigned int
LData::getGlobalNodeCount() const
{
    return d_global_node_count;
} // getGlobalNodeCount

inline unsigned int
LData::getLocalNodeCount() const
{
    return d_local_node_count;
} // getLocalNodeCount

inline unsigned int
LData::getGhostNodeCount() const
{
    return d_ghost_node_count;
} // getGhostNodeCount

inline unsigned int
LData::getDepth() const
{
    return d_depth;
} // getDepth

inline Vec
LData::getVec()
{
    restoreArrays();
    return d_global_vec;
} // getVec

inline boost::multi_array_ref<double, 1>*
LData::getArray()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_array) getArrayCommon();
    return &d_boost_array;
} // getArray

inline boost::multi_array_ref<double, 1>*
LData::getLocalFormArray()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_array) getArrayCommon();
    return &d_boost_local_array;
} // getLocalFormArray

inline boost::multi_array_ref<double, 1>*
LData::getGhostedLocalFormArray()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_ghosted_local_array) getGhostedLocalFormArrayCommon();
    return &d_boost_ghosted_local_array;
} // getGhostedLocalFormArray

inline boost::multi_array_ref<double, 2>*
LData::getVecArray()
{
    if (!d_array) getArrayCommon();
    return &d_boost_vec_array;
} // getVecArray

inline boost::multi_array_ref<double, 2>*
LData::getLocalFormVecArray()
{
    if (!d_array) getArrayCommon();
    return &d_boost_local_vec_array;
} // getLocalFormVecArray

inline boost::multi_array_ref<double, 2>*
LData::getGhostedLocalFormVecArray()
{
    if (!d_ghosted_local_array) getGhostedLocalFormArrayCommon();
    return &d_boost_vec_ghosted_local_array;
} // getGhostedLocalFormVecArray

inline void
LData::restoreArrays()
{
    int ierr;
    if (d_ghosted_local_array)
    {
        ierr = VecRestoreArray(d_ghosted_local_vec, &d_ghosted_local_array);
        IBTK_CHKERRQ(ierr);
        d_ghosted_local_array = nullptr;
    }
    if (d_ghosted_local_vec)
    {
        ierr = VecGhostRestoreLocalForm(d_global_vec, &d_ghosted_local_vec);
        IBTK_CHKERRQ(ierr);
        d_ghosted_local_vec = nullptr;
    }
    if (d_array)
    {
        ierr = VecRestoreArray(d_global_vec, &d_array);
        IBTK_CHKERRQ(ierr);
        d_array = nullptr;
    }
    return;
} // restoreArray

inline void
LData::beginGhostUpdate()
{
    const int ierr = VecGhostUpdateBegin(getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    return;
} // beginGhostUpdate

inline void
LData::endGhostUpdate()
{
    const int ierr = VecGhostUpdateEnd(getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    return;
} // endGhostUpdate

/////////////////////////////// PRIVATE //////////////////////////////////////

/**
 * Internal function for setting up a PETSc function of given size and
 * depth. See the documentation of LData for more information.
 */
inline Vec
setup_petsc_vector(const unsigned int num_local_nodes,
                   const unsigned int depth,
                   const std::vector<int>& nonlocal_petsc_indices)
{
    int ierr;
    Vec global_vec;
    if (depth == 1)
    {
        ierr = VecCreateGhost(PETSC_COMM_WORLD,
                              num_local_nodes,
                              PETSC_DECIDE,
                              static_cast<PetscInt>(nonlocal_petsc_indices.size()),
                              nonlocal_petsc_indices.data(),
                              &global_vec);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        ierr = VecCreateGhostBlock(PETSC_COMM_WORLD,
                                   depth,
                                   depth * num_local_nodes,
                                   PETSC_DECIDE,
                                   static_cast<PetscInt>(nonlocal_petsc_indices.size()),
                                   nonlocal_petsc_indices.data(),
                                   &global_vec);
        IBTK_CHKERRQ(ierr);
    }
    return global_vec;
}

template <std::size_t dim>
void
LData::destroy_ref(boost::multi_array_ref<double, dim>& ref)
{
    // destructor calls have difficult name lookup rules
    using T = typename boost::multi_array_ref<double, dim>;
    ref.~T();
}

inline void
LData::getArrayCommon()
{
    if (!d_array)
    {
        int ierr = VecGetArray(d_global_vec, &d_array);
        IBTK_CHKERRQ(ierr);
        int ilower, iupper;
        ierr = VecGetOwnershipRange(d_global_vec, &ilower, &iupper);
        IBTK_CHKERRQ(ierr);
        using range = boost::multi_array_types::extent_range;
        if (d_depth == 1)
        {
            using array_type = boost::multi_array<double, 1>;
            array_type::extent_gen extents;
            destroy_ref(d_boost_array);
            new (&d_boost_array) boost::multi_array_ref<double, 1>(d_array, extents[range(ilower, iupper)]);
            destroy_ref(d_boost_local_array);
            new (&d_boost_local_array) boost::multi_array_ref<double, 1>(d_array, extents[iupper - ilower]);
        }
        using array_type = boost::multi_array<double, 2>;
        array_type::extent_gen extents;
        // There is no way to reinitialize a boost array with a new data
        // pointer, so destroy the current arrays and build new ones in their
        // place:
        destroy_ref(d_boost_vec_array);
        new (&d_boost_vec_array)
            boost::multi_array_ref<double, 2>(d_array, extents[range(ilower / d_depth, iupper / d_depth)][d_depth]);
        destroy_ref(d_boost_local_vec_array);
        new (&d_boost_local_vec_array)
            boost::multi_array_ref<double, 2>(d_array, extents[(iupper - ilower) / d_depth][d_depth]);
    }
    return;
} // getArrayCommon

inline void
LData::getGhostedLocalFormArrayCommon()
{
    if (!d_ghosted_local_vec)
    {
        int ierr = VecGhostGetLocalForm(d_global_vec, &d_ghosted_local_vec);
        IBTK_CHKERRQ(ierr);
    }
    if (!d_ghosted_local_array)
    {
        int ierr = VecGetArray(d_ghosted_local_vec, &d_ghosted_local_array);
        IBTK_CHKERRQ(ierr);
        int ilower, iupper;
        ierr = VecGetOwnershipRange(d_ghosted_local_vec, &ilower, &iupper);
        IBTK_CHKERRQ(ierr);
        if (d_depth == 1)
        {
            destroy_ref(d_boost_ghosted_local_array);
            new (&d_boost_ghosted_local_array)
                boost::multi_array_ref<double, 1>(d_ghosted_local_array, boost::extents[iupper - ilower]);
        }
        destroy_ref(d_boost_vec_ghosted_local_array);
        new (&d_boost_vec_ghosted_local_array) boost::multi_array_ref<double, 2>(
            d_ghosted_local_array, boost::extents[(iupper - ilower) / d_depth][d_depth]);
    }
    return;
} // getGhostedLocalFormArrayCommon

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LData_inl_h
