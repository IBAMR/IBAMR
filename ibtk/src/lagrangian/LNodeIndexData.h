// Filename: LNodeIndexData.h
// Created on 04 Jun 2007 by Boyce Griffith
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

#ifndef included_LNodeIndexData
#define included_LNodeIndexData

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LNodeIndexSet.h>
#include <ibtk/IndexUtilities.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CellGeometry.h>
#include <CellIterator.h>
#include <IndexData.h>
#include <IntVector.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class LDataManager;
class LNodeIndexIterator;
}// namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LNodeIndexData is a specialization of the templated class
 * SAMRAI::pdat::IndexData that provides access to the Lagrangian and <A
 * HREF="http://www-unix.mcs.anl.gov/petsc">PETSc</A> index information and any
 * data associated with the Lagrangian nodes in the interior and ghost cell
 * region of a Cartesian grid patch.
 *
 * \see SAMRAI::pdat::IndexData
 */
class LNodeIndexData
    : public SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet,SAMRAI::pdat::CellGeometry<NDIM> >
{
public:
    friend class LDataManager;

    /*!
     * The cell iterator iterates over the elements of a cell centered box
     * geometry.
     */
    typedef SAMRAI::pdat::CellIterator<NDIM> CellIterator;

    /*!
     * The LNodeIndexSet iterator iterates over the LNodeIndexSet elements
     * within the IndexData patch data object.
     */
    typedef SAMRAI::pdat::IndexIterator<NDIM,LNodeIndexSet,SAMRAI::pdat::CellGeometry<NDIM> > LNodeIndexSetPatchIterator;

    /*!
     * The LNodeIndex iterator iterates over the LNodeIndex elements located
     * within a cell centered box geometry.
     */
    typedef IBTK::LNodeIndexIterator LNodeIndexIterator;

    /*!
     * Return an iterator to the first LNodeIndex in the specified region of
     * index space.
     */
    LNodeIndexIterator
    lnode_index_begin(
        const SAMRAI::hier::Box<NDIM>& box);

    /*!
     * Return an iterator pointing to the end of the collection of LNodeIndex
     * objects associated with the patch data object.
     */
    LNodeIndexIterator
    lnode_index_end();

    /*!
     * The constructor for an SAMRAI::pdat::IndexData<NDIM> object.  The box
     * describes the interior of the index space and the ghosts vector describes
     * the ghost nodes in each coordinate direction.
     */
    LNodeIndexData(
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * The virtual destructor for an LNodeIndexData object.
     */
    virtual
    ~LNodeIndexData();

    /*!
     * \return A constant reference to the set of local PETSc indices to
     * Lagrangian nodes which lie in the patch interior.
     *
     * \note This should be an ordered, contiguous set of indices.
     */
    const std::vector<int>&
    getInteriorLocalIndices() const;

    /*!
     * \return A constant reference to the set of local PETSc indices to
     * Lagrangian nodes which lie in ghost cell region of this
     * SAMRAI::hier::PatchData<NDIM> object.
     */
    const std::vector<int>&
    getGhostLocalIndices() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LNodeIndexData();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndexData(
        const LNodeIndexData& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndexData&
    operator=(
        const LNodeIndexData& that);

    std::vector<int> d_interior_local_indices;
    std::vector<int> d_ghost_local_indices;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/LNodeIndexData.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexData
