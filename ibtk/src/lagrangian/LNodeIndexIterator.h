// Filename: LNodeIndexIterator.h
// Created on 11 Dec 2009 by Boyce Griffith
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

#ifndef included_LNodeIndexIterator
#define included_LNodeIndexIterator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LNodeIndexSet.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CellGeometry.h>
#include <IndexData.h>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class LNodeIndexData;
}// namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LNodeIndexIterator is an iterator class which may be used to
 * iterate through LNodeIndex objects associated with a specified box in
 * cell-centered index space.
 */
class LNodeIndexIterator
    : public SAMRAI::tbox::DescribedClass
{
public:
    friend class LNodeIndexData;

    /*!
     * \brief Class constructor.
     */
    LNodeIndexIterator();

    /*!
     * \brief Class constructor.
     */
    LNodeIndexIterator(
        const LNodeIndexIterator& that);

    /*!
     * \brief Class destructor.
     */
    ~LNodeIndexIterator();

    /*!
     * \brief Assignment operator.
     */
    LNodeIndexIterator&
    operator=(
        const LNodeIndexIterator& that);

    /*!
     * \brief Test two iterators for equality.
     */
    bool
    operator==(
        const LNodeIndexIterator& that);

    /*!
     * \brief Test two iterators for inequality.
     */
    bool
    operator!=(
        const LNodeIndexIterator& that);

    /*!
     * \brief Prefix increment operator.
     */
    LNodeIndexIterator&
    operator++();

    /*!
     * \brief Postfix increment operator.
     */
    LNodeIndexIterator
    operator++(
        int);

    /*!
     * \brief Return a reference to the LNodeIndex referred to by the iterator.
     */
    LNodeIndex&
    operator*() const;

    /*!
     * \brief Return a reference to the LNodeIndex referred to by the iterator.
     */
    LNodeIndex&
    getLNodeIndex() const;

    /*!
     * \brief Return a const reference to the cell index referred to by the
     * iterator.
     */
    const SAMRAI::hier::Index<NDIM>&
    getCellIndex() const;

private:
    SAMRAI::hier::Box<NDIM> d_box;
    SAMRAI::pdat::IndexIterator<NDIM,LNodeIndexSet,SAMRAI::pdat::CellGeometry<NDIM> > d_index_it;
    LNodeIndexSet* d_node_set;
    LNodeIndexSet::iterator d_node_it;
};

}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/LNodeIndexIterator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexIterator
