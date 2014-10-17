// Filename: LSetDataIterator.h
// Created on 11 Dec 2009 by Boyce Griffith
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

#ifndef included_LSetDataIterator
#define included_LSetDataIterator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "Box.h"
#include "IndexData.h"
#include "ibtk/LSet.h"
#include "tbox/DescribedClass.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Index;
} // namespace hier
namespace pdat
{
template <int DIM>
class CellGeometry;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
template <class T>
class LSetData;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LSetDataIterator is an iterator class which may be used to
 * iterate through LSet objects associated with a specified box in cell-centered
 * index space.
 */
template <class T>
class LSetDataIterator : public SAMRAI::tbox::DescribedClass
{
public:
    friend class LSetData<T>;

    /*!
     * \brief Class constructor.
     */
    LSetDataIterator();

    /*!
     * \brief Class constructor.
     */
    LSetDataIterator(const LSetDataIterator& that);

    /*!
     * \brief Class destructor.
     */
    virtual ~LSetDataIterator();

    /*!
     * \brief Assignment operator.
     */
    LSetDataIterator<T>& operator=(const LSetDataIterator<T>& that);

    /*!
     * \brief Test two iterators for equality.
     */
    bool operator==(const LSetDataIterator<T>& that);

    /*!
     * \brief Test two iterators for inequality.
     */
    bool operator!=(const LSetDataIterator<T>& that);

    /*!
     * \brief Prefix increment operator.
     */
    LSetDataIterator<T>& operator++();

    /*!
     * \brief Postfix increment operator.
     */
    LSetDataIterator<T> operator++(int);

    /*!
     * \brief Return a reference to the Lagrangian data item referred to by the
     * iterator.
     */
    typename LSet<T>::value_type& operator*() const;

    /*!
     * \brief Return a reference to the Lagrangian data item referred to by the
     * iterator.
     */
    typename LSet<T>::value_type& getDataItem() const;

    /*!
     * \brief Return a const reference to the cell index referred to by the
     * iterator.
     */
    const SAMRAI::hier::Index<NDIM>& getCellIndex() const;

private:
    SAMRAI::hier::Box<NDIM> d_box;
    SAMRAI::pdat::IndexIterator<NDIM, LSet<T>, SAMRAI::pdat::CellGeometry<NDIM> > d_index_it;
    LSet<T>* d_node_set;
    typename LSet<T>::iterator d_node_it;
};

} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LSetDataIterator-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LSetDataIterator
