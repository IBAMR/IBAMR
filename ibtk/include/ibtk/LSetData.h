// Filename: LSetData.h
// Created on 04 Jun 2007 by Boyce Griffith
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

#ifndef included_LSetData
#define included_LSetData

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "Box.h"
#include "CellIterator.h"
#include "IndexData.h"
#include "IntVector.h"
#include "ibtk/LMarker.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSetDataIterator.h"

namespace IBTK
{
template <class T>
class LSet;
} // namespace IBTK
namespace SAMRAI
{
namespace pdat
{
template <int DIM>
class CellGeometry;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LSetData is a specialization of the templated class
 * SAMRAI::pdat::IndexData that provides access to Lagrangian objects that are
 * embedded in the a Cartesian grid patch.
 *
 * \see SAMRAI::pdat::IndexData
 */
template <class T>
class LSetData : public SAMRAI::pdat::IndexData<NDIM, LSet<T>, SAMRAI::pdat::CellGeometry<NDIM> >
{
public:
    /*!
     * This iterator iterates over the elements of a cell centered box geometry.
     */
    typedef SAMRAI::pdat::CellIterator<NDIM> CellIterator;

    /*!
     * This iterator iterates over the LSet elements within the IndexData patch
     * data object.
     */
    typedef SAMRAI::pdat::IndexIterator<NDIM, LSet<T>, SAMRAI::pdat::CellGeometry<NDIM> > SetIterator;

    /*!
     * This iterator iterates over the Lagrangian elements located within a cell
     * centered box geometry.
     */
    typedef IBTK::LSetDataIterator<T> DataIterator;

    /*!
     * Return an iterator to the first Lagrangian data object in the specified
     * region of index space.
     */
    DataIterator data_begin(const SAMRAI::hier::Box<NDIM>& box);

    /*!
     * Return an iterator pointing to the end of the collection of Lagrangian
     * data objects associated with the patch data object.
     */
    DataIterator data_end();

    /*!
     * The constructor for an SAMRAI::pdat::IndexData<NDIM> object.  The box
     * describes the interior of the index space and the ghosts vector describes
     * the ghost nodes in each coordinate direction.
     */
    LSetData(const SAMRAI::hier::Box<NDIM>& box, const SAMRAI::hier::IntVector<NDIM>& ghosts);

    /*!
     * The virtual destructor for an LSetData object.
     */
    virtual ~LSetData();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LSetData();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LSetData(const LSetData<T>& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LSetData& operator=(const LSetData<T>& that);
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LSetData-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LSetData
