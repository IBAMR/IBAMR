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

#ifndef included_IBTK_LSetData
#define included_IBTK_LSetData

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LMarker.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSetDataIterator.h"

#include "Box.h"
#include "CellIterator.h"
#include "IndexData.h"
#include "IntVector.h"

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
    using CellIterator = SAMRAI::pdat::CellIterator<NDIM>;

    /*!
     * This iterator iterates over the LSet elements within the IndexData patch
     * data object.
     */
    using SetIterator = SAMRAI::pdat::IndexIterator<NDIM, LSet<T>, SAMRAI::pdat::CellGeometry<NDIM> >;

    /*!
     * This iterator iterates over the Lagrangian elements located within a cell
     * centered box geometry.
     */
    using DataIterator = IBTK::LSetDataIterator<T>;

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
    LSetData(SAMRAI::hier::Box<NDIM> box, SAMRAI::hier::IntVector<NDIM> ghosts);

    /*!
     * The virtual destructor for an LSetData object.
     */
    virtual ~LSetData() = default;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LSetData() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LSetData(const LSetData<T>& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LSetData& operator=(const LSetData<T>& that) = delete;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LSetData-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LSetData
