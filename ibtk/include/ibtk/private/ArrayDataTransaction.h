// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#ifndef included_IBTK_ArrayDataTransaction
#define included_IBTK_ArrayDataTransaction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>
#include <tbox/Transaction.h>

#include <ArrayData.h>
#include <Box.h>
#include <BoxList.h>
#include <IntVector.h>
#include <VariableFillPattern.h>

#include <functional>
#include <iosfwd>
#include <vector>

namespace SAMRAI
{
namespace tbox
{
class AbstractStream;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Class ArrayDataTransaction is a concrete implementation of the abstract base class SAMRAI::tbox::Transaction.
 * It copies the values of a box of double precision array data from one process to another: the destination array
 * gets the values of the source array at the box shifted by -src_shift. Functions return the arrays when the
 * transaction runs, so that the transaction can copy data whose patch data index changes.
 */
class ArrayDataTransaction : public SAMRAI::tbox::Transaction
{
public:
    using ArrayAccessor = std::function<SAMRAI::pdat::ArrayData<NDIM, double>&()>;

    /*!
     * \brief Constructor. Only the process src_proc calls src_array, and only the process dst_proc calls dst_array.
     */
    ArrayDataTransaction(int src_proc,
                         int dst_proc,
                         const SAMRAI::hier::Box<NDIM>& dst_box,
                         const SAMRAI::hier::IntVector<NDIM>& src_shift,
                         int depth,
                         ArrayAccessor src_array,
                         ArrayAccessor dst_array);

    /*!
     * Return a boolean indicating whether this transaction can estimate the
     * size of an incoming message.
     */
    bool canEstimateIncomingMessageSize() override;

    /*!
     * Return the amount of buffer space needed for the incoming message.
     */
    int computeIncomingMessageSize() override;

    /*!
     * Return the buffer space needed for the outgoing message.
     */
    int computeOutgoingMessageSize() override;

    /*!
     * Return the sending processor for the communications transaction.
     */
    int getSourceProcessor() override;

    /*!
     * Return the receiving processor for the communications transaction.
     */
    int getDestinationProcessor() override;

    /*!
     * Pack the transaction data into the message stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream) override;

    /*!
     * Unpack the transaction data from the message stream.
     */
    void unpackStream(SAMRAI::tbox::AbstractStream& stream) override;

    /*!
     * Perform the local data copy for the transaction.
     */
    void copyLocalData() override;

    /*!
     * Print out transaction information.
     */
    void printClassData(std::ostream& stream) const override;

private:
    ArrayDataTransaction() = delete;
    ArrayDataTransaction(const ArrayDataTransaction& from) = delete;
    ArrayDataTransaction& operator=(const ArrayDataTransaction& that) = delete;

    const int d_src_proc, d_dst_proc;
    const SAMRAI::hier::Box<NDIM> d_dst_box;
    const SAMRAI::hier::IntVector<NDIM> d_src_shift;
    const int d_depth;
    const ArrayAccessor d_src_array, d_dst_array;
};

/*!
 * \brief The indices of each array of the data with index data_idx (one for cell data, NDIM for side data) that a
 * refine schedule to level ln with the fill pattern fill_pattern (null for the default pattern) can set on a patch
 * with the box patch_box, including the interior.
 */
/*!
 * \brief The properties of a fill pattern that determine fill_pattern_indices() for data with a given ghost width.
 */
struct FillPatternShape
{
    enum class Kind
    {
        FULL,
        CELL_NO_CORNERS,
        SIDE_NO_CORNERS,
        OTHER
    };

    //! The kind of pattern; FULL for null and for SAMRAI::xfer::BoxGeometryFillPattern.
    Kind kind = Kind::FULL;
    //! The fill width: the ghost width, limited by the stencil width of a NoCorners pattern.
    SAMRAI::hier::IntVector<NDIM> width = SAMRAI::hier::IntVector<NDIM>(0);
    //! Whether the indices include the interior of the patch.
    bool overwrite_interior = true;
    //! The pattern itself if its kind is OTHER, otherwise null.
    const SAMRAI::xfer::VariableFillPattern<NDIM>* other = nullptr;

    bool operator==(const FillPatternShape& that) const;
    bool operator!=(const FillPatternShape& that) const;
};

/*!
 * \brief The shape of the fill pattern fill_pattern (null for the default pattern) for data with the ghost width
 * ghost_width. Patterns are classified by exact type, so a subclass of a NoCorners pattern has the kind OTHER.
 */
FillPatternShape fill_pattern_shape(SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern,
                                    const SAMRAI::hier::IntVector<NDIM>& ghost_width);

std::vector<SAMRAI::hier::BoxList<NDIM>>
fill_pattern_indices(int data_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern,
                     const SAMRAI::hier::Box<NDIM>& patch_box,
                     int ln);
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_ArrayDataTransaction
