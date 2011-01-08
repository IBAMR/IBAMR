// Filename: SideSynchCopyTransaction.h
// Created on 17 Dec 2009 by Boyce Griffith
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

#ifndef included_SideSynchCopyTransaction
#define included_SideSynchCopyTransaction

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <BoxOverlap.h>
#include <PatchLevel.h>
#include <RefineClasses.h>
#include <SideData.h>
#include <tbox/Pointer.h>
#include <tbox/Transaction.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class SideSynchCopyTransaction represents a single side data
 * communication transaction between two processors or a local data copy for
 * refine schedules.
 *
 * The difference between class SideSynchCopyTransaction and the standard class
 * SAMRAI::xfer::RefineCopyTransaction is that class SideSynchCopyTransaction
 * synchronizes values shared between two distinct patches.
 *
 * Note that to there is an implicit hand-shaking between objects of this class
 * and the SAMRAI::xfer::RefineSchedule<NDIM> object that constructs them.
 * Following the refine schedule implementation, the source patch data index for
 * a transaction always refers to the source data and the destination patch data
 * index for a transaction is always the scratch data, all as defined in the
 * SAMRAI::xfer::RefineClasses<NDIM> class.
 *
 * \see SAMRAI::xfer::RefineSchedule
 * \see SAMRAI::xfer::RefineClasses
 * \see SAMRAI::tbox::Schedule
 * \see SAMRAI::tbox::Transaction
 */
class SideSynchCopyTransaction
    : public SAMRAI::tbox::Transaction
{
public:
    /*!
     * Static member function to set the array of refine class data items that
     * is shared by all object instances of this sum transaction class during
     * data transfers.  The array must be set before any transactions are
     * executed.  The array is set in the RefineSchedule<NDIM> class.
     */
    static void
    setRefineItems(
        const SAMRAI::xfer::RefineClasses<NDIM>::Data** refine_items,
        int num_refine_items);

    /*!
     * Static member function to unset the array of refine class data items that
     * is shared by all object instances of this sum transaction class during
     * data transfers.  The unset function is used to prevent erroneous
     * execution of different schedules.  The array is unset in the
     * RefineSchedule<NDIM> class.
     */
    static void
    unsetRefineItems();

    /*!
     * \brief Class constructor.
     *
     * Construct a transaction with the specified source and destination levels,
     * patches, and patch data components found in the refine class item with
     * the given id owned by the calling refine schedule.  In general, this
     * constructor is called by a SAMRAI::xfer::RefineSchedule<NDIM> object for
     * each data transaction that must occur.  This transaction will be
     * responsible for one of the following: (1) a local data copy and sum, or
     * (2) packing a message stream with source patch data, or (3) unpacking and
     * summing destination patch data from a message stream.
     *
     * \param dst_level        SAMRAI::tbox::Pointer to destination patch level.
     * \param src_level        SAMRAI::tbox::Pointer to source patch level.
     * \param overlap          SAMRAI::tbox::Pointer to overlap region between patches.
     * \param dst_patch        Integer index of destination patch in destination patch level.
     * \param src_patch        Integer index of source patch in source patch level.
     * \param refine_item_id   Integer id of refine data item owned by refine schedule.
     *
     * When assertion checking is active, an assertion will result if any of the
     * pointer arguments is null, or if any of the integer arguments is invalid
     * (i.e., < 0).
     */
    SideSynchCopyTransaction(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > dst_level,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > src_level,
        SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM> > overlap,
        int dst_patch,
        int src_patch,
        int refine_item_id);

    /*!
     * \brief The virtual destructor for the copy transaction releases all
     * memory associated with the transaction.
     */
    virtual
    ~SideSynchCopyTransaction();

    /*!
     * \brief Return a boolean indicating whether this transaction can estimate
     * the size of an incoming message.
     *
     * If this evaluates to false, then a different communication protocol kicks
     * in and the message size is transmitted between sides.
     */
    virtual bool
    canEstimateIncomingMessageSize();

    /*!
     * \brief Return the integer buffer space (in bytes) needed for the incoming
     * message.
     *
     * This routine is only called if the transaction can estimate the size of
     * the incoming message.
     *
     * \see canEstimateIncomingMessageSize()
     */
    virtual int
    computeIncomingMessageSize();

    /*!
     * \brief Return the integer buffer space (in bytes) needed for the outgoing
     * message.
     */
    virtual int
    computeOutgoingMessageSize();

    /*!
     * \brief Return the sending processor number for the communications
     * transaction.
     */
    virtual int
    getSourceProcessor();

    /*!
     * \brief Return the receiving processor number for the communications
     * transaction.
     */
    virtual int
    getDestinationProcessor();

    /*!
     * \brief Pack the transaction data into the message stream.
     */
    virtual void
    packStream(
        SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Unpack the transaction data from the message stream.
     */
    virtual void
    unpackStream(
        SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Perform the local data copy for the transaction.
     */
    virtual void
    copyLocalData();

    /*!
     * \brief Print out transaction information.
     */
    virtual void
    printClassData(
        std::ostream& stream) const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SideSynchCopyTransaction();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SideSynchCopyTransaction(
        const SideSynchCopyTransaction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    void
    operator=(
        const SideSynchCopyTransaction& that);

    /*!
     * \brief Synchronize values along the boundary of the destination patch
     * data object.
     */
    void
    synchronizeData(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > dst_data,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > src_data);

    static const SAMRAI::xfer::RefineClasses<NDIM>::Data** s_refine_items;
    static int s_num_refine_items;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > d_dst_level;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > d_src_level;
    SAMRAI::tbox::Pointer<SAMRAI::hier::BoxOverlap<NDIM> > d_overlap;
    int d_dst_patch, d_src_patch;
    int d_refine_item_id;
    int d_incoming_bytes, d_outgoing_bytes;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/SideSynchCopyTransaction.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SideSynchCopyTransaction
