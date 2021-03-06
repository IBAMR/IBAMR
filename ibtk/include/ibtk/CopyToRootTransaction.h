// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_CopyToRootTransaction
#define included_IBTK_CopyToRootTransaction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "IntVector.h"
#include "PatchData.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Transaction.h"

#include <iosfwd>
#include <string>

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
/*!
 * \brief Class CopyToRootTransaction is a concrete implementation of the
 * abstract base class SAMRAI::tbox::Transaction.  It is used to communicate
 * distributed patch data to a unified patch data object on a root MPI process.
 *
 * \note This class is designed to be used with uniform grid data only.
 */
class CopyToRootTransaction : public SAMRAI::tbox::Transaction
{
public:
    /*!
     * \brief Constructor
     */
    CopyToRootTransaction(int src_proc,
                          int dst_proc,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                          int src_patch_data_idx,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > dst_patch_data);

    /*!
     * \brief Destructor
     */
    ~CopyToRootTransaction() = default;

    /*!
     * Return a pointer to the data on the root process.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > getRootPatchData() const;

    /*!
     * Return a boolean indicating whether this transaction can estimate the
     * size of an incoming message.
     */
    bool canEstimateIncomingMessageSize() override;

    /*!
     * Return the amount of buffer space needed for the incoming message.
     * This routine is only called if the transaction can estimate the
     * size of the incoming message.
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
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CopyToRootTransaction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CopyToRootTransaction(const CopyToRootTransaction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CopyToRootTransaction& operator=(const CopyToRootTransaction& that) = delete;

    const int d_src_proc, d_dst_proc;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > d_patch_level;
    const int d_src_patch_data_idx;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > d_dst_patch_data;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CopyToRootTransaction
