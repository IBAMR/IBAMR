// Filename: CopyToRootTransaction.h
// Created on 04 May 2011 by Boyce Griffith
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

#ifndef included_CopyToRootTransaction
#define included_CopyToRootTransaction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <iosfwd>

#include "IntVector.h"
#include "PatchData.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Transaction.h"

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
    ~CopyToRootTransaction();

    /*!
     * Return a pointer to the data on the root process.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > getRootPatchData() const;

    /*!
     * Return a boolean indicating whether this transaction can estimate the
     * size of an incoming message.
     */
    bool canEstimateIncomingMessageSize();

    /*!
     * Return the amount of buffer space needed for the incoming message.
     * This routine is only called if the transaction can estimate the
     * size of the incoming message.
     */
    int computeIncomingMessageSize();

    /*!
     * Return the buffer space needed for the outgoing message.
     */
    int computeOutgoingMessageSize();

    /*!
     * Return the sending processor for the communications transaction.
     */
    int getSourceProcessor();

    /*!
     * Return the receiving processor for the communications transaction.
     */
    int getDestinationProcessor();

    /*!
     * Pack the transaction data into the message stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream);

    /*!
     * Unpack the transaction data from the message stream.
     */
    void unpackStream(SAMRAI::tbox::AbstractStream& stream);

    /*!
     * Perform the local data copy for the transaction.
     */
    void copyLocalData();

    /*!
     * Print out transaction information.
     */
    void printClassData(std::ostream& stream) const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CopyToRootTransaction();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CopyToRootTransaction(const CopyToRootTransaction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CopyToRootTransaction& operator=(const CopyToRootTransaction& that);

    const int d_src_proc, d_dst_proc;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > d_patch_level;
    const int d_src_patch_data_idx;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > d_dst_patch_data;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CopyToRootTransaction
