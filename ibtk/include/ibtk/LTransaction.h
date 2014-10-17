// Filename: LTransaction.h
// Created on 03 Mar 2010 by Boyce Griffith
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

#ifndef included_LTransaction
#define included_LTransaction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <iosfwd>
#include <vector>

#include "ibtk/LMarker.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSet.h"
#include "ibtk/ibtk_utilities.h"
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
 * \brief Class LTransaction represents a communication transaction between two
 * processors or a local data copy for communicating or copying Lagrangian
 * objects.
 *
 * \see SAMRAI::tbox::Schedule
 * \see SAMRAI::tbox::Transaction
 */
template <class T>
class LTransaction : public SAMRAI::tbox::Transaction
{
public:
    /*!
     * \brief Struct LTransaction::LTransactionComponent encapsulates the
     * individual data items that are communicated via class LTransaction.
     */
    class LTransactionComponent
    {
    public:
        /*!
         * \brief Default constructor.
         */
        inline LTransactionComponent(const typename LSet<T>::value_type& item = NULL, const Point& posn = Point::Zero())
            : item(item), posn(posn)
        {
            // intentionally blank
            return;
        } // LTransactionComponent

        /*!
         * \brief Copy constructor.
         *
         * \param from The value to copy to this object.
         */
        inline LTransactionComponent(const LTransactionComponent& from) : item(from.item), posn(from.posn)
        {
            // intentionally blank
            return;
        } // LTransactionComponent

        /*!
         * \brief Assignment operator.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        inline LTransactionComponent& operator=(const LTransactionComponent& that)
        {
            if (this != &that)
            {
                item = that.item;
                posn = that.posn;
            }
            return *this;
        } // operator=

        /*!
         * \brief Destructor.
         */
        inline ~LTransactionComponent()
        {
            // intentionally blank
            return;
        } // ~LTransactionComponent

        // Data.
        typename LSet<T>::value_type item;
        Point posn;
    };

    /*!
     * \brief Class constructor.
     */
    LTransaction(int src_proc, int dst_proc);

    /*!
     * \brief Class constructor.
     */
    LTransaction(int src_proc, int dst_proc, const std::vector<LTransactionComponent>& src_item_set);

    /*!
     * \brief The virtual destructor for the copy transaction releases all
     * memory associated with the transaction.
     */
    virtual ~LTransaction();

    /*!
     * \brief Return a constant reference to the source data.
     */
    inline const std::vector<LTransactionComponent>& getSourceData() const
    {
        return d_src_item_set;
    }

    /*!
     * \brief Return a constant reference to the destination data.
     */
    inline const std::vector<LTransactionComponent>& getDestinationData() const
    {
        return d_dst_item_set;
    }

    /*!
     * \brief Return a boolean indicating whether this transaction can estimate
     * the size of an incoming message.
     *
     * If this evaluates to false, then a different communication protocol kicks
     * in and the message size is transmitted between sides.
     */
    virtual bool canEstimateIncomingMessageSize();

    /*!
     * \brief Return the integer buffer space (in bytes) needed for the incoming
     * message.
     *
     * This routine is only called if the transaction can estimate the size of
     * the incoming message.
     *
     * \see canEstimateIncomingMessageSize()
     */
    virtual int computeIncomingMessageSize();

    /*!
     * \brief Return the integer buffer space (in bytes) needed for the outgoing
     * message.
     */
    virtual int computeOutgoingMessageSize();

    /*!
     * \brief Return the sending processor number for the communications
     * transaction.
     */
    virtual int getSourceProcessor();

    /*!
     * \brief Return the receiving processor number for the communications
     * transaction.
     */
    virtual int getDestinationProcessor();

    /*!
     * \brief Pack the transaction data into the message stream.
     */
    virtual void packStream(SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Unpack the transaction data from the message stream.
     */
    virtual void unpackStream(SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Perform the local data copy for the transaction.
     */
    virtual void copyLocalData();

    /*!
     * \brief Print out transaction information.
     */
    virtual void printClassData(std::ostream& stream) const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LTransaction();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LTransaction(const LTransaction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    void operator=(const LTransaction& that);

    std::vector<LTransactionComponent> d_src_item_set;
    int d_src_proc;
    int d_outgoing_bytes;

    std::vector<LTransactionComponent> d_dst_item_set;
    int d_dst_proc;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LTransaction
