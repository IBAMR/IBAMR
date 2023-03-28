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

#ifndef included_IBTK_LTransaction
#define included_IBTK_LTransaction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LMarker.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSet.h"
#include "ibtk/ibtk_utilities.h"

#include "tbox/Transaction.h"

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
        inline LTransactionComponent(const typename LSet<T>::value_type& item = nullptr,
                                     const Point& posn = Point::Zero())
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
    LTransaction(int src_proc, int dst_proc, std::vector<LTransactionComponent> src_item_set);

    /*!
     * \brief The virtual destructor for the copy transaction releases all
     * memory associated with the transaction.
     */
    virtual ~LTransaction() = default;

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
    virtual bool canEstimateIncomingMessageSize() override;

    /*!
     * \brief Return the integer buffer space (in bytes) needed for the incoming
     * message.
     *
     * This routine is only called if the transaction can estimate the size of
     * the incoming message.
     *
     * \see canEstimateIncomingMessageSize()
     */
    virtual int computeIncomingMessageSize() override;

    /*!
     * \brief Return the integer buffer space (in bytes) needed for the outgoing
     * message.
     */
    virtual int computeOutgoingMessageSize() override;

    /*!
     * \brief Return the sending processor number for the communications
     * transaction.
     */
    virtual int getSourceProcessor() override;

    /*!
     * \brief Return the receiving processor number for the communications
     * transaction.
     */
    virtual int getDestinationProcessor() override;

    /*!
     * \brief Pack the transaction data into the message stream.
     */
    virtual void packStream(SAMRAI::tbox::AbstractStream& stream) override;

    /*!
     * \brief Unpack the transaction data from the message stream.
     */
    virtual void unpackStream(SAMRAI::tbox::AbstractStream& stream) override;

    /*!
     * \brief Perform the local data copy for the transaction.
     */
    virtual void copyLocalData() override;

    /*!
     * \brief Print out transaction information.
     */
    virtual void printClassData(std::ostream& stream) const override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LTransaction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LTransaction(const LTransaction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    void operator=(const LTransaction& that) = delete;

    std::vector<LTransactionComponent> d_src_item_set;
    int d_src_proc;
    int d_outgoing_bytes = 0;

    std::vector<LTransactionComponent> d_dst_item_set;
    int d_dst_proc;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LTransaction
