// ---------------------------------------------------------------------
//
// Copyright (c) 2007 - 2022 by the IBAMR developers
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

#ifndef included_IBAMR_IBInstrumentationSpec
#define included_IBAMR_IBInstrumentationSpec

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/Streamable.h"
#include "ibtk/StreamableFactory.h"

#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class IntVector;
} // namespace hier
namespace tbox
{
class AbstractStream;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBInstrumentationSpec encapsulates the data required to
 * initialize distributed internal flow meters and pressure gauges.
 */
class IBInstrumentationSpec : public IBTK::Streamable
{
public:
    /*!
     * \brief Register this class and its factory class with the singleton
     * IBTK::StreamableManager object.  This method must be called before any
     * IBInstrumentationSpec objects are created.
     *
     * \note This method is collective on all MPI processes.  This is done to
     * ensure that all processes employ the same class ID for the
     * IBInstrumentationSpec class.
     */
    static void registerWithStreamableManager();

    /*!
     * \brief Returns a boolean indicating whether the class has been registered
     * with the singleton IBTK::StreamableManager object.
     */
    static bool getIsRegisteredWithStreamableManager();

    /*!
     * The unique class ID for this object type assigned by the
     * IBTK::StreamableManager.
     */
    static int STREAMABLE_CLASS_ID;

    /*!
     * \brief Set the names of the flow meters and pressure gauges.
     */
    static void setInstrumentNames(const std::vector<std::string>& names);

    /*!
     * \brief Get the names of the flow meters and pressure gauges.
     */
    static const std::vector<std::string>& getInstrumentNames();

    /*!
     * \brief Default constructor.
     */
    IBInstrumentationSpec(int master_idx = IBTK::invalid_index,
                          int meter_idx = IBTK::invalid_index,
                          int node_idx = IBTK::invalid_index);

    /*!
     * \brief Destructor.
     */
    ~IBInstrumentationSpec();

    /*!
     * \return A const reference to the master node index.
     */
    const int& getMasterNodeIndex() const;

    /*!
     * \return A non-const reference to the master node index.
     */
    int& getMasterNodeIndex();

    /*!
     * \return A const reference to the meter index associated with the master
     * node.
     */
    const int& getMeterIndex() const;

    /*!
     * \return A non-const reference to the meter index associated with the
     * master node.
     */
    int& getMeterIndex();

    /*!
     * \return A const reference to the node index associated with the master
     * node.
     */
    const int& getNodeIndex() const;

    /*!
     * \return A non-const reference to the node index associated with the master
     * node.
     */
    int& getNodeIndex();

    /*!
     * \brief Return the unique identifier used to specify the
     * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
     * extract Streamable objects from data streams.
     */
    int getStreamableClassID() const override;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t getDataStreamSize() const override;

    /*!
     * \brief Pack data into the output stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream) override;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBInstrumentationSpec(const IBInstrumentationSpec& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBInstrumentationSpec& operator=(const IBInstrumentationSpec& that) = delete;

    /*!
     * The names of the instrument names.
     */
    static std::vector<std::string> s_instrument_names;

    /*!
     * Data required to define the instrument.
     */
    int d_master_idx, d_meter_idx, d_node_idx;

    /*!
     * \brief A factory class to rebuild IBInstrumentationSpec objects from
     * SAMRAI::tbox::AbstractStream data streams.
     */
    class Factory : public IBTK::StreamableFactory
    {
    public:
        /*!
         * \brief Destructor.
         */
        ~Factory() = default;

        /*!
         * \brief Return the unique identifier used to specify the
         * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
         * extract IBInstrumentationSpec objects from data streams.
         */
        int getStreamableClassID() const override;

        /*!
         * \brief Set the unique identifier used to specify the
         * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
         * extract IBInstrumentationSpec objects from data streams.
         */
        void setStreamableClassID(int class_id) override;

        /*!
         * \brief Build an IBInstrumentationSpec object by unpacking data from
         * the data stream.
         */
        SAMRAI::tbox::Pointer<IBTK::Streamable> unpackStream(SAMRAI::tbox::AbstractStream& stream,
                                                             const SAMRAI::hier::IntVector<NDIM>& offset) override;

    private:
        /*!
         * \brief Default constructor.
         */
        Factory();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        Factory(const Factory& from) = delete;

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        Factory& operator=(const Factory& that) = delete;

        friend class IBInstrumentationSpec;
    };
    using IBInstrumentationSpecFactory = IBInstrumentationSpec::Factory;
};
} // namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibamr/private/IBInstrumentationSpec-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBInstrumentationSpec
