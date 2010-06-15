#ifndef included_IBInstrumentationSpec
#define included_IBInstrumentationSpec

// Filename: IBInstrumentationSpec.h
// Last modified: <15.Jun.2010 15:40:39 griffith@boyce-griffiths-mac-pro.local>
// Created on 11 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/Stashable.h>

// SAMRAI INCLUDES
#include <tbox/AbstractStream.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBInstrumentationSpec encapsulates the data required to
 * initialize distributed internal flow meters and pressure gauges.
 */
class IBInstrumentationSpec
    : public IBTK::Stashable
{
public:
    /*!
     * \brief Register this class and its factory class with the singleton
     * IBTK::StashableManager object.  This method must be called before any
     * IBInstrumentationSpec objects are created.
     *
     * \note This method is collective on all MPI processes.  This is done to
     * ensure that all processes employ the same stashable ID for the
     * IBInstrumentationSpec class.
     */
    static void
    registerWithStashableManager();

    /*!
     * \brief Returns a boolean indicating whether the class has been registered
     * with the singleton IBTK::StashableManager object.
     */
    static bool
    getIsRegisteredWithStashableManager();

    /*!
     * \brief Set the names of the flow meters and pressure gauges.
     */
    static void
    setInstrumentNames(
        const std::vector<std::string>& names);

    /*!
     * \brief Get the names of the flow meters and pressure gauges.
     */
    static const std::vector<std::string>&
    getInstrumentNames();

    /*!
     * \brief Default constructor.
     */
    IBInstrumentationSpec(
        const int master_idx=-1,
        const int meter_idx=-1,
        const int node_idx=-1);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBInstrumentationSpec();

    /*!
     * \return A const reference to the master node index.
     */
    const int&
    getMasterNodeIndex() const;

    /*!
     * \return A non-const reference to the master node index.
     */
    int&
    getMasterNodeIndex();

    /*!
     * \return A const reference to the meter index associated with the master
     * node.
     */
    const int&
    getMeterIndex() const;

    /*!
     * \return A non-const reference to the meter index associated with the
     * master node.
     */
    int&
    getMeterIndex();

    /*!
     * \return A const reference to the node index associated with the master
     * node.
     */
    const int&
    getNodeIndex() const;

    /*!
     * \return A non-const reference to the node index associated with the master
     * node.
     */
    int&
    getNodeIndex();

    /*!
     * \brief Return the unique identifier used to specify the IBTK::StashableFactory
     * object used by the IBTK::StashableManager to extract Stashable objects from
     * data streams.
     */
    virtual int
    getStashableID() const;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    virtual size_t
    getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    virtual void
    packStream(
        SAMRAI::tbox::AbstractStream& stream);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBInstrumentationSpec(
        const IBInstrumentationSpec& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBInstrumentationSpec&
    operator=(
        const IBInstrumentationSpec& that);

    /*!
     * Indicates whether the factory has been registered with the
     * IBTK::StashableManager.
     */
    static bool s_registered_factory;

    /*!
     * The stashable ID for this object type.
     */
    static int s_stashable_id;

    /*!
     * The names of the instrument names.
     */
    static std::vector<std::string> s_instrument_names;

    /*!
     * Data required to define the instrument.
     */
    int d_master_idx, d_meter_idx, d_node_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include "IBInstrumentationSpec.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBInstrumentationSpec
