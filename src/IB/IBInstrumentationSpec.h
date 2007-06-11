#ifndef included_IBInstrumentationSpec
#define included_IBInstrumentationSpec

// Filename: IBInstrumentationSpec.h
// Last modified: <11.Jun.2007 17:05:06 griffith@box221.cims.nyu.edu>
// Created on 11 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/Stashable.h>

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
    : public Stashable
{
public:
    /*!
     * \brief Register this class and its factory class with the singleton
     * StashableManager object.  This method must be called before any
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
     * with the singleton StashableManager object.
     */
    static bool
    getIsRegisteredWithStashableManager();

    /*!
     * \brief Default constructor.
     */
    IBInstrumentationSpec(
        const int master_idx=-1,
        const int flow_meter_idx=-1,
        const int pressure_gauge_idx=-1);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBInstrumentationSpec();

    /*!
     * \return A const refernce to the master node index.
     */
    const int&
    getMasterNodeIndex() const;

    /*!
     * \return A non-const reference to the master node index.
     */
    int&
    getMasterNodeIndex();

    /*!
     * \return A const refrence to the flow meter index associated with the
     * master node.
     */
    const int&
    getFlowMeterIndex() const;

    /*!
     * \return A non-const refrence to the flow meter index associated with the
     * master node.
     */
    int&
    getFlowMeterIndex();

    /*!
     * \return A const refrence to the pressure gauge index associated with the
     * master node.
     */
    const int&
    getPressureGaugeIndex() const;

    /*!
     * \return A non-const refrence to the pressure gauge index associated with
     * the master node.
     */
    int&
    getPressureGaugeIndex();

    /*!
     * \brief Return the unique identifier used to specify the StashableFactory
     * object used by the StashableManager to extract Stashable objects from
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
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBInstrumentationSpec();

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
     * StashableManager.
     */
    static bool s_registered_factory;

    /*!
     * The stashable ID for this object type.
     */
    static int s_stashable_id;

    /*!
     * Data required to define the instrument.
     */
    int d_master_idx, d_flow_meter_idx, d_pressure_gauge_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include "IBInstrumentationSpec.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBInstrumentationSpec
