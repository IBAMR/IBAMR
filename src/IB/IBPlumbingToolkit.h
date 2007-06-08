#ifndef included_IBPlumbingToolkit
#define included_IBPlumbingToolkit

// Filename: IBPlumbingToolkit.h
// Last modified: <12.May.2007 13:57:51 boyce@trasnaform2.local>
// Created on 12 May 2007 by Boyce Griffith (boyce@trasnaform2.local)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES

// SAMRAI INCLUDES
#include <tbox/Database.h>
#include <tbox/DescribedClass.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBPlumbingToolkit is a utility class that provides support for
 * flow meters and pressure taps.
 */
class IBPlumbingToolkit
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBPlumbingToolkit(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBPlumbingToolkit();

    /*!
     * \brief Initialize hierarchy-dependent data.
     */
    void
    initializeHierarchyData();

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBPlumbingToolkit(
        const IBPlumbingToolkit& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBPlumbingToolkit&
    operator=(
        const IBPlumbingToolkit& that);

    /*!
     * \brief Initialize the data associated with a specific flow meter.
     */
    void
    initializeFlowMeter(
        const int meter_number);

    /*!
     * \brief Initialize the data associated with a specific pressure tap.
     */
    void
    initializePressureTap(
        const int tap_number);

    /*!
     * \brief Hierarchy independent data.
     */
    std::string d_interp_type;
    int d_num_flow_meters, d_num_pressure_taps;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBPlumbingToolkit.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBPlumbingToolkit
