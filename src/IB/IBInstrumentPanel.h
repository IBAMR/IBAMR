#ifndef included_IBPlumbingToolkit
#define included_IBPlumbingToolkit

// Filename: IBPlumbingToolkit.h
// Last modified: <11.Jun.2007 16:56:38 griffith@box221.cims.nyu.edu>
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
 * \brief Class IBPlumbingToolkit provides support for flow meters and pressure
 * gauges.
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
    initializeHierarchyData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const double init_data_time,
        const bool initial_time,
        LDataManager* const lag_manager);

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
     * \brief Initialize the data associated with a specific pressure gauge.
     */
    void
    initializePressureGauge(
        const int gauge_number);

    /*!
     * \brief Hierarchy independent data.
     */
    std::string d_interp_type;
    int d_num_flow_meters, d_num_pressure_gauges;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBPlumbingToolkit.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBPlumbingToolkit
