#ifndef included_IBInstrumentPanel
#define included_IBInstrumentPanel

// Filename: IBInstrumentPanel.h
// Last modified: <11.Jun.2007 19:37:48 griffith@box221.cims.nyu.edu>
// Created on 12 May 2007 by Boyce Griffith (boyce@trasnaform2.local)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LDataManager.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <tbox/Database.h>
#include <tbox/DescribedClass.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBInstrumentPanel provides support for flow meters and pressure
 * gauges.
 */
class IBInstrumentPanel
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBInstrumentPanel(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBInstrumentPanel();

    /*!
     * \brief Compute the positions of all of the flow meters and pressure
     * gauges.
     */
    void
    computeMeterPositions(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const double init_data_time,
        const bool initial_time,
        LDataManager* const lag_manager);

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
    IBInstrumentPanel(
        const IBInstrumentPanel& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBInstrumentPanel&
    operator=(
        const IBInstrumentPanel& that);

    /*!
     * \brief Hierarchy independent data.
     */
    int d_interp_level;
    int d_num_meters;
    std::vector<int> d_num_meter_nodes;

    /*!
     * \brief Hierarchy and configuration dependent data.
     */
    std::vector<std::vector<double> > d_X_perimeter, d_X_centroid;
    std::vector<std::vector<std::vector<double> > > d_X_web;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBInstrumentPanel.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBInstrumentPanel
