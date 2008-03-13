#ifndef included_FeedbackForcer
#define included_FeedbackForcer

// Filename: FeedbackForcer.h
// Last modified: <19.Oct.2007 00:41:02 griffith@box221.cims.nyu.edu>
// Created on 19 Oct 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <PatchHierarchy.h>

// NAMESPACE
using namespace SAMRAI;
using namespace IBTK;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class FeedbackForcer is an implementation of the strategy class
 * SetDataStrategy that is used to specify velocity boundary conditions via a
 * feedback forcing (i.e., penalty) method.
 */
class FeedbackForcer
    : public SetDataStrategy
{
public:
    /*!
     * \brief Constructor
     */
    FeedbackForcer(
        const string& object_name,
        const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry);

    /*!
     * \brief Destructor.
     */
    virtual
    ~FeedbackForcer();

    /*!
     * Velocity patch data descriptor index.
     */
    int d_U_data_idx;

    /*!
     * Feedback forcing parameter.
     */
    double d_kappa;

    /*!
     * \name Implementation of SetDataStrategy interface.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete SetDataStrategy object is
     * time-dependent.
     */
    virtual bool
    isTimeDependent() const { return true; }

    /*!
     * \brief Set data on the specified patch interior.
     */
    virtual void
    setDataOnPatch(
        const int data_idx,
        tbox::Pointer<hier::Variable<NDIM> > var,
        hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FeedbackForcer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FeedbackForcer(
        const FeedbackForcer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FeedbackForcer&
    operator=(
        const FeedbackForcer& that);

    /*
     * The object name.
     */
    string d_object_name;

    /*
     * The grid geometry object.
     */
    const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <FeedbackForcer.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FeedbackForcer
