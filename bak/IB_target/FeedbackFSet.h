#ifndef included_FeedbackFSet
#define included_FeedbackFSet

// Filename: FeedbackFSet.h
// Last modified: <16.Feb.2007 16:32:58 griffith@box221.cims.nyu.edu>
// Created on 20 Nov 2006 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <CartesianGridGeometry.h>
#include <GridGeometry.h>
#include <VariableContext.h>
#include <tbox/Array.h>
#include <tbox/Database.h>

// NAMESPACE
using namespace SAMRAI;
using namespace STOOLS;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class to specify a feedback force that attempts to ensure
 * that the inflow velocity is U = (1,0).
 */
class FeedbackFSet
    : public SetDataStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    FeedbackFSet(
        const string& object_name,
        tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
        tbox::Pointer<tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~FeedbackFSet();

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    virtual bool isTimeDependent() const { return true; }

    /*!
     * Set the data on the patch interior to some values.
     */
    virtual void setDataOnPatch(
        const int data_idx,
        tbox::Pointer<hier::Variable<NDIM> > var,
        hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);

    /*!
     * The fluid velocity variable.
     */
    tbox::Pointer<pdat::CellVariable<NDIM,double> > d_U_var;

    /*!
     * The fluid velocity variable context.
     */
    tbox::Pointer<hier::VariableContext> d_U_context;
protected:

private:
    /*!
     * \brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    FeedbackFSet();

    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    FeedbackFSet(
        const FeedbackFSet& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FeedbackFSet& operator=(
        const FeedbackFSet& that);

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.
     */
    string d_object_name;

    /*
     * The grid geometry.
     */
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*
     * Feedback force parameter.
     */
    double d_kappa, d_width0, d_width1;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "FeedbackFSet.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FeedbackFSet
