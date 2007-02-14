#ifndef included_USet
#define included_USet

// Filename: USet.h
// Last modified: <13.Feb.2007 15:16:52 boyce@bigboy.nyconnect.com>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <GridGeometry.h>
#include <tbox/Array.h>
#include <tbox/Database.h>

// NAMESPACE
using namespace SAMRAI;
using namespace STOOLS;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the advection velocity u.
 */
class USet
    : public SetDataStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    USet(
        const string& object_name,
        tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
        tbox::Pointer<tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~USet();

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    virtual bool isTimeDependent() const { return false; }

    /*!
     * Set the data on the patch interior to some initial values.
     */
    virtual void setDataOnPatch(
        const int data_idx,
        tbox::Pointer<hier::Variable<NDIM> > var,
        hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    USet();

    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    USet(
        const USet& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    USet& operator=(
        const USet& that);
    /*!
     * Read input values, indicated above, from given database.
     */
    void getFromInput(
        tbox::Pointer<tbox::Database> db);

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
     * The center of the initial data.
     */
    tbox::Array<double> d_X;

    /*
     * The initialization type.
     */
    string d_init_type;

    /*
     * Parameters for uniform constant velocity.
     */
    tbox::Array<double> d_uniform_u;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "USet.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_USet
