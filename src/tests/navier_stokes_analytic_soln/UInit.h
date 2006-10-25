#ifndef included_UInit
#define included_UInit

// Filename: UInit.h
// Last modified: <04.Oct.2006 13:34:56 boyce@boyce-griffiths-powerbook-g4-15.local>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <tbox/Database.h>
#include <GridGeometry.h>
#include <Patch.h>
#include <tbox/Pointer.h>
#include <Variable.h>
#include <tbox/Array.h>

// NAMESPACE
using namespace IBAMR;
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the fluid velocity U.
 */
class UInit
    : public SetDataStrategy
{
public:
    /*!
     * \brief Default constructor.
     */
    UInit(
        const string& object_name,
        tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
        tbox::Pointer<tbox::Database> input_db,
        const double nu);

    /*!
     * \brief Destructor.
     */
    ~UInit();

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    bool isTimeDependent() const { return true; }

    /*!
     * Set the data on the patch interior to the exact answer.
     */
    void setDataOnPatch(
        const int data_idx,
        tbox::Pointer<hier::Variable<NDIM> > var,
        hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);

protected:

private:
    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    UInit(
        const UInit& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    UInit& operator=(
        const UInit& that);

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
     * The viscosity.
     */
    double d_nu;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "UInit.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_UInit
