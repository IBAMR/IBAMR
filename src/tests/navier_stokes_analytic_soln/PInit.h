#ifndef included_PInit
#define included_PInit

// Filename: PInit.h
// Last modified: <26.Sep.2006 11:04:52 boyce@boyce-griffiths-powerbook-g4-15.local>
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

// C++ STDLIB INCLUDES
#include <string>

// NAMESPACE
using namespace IBAMR;
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the hydrostatic pressure P
 * (used for convergence analysis only).
 */
class PInit
    : public SetDataStrategy
{
public:
    /*!
     * \brief Default constructor.
     */
    PInit(
        const string& object_name,
        tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
        tbox::Pointer<tbox::Database> input_db,
        const double nu);

    /*!
     * \brief Destructor.
     */
    ~PInit();

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
    PInit(
        const PInit& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PInit& operator=(
        const PInit& that);

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

//#include "PInit.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PInit
