#ifndef included_GravitationalBodyForce
#define included_GravitationalBodyForce

// Filename: GravitationalBodyForce.h
// Last modified: <12.Mar.2008 23:23:54 griffith@box221.cims.nyu.edu>
// Created on 03 May 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <tbox/Database.h>

// C++ STDLIB INCLUDES
#include <vector>

// NAMESPACE
using namespace IBTK;
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class to specify a uniform body force that arises from graviational
 * forces acting on an incompressible fluid with uniform mass density.
 */
class GravitationalBodyForce
    : public SetDataStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    GravitationalBodyForce(
        const string& object_name,
        tbox::Pointer<tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual
    ~GravitationalBodyForce();

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    virtual bool
    isTimeDependent() const
        { return true; }

    /*!
     * Set the data on the patch interior to some values.
     */
    virtual void
    setDataOnPatch(
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
     * \note This constructor is not implemented and should not be used.
     */
    GravitationalBodyForce();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    GravitationalBodyForce(
        const GravitationalBodyForce& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    GravitationalBodyForce&
    operator=(
        const GravitationalBodyForce& that);

    /*!
     * Read input values, indicated above, from given database.
     */
    void
    getFromInput(
        tbox::Pointer<tbox::Database> db);

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    string d_object_name;

    /*
     * The gravitational force vector.
     */
    vector<double> d_gravitational_force;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "GravitationalBodyForce.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_GravitationalBodyForce
