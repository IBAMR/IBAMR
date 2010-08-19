// Filename: GravitationalBodyForce.h
// Created on 03 May 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef included_GravitationalBodyForce
#define included_GravitationalBodyForce

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>

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
 * \brief Class to specify a uniform body force that arises from gravitational
 * forces acting on an incompressible fluid with uniform mass density.
 */
class GravitationalBodyForce
    : public CartGridFunction
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
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    virtual bool
    isTimeDependent() const
        { return true; }

    /*!
     * Evaluate the function on the patch interior.
     */
    virtual void
    setDataOnPatch(
        const int data_idx,
        tbox::Pointer<hier::Variable<NDIM> > var,
        tbox::Pointer<hier::Patch<NDIM> > patch,
        const double data_time,
        const bool initial_time=false,
        tbox::Pointer<hier::PatchLevel<NDIM> > level=tbox::Pointer<hier::PatchLevel<NDIM> >(NULL));

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
