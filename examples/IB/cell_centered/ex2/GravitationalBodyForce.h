// Filename: GravitationalBodyForce.h
// Created on 03 May 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_GravitationalBodyForce
#define included_GravitationalBodyForce

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>

// SAMRAI INCLUDES
#include <tbox/Database.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec.h>

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
    ~GravitationalBodyForce();

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool
    isTimeDependent() const
        { return true; }

    /*!
     * Evaluate the function on the patch interior.
     */
    void
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
    blitz::TinyVector<double,NDIM> d_gravitational_force;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "GravitationalBodyForce.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_GravitationalBodyForce
