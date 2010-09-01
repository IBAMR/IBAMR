// Filename: VelocityBcCoefs.h
// Created on 18 Dec 2007 by Boyce Griffith
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

#ifndef included_VelocityBcCoefs
#define included_VelocityBcCoefs

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <RobinBcCoefStrategy.h>

// NAMESPACE
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class VelocityBcCoefs is an implementation of the strategy class
 * solv::RobinBcCoefStrategy that is used to specify velocity boundary
 * conditions.
 */
class VelocityBcCoefs
    : public solv::RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor
     */
    VelocityBcCoefs(
        const string& object_name,
        const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry);

    /*!
     * \brief Destructor.
     */
    virtual
    ~VelocityBcCoefs();

    /*!
     * \name Implementation of solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     *
     * \see solv::RobinBcCoefStrategy::setBcCoefs()
     *
     * \param acoef_data  Boundary coefficient data.
     *        The array will have been defined to include index range
     *        for corresponding to the boundary box \a bdry_box and
     *        appropriate for the alignment of the given variable.  If
     *        this is a null pointer, then the calling function is not
     *        interested in a, and you can disregard it.
     * \param bcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it
     *        is to be filled with the b coefficient.
     * \param gcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it
     *        is to be filled with the g coefficient.
     * \param variable    Variable to set the coefficients for.
     *        If implemented for multiple variables, this parameter
     *        can be used to determine which variable's coefficients
     *        are being sought.
     * \param patch       Patch requiring bc coefficients.
     * \param bdry_box    Boundary box showing where on the boundary the coefficient data is needed.
     * \param fill_time   Solution time corresponding to filling, for use when coefficients are time-dependent.
     */
    virtual void
    setBcCoefs(
        tbox::Pointer<pdat::ArrayData<NDIM,double> >& acoef_data,
        tbox::Pointer<pdat::ArrayData<NDIM,double> >& bcoef_data,
        tbox::Pointer<pdat::ArrayData<NDIM,double> >& gcoef_data,
        const tbox::Pointer<hier::Variable<NDIM> >& variable,
        const hier::Patch<NDIM>& patch,
        const hier::BoundaryBox<NDIM>& bdry_box,
        double fill_time=0.0) const;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     *
     * The "extension" used here is the number of cells that a boundary box
     * extends past the patch in the direction parallel to the boundary.
     *
     * Note that the inability to fill the sufficient number of cells past the
     * edge or corner of the patch may preclude the child class from being used
     * in data refinement operations that require the extra data, such as linear
     * refinement.
     *
     * The boundary box that setBcCoefs() is required to fill should not extend
     * past the limits returned by this function.
     */
    virtual hier::IntVector<NDIM>
    numberOfExtensionsFillable() const;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VelocityBcCoefs();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VelocityBcCoefs(
        const VelocityBcCoefs& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VelocityBcCoefs&
    operator=(
        const VelocityBcCoefs& that);

    /*
     * The object name.
     */
    const string d_object_name;

    /*
     * The grid geometry object.
     */
    const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <VelocityBcCoefs.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VelocityBcCoefs
