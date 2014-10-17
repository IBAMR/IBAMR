// Filename: INSIntermediateVelocityBcCoef.h
// Created on 24 Jul 2008 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_INSIntermediateVelocityBcCoef
#define included_INSIntermediateVelocityBcCoef

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "IntVector.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class Patch;
template <int DIM>
class Variable;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class ArrayData;
} // namespace pdat
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSIntermediateVelocityBcCoef is a concrete
 * SAMRAI::solv::RobinBcCoefStrategy that is used to specify boundary conditions
 * for the intermediate velocity solve embedded in a projection method or in a
 * block preconditioner for the incompressible Navier-Stokes equations.
 *
 * This class interprets pure Dirichlet boundary conditions as prescribed
 * velocity boundary conditions, whereas pure Neumann boundary conditions are
 * interpreted as prescribed traction (stress) boundary conditions.  These are
 * translated into Dirichlet and Neumann boundary conditions, respectively, for
 * the intermediate velocity.
 */
class INSIntermediateVelocityBcCoef : public IBTK::ExtendedRobinBcCoefStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param comp_idx        Component of the velocity which this boundary condition
     *specification
     *is to operate on
     * \param bc_coefs        IBTK::Vector of boundary condition specification objects
     * \param homogeneous_bc  Whether to employ homogeneous (as opposed to inhomogeneous)
     *boundary
     *conditions
     *
     * \note Precisely NDIM boundary condition objects must be provided to the
     * class constructor.
     */
    INSIntermediateVelocityBcCoef(int comp_idx,
                                  const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                  bool homogeneous_bc = false);

    /*!
     * \brief Destructor.
     */
    ~INSIntermediateVelocityBcCoef();

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \param bc_coefs  IBTK::Vector of boundary condition specification objects
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void setSolutionTime(double solution_time);

    /*!
     * \brief Set the current time interval.
     */
    void setTimeInterval(double current_time, double new_time);

    /*!
     * \name Extended SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Set the target data index.
     */
    void setTargetPatchDataIndex(int target_idx);

    /*!
     * \brief Clear the target data index.
     */
    void clearTargetPatchDataIndex();

    /*!
     * \brief Set whether the class is filling homogeneous or inhomogeneous
     * boundary conditions.
     */
    void setHomogeneousBc(bool homogeneous_bc);

    //\}

    /*!
     * \name Implementation of SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     *
     * \note In the original SAMRAI::solv::RobinBcCoefStrategy interface, it was
     * assumed that \f$ b = (1-a) \f$.  In the new interface, \f$a\f$ and
     * \f$b\f$ are independent.
     *
     * \see SAMRAI::solv::RobinBcCoefStrategy::setBcCoefs()
     *
     * \param acoef_data  Boundary coefficient data.
     *        The array will have been defined to include index range for
     *        corresponding to the boundary box \a bdry_box and appropriate for
     *        the alignment of the given variable.  If this is a null pointer,
     *        then the calling function is not interested in a, and you can
     *        disregard it.
     * \param bcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it is to be
     *        filled with the b coefficient.
     * \param gcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it is to be
     *        filled with the g coefficient.
     * \param variable    Variable to set the coefficients for.
     *        If implemented for multiple variables, this parameter can be used
     *        to determine which variable's coefficients are being sought.
     * \param patch       Patch requiring bc coefficients.
     * \param bdry_box    Boundary box showing where on the boundary the coefficient data is
     *needed.
     * \param fill_time   Solution time corresponding to filling, for use when coefficients are
     *time-dependent.
     */
    void setBcCoefs(SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& acoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& bcoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& gcoef_data,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
                    const SAMRAI::hier::Patch<NDIM>& patch,
                    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const;

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
    SAMRAI::hier::IntVector<NDIM> numberOfExtensionsFillable() const;

    //\}

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSIntermediateVelocityBcCoef();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSIntermediateVelocityBcCoef(const INSIntermediateVelocityBcCoef& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSIntermediateVelocityBcCoef& operator=(const INSIntermediateVelocityBcCoef& that);

    /*
     * Component of the velocity which this boundary condition specification is
     * to operate on.
     */
    const int d_comp_idx;

    /*
     * The boundary condition specification objects for the updated velocity.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSIntermediateVelocityBcCoef
