// Filename: PoissonUtilities.h
// Created on 24 Aug 2010 by Boyce Griffith
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

#ifndef included_PoissonUtilities
#define included_PoissonUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "PoissonSpecifications.h"
#include "BoundaryBox.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Index;
template <int DIM>
class Patch;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class CellData;
template <int DIM, class TYPE>
class SideData;
} // namespace pdat
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PoissonUtilities provides utility functions for constructing
 * Poisson solvers.
 */
class PoissonUtilities
{
public:
    /*!
     * Compute the matrix coefficients corresponding to a cell-centered
     * discretization of the Laplacian.
     */
    static void computeMatrixCoefficients(SAMRAI::pdat::CellData<NDIM, double>& matrix_coefficients,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                          const std::vector<SAMRAI::hier::Index<NDIM> >& stencil,
                                          const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                          SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
                                          double data_time);

    /*!
     * Compute the matrix coefficients corresponding to a cell-centered
     * discretization of the Laplacian.
     */
    static void computeMatrixCoefficients(SAMRAI::pdat::CellData<NDIM, double>& matrix_coefficients,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                          const std::vector<SAMRAI::hier::Index<NDIM> >& stencil,
                                          const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                          const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                          double data_time);

    /*!
     * Compute the matrix coefficients corresponding to a side-centered
     * discretization of the Laplacian.
     */
    static void computeMatrixCoefficients(SAMRAI::pdat::SideData<NDIM, double>& matrix_coefficients,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                          const std::vector<SAMRAI::hier::Index<NDIM> >& stencil,
                                          const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                          const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                          double data_time);

    /*!
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a cell-centered discretization of the
     * Laplacian.
     */
    static void adjustRHSAtPhysicalBoundary(SAMRAI::pdat::CellData<NDIM, double>& rhs_data,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                            const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
                                            double data_time,
                                            bool homogeneous_bc);

    /*!
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a cell-centered discretization of the
     * Laplacian.
     */
    static void adjustRHSAtPhysicalBoundary(SAMRAI::pdat::CellData<NDIM, double>& rhs_data,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                            const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                            double data_time,
                                            bool homogeneous_bc);

    /*!
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a side-centered discretization of the
     * Laplacian.
     */
    static void adjustRHSAtPhysicalBoundary(SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                            const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                            double data_time,
                                            bool homogeneous_bc);

    /*!
     * Modify the right-hand side entries to account for coarse-fine interface boundary conditions corresponding to a
     * cell-centered discretization of the Laplacian.
     *
     * \note This function simply uses ghost cell values in sol_data to provide Dirichlet boundary values at coarse-fine
     * interfaces.  A more complete implementation would employ the interpolation stencil used at coarse-fine interfaces
     * to modify both the matrix coefficients and RHS values at coarse-fine interfaces.
     */
    static void
    adjustRHSAtCoarseFineBoundary(SAMRAI::pdat::CellData<NDIM, double>& rhs_data,
                                  const SAMRAI::pdat::CellData<NDIM, double>& sol_data,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                  const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                  const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >& type1_cf_bdry);

    /*!
     * Modify the right-hand side entries to account for coarse-fine interface boundary conditions corresponding to a
     * side-centered discretization of the Laplacian.
     *
     * \note This function simply uses ghost cell values in sol_data to provide Dirichlet boundary values at coarse-fine
     * interfaces.  A more complete implementation would employ the interpolation stencil used at coarse-fine interfaces
     * to modify both the matrix coefficients and RHS values at coarse-fine interfaces.
     */
    static void
    adjustRHSAtCoarseFineBoundary(SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
                                  const SAMRAI::pdat::SideData<NDIM, double>& sol_data,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                  const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                  const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >& type1_cf_bdry);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PoissonUtilities();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonUtilities(const PoissonUtilities& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonUtilities& operator=(const PoissonUtilities& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PoissonUtilities
