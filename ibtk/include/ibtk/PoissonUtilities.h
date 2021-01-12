// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_PoissonUtilities
#define included_IBTK_PoissonUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_enums.h"

#include "BoundaryBox.h"
#include "PoissonSpecifications.h"
#include "tbox/Pointer.h"

#include <map>
#include <vector>

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
     * Compute the matrix coefficients corresponding to a side-centered
     * discretization of the divergence of the viscous stress tensor.
     *
     * \note The scaling factors of \f$ C \f$ and \f$ D \f$ variables in
     * the PoissonSpecification object are passed separately and are denoted
     * by \f$ \beta \f$ and \f$ \alpha \f$, respectively.
     */
    static void computeVCSCViscousOpMatrixCoefficients(
        SAMRAI::pdat::SideData<NDIM, double>& matrix_coefficients,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        const std::vector<std::map<SAMRAI::hier::Index<NDIM>, int, IndexFortranOrder> >& stencil_map_vec,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        double alpha,
        double beta,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
        double data_time,
        VCInterpType mu_interp_type = VC_HARMONIC_INTERP);

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
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a side-centered discretization of the
     * variable-coefficient viscous operator.
     *
     * \note The scaling factors of \f$ D \f$ variable in the PoissonSpecification object
     * is passed separately and is denoted \f$ \alpha \f$.
     */
    static void
    adjustVCSCViscousOpRHSAtPhysicalBoundary(SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
                                             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                             const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                             double alpha,
                                             const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                             double data_time,
                                             bool homogeneous_bc,
                                             VCInterpType mu_interp_type = VC_HARMONIC_INTERP);

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

    /*!
     * Modify the right-hand side entries to account for coarse-fine interface boundary conditions corresponding to a
     * side-centered discretization of the variable coefficient viscous operator.
     *
     * \note This function simply uses ghost cell values in sol_data to provide Dirichlet boundary values at coarse-fine
     * interfaces.  A more complete implementation would employ the interpolation stencil used at coarse-fine interfaces
     * to modify both the matrix coefficients and RHS values at coarse-fine interfaces.
     *
     * \note The scaling factors of \f$ D \f$ variable in the PoissonSpecification object
     * is passed separately and is denoted \f$ \alpha \f$.
     */
    static void adjustVCSCViscousOpRHSAtCoarseFineBoundary(
        SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
        const SAMRAI::pdat::SideData<NDIM, double>& sol_data,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        double alpha,
        const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >& type1_cf_bdry,
        VCInterpType mu_interp_type = VC_HARMONIC_INTERP);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PoissonUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonUtilities(const PoissonUtilities& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonUtilities& operator=(const PoissonUtilities& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PoissonUtilities
