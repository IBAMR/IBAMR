// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/IndexUtilities.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Pointer.h>

#include <BoundaryBox.h>
#include <PoissonSpecifications.h>

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
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                          const std::vector<SAMRAI::hier::Index<NDIM>>& stencil,
                                          const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                          SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
                                          double data_time);

    /*!
     * Compute the matrix coefficients corresponding to a cell-centered
     * discretization of the Laplacian.
     */
    static void computeMatrixCoefficients(SAMRAI::pdat::CellData<NDIM, double>& matrix_coefficients,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                          const std::vector<SAMRAI::hier::Index<NDIM>>& stencil,
                                          const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                          const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                          double data_time);

    /*!
     * Compute the matrix coefficients corresponding to a side-centered
     * discretization of the Laplacian.
     */
    static void computeMatrixCoefficients(SAMRAI::pdat::SideData<NDIM, double>& matrix_coefficients,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                          const std::vector<SAMRAI::hier::Index<NDIM>>& stencil,
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
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
        const std::vector<std::map<SAMRAI::hier::Index<NDIM>, int, IndexFortranOrder>>& stencil_map_vec,
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
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
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
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                            const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                            double data_time,
                                            bool homogeneous_bc);

    /*!
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a side-centered discretization of the
     * Laplacian.
     *
     * The entry of a side that lies on a physical boundary, normal to its component, at which the Robin coefficient b
     * is exactly zero (a Dirichlet condition) is not modified: the matrix row of that side is an identity row, and
     * the entry holds the boundary value.
     */
    static void adjustRHSAtPhysicalBoundary(SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                            const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                            double data_time,
                                            bool homogeneous_bc);

    /*!
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a side-centered discretization of the
     * variable-coefficient viscous operator.
     *
     * The entry of a side that lies on a physical boundary, normal to its component, at which the Robin coefficient b
     * is exactly zero (a Dirichlet condition) is not modified: the matrix row of that side is an identity row, and
     * the entry holds the boundary value.
     *
     * \note The scaling factors of \f$ D \f$ variable in the PoissonSpecification object
     * is passed separately and is denoted \f$ \alpha \f$.
     */
    static void
    adjustVCSCViscousOpRHSAtPhysicalBoundary(SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
                                             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
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
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                  const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                  const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& type1_cf_bdry);

    /*!
     * Modify the right-hand side entries to account for coarse-fine interface boundary conditions corresponding to a
     * side-centered discretization of the Laplacian.
     *
     * The ghost values in sol_data take the place of the degrees of freedom that the matrix of a level solver does not
     * contain: for each side in the side box of the patch and each of its stencil neighbors, the contribution of the
     * neighbor is subtracted from the right-hand side if no patch of the level has the neighbor as a side. The
     * contributions of neighbors that lie on the other side of a physical boundary are not included. This function only
     * modifies rhs_data in the side box of the patch and only reads ghost values from sol_data, so rhs_data does not
     * need ghost cells.
     *
     * The boundary boxes of codimension 1 do not tell this function which ghost cells that touch the patch only along
     * an edge or a corner belong to the level. It takes such a ghost cell to belong to the level if either of the two
     * ghost cells that it neighbors toward the patch does. This is exact unless two patches of the level touch only at
     * an isolated corner or edge.
     *
     * \note This function cannot tell which rows of the matrix are identity rows, which are the rows of the sides that
     * lie on a physical boundary at which a Dirichlet condition is imposed on the component normal to the boundary. It
     * modifies the right-hand side entries of those sides as well. On a patch level with physical boundary conditions,
     * use the version of this function that takes the boundary boxes of codimension 2 and the boundary condition
     * coefficients.
     *
     * \note This function simply uses ghost cell values in sol_data to provide Dirichlet boundary values at coarse-fine
     * interfaces.  A more complete implementation would employ the interpolation stencil used at coarse-fine interfaces
     * to modify both the matrix coefficients and RHS values at coarse-fine interfaces.
     */
    static void
    adjustRHSAtCoarseFineBoundary(SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
                                  const SAMRAI::pdat::SideData<NDIM, double>& sol_data,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                  const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                  const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& type1_cf_bdry);

    /*!
     * Modify the right-hand side entries to account for coarse-fine interface boundary conditions corresponding to a
     * side-centered discretization of the Laplacian, as in the version of this function that takes only the boundary
     * boxes of codimension 1, but exact for any layout of the patches of the level.
     *
     * The boundary boxes of codimension 2 of the patch (the corner cells in 2D and the edge cells in 3D) determine
     * which ghost cells that touch the patch only along an edge or a corner belong to the level. The boundary
     * condition coefficients, one object for each component, determine which sides on a physical boundary have an
     * identity row in the matrix: this function does not modify the right-hand side entries of the sides that lie on
     * a physical boundary normal to their component and at which the coefficient b is exactly zero, because the value
     * of such an entry is the Dirichlet value. The coefficients are evaluated in the same way as in
     * adjustRHSAtPhysicalBoundary(), and only for the physical boundaries that contain a side whose entry would
     * otherwise be modified.
     */
    static void adjustRHSAtCoarseFineBoundary(SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
                                              const SAMRAI::pdat::SideData<NDIM, double>& sol_data,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                              const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                              const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& type1_cf_bdry,
                                              const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& type2_cf_bdry,
                                              const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                              double data_time,
                                              bool homogeneous_bc);

    /*!
     * Modify the right-hand side entries to account for coarse-fine interface boundary conditions corresponding to a
     * side-centered discretization of the variable coefficient viscous operator.
     *
     * This function corrects every stencil entry of a side that reaches across a boundary box of codimension 1 of the
     * patch. It is exact only if no two patches of the level touch only at a corner or an edge, the refined region has
     * no concave corner, and no coarse-fine boundary meets a physical boundary at which a Dirichlet condition is
     * imposed: otherwise it corrects couplings that the matrix of the level keeps, and it modifies the right-hand side
     * entries of identity rows. Use the version of this function that takes the boundary boxes of codimension 2 and
     * the boundary condition coefficients on a patch level with such features.
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
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        double alpha,
        const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& type1_cf_bdry,
        VCInterpType mu_interp_type = VC_HARMONIC_INTERP);

    /*!
     * Modify the right-hand side entries to account for coarse-fine interface boundary conditions corresponding to a
     * side-centered discretization of the variable coefficient viscous operator, as in the version of this function
     * that takes only the boundary boxes of codimension 1, but exact for any layout of the patches of the level.
     *
     * The ghost values in sol_data take the place of the degrees of freedom that the matrix of a level solver does not
     * contain: for each side in the side box of the patch and each entry of its stencil, including the entries that
     * couple the components of the velocity, the term of the entry is subtracted from the right-hand side if no patch
     * of the level has the side of the entry as a side. The coefficients of the terms are the coefficients of the
     * matrix of the level solver, so the terms of the sides that lie outside a physical boundary are not included, and
     * the right-hand side entry of an identity row is not modified. This function only modifies rhs_data in the side
     * box of the patch and only reads ghost values from sol_data, so rhs_data does not need ghost cells.
     *
     * The boundary boxes of codimension 2 of the patch (the corner cells in 2D and the edge cells in 3D) determine
     * which ghost cells that touch the patch only along an edge or a corner belong to the level. The boundary
     * condition coefficients are those of the matrix, and are evaluated at data_time.
     *
     * \note The scaling factors of \f$ D \f$ variable in the PoissonSpecification object
     * is passed separately and is denoted \f$ \alpha \f$.
     */
    static void adjustVCSCViscousOpRHSAtCoarseFineBoundary(
        SAMRAI::pdat::SideData<NDIM, double>& rhs_data,
        const SAMRAI::pdat::SideData<NDIM, double>& sol_data,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        double alpha,
        const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& type1_cf_bdry,
        const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>>& type2_cf_bdry,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
        double data_time,
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

#endif // #ifndef included_IBTK_PoissonUtilities
