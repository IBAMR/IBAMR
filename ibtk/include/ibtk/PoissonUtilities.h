// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIArray.h"
#include "SAMRAIBoundaryBox.h"
#include "SAMRAICellData.h"
#include "SAMRAIIndex.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPointer.h"
#include "SAMRAIPoissonSpecifications.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAISideData.h"

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
    static void computeMatrixCoefficients(SAMRAICellData<double>& matrix_coefficients,
                                          SAMRAIPointer<SAMRAIPatch> patch,
                                          const std::vector<SAMRAIIndex>& stencil,
                                          const SAMRAIPoissonSpecifications& poisson_spec,
                                          SAMRAIRobinBcCoefStrategy* bc_coef,
                                          double data_time);

    /*!
     * Compute the matrix coefficients corresponding to a cell-centered
     * discretization of the Laplacian.
     */
    static void computeMatrixCoefficients(SAMRAICellData<double>& matrix_coefficients,
                                          SAMRAIPointer<SAMRAIPatch> patch,
                                          const std::vector<SAMRAIIndex>& stencil,
                                          const SAMRAIPoissonSpecifications& poisson_spec,
                                          const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
                                          double data_time);

    /*!
     * Compute the matrix coefficients corresponding to a side-centered
     * discretization of the Laplacian.
     */
    static void computeMatrixCoefficients(SAMRAISideData<double>& matrix_coefficients,
                                          SAMRAIPointer<SAMRAIPatch> patch,
                                          const std::vector<SAMRAIIndex>& stencil,
                                          const SAMRAIPoissonSpecifications& poisson_spec,
                                          const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
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
        SAMRAISideData<double>& matrix_coefficients,
        SAMRAIPointer<SAMRAIPatch> patch,
        const std::vector<std::map<SAMRAIIndex, int, IndexFortranOrder>>& stencil_map_vec,
        const SAMRAIPoissonSpecifications& poisson_spec,
        double alpha,
        double beta,
        const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
        double data_time,
        VCInterpType mu_interp_type = VC_HARMONIC_INTERP);

    /*!
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a cell-centered discretization of the
     * Laplacian.
     */
    static void adjustRHSAtPhysicalBoundary(SAMRAICellData<double>& rhs_data,
                                            SAMRAIPointer<SAMRAIPatch> patch,
                                            const SAMRAIPoissonSpecifications& poisson_spec,
                                            SAMRAIRobinBcCoefStrategy* bc_coef,
                                            double data_time,
                                            bool homogeneous_bc);

    /*!
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a cell-centered discretization of the
     * Laplacian.
     */
    static void adjustRHSAtPhysicalBoundary(SAMRAICellData<double>& rhs_data,
                                            SAMRAIPointer<SAMRAIPatch> patch,
                                            const SAMRAIPoissonSpecifications& poisson_spec,
                                            const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
                                            double data_time,
                                            bool homogeneous_bc);

    /*!
     * Modify the right-hand side entries to account for physical boundary
     * conditions corresponding to a side-centered discretization of the
     * Laplacian.
     */
    static void adjustRHSAtPhysicalBoundary(SAMRAISideData<double>& rhs_data,
                                            SAMRAIPointer<SAMRAIPatch> patch,
                                            const SAMRAIPoissonSpecifications& poisson_spec,
                                            const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
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
    static void adjustVCSCViscousOpRHSAtPhysicalBoundary(SAMRAISideData<double>& rhs_data,
                                                         SAMRAIPointer<SAMRAIPatch> patch,
                                                         const SAMRAIPoissonSpecifications& poisson_spec,
                                                         double alpha,
                                                         const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
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
    static void adjustRHSAtCoarseFineBoundary(SAMRAICellData<double>& rhs_data,
                                              const SAMRAICellData<double>& sol_data,
                                              SAMRAIPointer<SAMRAIPatch> patch,
                                              const SAMRAIPoissonSpecifications& poisson_spec,
                                              const SAMRAIArray<SAMRAIBoundaryBox>& type1_cf_bdry);

    /*!
     * Modify the right-hand side entries to account for coarse-fine interface boundary conditions corresponding to a
     * side-centered discretization of the Laplacian.
     *
     * \note This function simply uses ghost cell values in sol_data to provide Dirichlet boundary values at coarse-fine
     * interfaces.  A more complete implementation would employ the interpolation stencil used at coarse-fine interfaces
     * to modify both the matrix coefficients and RHS values at coarse-fine interfaces.
     */
    static void adjustRHSAtCoarseFineBoundary(SAMRAISideData<double>& rhs_data,
                                              const SAMRAISideData<double>& sol_data,
                                              SAMRAIPointer<SAMRAIPatch> patch,
                                              const SAMRAIPoissonSpecifications& poisson_spec,
                                              const SAMRAIArray<SAMRAIBoundaryBox>& type1_cf_bdry);

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
    static void adjustVCSCViscousOpRHSAtCoarseFineBoundary(SAMRAISideData<double>& rhs_data,
                                                           const SAMRAISideData<double>& sol_data,
                                                           SAMRAIPointer<SAMRAIPatch> patch,
                                                           const SAMRAIPoissonSpecifications& poisson_spec,
                                                           double alpha,
                                                           const SAMRAIArray<SAMRAIBoundaryBox>& type1_cf_bdry,
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
