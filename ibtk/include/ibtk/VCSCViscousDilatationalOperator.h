// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2023 by the IBAMR developers
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

#ifndef included_IBTK_VCSCViscousDilatationalOperator
#define included_IBTK_VCSCViscousDilatationalOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/SCLaplaceOperator.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class VCSCViscousDilatationalOperator is a subclass of SCLaplaceOperator
 * which implements a globally second-order accurate side-centered finite
 * difference discretization of a vector elliptic operator of the form
 * \f$ L = \beta C I +  \nabla \cdot \mu ( (\nabla u) + (\nabla u)^T ) + \nabla (\lambda \nabla \cdot u) = \f$.
 *
 * Here \f$ u \f$ and \f$ C \f$ are vector valued side-centered fields,
 * \f$ \mu \f$ is a node-(2D) or edge-(3D) centered scalar field, and \f$ \lambda\f$
 * is a cell-centered field.
 *
 * The scaling factors of \f$ C \f$ and \f$ \mu \f$ variables are passed separately
 * and are denoted by \f$ \beta \f$ and \f$ \alpha \f$, respectively.
 */
class VCSCViscousDilatationalOperator : public SCLaplaceOperator
{
public:
    /*!
     * \brief Constructor for class VCSCViscousDilatationalOperator initializes the operator
     * coefficients and boundary conditions to default values.
     */
    VCSCViscousDilatationalOperator(std::string object_name, bool homogeneous_bc = true);

    /*!
     * \brief Destructor.
     */
    ~VCSCViscousDilatationalOperator();

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Compute y=Ax.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * Thus, the user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must
     * have enough ghost cells for computation of Ax).
     *
     * \note In general, the vectors x and y \em cannot be the same.
     *
     * Upon return from this function, the y vector will contain the result of
     * the application of A to x.
     *
     * initializeOperatorState must be called prior to any calls to
     * applyOperator.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Compute hierarchy-dependent data required for computing y=Ax (and
     * y=A'x).
     *
     * \param in input vector
     * \param out output vector
     *
     * \see KrylovLinearSolver::initializeSolverState
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Remove all hierarchy-dependent data computed by
     * initializeOperatorState().
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() even if the state is already
     * deallocated.
     *
     * \see initializeOperatorState
     * \see KrylovLinearSolver::deallocateSolverState
     */
    void deallocateOperatorState() override;

    //\}

    /*!
     * \brief Set the interpolation type to be used in computing the
     * variable coefficient viscous Laplacian.
     */
    void setDPatchDataInterpolationType(IBTK::VCInterpType D_interp_type);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VCSCViscousDilatationalOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCSCViscousDilatationalOperator(const VCSCViscousDilatationalOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCSCViscousDilatationalOperator& operator=(const VCSCViscousDilatationalOperator& that) = delete;

    /*
     * The interpolation type to be used in computing the variable coefficient viscous Laplacian.
     */
    IBTK::VCInterpType d_D_interp_type;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_VCSCViscousDilatationalOperator
