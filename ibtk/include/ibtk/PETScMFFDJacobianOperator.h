// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
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

#ifndef included_IBTK_PETScMFFDJacobianOperator
#define included_IBTK_PETScMFFDJacobianOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/GeneralOperator.h>
#include <ibtk/JacobianOperator.h>
#include <ibtk/PETScNewtonKrylovSolver.h>

#include <tbox/Pointer.h>

#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <IntVector.h>
#include <SAMRAIVectorReal.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScMFFDJacobianOperator provides a method for computing
 * Jacobian-vector products, i.e., \f$ F'[x]v \f$, via a matrix-free
 * finite-difference approach.
 */
class PETScMFFDJacobianOperator : public JacobianOperator
{
public:
    /*!
     * \brief Constructor.
     */
    PETScMFFDJacobianOperator(std::string object_name, std::string options_prefix = "");

    /*!
     * \brief Deallocate operator state and destroy the operator.
     */
    ~PETScMFFDJacobianOperator();

    /*!
     * \brief Set the operator to use in computing approximations to
     * Jacobian-vector products.
     *
     * The caller configures and initializes F before formJacobian(), which
     * evaluates it immediately in standalone use. This operator neither
     * initializes F nor forwards its own time settings to F.
     */
    void setOperator(SAMRAI::tbox::Pointer<GeneralOperator> F);

    /*!
     * \brief Set the PETScNewtonKrylov solver using this object to compute
     * Jacobian-vector products.
     */
    void setNewtonKrylovSolver(SAMRAI::tbox::Pointer<PETScNewtonKrylovSolver> nonlinear_solver);

    /*!
     * \name General Jacobian functionality.
     */
    //\{

    /*!
     * \brief Set the base state required for evaluating \f$F'[u]\f$.
     *
     * In standalone use, requires initialized operator state, copies u, and
     * caches F(u). F must already be ready to evaluate (see setOperator()).
     * Form this pair again after changing F or its time-dependent/shared state.
     * When associated with PETScNewtonKrylovSolver, uses that solver's current
     * solution and function vectors instead of an independent base from u.
     *
     * \param u standalone value where the Jacobian is to be evaluated
     */
    void formJacobian(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u) override;

    /*!
     * \brief Return the vector where the Jacobian is evaluated.
     *
     * In standalone use, the operator-owned vector is allocated at initialization,
     * but its contents are meaningful only after formJacobian(). Later formation
     * overwrites it; deallocation frees its components and resets the returned
     * Pointer to nullptr. Retaining a Pointer does not preserve those components.
     * Mutating the vector does not refresh the cached function value.
     * With an associated PETScNewtonKrylovSolver, returns its current solution;
     * the solver must have initialized its SNES state.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> getBaseVector() const override;

    //\}

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
     * The user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must have enough
     *   ghost cells for computation of Ax).
     *
     * In general, the vectors x and y \em cannot be the same.
     *
     * \note The operator MUST be initialized prior to calling apply.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Compute hierarchy dependent data required for computing y=Ax and
     * z=Ax+y.
     *
     * The vector arguments for apply(), applyAdd(), etc, need not match those
     * for initializeOperatorState().  However, there must be a certain degree
     * of similarity, including
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the input and output vectors
     *
     * \note It is generally necessary to reinitialize the operator state when
     * the hierarchy configuration changes.
     *
     * It is safe to call initializeOperatorState() when the state is already
     * initialized.  In this case, the operator state is first deallocated and
     * then reinitialized.
     *
     * Conditions on arguments:
     * - input and output vectors must have same hierarchy
     * - input and output vectors must have same structure, depth, etc.
     *
     * Call deallocateOperatorState() to remove any data allocated by this
     * method.
     *
     * \see deallocateOperatorState
     *
     * \param in input vector
     * \param out output vector
     *
     * Allocates this operator's matrix and vector storage, but does not
     * initialize F or form a base/function pair. See setOperator() and
     * formJacobian().
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * Remove all hierarchy dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() when the operator state is
     * already deallocated.
     *
     * \see initializeOperatorState
     *
     * Releases this operator's matrix and vector storage, not F's state.
     */
    void deallocateOperatorState() override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScMFFDJacobianOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScMFFDJacobianOperator(const PETScMFFDJacobianOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScMFFDJacobianOperator& operator=(const PETScMFFDJacobianOperator& that) = delete;

    static PetscErrorCode FormFunction_SAMRAI(void* p_ctx, Vec x, Vec f);

    SAMRAI::tbox::Pointer<GeneralOperator> d_F;
    SAMRAI::tbox::Pointer<PETScNewtonKrylovSolver> d_nonlinear_solver;
    Mat d_petsc_jac = nullptr;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_op_u, d_op_f_base, d_op_x, d_op_y;
    Vec d_petsc_u = nullptr, d_petsc_f_base = nullptr, d_petsc_x = nullptr, d_petsc_y = nullptr;
    std::string d_options_prefix;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PETScMFFDJacobianOperator
