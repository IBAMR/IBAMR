// Filename: LinearSolver.h
// Created on 08 Sep 2003 by Boyce Griffith
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

#ifndef included_LinearSolver
#define included_LinearSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/GeneralSolver.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LinearSolver provides an abstract interface for the
 * implementation of solvers for linear problems of the form \f$Ax=b\f$.
 */
class LinearSolver
    : public GeneralSolver
{
public:
    /*!
     * \brief Constructor.
     */
    LinearSolver(
        const std::string& object_name,
        bool homogeneous_bc=false);

    /*!
     * \brief Empty virtual destructor.
     */
    virtual
    ~LinearSolver();

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Set the nullspace of the linear system.
     *
     * Implementations can require the nullspace basis vectors to be orthogonal
     * but should not assume the basis vectors to be orthonormal.  If the basis
     * vectors are not orthonormal, the solver may normalize them in place.
     */
    virtual void
    setNullspace(
        bool contains_constant_vector,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >& nullspace_basis_vecs=std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >());

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    virtual void
    setInitialGuessNonzero(
        bool initial_guess_nonzero=true);

    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool
    getInitialGuessNonzero() const;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    virtual void
    setMaxIterations(
        int max_iterations);

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     */
    virtual int
    getMaxIterations() const;

    /*!
     * \brief Set the absolute residual tolerance for convergence.
     */
    virtual void
    setAbsoluteTolerance(
        double abs_residual_tol);

    /*!
     * \brief Get the absolute residual tolerance for convergence.
     */
    virtual double
    getAbsoluteTolerance() const;

    /*!
     * \brief Set the relative residual tolerance for convergence.
     */
    virtual void
    setRelativeTolerance(
        double rel_residual_tol);

    /*!
     * \brief Get the relative residual tolerance for convergence.
     */
    virtual double
    getRelativeTolerance() const;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    virtual int
    getNumIterations() const;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double
    getResidualNorm() const;

    //\}

protected:
    // Nullspace data.
    bool d_nullspace_contains_constant_vector;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > > d_nullspace_basis_vecs;

    // Solver parameters.
    bool d_initial_guess_nonzero;
    double d_rel_residual_tol;
    double d_abs_residual_tol;
    int d_max_iterations;
    int d_current_its;
    double d_current_residual_norm;

private:
    /*!
     * \brief Default constructor.
     *
     * \param from The value to copy to this object.
     */
    LinearSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LinearSolver(
        const LinearSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LinearSolver&
    operator=(
        const LinearSolver& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/LinearSolver.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LinearSolver
