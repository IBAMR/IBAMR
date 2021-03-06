// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_LinearSolver
#define included_IBTK_LinearSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/GeneralSolver.h"

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include <iosfwd>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LinearSolver provides an abstract interface for the
 * implementation of solvers for linear problems of the form \f$Ax=b\f$.
 */
class LinearSolver : public virtual GeneralSolver
{
public:
    /*!
     * \brief Constructor.
     */
    LinearSolver();

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~LinearSolver();

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
    virtual void setNullspace(
        bool nullspace_contains_constant_vec,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > >& nullspace_basis_vecs =
            std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > >());

    /*!
     * \brief Get whether the nullspace of the linear system contains th
     * constant vector.
     */
    virtual bool getNullspaceContainsConstantVector() const;

    /*!
     * \brief Get the basis vectors for the nullspace of the linear system.
     */
    virtual const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > >&
    getNullspaceBasisVectors() const;

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    virtual void setInitialGuessNonzero(bool initial_guess_nonzero = true);

    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool getInitialGuessNonzero() const;

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Print class data to stream.
     */
    virtual void printClassData(std::ostream& stream) override;

    //\}

protected:
    // Solver parameters.
    bool d_initial_guess_nonzero = true;

    // Nullspace data.
    bool d_nullspace_contains_constant_vec = false;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_nullspace_basis_vecs;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LinearSolver(const LinearSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LinearSolver& operator=(const LinearSolver& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LinearSolver
