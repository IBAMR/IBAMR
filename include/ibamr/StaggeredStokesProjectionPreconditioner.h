// Filename: StaggeredStokesProjectionPreconditioner.h
// Created on 29 Mar 2008 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBAMR_StaggeredStokesProjectionPreconditioner
#define included_IBAMR_StaggeredStokesProjectionPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "CellVariable.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesProjectionPreconditioner is a concrete StokesSolver
 * that implements a staggered grid (MAC) projection solver for the
 * incompressible Stokes operator.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class StaggeredStokesProjectionPreconditioner : public StaggeredStokesBlockPreconditioner
{
public:
    /*!
     * \brief Class constructor
     */
    StaggeredStokesProjectionPreconditioner(const std::string& object_name,
                                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                            const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesProjectionPreconditioner();

    /*!
     * \brief Static function to construct a
     * StaggeredStokesProjectionPreconditioner.
     */
    static SAMRAI::tbox::Pointer<StaggeredStokesSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        return new StaggeredStokesProjectionPreconditioner(object_name, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Compute the action of the preconditioner.
     */
    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a b must have same patch hierarchy
     * - vectors \a x and \a b must have same structure, depth, etc.
     *
     * \note It is safe to call initializeSolverState() when the solver state is
     * already initialized.
     *
     * \see deallocateSolverState
     *
     * \note A default implementation is provided which does nothing.
     */
    void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     *
     * \note A default implementation is provided which does nothing.
     */
    void deallocateSolverState();

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    void setInitialGuessNonzero(bool initial_guess_nonzero = true);

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    void setMaxIterations(int max_iterations);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesProjectionPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesProjectionPreconditioner(const StaggeredStokesProjectionPreconditioner& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesProjectionPreconditioner& operator=(const StaggeredStokesProjectionPreconditioner& that) = delete;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_Phi_bdry_fill_op, d_no_fill_op;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Phi_var, d_F_Phi_var;
    int d_Phi_scratch_idx, d_F_Phi_idx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_StaggeredStokesProjectionPreconditioner
