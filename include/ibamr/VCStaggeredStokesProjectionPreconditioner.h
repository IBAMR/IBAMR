// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_VCStaggeredStokesProjectionPreconditioner
#define included_IBAMR_VCStaggeredStokesProjectionPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesSolver.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"

#include "CellVariable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>

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
 * \brief Class VCStaggeredStokesProjectionPreconditioner is a concrete StokesSolver
 * that implements a staggered grid (MAC) projection solver for the
 * incompressible Stokes operator with variable coefficients.
 *
 * References
 * Cai et al, <A HREF="https://arxiv.org/pdf/1308.4605.pdf">Efficient Variable-Coefficient
 * Finite-Volume Stokes Solvers</A>
 *
 *
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class VCStaggeredStokesProjectionPreconditioner : public StaggeredStokesBlockPreconditioner
{
public:
    /*!
     * \brief Class constructor
     */
    VCStaggeredStokesProjectionPreconditioner(const std::string& object_name,
                                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                              const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~VCStaggeredStokesProjectionPreconditioner();

    /*!
     * \brief Static function to construct a
     * VCStaggeredStokesProjectionPreconditioner.
     */
    static SAMRAI::tbox::Pointer<StaggeredStokesSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        return new VCStaggeredStokesProjectionPreconditioner(object_name, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Compute the action of the preconditioner.
     */
    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                     SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

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
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

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
    void deallocateSolverState() override;

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    void setInitialGuessNonzero(bool initial_guess_nonzero = true) override;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    void setMaxIterations(int max_iterations) override;

    /*!
     * \brief Set the cell centered quantity representing the viscous velocity coefficient
     *
     * \note This quantity is needed on cell centers for the pressure update
     */
    void setVelocityCellCenteredDCoefficient(int velocity_D_cc_idx);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VCStaggeredStokesProjectionPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCStaggeredStokesProjectionPreconditioner(const VCStaggeredStokesProjectionPreconditioner& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCStaggeredStokesProjectionPreconditioner&
    operator=(const VCStaggeredStokesProjectionPreconditioner& that) = delete;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_Phi_bdry_fill_op, d_no_fill_op;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Phi_var, d_F_Phi_var;
    int d_Phi_scratch_idx, d_F_Phi_idx, d_velocity_D_cc_idx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_VCStaggeredStokesProjectionPreconditioner
