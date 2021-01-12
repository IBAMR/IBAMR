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

#ifndef included_IBAMR_StaggeredStokesBlockFactorizationPreconditioner
#define included_IBAMR_StaggeredStokesBlockFactorizationPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesSolver.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"

#include "CellVariable.h"
#include "SideVariable.h"
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
 * \brief Class StaggeredStokesBlockFactorizationPreconditioner is a concrete
 * StaggeredStokesBlockPreconditioner that implements a staggered grid (MAC)
 * block factorization (approximate Schur complement) preconditioner for the
 * incompressible Stokes equations.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class StaggeredStokesBlockFactorizationPreconditioner : public StaggeredStokesBlockPreconditioner
{
public:
    enum FactorizationType
    {
        LOWER_TRIANGULAR,
        UPPER_TRIANGULAR,
        SYMMETRIC,
        DIAGONAL
    };

    /*!
     * \brief Class constructor
     */
    StaggeredStokesBlockFactorizationPreconditioner(const std::string& object_name,
                                                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                    const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesBlockFactorizationPreconditioner();

    /*!
     * \brief Static function to construct a
     * StaggeredStokesBlockFactorizationPreconditioner.
     */
    static SAMRAI::tbox::Pointer<StaggeredStokesSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        return new StaggeredStokesBlockFactorizationPreconditioner(object_name, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * Set the factorization type.
     */
    void setFactorizationType(FactorizationType factorization_type);

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

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesBlockFactorizationPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesBlockFactorizationPreconditioner(const StaggeredStokesBlockFactorizationPreconditioner& from) =
        delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesBlockFactorizationPreconditioner&
    operator=(const StaggeredStokesBlockFactorizationPreconditioner& that) = delete;

    /*!
     * \brief Solve the pressure subsystem.
     */
    void solvePressureSubsystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b,
                                bool initial_guess_nonzero);

    /*!
     * \brief Solve the velocity subsystem.
     */
    void solveVelocitySubsystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b,
                                bool initial_guess_nonzero);

    // Solver configuration
    FactorizationType d_factorization_type = LOWER_TRIANGULAR;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_P_bdry_fill_op, d_no_fill_op;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_var;
    int d_F_U_mod_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_P_var;
    int d_P_scratch_idx = IBTK::invalid_index, d_F_P_mod_idx = IBTK::invalid_index;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_StaggeredStokesBlockFactorizationPreconditioner
