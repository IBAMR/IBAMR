// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2021 by the IBAMR developers
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

#ifndef included_IBTK_SCPoissonHypreLevelSolver
#define included_IBTK_SCPoissonHypreLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CoarseFineBoundary.h"
#include "Index.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"
IBTK_ENABLE_EXTRA_WARNINGS

#include <string>
#include <vector>

namespace SAMRAI
{
namespace pdat
{
template <int DIM, class TYPE>
class SideData;
} // namespace pdat
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class SCPoissonHypreLevelSolver is a concrete LinearSolver for solving
 * elliptic equations of the form \f$ \mbox{$L u$} = \mbox{$(C I + \nabla \cdot
 * D \nabla) u$} = f \f$ on a \em single SAMRAI::hier::PatchLevel using <A
 * HREF="https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html">hypre</A>.
 *
 * This solver class uses the \em hypre library to solve linear equations of the
 * form \f$ (C I + \nabla \cdot D \nabla ) u = f \f$, where \f$C\f$ and \f$D\f$
 * are constants, and \f$u\f$ and \f$f\f$ are side-centered arrays.  The
 * discretization is second-order accurate.
 *
 * Robin boundary conditions may be specified through the interface class
 * SAMRAI::solv::RobinBcCoefStrategy, but boundary conditions must be of either
 * Dirichlet or Neumann type.  In particular, mixed boundary conditions are \em
 * not presently supported.
 *
 * The user must perform the following steps to use class
 * SCPoissonHypreLevelSolver:
 *
 * -# Create a SCPoissonHypreLevelSolver object.
 * -# Set the problem specification via setPoissonSpecifications(),
 *    setPhysicalBcCoef(), and setHomogeneousBc().
 * -# Initialize SCPoissonHypreLevelSolver object using the function
 *    initializeSolverState().
 * -# Solve the linear system using the member function solveSystem(), passing
 *    in SAMRAI::solv::SAMRAIVectorReal objects corresponding to \f$u\f$ and
 *    \f$f\f$.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 enable_logging = FALSE         // see setLoggingEnabled()
 solver_type = "Split"          // choices are: "Split", "SysPFMG", "PCG", "GMRES", "FlexGMRES"
 ,
 "LGMRES", "BiCGSTAB"
 precond_type = "none"          // choices are: "Split", "SysPFMG"
 split_solver_type = "PFMG"     // choices are: "PFMG", "SMG", "Jacobi"
 max_iterations = 25            // see setMaxIterations()
 abs_residual_tol = 1.e-50      // see setAbsoluteTolerance() (only used by hypre Krylov
 solvers)
 rel_residual_tol = 1.0e-5      // see setRelativeTolerance()
 initial_guess_nonzero = FALSE  // see setInitialGuessNonzero()
 rel_change = 0                 // see hypre User's Manual (only used by SysPFMG or PCG solver)
 num_pre_relax_steps = 1        // number of pre-sweeps (only used by SysPFMG solver)
 num_post_relax_steps = 1       // number of post-sweeps (only used by SysPFMG solver)
 relax_type = 1                 // see hypre User's Manual (only used by SysPFMG solver or
 preconditioner)
 skip_relax = 1                 // see hypre User's Manual (only used by SysPFMG solver or
 preconditioner)
 two_norm = 1                   // see hypre User's Manual (only used by PCG solver)
 \endverbatim
 *
 * \em hypre is developed in the Center for Applied Scientific Computing (CASC)
 * at Lawrence Livermore National Laboratory (LLNL).  For more information about
 * \em hypre, see <A
 *
 HREF="https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html">https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html</A>.
 */
class SCPoissonHypreLevelSolver : public LinearSolver, public PoissonSolver
{
public:
    /*!
     * \brief Constructor.
     */
    SCPoissonHypreLevelSolver(const std::string& object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                              const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~SCPoissonHypreLevelSolver();

    /*!
     * \brief Static function to construct a SCPoissonHypreLevelSolver.
     */
    static SAMRAI::tbox::Pointer<PoissonSolver> allocate_solver(const std::string& object_name,
                                                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                                const std::string& default_options_prefix)
    {
        return new SCPoissonHypreLevelSolver(object_name, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Solve the linear system of equations \f$Ax=b\f$ for \f$x\f$.
     *
     * Before calling solveSystem(), the form of the solution \a x and
     * right-hand-side \a b vectors must be set properly by the user on all
     * patch interiors on the specified range of levels in the patch hierarchy.
     * The user is responsible for all data management for the quantities
     * associated with the solution and right-hand-side vectors.  In particular,
     * patch data in these vectors must be allocated prior to calling this
     * method.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a b must have same patch hierarchy
     * - vectors \a x and \a b must have same structure, depth, etc.
     *
     * \note The vector arguments for solveSystem() need not match those for
     * initializeSolverState().  However, there must be a certain degree of
     * similarity, including:\par
     * - hierarchy configuration (hierarchy pointer and range of levels)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the solution \a x and right-hand-side \a b
     *   vectors
     *
     * \note The solver need not be initialized prior to calling solveSystem();
     * however, see initializeSolverState() and deallocateSolverState() for
     * opportunities to save overhead when performing multiple consecutive
     * solves.
     *
     * \see initializeSolverState
     * \see deallocateSolverState
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                     SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     *
     * By default, the solveSystem() method computes some required hierarchy
     * dependent data before solving and removes that data after the solve.  For
     * multiple solves that use the same hierarchy configuration, it is more
     * efficient to:
     *
     * -# initialize the hierarchy-dependent data required by the solver via
     *    initializeSolverState(),
     * -# solve the system one or more times via solveSystem(), and
     * -# remove the hierarchy-dependent data via deallocateSolverState().
     *
     * Note that it is generally necessary to reinitialize the solver state when
     * the hierarchy configuration changes.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a b must have same patch hierarchy
     * - vectors \a x and \a b must have same structure, depth, etc.
     *
     * \note The vector arguments for solveSystem() need not match those for
     * initializeSolverState().  However, there must be a certain degree of
     * similarity, including:\par
     * - hierarchy configuration (hierarchy pointer and range of levels)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the solution \a x and right-hand-side \a b
     *   vectors
     *
     * \note It is safe to call initializeSolverState() when the state is
     * already initialized.  In this case, the solver state is first deallocated
     * and then reinitialized.
     *
     * \see deallocateSolverState
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
     */
    void deallocateSolverState() override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SCPoissonHypreLevelSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SCPoissonHypreLevelSolver(const SCPoissonHypreLevelSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SCPoissonHypreLevelSolver& operator=(const SCPoissonHypreLevelSolver& that) = delete;

    /*!
     * \brief Functions to allocate, initialize, access, and deallocate hypre
     * data structures.
     */
    void allocateHypreData();
    void setMatrixCoefficients();
    void setupHypreSolver();
    bool solveSystem(int x_idx, int b_idx);
    void destroyHypreSolver();
    void deallocateHypreData();

    /*!
     * \brief Associated hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    /*!
     * \brief Associated patch level and C-F boundary (for level numbers > 0).
     */
    int d_level_num = IBTK::invalid_level_number;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > d_level;
    SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM> > d_cf_boundary;

    /*!
     * \name hypre objects.
     */
    //\{
    static const int PART = 0;
    static const int NPARTS = 1;
    static const int NVARS = NDIM;
    static const int X_VAR = 0;
    static const int Y_VAR = 1;
    static const int Z_VAR = 2;

    HYPRE_SStructGrid d_grid = nullptr;
    HYPRE_SStructStencil d_stencil[NVARS];
    HYPRE_SStructGraph d_graph = nullptr;
    HYPRE_SStructMatrix d_matrix = nullptr;
    HYPRE_SStructVector d_rhs_vec = nullptr, d_sol_vec = nullptr;
    HYPRE_SStructSolver d_solver = nullptr, d_precond = nullptr;
    std::vector<SAMRAI::hier::Index<NDIM> > d_stencil_offsets;

    std::string d_solver_type = "Split", d_precond_type = "none", d_split_solver_type = "PFMG";
    int d_rel_change = 0;
    int d_num_pre_relax_steps = 1, d_num_post_relax_steps = 1;
    int d_relax_type;
    int d_skip_relax = 1;
    int d_two_norm = 1;
    //\}
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_SCPoissonHypreLevelSolver
