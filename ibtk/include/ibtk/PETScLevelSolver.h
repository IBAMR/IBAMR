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

#ifndef included_IBTK_PETScLevelSolver
#define included_IBTK_PETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/DOFCoverage.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/PETScLevelSolverSubdomainSolver.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include <CoarseFineBoundary.h>
#include <IntVector.h>
#include <PatchHierarchy.h>
#include <SAMRAIVectorReal.h>

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScLevelSolver is an abstract LinearSolver for solving systems
 * of linear equations on a \em single SAMRAI::hier::PatchLevel using <A
 * HREF="http://www.mcs.anl.gov/petsc/petsc-as">PETSc</A>.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 options_prefix = ""                       // see setOptionsPrefix()
 ksp_type = "gmres"                        // see setKSPType()
 pc_type = "ilu"                           // the PETSc preconditioner type
 initial_guess_nonzero = TRUE              // see setInitialGuessNonzero()
 rel_residual_tol = 1.0e-5                 // see setRelativeTolerance()
 abs_residual_tol = 1.0e-50                // see setAbsoluteTolerance()
 max_iterations = 10000                    // see setMaxIterations()
 enable_logging = FALSE                    // see setLoggingEnabled()
 subdomain_box_size = 2, 2                 // the size of the ASM subdomains, one entry per direction
 subdomain_overlap_size = 1, 1             // the overlap of the ASM subdomains, one entry per direction
 shell_pc_type = "additive"                // no default; see "Shell preconditioners" below
 subdomain_solver = "petsc"                // see "Subdomain solvers" below
 check_subdomain_coverage = FALSE          // TRUE by default in debug builds
 \endverbatim
 *
 * <b>Shell preconditioners</b>
 *
 * With pc_type = "shell", this class applies a Schwarz preconditioner itself, on the subdomains from
 * generateASMSubdomains(). shell_pc_type chooses the method and must be set, even when pc_type =
 * "shell" comes from the PETSc options:
 *
 * - "additive": restricted additive Schwarz. Every subdomain is solved with the same right-hand side,
 *   and each keeps its solution only on its nonoverlapping subset, so those subsets must partition the
 *   DOFs (checked when check_subdomain_coverage is TRUE).
 * - "multiplicative": multiplicative Schwarz. The subdomains are solved one after another, each with
 *   the residual left by the previous solves; ranks work through their own subdomains in parallel.
 *
 * <b>Subdomain solvers</b>
 *
 * subdomain_solver chooses how each subdomain problem is solved: "petsc" (default) or "blas-lapack";
 * see make_petsc_subdomain_solver() and make_blas_lapack_subdomain_solver() for their settings.
 * setSubdomainSolver() or a derived class can supply others.
 *
 * PETSc is developed at the Argonne National Laboratory Mathematics and
 * Computer Science Division.  For more information about \em PETSc, see <A
 * HREF="http://www.mcs.anl.gov/petsc">http://www.mcs.anl.gov/petsc</A>.
 */
class PETScLevelSolver : public LinearSolver
{
public:
    /*!
     * \brief A function that creates a subdomain solver, with its settings read from the input database of the
     * level solver, which may be null.
     *
     * The subdomain solver owns whatever it needs after the function returns, or borrows only objects that outlive
     * it: the factory itself is not retained.
     */
    using SubdomainSolverFactory =
        std::function<PETScLevelSolverSubdomainSolver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>)>;

    /*!
     * \brief Names and factories of subdomain solvers that subdomain_solver can select in addition to the built-in
     * ones. Names are compared without regard to case.
     */
    using SubdomainSolverFactories = std::vector<std::pair<std::string, SubdomainSolverFactory>>;

    /*!
     * \brief Default constructor.
     */
    PETScLevelSolver();

    /*!
     * \brief Destructor.
     */
    ~PETScLevelSolver();

    /*!
     * \brief Set the KSP type.
     */
    void setKSPType(const std::string& ksp_type);

    /*!
     * \brief Set the options prefix used by this PETSc solver object.
     */
    void setOptionsPrefix(const std::string& options_prefix);

    /*!
     * \brief Use subdomain_solver for the local problems of the shell preconditioners,
     * in place of the built-in subdomain solver that uses PETSc.
     *
     * The solver takes ownership of subdomain_solver, which must not be empty, and
     * retains it across reinitialization of the solver state. It is initialized
     * and deallocated only when pc_type = "shell". Call this before initializing
     * the level solver, or after calling its deallocateSolverState().
     */
    void setSubdomainSolver(PETScLevelSolverSubdomainSolver subdomain_solver);

    /*!
     * \brief Get the PETSc KSP object.
     */
    const KSP& getPETScKSP() const;

    /*!
     * \brief Get ASM subdomains.
     */
    void getASMSubdomains(std::vector<IS>** nonoverlapping_subdomains, std::vector<IS>** overlapping_subdomains);

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Set the nullspace of the linear system.
     */
    void setNullSpace(
        bool contains_constant_vec,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>>& nullspace_basis_vecs =
            std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>>()) override;

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
     * \note Subclasses of class PETScLevelSolver should \em not override this
     * method.  Instead, they should override the protected method
     * initializeSolverStateSpecialized().
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
     * \note Subclasses of class PETScLevelSolver should \em not override this
     * method.  Instead, they should override the protected method
     * deallocatedSolverStateSpecialized().
     *
     * \see initializeSolverState
     */
    void deallocateSolverState() override;

    //\}

protected:
    /*!
     * \brief Basic initialization.
     *
     * Reads the settings from input_db and creates the subdomain solver that subdomain_solver names: a built-in one
     * or the one of the entry of subdomain_solver_factories with that name. It is an error for a name to be empty, to
     * be that of a built-in subdomain solver, or to be supplied twice, for a function to be empty or to return an
     * empty subdomain solver, and for subdomain_solver to name none of these.
     */
    void init(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
              const std::string& default_options_prefix,
              const SubdomainSolverFactories& subdomain_solver_factories = {});

    /*!
     * \brief Generate IS/subdomains for Schwarz type preconditioners.
     *
     * The subdomains of this rank are listed in order, and overlap_is[i] and nonoverlap_is[i]
     * describe the same subdomain with global DOF indices. Each overlapping set contains the
     * DOFs of its subdomain, and each nonoverlapping set is a subset of the overlapping set of
     * the same subdomain. The additive shell preconditioner and the restricted ASM
     * preconditioner also need the nonoverlapping sets of this rank to partition the DOFs that it
     * owns. Initialization reports a violation of these requirements. A rank may have no
     * subdomains.
     */
    virtual void generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                                       std::vector<std::set<int>>& nonoverlap_is);

    /*!
     * \brief Generate IS/subdomains for fieldsplit type preconditioners.
     */
    virtual void generateFieldSplitSubdomains(std::vector<std::string>& field_names,
                                              std::vector<std::set<int>>& field_is);

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     */
    virtual void initializeSolverStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                                  const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) = 0;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverStateSpecialized().
     */
    virtual void deallocateSolverStateSpecialized() = 0;

    /*!
     * \brief Copy a generic vector to the PETSc representation.
     */
    virtual void copyToPETScVec(Vec& petsc_x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) = 0;

    /*!
     * \brief Copy a generic vector from the PETSc representation.
     */
    virtual void copyFromPETScVec(Vec& petsc_x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) = 0;

    /*!
     * \brief Copy solution and right-hand-side data to the PETSc
     * representation, including any modifications to account for boundary
     * conditions.
     */
    virtual void setupKSPVecs(Vec& petsc_x,
                              Vec& petsc_b,
                              SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                              SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) = 0;

    /*!
     * \brief Setup the solver nullspace (if any).
     */
    virtual void setupNullSpace();

    /*!
     * \brief Associated hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    /*!
     * \brief Associated patch level and C-F boundary (for level numbers > 0).
     */
    int d_level_num = IBTK::invalid_level_number;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> d_level;
    SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM>> d_cf_boundary;

    /*!
     * \brief Scratch data.
     */
    SAMRAIDataCache d_cached_eulerian_data;

    /*!
     * \name Solver settings.
     */
    //\{
    //! How a shell preconditioner composes the corrections of its subdomains.
    enum class ShellComposition
    {
        ADDITIVE,
        MULTIPLICATIVE
    };
    std::string d_ksp_type = KSPGMRES, d_pc_type = PCILU;
    std::string d_options_prefix;
    //! Set from shell_pc_type; required only when a shell preconditioner is selected, possibly through PETSc options.
    std::optional<ShellComposition> d_shell_composition;
    //\}

    /*!
     * \name PETSc objects.
     */
    //\{
    KSP d_petsc_ksp = nullptr;
    Mat d_petsc_mat = nullptr, d_petsc_pc = nullptr;
    MatNullSpace d_petsc_nullsp = nullptr;
    Vec d_petsc_x = nullptr, d_petsc_b = nullptr;
    //\}

    /*!
     * \name ASM subdomains and the storage of the shell preconditioners.
     */
    //\{
    SAMRAI::hier::IntVector<NDIM> d_box_size, d_overlap_size;
    int d_n_local_subdomains = 0;
    std::vector<IS> d_overlap_is, d_nonoverlap_is;
    Mat* d_sub_mat = nullptr;

    /*!
     * The right-hand sides and solutions of the local problems of the subdomains of this rank, packed
     * in order into sequential vectors: subdomain i occupies the entries from d_subdomain_offsets[i] up
     * to d_subdomain_offsets[i + 1]. One scatter gathers all of the right-hand sides. The vectors are
     * host memory (VECSEQ), and the vectors that view the entries of one subdomain share their arrays
     * for as long as the packed vectors exist.
     */
    Vec d_subdomain_rhs = nullptr, d_subdomain_solution = nullptr;
    VecScatter d_restriction = nullptr;
    //! Whether gathering the right-hand sides needs the scatter, or only the local array at these indices.
    bool d_restriction_communicates = true;
    std::vector<PetscInt> d_gather_indices;
    std::vector<PetscInt> d_subdomain_offsets;

    /*!
     * The entries of the packed solutions that a subdomain writes to the output, from
     * d_write_offsets[i] up to d_write_offsets[i + 1] in d_write_sources, the packed positions of the
     * nonoverlapping DOFs, and d_write_targets, their local indices in the output.
     */
    std::vector<PetscInt> d_write_offsets, d_write_sources, d_write_targets;
    //\}

    /*!
     * \name Multiplicative shell preconditioner.
     *
     * The multiplicative shell visits stages 0, ..., d_n_stages - 1, the largest number of subdomains on any
     * rank, and uses subdomain i of this rank at stage i, if there is one. For subdomain i,
     * d_residual_matrices[i] is the matrix -A(O_i, C_i) of its rows and the columns C_i that they couple to.
     * For stage s, d_halo_vectors[s] holds the current output on C_i for the subdomain i of the stage, and
     * the scatters of the stage gather those values and add the corrections from O_i to the output.
     */
    //\{
    int d_n_stages = 0;
    //! Vectors that share the storage of each subdomain's packed right-hand side and solution.
    std::vector<Vec> d_subdomain_rhs_views, d_subdomain_solution_views;
    std::vector<Mat> d_residual_matrices;
    std::vector<Vec> d_halo_vectors;
    std::vector<VecScatter> d_halo_scatters, d_correction_scatters;
    //! Whether each stage needs its scatter, or works on the local indices of the output, which are also kept.
    std::vector<bool> d_halo_communicates, d_correction_communicates;
    std::vector<std::vector<PetscInt>> d_halo_local_indices, d_correction_local_indices;
    //\}

    /*!
     * \name Field split preconditioner.
     */
    //\{
    std::vector<std::string> d_field_name;
    std::vector<IS> d_field_is;
    //\}

private:
    //! The subdomain_solver input.
    std::string d_subdomain_solver_type = "petsc";
    //! The subdomain solver of the shell preconditioners.
    std::optional<PETScLevelSolverSubdomainSolver> d_subdomain_solver;
    //! Whether d_subdomain_solver is initialized for the current solver state.
    bool d_subdomain_solver_initialized = false;
    /*!
     * Whether initializeSolverState() converted std::set<int> subdomains from
     * generateASMSubdomains() into d_overlap_is and d_nonoverlap_is itself, so that
     * deallocateSolverState() should destroy and clear them. Subclasses that construct
     * PETSc index sets directly instead manage their own regeneration.
     */
    bool d_generated_subdomain_is = false;
    //! Whether initialization checks that the subdomains cover the DOFs as the preconditioner requires.
    bool d_check_subdomain_coverage = default_check_dof_coverage();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is deleted.
     *
     * \param from The value to copy to this object.
     */
    PETScLevelSolver(const PETScLevelSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is deleted.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScLevelSolver& operator=(const PETScLevelSolver& that) = delete;

    /*!
     * \brief Gather the right-hand sides of all subdomains from x into the packed vector.
     */
    PetscErrorCode gatherSubdomainRhs(Vec x) const;

    /*!
     * \brief Set up the stages of the multiplicative shell from the rows of the operator on each
     * overlapping subdomain, with every column.
     */
    void initializeMultiplicativeShell(Mat* rows);

    /*!
     * \brief Write the nonoverlapping parts of the packed solutions of subdomains first, ..., last - 1
     * to y.
     */
    PetscErrorCode writeSubdomainSolutions(int first, int last, Vec y) const;

    /*!
     * \brief Apply the additive shell preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_Additive(PC pc, Vec x, Vec y);

    /*!
     * \brief Apply the multiplicative shell preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_Multiplicative(PC pc, Vec x, Vec y);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PETScLevelSolver
