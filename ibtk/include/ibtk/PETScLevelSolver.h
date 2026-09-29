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

 options_prefix = "..."                    // the default is chosen by the derived class; see setOptionsPrefix()
 ksp_type = "gmres"                        // see setKSPType()
 pc_type = "ilu"                           // the PETSc preconditioner type
 initial_guess_nonzero = TRUE              // see setInitialGuessNonzero()
 rel_residual_tol = 1.0e-5                 // see setRelativeTolerance()
 abs_residual_tol = 1.0e-50                // see setAbsoluteTolerance()
 max_iterations = 10000                    // see setMaxIterations()
 enable_logging = FALSE                    // see setLoggingEnabled()
 subdomain_box_size = 2, 2                 // the size of the ASM subdomains, one entry per direction
 subdomain_overlap_size = 1, 1             // the overlap of the ASM subdomains, one entry per direction
 check_subdomain_coverage = FALSE          // whether to check that the subdomains cover the DOFs; TRUE in debug builds
 subdomain_relaxation {                    // see "Subdomain relaxation" below
    composition = "MULTIPLICATIVE"         // no default: "ADDITIVE" or "MULTIPLICATIVE"
    grouping = "RANK"                      // no default: "RANK" or "SAMRAI_PATCH"; only for MULTIPLICATIVE
    output = "FULL"                        // no default: "FULL" or "OWNED"
    traversal = "FORWARD"                  // "FORWARD" (default), "REVERSE", or "SYMMETRIC"
    subdomain_solver {                     // optional; see "Subdomain solvers" below
       type = "petsc"                      // the default
    }
 }
 \endverbatim
 *
 * <b>Subdomain relaxation</b>
 *
 * With pc_type = "shell", IBAMR configures subdomain relaxation using the subdomain_relaxation database.
 * A subdomain is a set of DOFs on which a local problem is solved; generateASMSubdomains() forms the
 * subdomains of each rank. With A the level operator, r the input of the preconditioner, and P_i the DOFs of
 * subdomain i, a solve computes the correction d_i = A(P_i, P_i)^{-1} q_i:
 *
 * - composition = "ADDITIVE": each subdomain is solved independently, with q_i = r(P_i).
 * - composition = "MULTIPLICATIVE": the subdomains of each group are solved one after another. The
 *   correction z of the group starts at zero, and each solve uses the residual that the previous solves of
 *   the group leave, q_i = r(P_i) - A(P_i, :) z, and then adds its correction, z(P_i) += d_i. A group does
 *   not use the corrections of other groups. grouping = "RANK" makes the subdomains of each rank one group,
 *   and grouping = "SAMRAI_PATCH" forms one group for each SAMRAI patch (a patch of the grid) with
 *   generateSubdomainGroups(), which only some derived classes provide. A subdomain may be in several
 *   groups.
 *   traversal sets the order in which a group visits its subdomains: FORWARD (default), REVERSE, or
 *   SYMMETRIC, which visits them forward and then backward without repeating the last one. ADDITIVE
 *   composition accepts only FORWARD.
 *
 * output chooses how the corrections make up the result:
 *
 * - output = "FULL": the result is the sum of the corrections of the subdomains (ADDITIVE) or of the groups
 *   (MULTIPLICATIVE), including their overlapping entries and the entries of other ranks.
 * - output = "OWNED": each DOF takes its value from the one correction that owns it. A subdomain owns the
 *   nonoverlapping set that generateASMSubdomains() gives it, and a group owns those of the subdomains that
 *   it owns: all of its subdomains with RANK grouping, or those that generateSubdomainGroups() assigns to
 *   it with SAMRAI_PATCH grouping. The
 *   owned DOFs must partition the DOFs of each rank, which is always checked. ADDITIVE composition with
 *   OWNED output is restricted additive Schwarz.
 *
 * composition and output must be set, and grouping as well for MULTIPLICATIVE composition, whenever the
 * preconditioner is a shell, including when pc_type = "shell" comes from the PETSc options. For example,
 * \verbatim
 subdomain_relaxation {
    composition = "MULTIPLICATIVE"
    grouping = "RANK"
    output = "FULL"
 }
 \endverbatim
 *
 * <b>Subdomain solvers</b>
 *
 * The optional database subdomain_relaxation.subdomain_solver selects the subdomain solver by its type
 * and holds its settings. The type is "petsc" (default), "blas-lapack", "eigen", or
 * "eigen-pseudoinverse"; see make_petsc_subdomain_solver(), make_blas_lapack_subdomain_solver(),
 * make_eigen_subdomain_solver(), and make_eigen_pseudoinverse_subdomain_solver() for their settings. The
 * type can also name one of the SubdomainSolverFactories that a derived class or its user supplies, and
 * setSubdomainSolver() supplies another subdomain solver.
 *
 * PETSc is developed at the Argonne National Laboratory Mathematics and
 * Computer Science Division.  For more information about \em PETSc, see <A
 * HREF="http://www.mcs.anl.gov/petsc">http://www.mcs.anl.gov/petsc</A>.
 */
class PETScLevelSolver : public LinearSolver
{
public:
    /*!
     * \brief A function that creates a subdomain solver from the subdomain_solver database, which may be null.
     *
     * The function validates the settings of the database, which include the type that selected it and may include
     * petsc_settings databases. The level solver keeps a copy of the selected function until it first sets up a
     * shell preconditioner. The subdomain solver owns whatever it needs after the function returns, or borrows only
     * objects that outlive it.
     */
    using SubdomainSolverFactory =
        std::function<PETScLevelSolverSubdomainSolver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>)>;

    /*!
     * \brief Names and factories of subdomain solvers that the type of the subdomain_solver database can select
     * in addition to the built-in ones. Names are compared without regard to case.
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
     * \brief Use subdomain_solver for the local problems of subdomain relaxation, in place
     * of the one that the subdomain_solver database selects.
     *
     * The solver takes ownership of subdomain_solver, which must not be empty, and
     * retains it across reinitialization of the solver state. It takes precedence
     * over the subdomain_solver database. It is initialized
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
     * Reads the settings from input_db and selects the factory of the subdomain solver that the type of the
     * subdomain_solver database names: a built-in one or the entry of subdomain_solver_factories with that name.
     * It is an error for a name to be empty, to be that of a built-in subdomain solver, or to be supplied twice, for
     * a function to be empty, and for the type to name none of these. The factory is called only when a shell
     * preconditioner is set up without a subdomain solver, and it is an error for it to return an empty one.
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
     * the same subdomain. Subdomain relaxation with OWNED output and the restricted ASM
     * preconditioner also need the nonoverlapping sets of this rank to partition the DOFs that it
     * owns. Initialization reports a violation of these requirements. A rank may have no
     * subdomains.
     */
    virtual void generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                                       std::vector<std::set<int>>& nonoverlap_is);

    /*!
     * \brief Validate the preconditioner type in effect for this initialization.
     *
     * Called by initializeSolverState() after the PETSc options database has been applied, so that
     * d_pc_type is the type that will be used, including a command-line override of the type from
     * the input database. The default accepts every type; a derived class that requires a particular
     * type reports a violation.
     */
    virtual void validatePreconditionerType()
    {
    }

    /*!
     * \brief Generate the groups of multiplicative subdomain relaxation with SAMRAI_PATCH grouping.
     *
     * group_subdomains[g] lists the subdomains of this rank that group g solves, in the order of a FORWARD
     * traversal; a subdomain may be in several groups. owning_groups[i] is the group that owns subdomain i and
     * writes its nonoverlapping DOFs with OWNED output, which must be one of the groups that solve it.
     * Initialization reports a violation. The default reports that SAMRAI_PATCH grouping is not supported.
     */
    virtual void generateSubdomainGroups(std::vector<std::vector<int>>& group_subdomains,
                                         std::vector<int>& owning_groups);

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
    //! How subdomain relaxation composes the solves of its subdomains.
    enum class SubdomainComposition
    {
        ADDITIVE,
        MULTIPLICATIVE
    };
    //! How multiplicative subdomain relaxation groups the subdomains.
    enum class SubdomainGrouping
    {
        RANK,
        SAMRAI_PATCH
    };
    //! The order in which multiplicative subdomain relaxation visits the subdomains of a group.
    enum class SubdomainTraversal
    {
        FORWARD,
        REVERSE,
        SYMMETRIC
    };
    //! Which entries of the corrections subdomain relaxation adds to its result.
    enum class SubdomainOutput
    {
        FULL,
        OWNED
    };
    //! d_pc_type is the preconditioner type in effect, which the PETSc options may override at initialization.
    std::string d_ksp_type = KSPGMRES, d_pc_type = PCILU;
    std::string d_options_prefix;
    //! Set from the subdomain_relaxation database; required only when the preconditioner is a shell, possibly
    //! selected through PETSc options.
    std::optional<SubdomainComposition> d_subdomain_composition;
    std::optional<SubdomainGrouping> d_subdomain_grouping;
    std::optional<SubdomainOutput> d_subdomain_output;
    SubdomainTraversal d_subdomain_traversal = SubdomainTraversal::FORWARD;
    //! Whether initialization checks that the subdomains cover the DOFs.
    bool d_check_subdomain_coverage = default_check_dof_coverage();
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
     * \name ASM subdomains and the storage of subdomain relaxation.
     */
    //\{
    SAMRAI::hier::IntVector<NDIM> d_box_size, d_overlap_size;
    int d_n_local_subdomains = 0;
    std::vector<IS> d_overlap_is, d_nonoverlap_is;
    Mat* d_sub_mat = nullptr;

    /*!
     * The right-hand sides and solutions of the local problems of the subdomains of this rank, packed
     * in order into sequential vectors: subdomain i occupies the entries from d_subdomain_offsets[i] up
     * to d_subdomain_offsets[i + 1]. One scatter gathers all of the right-hand sides.
     */
    Vec d_subdomain_rhs = nullptr, d_subdomain_solution = nullptr;
    VecScatter d_restriction = nullptr;
    //! Whether gathering the right-hand sides needs the scatter, or only the local array at these indices.
    bool d_restriction_communicates = true;
    std::vector<PetscInt> d_gather_indices;
    std::vector<PetscInt> d_subdomain_offsets;

    /*!
     * For ADDITIVE composition with OWNED output, the entries of the packed solutions that a subdomain
     * writes to the output, from d_write_offsets[i] up to d_write_offsets[i + 1] in d_write_sources, the
     * packed positions of the nonoverlapping DOFs, and d_write_targets, their local indices in the output.
     */
    std::vector<PetscInt> d_write_offsets, d_write_sources, d_write_targets;
    //\}

    /*!
     * \name Subdomain residuals of multiplicative subdomain relaxation.
     *
     * For subdomain i, d_residual_matrices[i] is -A(P_i, C_i), in which C_i is the sorted list of the
     * columns that the rows of P_i couple to, and d_halo_vectors[i] holds the current correction at C_i.
     * The residual of subdomain i is formed in its view of the packed residuals from its view of the packed
     * right-hand sides.
     */
    //\{
    Vec d_subdomain_residual = nullptr;
    std::vector<Vec> d_subdomain_rhs_views, d_subdomain_residual_views;
    std::vector<Mat> d_residual_matrices;
    std::vector<Vec> d_halo_vectors;
    //\}

    /*!
     * \name Field split preconditioner.
     */
    //\{
    std::vector<std::string> d_field_name;
    std::vector<IS> d_field_is;
    //\}

private:
    //! The preconditioner type that the input selects, which each initialization applies before the PETSc options.
    std::string d_selected_pc_type = PCILU;
    //! The subdomain_relaxation.subdomain_solver database, if there is one.
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_subdomain_solver_db;
    //! The type of the subdomain solver and the factory that creates it from d_subdomain_solver_db.
    std::string d_subdomain_solver_type = "petsc";
    SubdomainSolverFactory d_subdomain_solver_factory;
    //! The subdomain solver of subdomain relaxation.
    std::optional<PETScLevelSolverSubdomainSolver> d_subdomain_solver;
    //! Whether d_subdomain_solver is initialized for the current solver state.
    bool d_subdomain_solver_initialized = false;

    /*!
     * \name Groups of multiplicative subdomain relaxation.
     *
     * The corrections are stored at the positions of their DOFs in the sorted list of the distinct DOFs of the
     * subdomains of this rank. Entry k of the packed vectors has the DOF at position d_packed_positions[k], and
     * d_position_slots[p] is the first entry of the packed vectors with the DOF at position p. Column k of the
     * residual matrix of subdomain i is at position d_halo_positions[d_halo_offsets[i] + k], or -1 if it is not
     * in the list.
     *
     * Group g visits the subdomains d_group_visits[j] for j from d_group_visit_offsets[g] up to
     * d_group_visit_offsets[g + 1], and its correction is zero except at the positions
     * d_group_support[j] for j from d_group_support_offsets[g] up to d_group_support_offsets[g + 1]. With OWNED
     * output, the group writes the entries d_group_output_sources[j] of its correction to the local entries
     * d_group_output_targets[j] of the output, for j from d_group_output_offsets[g] up to
     * d_group_output_offsets[g + 1]. d_group_correction is the correction of the group being applied, and
     * d_full_correction sums the corrections of the groups for FULL output.
     */
    //\{
    std::vector<PetscInt> d_packed_positions, d_position_slots, d_halo_offsets, d_halo_positions;
    std::vector<int> d_group_visit_offsets, d_group_visits;
    std::vector<PetscInt> d_group_support_offsets, d_group_support, d_group_output_offsets, d_group_output_sources,
        d_group_output_targets;
    std::vector<PetscScalar> d_group_correction, d_full_correction;
    //\}

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
     * \brief Set up the output of ADDITIVE composition with OWNED output, for the DOFs of this rank from n_lo up
     * to n_hi.
     */
    void initializeOwnedSubdomainOutput(PetscInt n_lo, PetscInt n_hi);

    /*!
     * \brief Set up the residuals of the subdomains of multiplicative subdomain relaxation from rows, the rows
     * A(P_i, :) of each subdomain with all of the columns, and dofs, the sorted list of the distinct DOFs of the
     * subdomains of this rank.
     */
    void initializeSubdomainResiduals(Mat* rows, const std::vector<PetscInt>& dofs);

    /*!
     * \brief Set up the groups of multiplicative subdomain relaxation from gathered_indices, the DOFs of the
     * entries of the packed vectors, and dofs, the sorted list of the distinct DOFs of the subdomains of this rank,
     * for the DOFs of this rank from n_lo up to n_hi.
     */
    void initializeSubdomainGroups(const std::vector<PetscInt>& gathered_indices,
                                   const std::vector<PetscInt>& dofs,
                                   PetscInt n_lo,
                                   PetscInt n_hi);

    /*!
     * \brief Return the positions, in the order in which a group visits them, for a traversal of its n
     * subdomains: FORWARD visits 0, ..., n - 1, REVERSE visits n - 1, ..., 0, and SYMMETRIC visits
     * 0, ..., n - 1, n - 2, ..., 0, the last subdomain once.
     */
    static std::vector<int> subdomainVisitOrder(SubdomainTraversal traversal, int n);

    /*!
     * \brief Gather the right-hand sides of all subdomains from x into the packed vector.
     */
    PetscErrorCode gatherSubdomainRhs(Vec x) const;

    /*!
     * \brief Write the nonoverlapping parts of the packed solutions of subdomains first, ..., last - 1
     * to y.
     */
    PetscErrorCode writeSubdomainSolutions(int first, int last, Vec y) const;

    /*!
     * \brief Add every entry of the packed solutions to y at its DOF, including the DOFs of other ranks.
     */
    PetscErrorCode addPackedSolutions(Vec y) const;

    /*!
     * \brief Apply ADDITIVE subdomain relaxation to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_Additive(PC pc, Vec x, Vec y);

    /*!
     * \brief Apply MULTIPLICATIVE subdomain relaxation to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_Multiplicative(PC pc, Vec x, Vec y);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PETScLevelSolver
