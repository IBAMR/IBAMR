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

#ifndef included_IBAMR_StaggeredStokesPETScLevelSolver
#define included_IBAMR_StaggeredStokesPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/StaggeredStokesSolver.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/PETScLevelSolver.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petscmat.h>
#include <petscvec.h>

#include <CellVariable.h>
#include <IntVector.h>
#include <RefineSchedule.h>
#include <SideVariable.h>
#include <VariableContext.h>

#include <memory>
#include <set>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
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
 * \brief Class StaggeredStokesPETScLevelSolver is a concrete PETScLevelSolver
 * for a staggered-grid (MAC) discretization of the incompressible Stokes
 * equations.
 *
 * The input asm_subdomain_construction_mode selects how the subdomains of the ASM
 * and shell preconditioners are constructed:
 *
 * - GEOMETRICAL (default): subdomains defined by the geometry of the level.
 * - COUPLING_AWARE: subdomains built from the couplings stored in a matrix, on one MPI rank only.
 *
 * Coupling-aware construction uses the following terms. The standard Vanka patch
 * of a cell consists of its pressure DOF and the 2*NDIM velocity DOFs on its
 * faces, which are the DOFs of the cell's discrete divergence equation. A seed
 * is a DOF from which one subdomain is generated. Two velocity DOFs are coupled
 * if the matrix entry between them is larger than a threshold set by
 * relative_zero_tol. Expanding a seed adds the velocity DOFs coupled to it. The
 * closure policy selects the cells whose
 * standard Vanka patches are joined into the subdomain: RELAXED joins every cell
 * incident to an expanded velocity DOF (the closure of Gruninger and Griffith,
 * arXiv:2608.14310), and STRICT joins only cells whose
 * complete velocity stencil is in the expanded set. The construction contract is
 * given by
 * StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(),
 * which reads the couplings from the assembled level operator, including any augmentation.
 * Vanka smoothing needs every DOF to be in a subdomain. If check_subdomain_coverage is
 * TRUE, which is its default in debug builds only, initialization checks this with
 * IBTK::check_dof_coverage(); a seed_stride above 1 can leave DOFs out.
 *
 * The coupling-aware settings are read at construction. Each input key has the
 * prefix coupling_aware_asm_:
 *
 * - seed_axis (default 0): the velocity component whose DOFs are the seeds.
 * - seed_stride (default 1): every seed_stride-th seed is used.
 * - seed_traversal_order (default I_J in 2D and I_J_K in 3D): the order in which
 *   the seeds are visited, from the slowest to the fastest varying logical
 *   coordinate. J_I is the alternative in 2D, and J_K_I and K_I_J are the
 *   alternatives in 3D.
 * - closure_policy (default RELAXED): RELAXED or STRICT.
 * - relative_zero_tol (default 1.0e-14): the relative threshold below which
 *   couplings are ignored, as defined by the construction contract.
 *
 * \see INSStaggeredHierarchyIntegrator
 *
 * The subdomain solver "eigen-schur-complement" solves the local problems of the shell
 * preconditioners with the Schur complement of the pressure block, and uses the
 * velocity and pressure fields of the level. a00_solver_type and schur_solver_type
 * default to FULL_PIV_HOUSEHOLDER_QR; a00_solver_threshold and schur_solver_threshold
 * default to -1. The types and thresholds are those of IBTK::make_eigen_subdomain_solver().
 * The solve matrices of A00 and of the Schur complement are precomputed.
 */
class StaggeredStokesPETScLevelSolver : public IBTK::PETScLevelSolver, public StaggeredStokesSolver
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesPETScLevelSolver(const std::string& object_name,
                                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                    const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesPETScLevelSolver();

    /*!
     * \brief Static function to construct a StaggeredStokesPETScLevelSolver.
     */
    static SAMRAI::tbox::Pointer<StaggeredStokesSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        return new StaggeredStokesPETScLevelSolver(object_name, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \brief Set a full-level PETSc matrix instead of rediscretizing the Stokes operator.
     *
     * The assembled matrix must use the full coupled velocity-pressure numbering
     * and local row/column distribution defined by
     * StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices() for this
     * level, on PETSC_COMM_WORLD. Its boundary treatment must match the solver's
     * boundary configuration, and any configured nullspace must be a nullspace
     * of the supplied system.
     * Without augmentation the solver uses the exact supplied matrix handle,
     * and initialization attaches the configured nullspace, if any, to that
     * matrix; deallocateSolverState() removes it and restores the nullspace that
     * the matrix had before initialization, if any.
     * Passing nullptr restores rediscretization without clearing an augmentation.
     *
     * The solver retains a PETSc reference, so the caller may release its own
     * reference after this call. The matrix must not be modified or reassembled
     * through any alias while installed. The retained reference survives
     * deallocateSolverState() and is released on replacement, clearing with
     * nullptr, or destruction of the solver.
     * Set, replace, or clear this handle only while the solver is deallocated.
     */
    void setOperatorMat(Mat operator_mat);

    /*!
     * \brief Set a matrix contribution to add to the supplied or rediscretized level operator.
     *
     * A caller's contribution is typically a purely velocity-space term (for example a structural
     * coupling stiffness) added to the velocity-velocity block of the saddle-point system
     * [A_0 + A_aug, B; B^T, 0]; supplying it in compact velocity-only numbering saves such a caller
     * from embedding it into the full coupled system itself. The assembled contribution may use
     * either the full coupled numbering and distribution in setOperatorMat(), or that compact
     * velocity-only numbering.
     * For the latter, each rank owns a contiguous range with one row per locally
     * owned velocity DOF, ordered by increasing coupled global index within the
     * velocity field from StaggeredStokesPETScMatUtilities::constructPatchLevelFields().
     * Columns use the same global compact-to-coupled mapping. All contributions
     * must use PETSC_COMM_WORLD. The compact velocity matrix must support PETSc
     * row access, and the base matrix must support addition of the full-system
     * or embedded velocity contribution.
     *
     * Addition preserves both installed matrices' entries. Configure boundary
     * conditions and nullspaces consistently with the resulting system.
     * Passing nullptr clears only the augmentation. The installed-reference,
     * immutability, and deallocated-state requirements of setOperatorMat() apply
     * independently to this handle.
     */
    void setAugmentedOperatorMat(Mat augmented_operator_mat);

protected:
    /*!
     * \brief Require pc_type = asm or shell when the ASM subdomains are coupling-aware.
     */
    void validatePreconditionerType() override;

    /*!
     * \brief Generate IS/subdomains for Schwarz type preconditioners.
     */
    void generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                               std::vector<std::set<int>>& nonoverlap_is) override;

    /*!
     * \brief Generate IS/subdomains for fieldsplit type preconditioners.
     */
    void generateFieldSplitSubdomains(std::vector<std::string>& field_names,
                                      std::vector<std::set<int>>& field_is) override;

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     */
    void initializeSolverStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                          const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverStateSpecialized().
     */
    void deallocateSolverStateSpecialized() override;

    /*!
     * \brief Copy a generic vector to the PETSc representation.
     */
    void copyToPETScVec(Vec& petsc_x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) override;

    /*!
     * \brief Copy a generic vector from the PETSc representation.
     */
    void copyFromPETScVec(Vec& petsc_x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) override;

    /*!
     * \brief Copy solution and right-hand-side data to the PETSc
     * representation, including any modifications to account for boundary
     * conditions.
     */
    void setupKSPVecs(Vec& petsc_x,
                      Vec& petsc_b,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

private:
    /*!
     * \brief Return the subdomain solver "eigen-schur-complement", with its settings from input_db. It obtains the
     * velocity and pressure fields of the level each time the solver state is initialized.
     */
    IBTK::PETScLevelSolverSubdomainSolver
    makeEigenSchurComplementSubdomainSolver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesPETScLevelSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesPETScLevelSolver(const StaggeredStokesPETScLevelSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesPETScLevelSolver& operator=(const StaggeredStokesPETScLevelSolver& that) = delete;

    /*!
     * \name PETSc objects.
     */
    //\{

    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    std::vector<int> d_num_dofs_per_proc;
    int d_u_dof_index_idx = IBTK::invalid_index, d_p_dof_index_idx = IBTK::invalid_index;
    int d_u_nullspace_idx = IBTK::invalid_index, d_p_nullspace_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, int>> d_u_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_u_nullspace_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, int>> d_p_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_p_nullspace_var;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>> d_data_synch_sched, d_ghost_fill_sched;
    ASMSubdomainConstructionMode d_asm_mode = ASMSubdomainConstructionMode::GEOMETRICAL;
    int d_ca_seed_axis = 0, d_ca_seed_stride = 1;
#if (NDIM == 2)
    CouplingAwareASMSeedTraversalOrder d_ca_order = CouplingAwareASMSeedTraversalOrder::I_J;
#else
    CouplingAwareASMSeedTraversalOrder d_ca_order = CouplingAwareASMSeedTraversalOrder::I_J_K;
#endif
    CouplingAwareASMClosurePolicy d_ca_policy = CouplingAwareASMClosurePolicy::RELAXED;
    double d_ca_relative_zero_tol = 1.0e-14;
    // Owned references to the installed inputs, independent of solver state.
    Mat d_operator_mat = nullptr;
    Mat d_augmented_operator_mat = nullptr;
    //! The caller's nullspace of a supplied matrix that the solver uses directly, restored on deallocation.
    MatNullSpace d_operator_mat_prior_nullsp = nullptr;

    //\}
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesPETScLevelSolver
