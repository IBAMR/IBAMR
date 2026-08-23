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

#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesSolver.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/PETScLevelSolver.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petscao.h>
#include <petscis.h>
#include <petscmat.h>
#include <petscvec.h>

#include <CellVariable.h>
#include <IntVector.h>
#include <RefineSchedule.h>
#include <SideVariable.h>
#include <VariableContext.h>

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
 * \see INSStaggeredHierarchyIntegrator
 */
class StaggeredStokesPETScLevelSolver : public IBTK::PETScLevelSolver, public StaggeredStokesSolver
{
public:
    /*!
     * \brief Nonowning read-only view of an initialized level operator state.
     *
     * The PETSc matrix handle and DOF containers remain owned by the solver and
     * are valid only until deallocateSolverState() or reinitialization. PETSc's
     * C API does not provide const matrix handles: callers must not modify or
     * retain the matrix. DOF values use the full coupled velocity-pressure
     * global numbering; the referenced vectors contain the locally owned DOFs.
     */
    struct LiveOperatorStateView
    {
        bool initialized = false;
        Mat operator_mat = nullptr;
        const std::vector<int>* num_dofs_per_proc = nullptr;
        const std::vector<PetscInt>* locally_owned_velocity_dofs = nullptr;
        const std::vector<PetscInt>* locally_owned_pressure_dofs = nullptr;
        int level_number = IBTK::invalid_level_number;
        bool operator_was_provided = false;
        bool includes_augmented_operator = false;
        bool velocity_nullspace_declared = false;
        bool pressure_nullspace_declared = false;
        bool operator_nullspace_attached = false;
        bool zero_mean_pressure_correction = false;
    };

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
     * \brief Set a full-level PETSc operator matrix to use directly instead of
     * rediscretizing the Stokes operator in initializeSolverStateSpecialized().
     *
     * The provided matrix must live in the full coupled Stokes DOF space
     * (velocity+pressure). When no augmented operator is configured, the solver
     * uses this matrix as a nonowning alias until deallocateSolverState(); the
     * caller must retain, and must not modify or reassemble, it while the solver
     * is initialized. An augmented operator requires a private copy since its
     * entries are added during initialization.
     */
    void setOperatorMat(Mat operator_mat);

    /*!
     * \brief Set a matrix contribution to be added to the full level operator
     * matrix assembled in initializeSolverStateSpecialized().
     *
     * The augmentation may be either:
     * 1) a full coupled Stokes matrix contribution (velocity+pressure space),
     * or
     * 2) an A00 velocity-block contribution; this is embedded in the full
     *    coupled matrix by velocity DOF mapping.
     */
    void setAugmentedOperatorMat(Mat augmented_operator_mat);

    /*!
     * \brief Set the live Eulerian elasticity contribution used only for
     * pressure-cell-seeded CAV patch construction.
     *
     * The matrix must use the full coupled level's global velocity-pressure
     * numbering, with numerically zero pressure rows and columns. The solver
     * borrows it until deallocateSolverState() and does not copy or modify it.
     * Local patch matrices and residual updates continue to use the full
     * original level operator.
     */
    void setCouplingAwareASMConstructionMat(Mat construction_mat);

    /*!
     * \brief Return the ordered pressure seed DOFs used by the initialized
     * pressure-cell-seeded CAV construction.
     *
     * The returned reference is invalidated by deallocateSolverState() or
     * reinitialization. It is empty for the legacy velocity-seeded family.
     */
    const std::vector<int>& getCouplingAwareASMPressureSeedDOFs() const;

    /*!
     * \brief Return a coherent nonowning view of the current live operator,
     * block partition, ownership, nullspace, and pressure-gauge declarations.
     *
     * Before initialization and after deallocation this returns an empty view
     * with \c initialized equal to false. This function does not allocate,
     * extract matrix blocks, translate indices, or copy algebraic objects.
     */
    LiveOperatorStateView getLiveOperatorStateView() const;

    /*!
     * \brief Return whether \p dof is classified as a velocity DOF in the
     * current cached level state.
     */
    bool isVelocityDOF(int dof) const;

    /*!
     * \brief Return whether \p dof is classified as a pressure DOF in the
     * current cached level state.
     */
    bool isPressureDOF(int dof) const;

protected:
    /*!
     * \brief Generate IS/subdomains for Schwartz type preconditioners.
     */
    void generateASMSubdomains(std::vector<std::vector<int>>& overlap_dofs,
                               std::vector<std::vector<int>>& nonoverlap_dofs) override;

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

    /*!
     * \brief Postprocess shell preconditioner output.
     */
    void postprocessShellResult(Vec& y) override;

private:
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
    SAMRAI::hier::IntVector<NDIM> d_box_size = 2, d_overlap_size = 1;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, int>> d_u_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_u_nullspace_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, int>> d_p_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_p_nullspace_var;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>> d_data_synch_sched, d_ghost_fill_sched;
    IS d_velocity_field_is_local = nullptr;
    IS d_pressure_is_local = nullptr;
    AO d_velocity_field_ao = nullptr;
    ASMSubdomainConstructionMode d_asm_subdomain_construction_mode = ASMSubdomainConstructionMode::GEOMETRICAL;
    CouplingAwareASMPatchSeedType d_coupling_aware_asm_patch_seed_type =
        CouplingAwareASMPatchSeedType::VELOCITY_COMPONENT;
    int d_coupling_aware_asm_seed_axis = 0;
    int d_coupling_aware_asm_seed_stride = 1;
#if (NDIM == 2)
    CouplingAwareASMSeedTraversalOrder d_coupling_aware_asm_seed_traversal_order =
        CouplingAwareASMSeedTraversalOrder::I_J;
#else
    CouplingAwareASMSeedTraversalOrder d_coupling_aware_asm_seed_traversal_order =
        CouplingAwareASMSeedTraversalOrder::I_J_K;
#endif
    CouplingAwareASMClosurePolicy d_coupling_aware_asm_closure_policy = CouplingAwareASMClosurePolicy::RELAXED;
    StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData d_coupling_aware_asm_map_data;
    bool d_coupling_aware_asm_map_data_is_initialized = false;
    bool d_log_asm_subdomains = false;
    double d_coupling_aware_asm_relative_zero_tol = IBTK_RELATIVE_NUMERICAL_ZERO_TOL;
    Mat d_operator_mat = nullptr;
    Mat d_augmented_operator_mat = nullptr;
    Mat d_coupling_aware_asm_construction_mat = nullptr;
    std::vector<int> d_coupling_aware_asm_pressure_seed_dofs;
    std::vector<PetscInt> d_velocity_dofs;
    std::vector<PetscInt> d_pressure_dofs;

    //\}
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesPETScLevelSolver
