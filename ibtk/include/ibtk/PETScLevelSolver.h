// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Sparse>

#include <CoarseFineBoundary.h>
#include <IntVector.h>
#include <PatchHierarchy.h>
#include <SAMRAIVectorReal.h>

#include <memory>
#include <set>
#include <string>
#include <unordered_map>
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

 options_prefix = ""           // see setOptionsPrefix()
 ksp_type = "gmres"            // see setKSPType()
 initial_guess_nonzero = TRUE  // see setInitialGuessNonzero()
 rel_residual_tol = 1.0e-5     // see setRelativeTolerance()
 abs_residual_tol = 1.0e-50    // see setAbsoluteTolerance()
 max_iterations = 10000        // see setMaxIterations()
 enable_logging = FALSE        // see setLoggingEnabled()
 \endverbatim
 *
 * PETSc is developed at the Argonne National Laboratory Mathematics and
 * Computer Science Division.  For more information about \em PETSc, see <A
 * HREF="http://www.mcs.anl.gov/petsc">http://www.mcs.anl.gov/petsc</A>.
 */
class PETScLevelSolver : public LinearSolver
{
public:
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
     * \brief Get the PETSc KSP object.
     */
    const KSP& getPETScKSP() const;

    /*!
     * \brief Get the stored ASM-like subdomain description.
     *
     * The returned \p subdomain_dofs define the overlapping subdomains used by
     * ASM-like methods. The returned \p nonoverlap_subdomain_dofs define
     * nonoverlapping subsets of those subdomains whose union covers the local
     * domain.
     */
    void getASMSubdomains(std::vector<std::vector<int>>** nonoverlap_subdomain_dofs,
                          std::vector<std::vector<int>>** subdomain_dofs);

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
     */
    void init(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, const std::string& default_options_prefix);

    /*!
     * \brief Generate overlapping subdomains and nonoverlapping partition subsets
     * for ASM-like preconditioners.
     *
     * The \p subdomain_dofs argument stores the actual overlapping subdomains.
     * The \p nonoverlap_subdomain_dofs argument stores nonoverlapping subsets
     * of those subdomains whose union covers the local domain. Set-like input
     * containers keep construction simple for callers; this class normalizes
     * them into sorted vectors before storing or materializing PETSc objects.
     */
    virtual void generateASMSubdomains(std::vector<std::set<int>>& subdomain_dofs,
                                       std::vector<std::set<int>>& nonoverlap_subdomain_dofs);

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
     * \brief Optional postprocess hook for shell smoothers.
     *
     * Subclasses may override this method to apply problem-specific postprocess
     * operations after shell preconditioner application.
     */
    virtual void postprocessShellResult(Vec& /*y*/);

    /*!
     * \brief Copy a PETSc matrix into an Eigen sparse matrix.
     *
     * Subclasses can use this helper when experimenting with alternative
     * Eigen-based local solvers built from PETSc-owned operators.
     */
    static Eigen::SparseMatrix<double, Eigen::RowMajor> copyPETScMatToEigenSparse(Mat mat);

    class ConstPetscVecArrayMap
    {
    public:
        ConstPetscVecArrayMap(Vec vec, Eigen::Index size) : d_vec(vec), d_size(size)
        {
            const int ierr = VecGetArrayRead(d_vec, &d_array);
            IBTK_CHKERRQ(ierr);
        }

        ~ConstPetscVecArrayMap()
        {
            const int ierr = VecRestoreArrayRead(d_vec, &d_array);
            IBTK_CHKERRQ(ierr);
        }

        Eigen::Map<const Eigen::VectorXd> getMap() const
        {
            return Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double*>(d_array), d_size);
        }

    private:
        Vec d_vec = nullptr;
        Eigen::Index d_size = 0;
        const PetscScalar* d_array = nullptr;
    };

    class PetscVecArrayMap
    {
    public:
        PetscVecArrayMap(Vec vec, Eigen::Index size) : d_vec(vec), d_size(size)
        {
            const int ierr = VecGetArray(d_vec, &d_array);
            IBTK_CHKERRQ(ierr);
        }

        ~PetscVecArrayMap()
        {
            const int ierr = VecRestoreArray(d_vec, &d_array);
            IBTK_CHKERRQ(ierr);
        }

        Eigen::Map<Eigen::VectorXd> getMap() const
        {
            return Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*>(d_array), d_size);
        }

    private:
        Vec d_vec = nullptr;
        Eigen::Index d_size = 0;
        PetscScalar* d_array = nullptr;
    };

    enum class EigenSubdomainSolverType
    {
        CUSTOM,
        LLT,
        LDLT,
        PARTIAL_PIV_LU,
        FULL_PIV_LU,
        HOUSEHOLDER_QR,
        COL_PIV_HOUSEHOLDER_QR,
        COMPLETE_ORTHOGONAL_DECOMPOSITION,
        FULL_PIV_HOUSEHOLDER_QR,
        JACOBI_SVD,
        BDC_SVD
    };

    template <class Handler>
    void dispatchEigenSolverType(EigenSubdomainSolverType solver_type, Handler&& handler) const;

    EigenSubdomainSolverType parseEigenSubdomainSolverType(const std::string& type) const;

    bool usingCustomEigenSolveSubdomainSolver() const;

    bool usingCustomEigenPseudoinverseSubdomainSolver() const;

    const std::vector<int>& getEigenOverlapDofs(std::size_t subdomain_num) const;

    const std::vector<int>& getEigenUpdateDofs(std::size_t subdomain_num) const;

    const std::vector<int>& getEigenUpdateLocalPositions(std::size_t subdomain_num) const;

    const std::vector<int>& getEigenActiveResidualUpdateRows(std::size_t subdomain_num) const;

    const Eigen::SparseMatrix<double, Eigen::RowMajor>&
    getEigenActiveResidualUpdateMat(std::size_t subdomain_num) const;

    template <class InitializeSubdomainSolver>
    void initializeEigenShellDataWithLocalOperatorHook(InitializeSubdomainSolver initialize_subdomain_solver);

    void initializeBuiltinEigenSubdomainSolver(EigenSubdomainSolverType solver_type,
                                               const Eigen::MatrixXd& local_operator,
                                               std::size_t subdomain_num);

    Eigen::VectorXd solveBuiltinEigenSubdomainSystem(EigenSubdomainSolverType solver_type,
                                                     const Eigen::VectorXd& rhs,
                                                     std::size_t subdomain_num) const;

    virtual void initializeEigenSubdomainSolver(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num);

    virtual Eigen::VectorXd solveEigenSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const;

    virtual void initializeCustomEigenShellData();

    virtual void applyAdditiveCustomEigen(Vec x, Vec y);

    virtual void applyMultiplicativeCustomEigen(Vec x, Vec y);

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
     * \name PETSc objects.
     */
    //\{
    std::string d_ksp_type = KSPGMRES, d_pc_type = PCILU, d_shell_pc_type;
    std::string d_options_prefix;
    KSP d_petsc_ksp = nullptr;
    Mat d_petsc_mat = nullptr, d_petsc_pc = nullptr;
    MatNullSpace d_petsc_nullsp = nullptr;
    Vec d_petsc_x = nullptr, d_petsc_b = nullptr;
    //\}

    /*!
     * \name Field split preconditioner.
     */
    //\{
    std::vector<std::string> d_field_name;
    std::vector<IS> d_field_is;
    //\}

    /*!
     * \brief Overlapping subdomains and nonoverlapping partition subsets used
     * by ASM-like preconditioners.
     *
     * The overlapping subdomains define the subdomain solves. The partition
     * subsets define nonoverlapping subsets of those subdomains whose union
     * covers the local domain.
     */
    std::vector<std::vector<int>> d_subdomain_dofs, d_nonoverlap_subdomain_dofs;

    /*!
     * \brief PETSc index-set caches derived from \ref d_subdomain_dofs and
     * \ref d_nonoverlap_subdomain_dofs when PETSc needs explicit IS objects.
     */
    std::vector<IS> d_petsc_subdomain_is, d_petsc_nonoverlap_subdomain_is;

private:
    enum class PreconditionerType
    {
        OTHER,
        ASM,
        FIELDSPLIT,
        SHELL
    };

    enum class ShellSmootherComposition
    {
        ADDITIVE,
        MULTIPLICATIVE
    };

    enum class ShellSmootherBackend
    {
        PETSC,
        EIGEN,
        EIGEN_PSEUDOINVERSE
    };

    enum class ShellSmootherPartition
    {
        BASIC,
        RESTRICT
    };

    struct PetscShellSmootherData
    {
        Vec shell_r = nullptr;
        InsertMode prolongation_insert_mode = INSERT_VALUES;
        IS owned_residual_update_rows_is = nullptr;
        std::vector<IS> local_overlap_is, local_nonoverlap_is;
        std::vector<VecScatter> restriction, prolongation;
        std::vector<KSP> sub_ksp;
        Mat *sub_mat = nullptr, *active_residual_update_mat = nullptr;
        std::vector<Vec> sub_x, sub_y, active_residual_update_x, active_residual_update_y;
        std::vector<std::vector<int>> active_update_local_positions;
    };

    struct EigenSubdomainCommonCache
    {
        std::vector<int> overlap_dofs;
        std::vector<int> nonoverlap_dofs;
        std::vector<int> nonoverlap_local_positions;
        std::vector<int> update_dofs;
        std::vector<int> update_local_positions;
        std::vector<int> active_residual_update_rows;
        Eigen::SparseMatrix<double, Eigen::RowMajor> active_residual_update_mat;
        Eigen::MatrixXd local_pseudoinverse;
        Eigen::VectorXd rhs_workspace;
        Eigen::VectorXd delta_workspace;
        Eigen::VectorXd residual_input_workspace;
        Eigen::VectorXd residual_delta_workspace;
    };

    struct EigenSolveStorageBase
    {
        virtual ~EigenSolveStorageBase() = default;
    };

    template <class SolverType>
    struct EigenTypedSolveStorage : public EigenSolveStorageBase
    {
        std::vector<SolverType> solvers;
    };

    struct EigenShellSmootherData
    {
        Eigen::Index n_dofs = 0;
        std::vector<EigenSubdomainCommonCache> common_subdomains;
        std::unique_ptr<EigenSolveStorageBase> solve_storage;
    };

    template <class SolverType>
    struct EigenSolverTag
    {
        using type = SolverType;
    };

    template <class SolverType>
    EigenTypedSolveStorage<SolverType>& getEigenSolveStorage();

    template <class SolverType>
    const EigenTypedSolveStorage<SolverType>& getEigenSolveStorage() const;

    template <class SolverType>
    void initializeBuiltinEigenSolveStorage();

    template <class SolverType>
    void applyAdditiveEigenSolveImpl(Vec x, Vec y);

    template <class SolverType>
    void applyMultiplicativeEigenSolveImpl(Vec x, Vec y);

    struct ShellSmootherData
    {
        int n_local_subdomains = 0;
        std::unique_ptr<PetscShellSmootherData> petsc_data;
        std::unique_ptr<EigenShellSmootherData> eigen_data;
    };

    /*!
     * \brief Parse the configured shell smoother string into backend,
     * composition, and partition settings.
     */
    void configureShellSmootherType();

    std::string normalizeShellSmootherType(const std::string& type) const;
    PreconditionerType parsePreconditionerType(const std::string& type) const;
    ShellSmootherBackend parseShellSmootherBackend(const std::string& type) const;
    ShellSmootherComposition parseShellSmootherComposition(const std::string& type) const;
    ShellSmootherPartition parseShellSmootherPartition(const std::string& type,
                                                       ShellSmootherComposition composition) const;
    EigenSubdomainSolverType getEigenSolveSolverType() const;
    EigenSubdomainSolverType getEigenPseudoinverseSolverType() const;

    /*!
     * \brief Cache set-like subdomain descriptions in the stored vector form.
     */
    void cacheASMSubdomains(const std::vector<std::set<int>>& subdomain_dofs,
                            const std::vector<std::set<int>>& nonoverlap_subdomain_dofs);

    /*!
     * \brief Generate and cache the stored ASM-like subdomain description.
     */
    void cacheGeneratedASMSubdomains();

    /*!
     * \brief Configure a PETSc ASM preconditioner from the cached subdomain
     * description.
     */
    void configureASMPreconditioner(PC ksp_pc);

    /*!
     * \brief Configure a PETSc fieldsplit preconditioner.
     */
    void configureFieldSplitPreconditioner(PC ksp_pc);

    /*!
     * \brief Configure a shell preconditioner from the cached subdomain
     * description.
     */
    void configureShellPreconditioner(PC ksp_pc);

    /*!
     * \brief Initialize PETSc data shared by additive and multiplicative shell
     * smoothers.
     */
    void initializePetscShellData();

    /*!
     * \brief Deallocate PETSc data shared by additive and multiplicative shell
     * smoothers.
     */
    void deallocatePetscShellData();

    /*!
     * \brief Initialize Eigen data shared by additive and multiplicative shell
     * smoothers.
     */
    void initializeEigenShellData();

    /*!
     * \brief Deallocate Eigen data shared by additive and multiplicative shell
     * smoothers.
     */
    void deallocateEigenShellData();

    /*!
     * \brief Deallocate shell smoother data.
     */
    void deallocateShellData();

    /*!
     * \brief Configure the shell apply callback.
     */
    void configureShellApply(PC ksp_pc);

    /*!
     * \brief Begin accumulating a PETSc subdomain correction using the
     * configured partition mode.
     */
    void beginAccumulateCorrectionPetsc(int subdomain_num, Vec sub_y, Vec y);

    /*!
     * \brief Finish accumulating a PETSc subdomain correction using the
     * configured partition mode.
     */
    void endAccumulateCorrectionPetsc(int subdomain_num, Vec sub_y, Vec y);

    /*!
     * \brief Accumulate a PETSc subdomain correction using the configured
     * partition mode.
     */
    void accumulateCorrectionPetsc(int subdomain_num, Vec sub_y, Vec y);

    /*!
     * \brief Apply the PETSc-based additive shell preconditioner.
     */
    void applyAdditivePetsc(Vec x, Vec y);

    /*!
     * \brief Apply the PETSc-based multiplicative shell preconditioner.
     */
    void applyMultiplicativePetsc(Vec x, Vec y);

    /*!
     * \brief Update the cached PETSc shell residual after one multiplicative
     * subdomain correction.
     */
    void updateResidualPetsc(int subdomain_num, Vec sub_y, Vec residual);

    /*!
     * \brief Apply the Eigen-based additive shell preconditioner.
     */
    void applyAdditiveEigen(Vec x, Vec y);

    /*!
     * \brief Apply the Eigen-based multiplicative shell preconditioner.
     */
    void applyMultiplicativeEigen(Vec x, Vec y);

    Eigen::MatrixXd buildEigenSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator) const;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScLevelSolver(const PETScLevelSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScLevelSolver& operator=(const PETScLevelSolver& that) = delete;

    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_AdditivePetsc(PC pc, Vec x, Vec y);

    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_MultiplicativePetsc(PC pc, Vec x, Vec y);

    /*!
     * \brief Apply the Eigen-based additive shell preconditioner.
     */
    static PetscErrorCode PCApply_AdditiveEigen(PC pc, Vec x, Vec y);

    /*!
     * \brief Apply the Eigen-based multiplicative shell preconditioner.
     */
    static PetscErrorCode PCApply_MultiplicativeEigen(PC pc, Vec x, Vec y);

    PreconditionerType d_preconditioner_type = PreconditionerType::OTHER;
    ShellSmootherComposition d_shell_smoother_composition = ShellSmootherComposition::MULTIPLICATIVE;
    ShellSmootherBackend d_shell_smoother_backend = ShellSmootherBackend::PETSC;
    ShellSmootherPartition d_shell_smoother_partition = ShellSmootherPartition::BASIC;
    EigenSubdomainSolverType d_eigen_subdomain_solver_type = EigenSubdomainSolverType::PARTIAL_PIV_LU;
    double d_eigen_subdomain_solver_threshold = -1.0;
    EigenSubdomainSolverType d_eigen_subdomain_pseudoinverse_type = EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
    double d_eigen_subdomain_pseudoinverse_threshold = -1.0;
    ShellSmootherData d_shell_data;
};

template <class InitializeSubdomainSolver>
inline void
PETScLevelSolver::initializeEigenShellDataWithLocalOperatorHook(InitializeSubdomainSolver initialize_subdomain_solver)
{
    auto& shell = d_shell_data;
    if (!shell.eigen_data) shell.eigen_data = std::make_unique<EigenShellSmootherData>();
    auto& eigen = *shell.eigen_data;
    const bool use_multiplicative = d_shell_smoother_composition == ShellSmootherComposition::MULTIPLICATIVE;
    const Eigen::SparseMatrix<double, Eigen::RowMajor> eigen_level_mat = copyPETScMatToEigenSparse(d_petsc_mat);
    Eigen::SparseMatrix<double> eigen_level_mat_transpose;
    if (use_multiplicative) eigen_level_mat_transpose = eigen_level_mat.transpose();
    eigen.n_dofs = eigen_level_mat.rows();
    eigen.common_subdomains.clear();
    eigen.common_subdomains.resize(static_cast<std::size_t>(shell.n_local_subdomains));
    for (int subdomain_num = 0; subdomain_num < shell.n_local_subdomains; ++subdomain_num)
    {
        auto& cache = eigen.common_subdomains[static_cast<std::size_t>(subdomain_num)];
        cache.overlap_dofs = d_subdomain_dofs[static_cast<std::size_t>(subdomain_num)];
        const Eigen::Index overlap_size = static_cast<Eigen::Index>(cache.overlap_dofs.size());
        std::unordered_map<int, int> overlap_col_map;
        overlap_col_map.reserve(static_cast<std::size_t>(overlap_size));
        for (Eigen::Index local_col = 0; local_col < overlap_size; ++local_col)
        {
            overlap_col_map[cache.overlap_dofs[static_cast<std::size_t>(local_col)]] = static_cast<int>(local_col);
        }

        cache.nonoverlap_dofs = d_nonoverlap_subdomain_dofs[static_cast<std::size_t>(subdomain_num)];
        cache.nonoverlap_local_positions.resize(cache.nonoverlap_dofs.size());
        std::unordered_map<int, int> nonoverlap_col_map;
        nonoverlap_col_map.reserve(cache.nonoverlap_dofs.size());
        for (std::size_t local_col = 0; local_col < cache.nonoverlap_dofs.size(); ++local_col)
        {
            const int dof = cache.nonoverlap_dofs[local_col];
            const auto overlap_pos_it = overlap_col_map.find(dof);
#if !defined(NDEBUG)
            TBOX_ASSERT(overlap_pos_it != overlap_col_map.end());
#endif
            cache.nonoverlap_local_positions[local_col] = overlap_pos_it->second;
            nonoverlap_col_map[dof] = static_cast<int>(local_col);
        }

        Eigen::MatrixXd local_operator = Eigen::MatrixXd::Zero(overlap_size, overlap_size);
        for (Eigen::Index local_row = 0; local_row < overlap_size; ++local_row)
        {
            const int global_row = cache.overlap_dofs[static_cast<std::size_t>(local_row)];
            for (auto it = Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(eigen_level_mat, global_row); it;
                 ++it)
            {
                const auto col_it = overlap_col_map.find(static_cast<int>(it.col()));
                if (col_it != overlap_col_map.end()) local_operator(local_row, col_it->second) = it.value();
            }
        }

        initialize_subdomain_solver(local_operator, static_cast<std::size_t>(subdomain_num));

        switch (d_shell_smoother_partition)
        {
        case ShellSmootherPartition::BASIC:
            cache.update_dofs = cache.overlap_dofs;
            cache.update_local_positions.resize(cache.overlap_dofs.size());
            for (std::size_t local_pos = 0; local_pos < cache.update_local_positions.size(); ++local_pos)
            {
                cache.update_local_positions[local_pos] = static_cast<int>(local_pos);
            }
            break;
        case ShellSmootherPartition::RESTRICT:
            cache.update_dofs = cache.nonoverlap_dofs;
            cache.update_local_positions = cache.nonoverlap_local_positions;
            break;
        }
        cache.rhs_workspace.resize(overlap_size);
        cache.delta_workspace.resize(overlap_size);

        if (use_multiplicative)
        {
            const std::vector<int>& active_update_dofs = cache.update_dofs;
            const int active_num_cols = static_cast<int>(active_update_dofs.size());
            std::vector<Eigen::Triplet<double>> triplets;
            std::unordered_map<int, int> row_map;
            row_map.reserve(static_cast<std::size_t>(active_num_cols));
            for (int local_col = 0; local_col < active_num_cols; ++local_col)
            {
                const int global_col = active_update_dofs[static_cast<std::size_t>(local_col)];
                for (auto it = Eigen::SparseMatrix<double>::InnerIterator(eigen_level_mat_transpose, global_col); it;
                     ++it)
                {
                    const int row = static_cast<int>(it.row());
                    const auto row_it = row_map.find(row);
                    int local_row = -1;
                    if (row_it == row_map.end())
                    {
                        local_row = static_cast<int>(cache.active_residual_update_rows.size());
                        cache.active_residual_update_rows.push_back(row);
                        row_map.emplace(row, local_row);
                    }
                    else
                    {
                        local_row = row_it->second;
                    }
                    triplets.emplace_back(local_row, local_col, it.value());
                }
            }

            cache.active_residual_update_mat.resize(static_cast<int>(cache.active_residual_update_rows.size()),
                                                    active_num_cols);
            cache.active_residual_update_mat.setFromTriplets(triplets.begin(), triplets.end());
            cache.residual_input_workspace.resize(active_num_cols);
            cache.residual_delta_workspace.resize(static_cast<Eigen::Index>(cache.active_residual_update_rows.size()));
        }
    }
    return;
}

template <class Handler>
inline void
PETScLevelSolver::dispatchEigenSolverType(const EigenSubdomainSolverType solver_type, Handler&& handler) const
{
    switch (solver_type)
    {
    case EigenSubdomainSolverType::CUSTOM:
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::dispatchEigenSolverType()\n"
                                 << "CUSTOM is not a built-in Eigen subdomain solver type." << std::endl);
        break;
    case EigenSubdomainSolverType::LLT:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::LLT<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::LDLT:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::LDLT<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::PARTIAL_PIV_LU:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::PartialPivLU<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::FULL_PIV_LU:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::FullPivLU<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::HOUSEHOLDER_QR:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::HouseholderQR<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::COMPLETE_ORTHOGONAL_DECOMPOSITION:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::JACOBI_SVD:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::JacobiSVD<Eigen::MatrixXd>>{});
        return;
    case EigenSubdomainSolverType::BDC_SVD:
        std::forward<Handler>(handler)(EigenSolverTag<Eigen::BDCSVD<Eigen::MatrixXd>>{});
        return;
    }
    return;
}
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PETScLevelSolver
