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
#include <ibtk/private/PETScLevelSolverEigenShellBackendCommon.h>

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

#include <map>
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
class PETScLevelSolverEigenShellBackendBase;
class PETScLevelSolverPetscShellBackend;
class PETScLevelSolverEigenShellBackend;
class PETScLevelSolverEigenPseudoinverseShellBackend;
class PETScLevelSolverEigenReferenceShellBackend;
class PETScLevelSolver;
class PETScLevelSolverBackendContext
{
public:
    virtual ~PETScLevelSolverBackendContext() = default;

    virtual const std::string& getObjectNameForBackend() const = 0;
    virtual const std::string& getOptionsPrefixForBackend() const = 0;

    virtual bool isShellMultiplicativeForBackend() const = 0;
    virtual bool useRestrictPartitionForBackend() const = 0;

    virtual const std::vector<std::vector<int>>& getSubdomainDOFsForBackend() const = 0;
    virtual const std::vector<std::vector<int>>& getNonoverlapSubdomainDOFsForBackend() const = 0;

    virtual Mat getPETScMatForBackend() const = 0;
    virtual Vec getPETScXForBackend() const = 0;
    virtual Vec getPETScBForBackend() const = 0;

    virtual void postprocessShellResultForBackend(Vec y) = 0;
};

class PETScLevelSolverShellBackend
{
public:
    virtual ~PETScLevelSolverShellBackend() = default;

    virtual const std::string& getTypeKey() const = 0;

    virtual void configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) = 0;

    virtual const char* getPCNameSuffixAdditive() const = 0;

    virtual const char* getPCNameSuffixMultiplicative() const = 0;

    virtual void initialize() = 0;

    virtual void deallocate() = 0;

    virtual void applyAdditive(Vec x, Vec y) = 0;

    virtual void applyMultiplicative(Vec x, Vec y) = 0;
};

class PETScLevelSolverShellBackendManager
{
public:
    using ShellBackendMaker =
        std::unique_ptr<PETScLevelSolverShellBackend> (*)(PETScLevelSolver& solver,
                                                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    static PETScLevelSolverShellBackendManager* getManager();

    static void freeManager();

    std::unique_ptr<PETScLevelSolverShellBackend>
    allocateShellBackend(const std::string& type_key,
                         PETScLevelSolver& solver,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) const;

    void registerShellBackendFactoryFunction(const std::string& type_key, ShellBackendMaker backend_maker);

    std::vector<std::string> getRegisteredShellBackendTypes() const;

protected:
    PETScLevelSolverShellBackendManager();
    ~PETScLevelSolverShellBackendManager() = default;

private:
    PETScLevelSolverShellBackendManager(const PETScLevelSolverShellBackendManager& from) = delete;
    PETScLevelSolverShellBackendManager& operator=(const PETScLevelSolverShellBackendManager& that) = delete;

    static PETScLevelSolverShellBackendManager* s_shell_backend_manager_instance;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    std::map<std::string, ShellBackendMaker> d_shell_backend_maker_map;
};

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
class PETScLevelSolver : public LinearSolver, public PETScLevelSolverBackendContext
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

    const std::string& getObjectNameForBackend() const override;

    const std::string& getOptionsPrefixForBackend() const override;

    bool isShellMultiplicativeForBackend() const override;

    bool useRestrictPartitionForBackend() const override;

    const std::vector<std::vector<int>>& getSubdomainDOFsForBackend() const override;

    const std::vector<std::vector<int>>& getNonoverlapSubdomainDOFsForBackend() const override;

    Mat getPETScMatForBackend() const override;

    Vec getPETScXForBackend() const override;

    Vec getPETScBForBackend() const override;

    void postprocessShellResultForBackend(Vec y) override;

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
     * \brief Setup the solver nullspace (if any).
     */
    virtual void setupNullSpace();

    PETScLevelSolverShellBackend* getShellBackend(const std::string& type_key);
    const PETScLevelSolverShellBackend* getShellBackend(const std::string& type_key) const;

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

    enum class ShellSmootherPartition
    {
        BASIC,
        RESTRICT
    };

    /*!
     * \brief Parse the configured shell smoother string into backend,
     * composition, and partition settings.
     */
    void configureShellSmootherType();
    void loadShellBackends(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    std::string normalizeShellSmootherType(const std::string& type) const;
    std::string extractShellSmootherTypeKey(const std::string& type) const;
    PreconditionerType parsePreconditionerType(const std::string& type) const;
    std::string parseShellSmootherBackendKey(const std::string& type_key) const;
    ShellSmootherComposition parseShellSmootherComposition(const std::string& type) const;
    ShellSmootherPartition parseShellSmootherPartition(const std::string& type,
                                                       ShellSmootherComposition composition) const;
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
     * \brief Deallocate shell smoother data.
     */
    void deallocateShellData();

    /*!
     * \brief Configure the shell apply callback.
     */
    void configureShellApply(PC ksp_pc);

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
    static PetscErrorCode PCApply_AdditiveShell(PC pc, Vec x, Vec y);

    static PetscErrorCode PCApply_MultiplicativeShell(PC pc, Vec x, Vec y);

    PreconditionerType d_preconditioner_type = PreconditionerType::OTHER;
    ShellSmootherComposition d_shell_smoother_composition = ShellSmootherComposition::MULTIPLICATIVE;
    std::string d_shell_smoother_backend_key = "petsc";
    ShellSmootherPartition d_shell_smoother_partition = ShellSmootherPartition::BASIC;
    std::unordered_map<std::string, std::unique_ptr<PETScLevelSolverShellBackend>> d_shell_backends;
    PETScLevelSolverShellBackend* d_active_shell_backend = nullptr;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PETScLevelSolver
