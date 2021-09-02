// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_KrylovMobilitySolver
#define included_IBAMR_KrylovMobilitySolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/HierarchyGhostCellInterpolation.h"

#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include "petscksp.h"

#include <vector>

namespace IBAMR
{
class StaggeredStokesSolver;
class INSStaggeredHierarchyIntegrator;
class StaggeredStokesPhysicalBoundaryHelper;
class CIBStrategy;
} // namespace IBAMR
namespace IBTK
{
class PoissonSolver;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * We are trying to solve the problem
 *
 * \f$ Mx = [J L^{-1} S]x = b \f$; for \f$ x \f$.
 *
 * Here, \f$ M \f$ is the mobility matrix, \f$ J \f$ is the interpolation
 * operator, \f$ L \f$ is the Stokes operator, and \f$ S \f$ is the spreading
 * operator.
 *
 */
class KrylovMobilitySolver : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor for mobility solver that employs the
     * PETSc KSP solver framework.
     */
    KrylovMobilitySolver(std::string object_name,
                         SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                         SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         std::string default_options_prefix,
                         MPI_Comm petsc_comm = PETSC_COMM_WORLD);

    /*!
     * \brief Destructor.
     */
    virtual ~KrylovMobilitySolver();

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
     * \brief Return the Stokes solver used in the preconditioner of the
     * solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesSolver> getStokesSolver() const;

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& u_problem_coefs);

    /*!
     * \brief Set the StokesSpecifications object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesPhysicalBoundaryHelper> bc_helper);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a u_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a p_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param u_bc_coefs Vector of pointers to objects that can set the Robin
     * boundary condition coefficients for the velocity.
     *
     * \param p_bc_coef Pointer to object that can set the Robin boundary
     * condition coefficients for the pressure.
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_bc_coef);

    /*!
     * \brief Solve the linear system of equations \f$ Mx=b \f$ for \f$ x \f$.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveSystem(Vec x, Vec b);

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$ Mx = b \f$.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     */
    void initializeSolverState(Vec x, Vec b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     */
    void deallocateSolverState();

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void setSolutionTime(double solution_time);

    /*!
     * \brief Set the current time interval.
     */
    void setTimeInterval(double current_time, double new_time);

    /*!
     * \brief Set scale factor for interp operator.
     */
    void setInterpScale(const double scale_interp);

    /*!
     * \brief Set scale factor for spread operator.
     */
    void setSpreadScale(const double scale_spread);

    /*!
     * \brief Set scale factor for regularizing mobility matrix.
     */
    void setRegularizeMobilityScale(const double scale_reg_mob);

    /*!
     * \brief Set if the mean of the Lagrangian force is to be subtracted
     * from the Eulerian force variable.
     *
     * \note This operation is needed for certain situations like Stokes flow
     * with periodic BCs.
     */
    void setNormalizeSpreadForce(const bool normalize_force);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    KrylovMobilitySolver(const KrylovMobilitySolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    KrylovMobilitySolver& operator=(const KrylovMobilitySolver& that) = delete;

    /*!
     * \brief Get solver settings from the input file.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Initialize the Stokes solver needed in the mobility matrix.
     */
    void initializeStokesSolver(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& sol_vec,
                                const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs_vec);

    /*!
     * \brief Routine to setup KSP object.
     */
    void initializeKSP();

    /*!
     * \brief Reset the values of the convergence tolerances for the PETSc KSP
     * object.
     */
    void resetKSPOptions();

    /*!
     * \brief Reset the KSP operators to correspond to the supplied
     * LinearOperator.
     */
    void resetKSPOperators();

    /*!
     * \brief Reset the KSP PC to correspond to the supplied preconditioner.
     */
    void resetKSPPC();

    /*!
     *\brief Routine to destroy KSP object.
     */
    void destroyKSP();

    /*!
     * \name Static functions for use by PETSc KSP and MatShell objects.
     */
    //\{

    /*!
     * \brief Compute the matrix vector product \f$y=Ax\f$.
     */
    static PetscErrorCode MatVecMult_KMInv(Mat A, Vec x, Vec y);

    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_KMInv(PC pc, Vec x, Vec y);

    /*!
     * \brief Set KSP monitoring routine for the KSP.
     */
    static PetscErrorCode monitorKSP(KSP ksp, int it, PetscReal rnorm, void* mctx);

    //\}

    // Solver stuff
    std::string d_object_name, d_ksp_type = KSPGMRES, d_pc_type = "none";
    bool d_is_initialized = false;
    bool d_reinitializing_solver = false;
    Vec d_petsc_x = nullptr, d_petsc_b = nullptr;
    std::string d_options_prefix;
    MPI_Comm d_petsc_comm;
    KSP d_petsc_ksp = nullptr;
    Mat d_petsc_mat = nullptr;

    // Linear operator.
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, PetscScalar> > > d_samrai_temp;
    SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> d_ins_integrator;
    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;
    SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesSolver> d_LInv;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_velocity_solver, d_pressure_solver;

    // KSP options and settings.
    int d_max_iterations = 10000, d_current_iterations;
    double d_abs_residual_tol = 1.0e-50, d_rel_residual_tol = 1.0e-5;
    double d_current_residual_norm;
    bool d_initial_guess_nonzero = false;
    bool d_enable_logging = false;

    // Velocity BCs and cached communication operators for interpolation operation.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_u_bc_coefs;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_fill_pattern;
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill;

    // Nullspace vectors for LInv
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_nul_vecs;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_U_nul_vecs;

    /*!
     * This boolean value determines whether the pressure is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     */
    bool d_normalize_pressure = false;

    /*!
     * This boolean value determines whether the velocity is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     *
     * This parameter only affects the case in which rho=0 (i.e. the steady
     * Stokes equations).
     */
    bool d_normalize_velocity = false;

    // The current and new time.
    double d_current_time = std::numeric_limits<double>::signaling_NaN(),
           d_new_time = std::numeric_limits<double>::signaling_NaN();

    // Scaling parameters and force normalization of the problem.
    double d_scale_interp = 1.0, d_scale_spread = 1.0, d_reg_mob_factor = 0.0, d_normalize_spread_force;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_KrylovMobilitySolver
