// Filename: KrylovMobilityInverse.h
// Created on 28 Oct 2013 by Amneet Bhalla on Taylor@mech.northwestern.edu


#ifndef included_KrylovMobilityInverse
#define included_KrylovMobilityInverse

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "petscksp.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Database.h"
#include "SAMRAIVectorReal.h"
#include "RobinBcCoefStrategy.h"
#include "PoissonSpecifications.h"
#include "vector"

namespace IBAMR
{
class InterpOperator;
class SpreadOperator;
class StaggeredStokesSolver;
class INSStaggeredHierarchyIntegrator;
class StaggeredStokesPhysicalBoundaryHelper;
class DirectMobilityInverse;
class cIBMethod;
}
namespace IBTK
{
class PoissonSolver;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
  
/*! 
 * We are trying to solve the problem
 * 
 * \f$ Mx = [J inv(L) S]x = b \f$; for \f$ x \f$. 
 * 
 * Here, \f$ M \f$ is the mobility matrix, \f$ J \f$ is the interpolation 
 * operator, \f$ L \f$ is the Stokes solver, and \f$ S \f$ is the spreading 
 * operator.
 * 
 */ 

class KrylovMobilityInverse
    : public SAMRAI::tbox::DescribedClass
{
public:
  
    /*!
     * \brief Constructor for [J L^-1 S]^-1 solver that employs the
     * PETSc KSP solver framework.
     */
    KrylovMobilityInverse(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
        SAMRAI::tbox::Pointer<IBAMR::cIBMethod> cib_method,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::string& default_options_prefix,
        MPI_Comm petsc_comm=PETSC_COMM_WORLD);

    /*!
     * \brief Destructor.
     */
    ~KrylovMobilityInverse();
      
    /*!
     * \brief Set the KSP type.
     */
    void
    setKSPType(
        const std::string& ksp_type);

    /*!
     * \brief Set the options prefix used by this PETSc solver object.
     */
    void
    setOptionsPrefix(
        const std::string& options_prefix);

    /*!
     * \name Functions to access the underlying PETSc objects.
     */
    //\{

    /*!
     * \brief Get the PETSc KSP object.
     */
    const KSP&
    getPETScKSP() const;

   /*!
    * \brief Set the PoissonSpecifications object used to specify the
    * coefficients for the momentum equation in the incompressible Stokes
    * operator.
    */
    void
    setVelocityPoissonSpecifications(
        const SAMRAI::solv::PoissonSpecifications& u_problem_coefs);

    /*!
     * \brief Set the StokesSpecifications object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    void
    setPhysicalBoundaryHelper(
        SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesPhysicalBoundaryHelper> bc_helper);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param u_bc_coefs  IBTK::Vector of pointers to objects that can set the Robin boundary condition coefficients for the velocity
     * \param p_bc_coef   Pointer to object that can set the Robin boundary condition coefficients for the pressure
     */
    void
    setPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_bc_coef);


    /*!
     * \brief Set the linear operators used when solving \f$ Mx=b \f$.
     */
    void
    setLinearOperators(
	SAMRAI::tbox::Pointer<IBAMR::InterpOperator> J,
	SAMRAI::tbox::Pointer<IBAMR::SpreadOperator> S);

    /*!
     * \brief Set DirectMobility as a preconditioner for this solver.
     */
    void
    setPreconditioner(
        SAMRAI::tbox::Pointer<IBAMR::DirectMobilityInverse> pc);
    
    /*!
     * \brief Solve the linear system of equations \f$ Mx=b \f$ for \f$ x \f$.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool
    solveSystem(
        Vec x,
        Vec b);
    
    /*!
     * \brief Compute hierarchy dependent data required for solving \f$ Mx = b \f$.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     */
    void
    initializeSolverState(
        Vec x,
        Vec b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     */
    void
    deallocateSolverState();

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void
    setSolutionTime(
        double solution_time);

    /*!
     * \brief Set the current time interval.
     */
    void
    setTimeInterval(
        double current_time,
        double new_time);

    /*!
     * \brief Set scale factor for interp operator.
     */
    void
    setInterpScaleFactor(
        const double beta);
 
    /*!
     * \brief Set scale factor for spread operator.
     */
    void
    setSpreadScaleFactor(
        const double gamma);  

    /*!
     * \brief Set scale factor for regularizing mobility matrix.
     */
    void
    setRegularizeMobilityFactor(
        const double delta); 

private:
    
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    KrylovMobilityInverse();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    KrylovMobilityInverse(
        const KrylovMobilityInverse& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    KrylovMobilityInverse&
    operator=(
        const KrylovMobilityInverse& that);

    /*!
     * \brief Initialize the Stokes solver needed in the mobility matrix.
     */
    void
    initializeStokesSolver(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& sol_vec,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs_vec);
         
    /*!
     * \brief Report the KSPConvergedReason.
     */
    void
    reportKSPConvergedReason(
        const KSPConvergedReason& reason,
        std::ostream& os) const;

    /*!
     * \brief Routine to setup KSP object.
     */
    void
    initializeKSP();

    /*!
     * \brief Reset the values of the convergence tolerances for the PETSc KSP
     * object.
     */
    void
    resetKSPOptions();

    /*!
     * \brief Reset the KSP operators to correspond to the supplied
     * LinearOperator.
     */
    void
    resetKSPOperators();

    /*!
     * \brief Reset the KSP PC to correspond to the supplied preconditioner.
     */
    void
    resetKSPPC();

    /*!
     *\brief Routine to destroy KSP object.
     */
    void
    destroyKSP();

    /*!
     * \name Static functions for use by PETSc KSP and MatShell objects.
     */
    //\{

    /*!
     * \brief Compute the matrix vector product \f$y=Ax\f$.
     */
    static PetscErrorCode
    MatVecMult_KMInv(
        Mat A,
        Vec x,
        Vec y);

    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode
    PCApply_KMInv(
        PC pc,
        Vec x,
        Vec y);

    /*!
     * \brief Set KSP monitoring routine for the KSP.
     */
    static PetscErrorCode
    monitorKSP(
        KSP ksp,
        int it,
        PetscReal rnorm,
        void *mctx);

    //\}

    // Solver stuff
    std::string d_object_name, d_ksp_type, d_pc_type;
    bool d_is_initialized;
    bool d_reinitializing_solver;

    Vec d_petsc_x, d_petsc_b;
    
    std::string d_options_prefix;

    MPI_Comm     d_petsc_comm;
    KSP          d_petsc_ksp;
    Mat          d_petsc_mat;
    
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,PetscScalar> > > d_samrai_temp; 
    
    SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> d_ins_integrator;
    SAMRAI::tbox::Pointer<IBAMR::cIBMethod> d_cib_method;
    SAMRAI::tbox::Pointer<IBAMR::SpreadOperator> d_S;
    SAMRAI::tbox::Pointer<IBAMR::InterpOperator> d_J;
    SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesSolver> d_LInv;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_velocity_solver, d_pressure_solver;    
    
    int d_max_iterations, d_current_iterations;
    double d_abs_residual_tol, d_rel_residual_tol;
    double d_current_residual_norm;
    bool d_initial_guess_nonzero; 
    bool d_enable_logging;

    // Preconditioner for KrylovMobility
    SAMRAI::tbox::Pointer<IBAMR::DirectMobilityInverse> d_DMInv;

    //Nullspace vectors for LInv
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > > d_nul_vecs;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > > d_U_nul_vecs;

    /*!
     * This boolean value determines whether the pressure is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     */
    bool d_normalize_pressure;

    /*!
     * This boolean value determines whether the velocity is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     *
     * This parameter only affects the case in which rho=0 (i.e. the steady
     * Stokes equations).
     */
    bool d_normalize_velocity;

    // The current, new and time-interval of integration.
    double d_current_time, d_new_time, d_dt;

    // Scaling parameters of the problem.
    double d_scale_interp, d_scale_spread, d_reg_mob_factor;
};
}// namespace IBAMR


//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_KrylovMobilityInverse
