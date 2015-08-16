// Filename: FreeBodyMobilitySolver.h

#ifndef included_FreeBodyMobilitySolver
#define included_FreeBodyMobilitySolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "petsc-private/petscimpl.h"

#include "petscksp.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Database.h"
#include <vector>

namespace IBAMR
{
class CIBStrategy;
class CIBMobilitySolver;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
  
/*! 
 * We are trying to solve the problem
 * 
 * \f$ Nx = [T M^-1 T*]x = b \f$; for \f$ x \f$. 
 * 
 * Here, \f$ N \f$ is the body mobility matrix, \f$ M \f$ is the mobility matrix, 
 * \f$ J \f$ is the interpolation operator, \f$ L \f$ is the Stokes solver, 
 * and \f$ S \f$ is the spreading operator.
 * 
 */ 

class FreeBodyMobilitySolver
    : public SAMRAI::tbox::DescribedClass
{
public:
  
    /*!
     * \brief Constructor for [T M^-1 T*]^-1 solver that employs the
     * PETSc KSP solver framework.
     */
    FreeBodyMobilitySolver(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::string& default_options_prefix,
        SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy,
        MPI_Comm petsc_comm=PETSC_COMM_WORLD);

    /*!
     * \brief Destructor.
     */
    ~FreeBodyMobilitySolver();
      

    /*!
     * \brief Set  as a mobility solver.
     */
    void
    setMobilitySolver(
        SAMRAI::tbox::Pointer<IBAMR::CIBMobilitySolver> minv);
    
    /*!
     * \brief Solve the linear system of equations \f$ Nx = b \f$ for \f$ x \f$.
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
     * \brief Compute hierarchy dependent data required for solving \f$ Nx = b \f$.
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

private:
    
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FreeBodyMobilitySolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FreeBodyMobilitySolver(
        const FreeBodyMobilitySolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FreeBodyMobilitySolver&
    operator=(
        const FreeBodyMobilitySolver& that);
         
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
    MatVecMult_KBMInv(
        Mat A,
        Vec x,
        Vec y);

    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode
    PCApply_KBMInv(
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

    Vec d_petsc_b, d_petsc_temp_v, d_petsc_temp_f;
    
    std::string d_options_prefix;

    MPI_Comm     d_petsc_comm;
    KSP          d_petsc_ksp;
    Mat          d_petsc_mat;

    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;
    SAMRAI::tbox::Pointer<IBAMR::CIBMobilitySolver> d_MInv;
    
    int d_max_iterations, d_current_iterations;
    double d_abs_residual_tol, d_rel_residual_tol;
    double d_current_residual_norm;
    bool d_initial_guess_nonzero; 
    bool d_enable_logging;

    // Number of rigid bodies
    unsigned d_num_rigid_parts;

    // The current, new and time-interval of integration.
    double d_current_time, d_new_time, d_dt;
    
    //viscosity, effective body hydro radius for preconditioner
    double d_mu, d_body_effect_hrad;

};
}// namespace IBAMR


//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FreeBodyMobilitySolver
