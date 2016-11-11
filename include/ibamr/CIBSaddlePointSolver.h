// Filename: CIBSaddlePointSolver.h
// Created on 10 Nov 2014 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBAMR_CIBSaddlePointSolver
#define included_IBAMR_CIBSaddlePointSolver

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <string>
#include <vector>

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "petscksp.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace solv
{
class PoissonSpecifications;
template <int DIM, typename T>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
namespace IBAMR
{
class CIBStrategy;
class CIBStaggeredStokesOperator;
class StaggeredStokesPhysicalBoundaryHelper;
class INSStaggeredHierarchyIntegrator;
class StaggeredStokesSolver;
class CIBMobilitySolver;
} // namespace IBAMR
namespace IBTK
{
class LinearOperator;
class PoissonSolver;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CIBSaddlePointSolver solves for the fluid velocity
 * \f$ \vec{v} \f$, fluid pressure \f$ p \f$, and the Lagrange multiplier \f$ \vec{\lambda} \f$
 * maintaining the rigidity constraint.
 *
 * For free-moving (self-moving bodies)
 * it also solves for the rigid body translational and rotational velocities
 * \f$  \vec{U} \f$.
 *
 * - For free-moving case, the class solves a saddle-point problem of the
 * type:
 *
 * \f{eqnarray*}{
 *  A \vec{v} + \nabla p - \gamma S \vec{\lambda}   &=& \vec{g}  \\
 *  -\nabla \cdot \vec{v}                           &=&  (h = 0) \\
 *  -\beta J \vec{v} + \beta T^{*} \vec{U} -\beta \delta \vec{\lambda} &=&  (\vec{w} = -\beta \vec{U}_{\mbox{\scriptsize
 def}}) \\
 *  T \vec{\lambda} &=&  \vec{F}
 *  \f}

 *
 * - For imposed kinematics case, the class solves:
 *
 * \f{eqnarray*}{
 *   A \vec{v} + \nabla p - \gamma S \vec{\lambda} &=& \vec{g} \\
 *   -\nabla \cdot \vec{v}   &=& (h = 0) \\
 *   -\beta J \vec{v}   - \beta \delta \vec{\lambda}  &=& (w = -\beta (T^{*}\vec{U}+ \vec{U}_{\mbox{\scriptsize def}}) =
 -\beta \vec{U}_{\mbox{\scriptsize imposed}})
 * \f}
 *
 * Here, \f$ A \f$ is the Helmholtz operator, \f$ S \f$ is
 * the spreading operator, \f$ G \f$ is the gradient operator, \f$ J \f$
 * is the interpolation operator, and \f$ T \f$ is the rigid body operator
 * (with \f$ T^{*} \f$ denoting its adjoint), \f$ \vec{F} \f$ is the net external
 * force/torque on the body (excluding the hydrodynamic force/torque),
 * \f$ \vec{U}_{\mbox{\scriptsize def}} \f$ is the deformational velocity of the body,
 * and \f$ \vec{U}_{\mbox{\scriptsize imposed}} \f$ is
 * the velocity imposed on the body. \f$ \gamma \f$ and \f$ \beta \f$ are two scaling
 * parameters to make the system well-scaled. \f$ \delta \f$ is the regularization parameter
 * that regularizes the problem in presence of many IB/control points in a grid cell.
 *
 * Here, we employ the Krylov solver to solve the above saddle-point problem. We use
 * Schur complement preconditioner to precondition the iterative solver.
 */
class CIBSaddlePointSolver : public SAMRAI::tbox::DescribedClass
{
    //////////////////////////////////////////////////////////////////////////////
public:
    /*!
     * \brief Constructor for a saddle-point solver that employs the
     * PETSc KSP solver framework.
     */
    CIBSaddlePointSolver(const std::string& object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                         SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy,
                         const std::string& default_options_prefix,
                         MPI_Comm petsc_comm = PETSC_COMM_WORLD);

    /*!
     * \brief Destructor.
     */
    ~CIBSaddlePointSolver();

    /*!
     * \brief Set the KSP type.
     */
    void setKSPType(const std::string& ksp_type);

    /*!
     * \brief Set the options prefix used by this PETSc solver object.
     */
    void setOptionsPrefix(const std::string& options_prefix);

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& u_problem_coefs);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a u_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a p_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param u_bc_coefs  IBTK::Vector of pointers to objects that can set the
     * Robin boundary condition coefficients for the velocity.
     *
     * \param p_bc_coef   Pointer to object that can set the Robin boundary
     * condition coefficients for the pressure.
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_bc_coef);

    /*!
     * \brief Set the StokesSpecifications object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesPhysicalBoundaryHelper> bc_helper);

    /*!
     * \brief Return the linear operator for the saddle-point solver.
     */
    SAMRAI::tbox::Pointer<IBTK::LinearOperator> getA() const;

    /*!
     * \brief Return the Stokes solver used in the preconditioner of the
     * solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesSolver> getStokesSolver() const;

    //\{
    // Return the scaling factors used by scale the system of equations.
    /*!
     * \brief Return the interpolation scale factor.
     */
    inline double getInterpScale() const
    {
        return d_scale_interp;

    } // getInterpScale

    /*!
     * \brief Return the spreading scale factor.
     */
    inline double getSpreadScale() const
    {
        return d_scale_spread;

    } // getSpreadScale

    /*!
     * \brief Return the regularization scale factor.
     */
    inline double getRegularizationWeight() const
    {
        return d_reg_mob_factor;

    } // getRegularizationWeight

    /*!
     * \brief Return if the spread force is to be normalized.
     */
    inline bool getNormalizeSpreadForce() const
    {
        return d_normalize_spread_force;

    } // getNormalizeSpreadForce

    //\}

    /*!
     * \brief Solve the linear system of equations \f$Ax=b\f$ for \f$x\f$.
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveSystem(Vec x, Vec b);

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     *
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
     * \brief Get the mobility solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::CIBMobilitySolver> getCIBMobilitySolver() const;

    //////////////////////////////////////////////////////////////////////////////
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CIBSaddlePointSolver(const CIBSaddlePointSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CIBSaddlePointSolver& operator=(const CIBSaddlePointSolver& that);

    /*!
     * \brief Report the KSPConvergedReason.
     */
    void reportKSPConvergedReason(const KSPConvergedReason& reason, std::ostream& os) const;

    /*!
     * \brief Get options from input file.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     *\brief Routine to destroy KSP object.
     */
    void destroyKSP();

    /*!
     * \brief Initialize the Stokes solver needed in the preconditioning step.
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
     * \name Static functions for use by PETSc KSP and MatShell objects.
     */
    //\{

    /*!
     * \brief Compute the matrix vector product \f$y=Ax\f$.
     */
    static PetscErrorCode MatVecMult_SaddlePoint(Mat A, Vec x, Vec y);

    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_SaddlePoint(PC pc, Vec x, Vec y);

    /*!
     * \brief Set KSP monitoring routine for the KSP.
     */
    static PetscErrorCode monitorKSP(KSP ksp, int it, PetscReal rnorm, void* mctx);

    //\}

    // Solver stuff
    std::string d_object_name, d_ksp_type, d_pc_type;
    bool d_is_initialized;
    bool d_reinitializing_solver;

    Vec d_petsc_x, d_petsc_b;
    std::string d_options_prefix;

    MPI_Comm d_petsc_comm;
    KSP d_petsc_ksp;
    Mat d_petsc_mat;
    SAMRAI::tbox::Pointer<IBAMR::CIBStaggeredStokesOperator> d_A;

    int d_max_iterations, d_current_iterations;
    double d_abs_residual_tol, d_rel_residual_tol;
    double d_current_residual_norm;
    bool d_initial_guess_nonzero;
    bool d_enable_logging;

    // Preconditioner stuff
    SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> d_ins_integrator;
    SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesSolver> d_LInv;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_velocity_solver, d_pressure_solver;
    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;
    SAMRAI::tbox::Pointer<IBAMR::CIBMobilitySolver> d_mob_solver;

    // Book-keeping
    const unsigned int d_num_rigid_parts;

    // Scales used for interpolation and spreading operators.
    double d_scale_interp, d_scale_spread, d_reg_mob_factor;
    bool d_normalize_spread_force;

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
    double d_current_time, d_new_time;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_CIBSaddlePointSolver
