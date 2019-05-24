// Filename: KrylovFreeBodyMobilitySolver.h
// Created on 16 Aug 2015 by Amneet Bhalla and Bakytzhan Kallemov.
//
// Copyright (c) 2002-2017, Amneet Bhalla, Bakytzhan Kallemov, and
// Boyce Griffith.
//
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

#ifndef included_IBAMR_KrylovFreeBodyMobilitySolver
#define included_IBAMR_KrylovFreeBodyMobilitySolver

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <vector>

#include "petscksp.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace IBAMR
{
class CIBStrategy;
class CIBMobilitySolver;
class StokesSpecifications;
class DirectMobilitySolver;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * We are trying to solve the problem
 *
 * \f$ Nx = [T M^{-1} T^{*}]x = b \f$; for \f$ x \f$,
 *
 * in this class. Here, \f$ N \f$ is the body mobility matrix, \f$ M \f$ is the mobility matrix,
 * \f$ T \f$ is the rigid body operator, and \f$ T^{*} \f$ is the rigid body velocity
 * operator. We employ Krylov solver preconditioned by direct solver to solve the
 * body-mobility equation.
 *
 */

class KrylovFreeBodyMobilitySolver : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor for \f$ [T M^{-1} T^{*}]^{-1} \f$ solver that employs the
     * PETSc KSP solver framework.
     */
    KrylovFreeBodyMobilitySolver(std::string object_name,
                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                 std::string default_options_prefix,
                                 SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy,
                                 MPI_Comm petsc_comm = PETSC_COMM_WORLD);

    /*!
     * \brief Destructor.
     */
    virtual ~KrylovFreeBodyMobilitySolver();

    /*!
     * \brief Set the mobility solver for this class.
     */
    void setMobilitySolver(SAMRAI::tbox::Pointer<IBAMR::CIBMobilitySolver> mobility_solver);

    /*!
     * \brief Set scale for interp operator.
     */
    void setInterpScale(const double interp_scale);

    /*!
     * \brief Set scale for spread operator.
     */
    void setSpreadScale(const double spread_scale);

    /*!
     * \brief Set stokes specifications.
     */
    void setStokesSpecifications(const IBAMR::StokesSpecifications& stokes_spec);

    /*!
     * \brief Set the KSP type.
     */
    void setKSPType(const std::string& ksp_type);

    /*!
     * \brief Set the options prefix used by this PETSc solver object.
     */
    void setOptionsPrefix(const std::string& options_prefix);

    /*!
     * \name Functions to access the underlying PETSc objects.
     */
    //\{

    /*!
     * \brief Get the PETSc KSP object.
     */
    const KSP& getPETScKSP() const;

    //\}

    /*!
     * \brief Solve the linear system of equations \f$ Nx = b \f$ for \f$ x \f$.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveSystem(Vec x, Vec b);

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$ Nx = b \f$.
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

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    KrylovFreeBodyMobilitySolver(const KrylovFreeBodyMobilitySolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    KrylovFreeBodyMobilitySolver& operator=(const KrylovFreeBodyMobilitySolver& that) = delete;

    /*!
     * \brief Report the KSPConvergedReason.
     */
    void reportKSPConvergedReason(const KSPConvergedReason& reason, std::ostream& os) const;

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
    static PetscErrorCode MatVecMult_KFBMSolver(Mat A, Vec x, Vec y);

    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply_KFBMSolver(PC pc, Vec x, Vec y);

    /*!
     * \brief Set KSP monitoring routine for the KSP.
     */
    static PetscErrorCode monitorKSP(KSP ksp, int it, PetscReal rnorm, void* mctx);

    //\}

    // Solver stuff
    std::string d_object_name, d_ksp_type = KSPGMRES, d_pc_type = "shell";
    bool d_is_initialized = false;
    bool d_reinitializing_solver = false;
    Vec d_petsc_b = nullptr, d_petsc_temp_v = nullptr, d_petsc_temp_f = nullptr;
    std::string d_options_prefix;
    MPI_Comm d_petsc_comm;
    KSP d_petsc_ksp = nullptr;
    Mat d_petsc_mat = nullptr;

    int d_max_iterations = 10000, d_current_iterations;
    double d_abs_residual_tol = 1.0e-50, d_rel_residual_tol = 1.0e-5;
    double d_current_residual_norm;
    bool d_initial_guess_nonzero = true;
    bool d_enable_logging = false;

    // Pointers
    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;
    SAMRAI::tbox::Pointer<IBAMR::CIBMobilitySolver> d_mobility_solver;

    // The current, new, solution and time-interval of integration.
    double d_current_time = std::numeric_limits<double>::signaling_NaN(),
           d_new_time = std::numeric_limits<double>::signaling_NaN(),
           d_solution_time = std::numeric_limits<double>::signaling_NaN(),
           d_dt = std::numeric_limits<double>::signaling_NaN();

    // Interp and spread scale.
    double d_interp_scale, d_spread_scale;

    // System physical parameters.
    double d_mu = 1.0;
    double d_rho = 1.0;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_KrylovFreeBodyMobilitySolver
