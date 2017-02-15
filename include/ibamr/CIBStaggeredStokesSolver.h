// Filename: CIBStaggeredStokesSolver.h
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

#ifndef included_IBAMR_CIBStaggeredStokesSolver
#define included_IBAMR_CIBStaggeredStokesSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "IntVector.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
namespace IBAMR
{
class CIBSaddlePointSolver;
class CIBStrategy;
class INSStaggeredHierarchyIntegrator;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CIBStaggeredStokesSolver is an extension of
 * IBAMR::StaggeredStokesSolver class that solves for the Langrange multiplier
 * field \f$ \vec{\lambda} \f$ which imposes the rigidity constraint and solves
 * for rigid body velocity (if it is freely-moving) along with the fluid
 * velocity and pressure.
 */
class CIBStaggeredStokesSolver : public StaggeredStokesSolver
{
    //////////////////////////////////////////////////////////////////////////////
public:
    /*!
     * \brief Constructor.
     */
    CIBStaggeredStokesSolver(const std::string& object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                             SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy,
                             const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~CIBStaggeredStokesSolver();

    /*!
     * Initialize the solver before solving the system of equations.
     */
    virtual void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                       const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

    /*!
     * Deallocate the solver when hierarchy changes.
     */
    virtual void deallocateSolverState();

    /*!
     * \brief Solve the system of equations.
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise.
     */
    virtual bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                             SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    virtual void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& u_problem_coefs);

    /*!
     *\brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a u_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param u_bc_coefs Vector of pointers to objects that can set the
     * Robin boundary condition coefficients for the velocity.
     *
     * \param p_bc_coef Pointer to object that can set the Robin boundary condition
     * coefficients for the pressure.
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_bc_coef);

    /*!
     * \brief Set the StaggeredStokesPhysicalBoundaryHelper object to be used
     * in the StokesOperator.
     */
    virtual void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper);

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    virtual void setSolutionTime(double solution_time);

    /*!
     * \brief Set the current time interval.
     */
    virtual void setTimeInterval(double current_time, double new_time);

    /*!
     * \brief Get the saddle-point solver for the extended Stokes
     * system.
     */
    SAMRAI::tbox::Pointer<IBAMR::CIBSaddlePointSolver> getSaddlePointSolver() const;

    //////////////////////////////////////////////////////////////////////////////
protected:
    //////////////////////////////////////////////////////////////////////////////

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CIBStaggeredStokesSolver(const CIBStaggeredStokesSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CIBStaggeredStokesSolver& operator=(const CIBStaggeredStokesSolver& that);

    // Pointer to the constraint IB strategy.
    SAMRAI::tbox::Pointer<CIBStrategy> d_cib_strategy;

    // Number of structures.
    const unsigned int d_num_rigid_parts;

    // Pointer to saddle-point solver.
    SAMRAI::tbox::Pointer<IBAMR::CIBSaddlePointSolver> d_sp_solver;

    // Patch data to support delta function.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_wide_u_var, d_wide_f_var;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_wide_ctx;
    int d_wide_u_idx, d_wide_f_idx;

    // SVR for holding widened u/f.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_x_wide, d_b_wide;

    // Bools to control initialization and deallocation
    bool d_is_initialized, d_reinitializing_solver;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_CIBStaggeredStokesSolver
