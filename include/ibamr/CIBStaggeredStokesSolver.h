// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_CIBStaggeredStokesSolver
#define included_IBAMR_CIBStaggeredStokesSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/StaggeredStokesSolver.h>

#include <ibtk/ibtk_utilities.h>
#include <ibtk/samrai_compatibility_names.h>

#include <petscvec.h>

#include <SAMRAIDatabase.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPointer.h>
#include <SAMRAIPoissonSpecifications.h>
#include <SAMRAIRobinBcCoefStrategy.h>
#include <SAMRAISAMRAIVectorReal.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIVariableContext.h>

#include <string>

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
                             SAMRAIPointer<SAMRAIDatabase> input_db,
                             SAMRAIPointer<IBAMR::INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                             SAMRAIPointer<IBAMR::CIBStrategy> cib_strategy,
                             const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~CIBStaggeredStokesSolver();

    /*!
     * Initialize the solver before solving the system of equations.
     */
    virtual void initializeSolverState(const SAMRAISAMRAIVectorReal<double>& x,
                                       const SAMRAISAMRAIVectorReal<double>& b) override;

    /*!
     * Deallocate the solver when hierarchy changes.
     */
    virtual void deallocateSolverState() override;

    /*!
     * \brief Solve the system of equations.
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise.
     */
    virtual bool solveSystem(SAMRAISAMRAIVectorReal<double>& x, SAMRAISAMRAIVectorReal<double>& b) override;

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    virtual void setVelocityPoissonSpecifications(const SAMRAIPoissonSpecifications& u_problem_coefs) override;

    /*!
     *\brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a u_bc_coefs may be nullptr.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be nullptr; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param u_bc_coefs Vector of pointers to objects that can set the
     * Robin boundary condition coefficients for the velocity.
     *
     * \param p_bc_coef Pointer to object that can set the Robin boundary condition
     * coefficients for the pressure.
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& u_bc_coefs,
                                    SAMRAIRobinBcCoefStrategy* p_bc_coef) override;

    /*!
     * \brief Set the StaggeredStokesPhysicalBoundaryHelper object to be used
     * in the StokesOperator.
     */
    virtual void setPhysicalBoundaryHelper(SAMRAIPointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper) override;

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    virtual void setSolutionTime(double solution_time) override;

    /*!
     * \brief Set the current time interval.
     */
    virtual void setTimeInterval(double current_time, double new_time) override;

    /*!
     * \brief Get the saddle-point solver for the extended Stokes
     * system.
     */
    SAMRAIPointer<IBAMR::CIBSaddlePointSolver> getSaddlePointSolver() const;

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
    CIBStaggeredStokesSolver(const CIBStaggeredStokesSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CIBStaggeredStokesSolver& operator=(const CIBStaggeredStokesSolver& that) = delete;

    // Pointer to the constraint IB strategy.
    SAMRAIPointer<CIBStrategy> d_cib_strategy;

    // Number of structures.
    const unsigned int d_num_rigid_parts;

    // Pointer to saddle-point solver.
    SAMRAIPointer<IBAMR::CIBSaddlePointSolver> d_sp_solver;

    // Patch data to support delta function.
    SAMRAIPointer<SAMRAISideVariable<double>> d_wide_u_var, d_wide_f_var;
    SAMRAIPointer<SAMRAIVariableContext> d_wide_ctx;
    int d_wide_u_idx = IBTK::invalid_index, d_wide_f_idx = IBTK::invalid_index;

    // SVR for holding widened u/f.
    SAMRAIPointer<SAMRAISAMRAIVectorReal<double>> d_x_wide, d_b_wide;

    // Bools to control initialization and deallocation
    bool d_is_initialized = false, d_reinitializing_solver = false;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_CIBStaggeredStokesSolver
