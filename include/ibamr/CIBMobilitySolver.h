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

#ifndef included_IBAMR_CIBMobilitySolver
#define included_IBAMR_CIBMobilitySolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/CIBStrategy.h"
#include "ibamr/DirectMobilitySolver.h"
#include "ibamr/KrylovFreeBodyMobilitySolver.h"
#include "ibamr/KrylovMobilitySolver.h"
#include "ibamr/ibamr_enums.h"

#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include "petscksp.h"
#include "petscvec.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox

namespace solv
{
class PoissonSpecifications;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI
namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
class StaggeredStokesPhysicalBoundaryHelper;
class CIBStrategy;
class DirectMobilitySolver;
class KrylovMobilitySolver;
class KrylovFreeBodyMobilitySolver;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{
/*!
 * \brief Class CIBMobilitySolver solves for the constraint forces \f$ \vec{\lambda}\f$
 * and rigid body velocity \f$ \vec{U}\f$ of the structure(s).
 * Specifically, the class solves two types of matrix-equations
 *
 *  \f{eqnarray*}{
 *      \beta J L^{-1} S \gamma \vec{\lambda} &=& M \vec{\lambda} = \vec{w} \\
 *      T M T^* \vec{U} &=& \vec{F}.
 * \f}
 *
 * Here, \f$ J \f$ is the interpolation operator, \f$ S \f$ is the spreading
 * operator, \f$ L \f$ is the incompressible Stokes operator, \f$ T \f$ is the
 * rigid body operator, \f$ \vec{w} \f$ is the desired velocity at the nodes
 * of the structure(s), and \f$ \vec{F} \f$ is the net external force and torque
 * on the body.
 *
 * This class employs direct solver for the approximate mobility and body-mobility
 * sub-problems. The approximate mobility matrix is intended to be used in the
 * preconditioning step of the overall constraint solver. The overall preconditioner
 * is implemented in \see IBAMR::CIBSaddlePointSolver class. The class also supports
 * Krylov body mobility solver for bodies moving under external force and torque.
 */

class CIBMobilitySolver : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief The only constructor of this class.
     */
    CIBMobilitySolver(std::string object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      SAMRAI::tbox::Pointer<IBAMR::INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                      SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy);

    /*!
     * \brief Destructor for this class.
     */
    virtual ~CIBMobilitySolver();

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void setSolutionTime(const double solution_time);

    /*!
     * \brief Set the time interval of integration.
     */
    void setTimeInterval(const double current_time, const double new_time);

    /*!
     * \brief Set scale for interp operator.
     */
    void setInterpScale(const double interp_scale);

    /*!
     * \brief Set scale for spread operator.
     */
    void setSpreadScale(const double spread_scale);

    /*!
     * \brief Set scale for regularizing mobility matrix.
     */
    void setRegularizeMobilityScale(const double reg_mob_scale);

    /*!
     * \brief Set if the mean of the Lagrangian force is to be subtracted
     * from the Eulerian force variable.
     *
     * \note This operation is needed for certain situations like Stokes flow
     * with periodic BCs.
     */
    void setNormalizeSpreadForce(const bool normalize_force);

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
     * \brief Solves the mobility problem.
     *
     * \param x Vec storing the Lagrange multiplier.
     *
     * \param b Vec storing the desired velocity.
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveMobilitySystem(Vec x, Vec b);

    /*!
     * \brief Solves the mobility problem.
     *
     * \param x Vec storing the rigid body velocity
     *
     * \param b Vec storing the net external generalized force
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveBodyMobilitySystem(Vec x, Vec b);

    /*!
     * \brief Initialize the mobility solver.
     *
     * \param x Vec storing the Lagrange multiplier
     *
     * \param b Vec storing the desired velocity
     *
     */
    void initializeSolverState(Vec x, Vec b);

    /*!
     * \brief Deallocate the mobility solver.
     */
    void deallocateSolverState();

    /*!
     * \brief Get access to mobility solvers.
     * \note A null argument is simply skipped with no corresponding solver
     * return.
     */
    void getMobilitySolvers(IBAMR::KrylovMobilitySolver** km_solver = nullptr,
                            IBAMR::DirectMobilitySolver** dm_solver = nullptr,
                            IBAMR::KrylovFreeBodyMobilitySolver** fbm_solver = nullptr);

private:
    /*!
     * \brief Get various options from input db.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    // Name of this object.
    std::string d_object_name;

    // Number of rigid bodies
    unsigned d_num_rigid_parts;

    // Pointers.
    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;
    SAMRAI::tbox::Pointer<IBAMR::DirectMobilitySolver> d_direct_mob_solver;
    SAMRAI::tbox::Pointer<IBAMR::KrylovMobilitySolver> d_krylov_mob_solver;
    SAMRAI::tbox::Pointer<IBAMR::KrylovFreeBodyMobilitySolver> d_krylov_freebody_mob_solver;

    // Other parameters.
    double d_solution_time, d_current_time, d_new_time;
    bool d_is_initialized = false, d_reinitializing_solver, d_has_free_parts;
    double d_interp_scale = 1.0, d_spread_scale = 1.0, d_reg_mob_scale;

    // Type of mobility solver to be used
    enum MobilitySolverType
    {
        DIRECT,
        KRYLOV,
        UNKNOWN_MOBILITY_SOLVER_TYPE = -1
    };
    MobilitySolverType d_mobility_solver_type;

}; // CIBMobilitySolver

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_CIBMobilitySolver
