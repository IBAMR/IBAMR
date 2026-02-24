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

#ifndef included_IBAMR_StaggeredStokesSolver
#define included_IBAMR_StaggeredStokesSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/GeneralSolver.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAILocationIndexRobinBcCoefs.h>
#include <SAMRAIPointer.h>
#include <SAMRAIPoissonSpecifications.h>
#include <SAMRAIRobinBcCoefStrategy.h>

#include <vector>

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesSolver is an abstract base class for
 * staggered-grid Stokes solvers.
 */
class StaggeredStokesSolver : public virtual IBTK::GeneralSolver
{
public:
    /*!
     * \brief Deafult constructor.
     */
    StaggeredStokesSolver();

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesSolver() = default;

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    virtual void setVelocityPoissonSpecifications(const SAMRAIPoissonSpecifications& U_problem_coefs);

    /*!
     * \brief Set if velocity and pressure have nullspace.
     */
    virtual void setComponentsHaveNullSpace(const bool has_velocity_nullspace, const bool has_pressure_nullspace);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be nullptr.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be nullptr; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  IBTK::Vector of pointers to objects that can set the Robin boundary
     *condition coefficients for the velocity
     * \param P_bc_coef   Pointer to object that can set the Robin boundary condition
     *coefficients
     *for the pressure
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& U_bc_coefs,
                                    SAMRAIRobinBcCoefStrategy* P_bc_coef);

    /*!
     * \brief Set the StokesSpecifications object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    virtual void setPhysicalBoundaryHelper(SAMRAIPointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper);

protected:
    // Problem specification.
    SAMRAIPoissonSpecifications d_U_problem_coefs;
    SAMRAILocationIndexRobinBcCoefs d_default_U_bc_coef;
    std::vector<SAMRAIRobinBcCoefStrategy*> d_U_bc_coefs;
    SAMRAILocationIndexRobinBcCoefs d_default_P_bc_coef;
    SAMRAIRobinBcCoefStrategy* d_P_bc_coef;

    // Boundary condition helper object.
    SAMRAIPointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;

    // Null space info
    bool d_has_velocity_nullspace = false, d_has_pressure_nullspace = false;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesSolver(const StaggeredStokesSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesSolver& operator=(const StaggeredStokesSolver& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesSolver
