// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_KrylovLinearSolverStaggeredStokesSolverInterface
#define included_IBAMR_KrylovLinearSolverStaggeredStokesSolverInterface

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StaggeredStokesSolver.h"

#include "PoissonSpecifications.h"
#include "tbox/Pointer.h"

#include <vector>

namespace IBAMR
{
class StaggeredStokesPhysicalBoundaryHelper;
} // namespace IBAMR
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
 * \brief Class KrylovLinearSolverStaggeredStokesSolverInterface provides an
 * interface for KrylovLinearSolvers that are to be used as staggered Stokes
 * solvers.
 *
 * This class is intented to be used to create a (trivial) subclass of an
 * existing implementation of KrylovLinearSolver that also supports the
 * StaggeredStokesSolver interface.
 *
 * \see PETScKrylovStaggeredStokesSolver
 */
class KrylovLinearSolverStaggeredStokesSolverInterface : public StaggeredStokesSolver
{
public:
    /*!
     * Default constructor.
     */
    KrylovLinearSolverStaggeredStokesSolverInterface() = default;

    /*!
     * Destructor.
     */
    ~KrylovLinearSolverStaggeredStokesSolverInterface() = default;

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& U_problem_coefs) override;

    /*!
     * \brief Set if velocity and pressure have nullspace.
     */
    void setComponentsHaveNullspace(const bool has_velocity_nullspace, const bool has_pressure_nullspace) override;

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  IBTK::Vector of pointers to objects that can set the Robin boundary
     *condition coefficients for the velocity
     * \param P_bc_coef   Pointer to object that can set the Robin boundary condition
     *coefficients
     *for the pressure
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef) override;

    /*!
     * \brief Set the StokesSpecifications object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper) override;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    KrylovLinearSolverStaggeredStokesSolverInterface(const KrylovLinearSolverStaggeredStokesSolverInterface& from) =
        delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    KrylovLinearSolverStaggeredStokesSolverInterface&
    operator=(const KrylovLinearSolverStaggeredStokesSolverInterface& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_KrylovLinearSolverStaggeredStokesSolverInterface
