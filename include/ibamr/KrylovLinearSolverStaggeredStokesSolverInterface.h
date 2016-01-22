// Filename: KrylovLinearSolverStaggeredStokesSolverInterface.h
// Created on 16 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_KrylovLinearSolverStaggeredStokesSolverInterface
#define included_KrylovLinearSolverStaggeredStokesSolverInterface

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "PoissonSpecifications.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "tbox/Pointer.h"

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
    KrylovLinearSolverStaggeredStokesSolverInterface();

    /*!
     * Destructor.
     */
    ~KrylovLinearSolverStaggeredStokesSolverInterface();

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& U_problem_coefs);

    /*!
     * \brief Set if velocity and pressure have nullspace.
     */
    void setComponentsHaveNullspace(const bool has_velocity_nullspace, const bool has_pressure_nullspace);

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
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef);

    /*!
     * \brief Set the StokesSpecifications object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    KrylovLinearSolverStaggeredStokesSolverInterface(const KrylovLinearSolverStaggeredStokesSolverInterface& from);

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
    operator=(const KrylovLinearSolverStaggeredStokesSolverInterface& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_KrylovLinearSolverStaggeredStokesSolverInterface
