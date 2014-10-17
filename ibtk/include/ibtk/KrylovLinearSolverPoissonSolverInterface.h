// Filename: KrylovLinearSolverPoissonSolverInterface.h
// Created on 13 Aug 2012 by Boyce Griffith
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

#ifndef included_KrylovLinearSolverPoissonSolverInterface
#define included_KrylovLinearSolverPoissonSolverInterface

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "PoissonSpecifications.h"
#include "ibtk/PoissonSolver.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class KrylovLinearSolverPoissonSolverInterface provides an interface
 * for KrylovLinearSolvers that are to be used as Poisson solvers.
 *
 * This class is intented to be used to create a (trivial) subclass of an
 * existing implementation of KrylovLinearSolver that also supports the
 * PoissonSolver interface.
 *
 * \see PETScKrylovPoissonSolver
 */
class KrylovLinearSolverPoissonSolverInterface : public PoissonSolver
{
public:
    /*!
     * Default constructor.
     */
    KrylovLinearSolverPoissonSolverInterface();

    /*!
     * Destructor.
     */
    ~KrylovLinearSolverPoissonSolverInterface();

    /*!
     * \brief Set the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients for the scalar-valued or vector-valued Laplace operator.
     */
    void setPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& poisson_spec);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy object used to specify
     * physical boundary conditions.
     *
     * \note \a bc_coef may be NULL.  In this case, default boundary conditions
     * (as supplied to the class constructor) are employed.
     *
     * \param bc_coef  Pointer to an object that can set the Robin boundary condition
     *coefficients
     */
    void setPhysicalBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a bc_coefs may be NULL.  In this case,
     * default boundary conditions (as supplied to the class constructor) are
     * employed for that data depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition
     *coefficients
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    KrylovLinearSolverPoissonSolverInterface(const KrylovLinearSolverPoissonSolverInterface& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    KrylovLinearSolverPoissonSolverInterface& operator=(const KrylovLinearSolverPoissonSolverInterface& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_KrylovLinearSolverPoissonSolverInterface
