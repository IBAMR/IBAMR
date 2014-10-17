// Filename: LaplaceOperator.h
// Created on 07 Apr 2012 by Boyce Griffith
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

#ifndef included_LaplaceOperator
#define included_LaplaceOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "PoissonSpecifications.h"
#include "ibtk/LinearOperator.h"

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
 * \brief Class LaplaceOperator is an abstract base class for a Laplace-type
 * operators.
 */
class LaplaceOperator : public LinearOperator
{
public:
    /*!
     * \brief Constructor.
     */
    LaplaceOperator(const std::string& object_name, bool homogeneous_bc = false);

    /*!
     * \brief Destructor.
     */
    ~LaplaceOperator();

    /*!
     * \brief Set the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients for the scalar-valued or vector-valued Laplace operator.
     */
    virtual void setPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& poisson_spec);

    /*!
     * \brief Get the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients for the scalar-valued or vector-valued Laplace operator.
     */
    virtual const SAMRAI::solv::PoissonSpecifications& getPoissonSpecifications() const;

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
    virtual void setPhysicalBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef);

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
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Get the SAMRAI::solv::RobinBcCoefStrategy object(s) used to
     * specify physical boundary conditions.
     */
    virtual const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& getPhysicalBcCoefs() const;

protected:
    // Problem specification.
    SAMRAI::solv::PoissonSpecifications d_poisson_spec;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_default_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LaplaceOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LaplaceOperator(const LaplaceOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LaplaceOperator& operator=(const LaplaceOperator& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LaplaceOperator
