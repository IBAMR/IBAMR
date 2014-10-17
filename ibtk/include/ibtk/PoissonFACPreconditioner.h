// Filename: PoissonFACPreconditioner.h
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

#ifndef included_PoissonFACPreconditioner
#define included_PoissonFACPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "PoissonSpecifications.h"
#include "ibtk/FACPreconditioner.h"
#include "ibtk/PoissonSolver.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class PoissonFACPreconditionerStrategy;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PoissonFACPreconditioner is a FACPreconditioner that has been
 * specialized for Poisson problems.
*/
class PoissonFACPreconditioner : public FACPreconditioner, public PoissonSolver
{
public:
    /*!
     * Constructor.
     */
    PoissonFACPreconditioner(const std::string& object_name,
                             SAMRAI::tbox::Pointer<PoissonFACPreconditionerStrategy> fac_strategy,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             const std::string& default_options_prefix);

    /*!
     * Destructor.
     */
    ~PoissonFACPreconditioner();

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
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PoissonFACPreconditioner();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonFACPreconditioner(const PoissonFACPreconditioner& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonFACPreconditioner& operator=(const PoissonFACPreconditioner& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PoissonFACPreconditioner
