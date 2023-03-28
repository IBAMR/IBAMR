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

#ifndef included_IBTK_PoissonFACPreconditioner
#define included_IBTK_PoissonFACPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/FACPreconditioner.h"
#include "ibtk/PoissonSolver.h"

#include "PoissonSpecifications.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

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
                             std::string default_options_prefix);

    /*!
     * Destructor.
     */
    ~PoissonFACPreconditioner() = default;

    /*!
     * \brief Set the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients for the scalar-valued or vector-valued Laplace operator.
     */
    void setPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& poisson_spec) override;

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
    void setPhysicalBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef) override;

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
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PoissonFACPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonFACPreconditioner(const PoissonFACPreconditioner& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonFACPreconditioner& operator=(const PoissonFACPreconditioner& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PoissonFACPreconditioner
