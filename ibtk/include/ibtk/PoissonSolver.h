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

#ifndef included_IBTK_PoissonSolver
#define included_IBTK_PoissonSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/GeneralSolver.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIPointer.h"
#include "SAMRAIPoissonSpecifications.h"
#include "SAMRAIRobinBcCoefStrategy.h"

#include <memory>
#include <string>
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

namespace IBTK
{
/*!
 * \brief Class PoissonSolver is an abstract base class for Poisson solvers.
 */
class PoissonSolver : public virtual GeneralSolver
{
public:
    /*!
     * \brief Default constructor.
     */
    PoissonSolver() = default;

    /*!
     * \brief Default destructor.
     */
    ~PoissonSolver() = default;

    /*!
     * \brief Set the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients for the scalar-valued or vector-valued Laplace operator.
     */
    virtual void setPoissonSpecifications(const SAMRAIPoissonSpecifications& poisson_spec);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy object used to specify
     * physical boundary conditions.
     *
     * \note \a bc_coef may be nullptr.  In this case, default boundary conditions
     * (as supplied to the class constructor) are employed.
     *
     * \param bc_coef  Pointer to an object that can set the Robin boundary condition
     *coefficients
     */
    virtual void setPhysicalBcCoef(SAMRAIRobinBcCoefStrategy* bc_coef);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a bc_coefs may be nullptr.  In this case,
     * default boundary conditions (as supplied to the class constructor) are
     * employed for that data depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition
     *coefficients
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs);

protected:
    // Specialized solver initialization.
    void initSpecialized(const std::string& object_name, bool homogeneous_bc) override;

    // Problem specification.
    SAMRAIPoissonSpecifications d_poisson_spec = SAMRAIPoissonSpecifications("");
    std::unique_ptr<SAMRAIRobinBcCoefStrategy> d_default_bc_coef;
    std::vector<SAMRAIRobinBcCoefStrategy*> d_bc_coefs;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonSolver(const PoissonSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonSolver& operator=(const PoissonSolver& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PoissonSolver
