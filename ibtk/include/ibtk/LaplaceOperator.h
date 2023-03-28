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

#ifndef included_IBTK_LaplaceOperator
#define included_IBTK_LaplaceOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LinearOperator.h"

#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"

#include <memory>
#include <string>
#include <vector>

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
    LaplaceOperator(std::string object_name, bool homogeneous_bc = false);

    /*!
     * \brief Default destructor.
     */
    ~LaplaceOperator() = default;

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
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM> > d_default_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LaplaceOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LaplaceOperator(const LaplaceOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LaplaceOperator& operator=(const LaplaceOperator& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LaplaceOperator
