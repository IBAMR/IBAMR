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

#ifndef included_IBAMR_StokesSpecifications
#define included_IBAMR_StokesSpecifications

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "tbox/DescribedClass.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StokesSpecifications is a lightweight utility class that is used to
 * specify the physical parameters of the incompressible Navier-Stokes
 * equations.
 */
class StokesSpecifications : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    inline StokesSpecifications(double rho = 0.0, double mu = 0.0, double lambda = 0.0)
        : d_rho(rho), d_mu(mu), d_lambda(lambda)
    {
        // intentionally blank
        return;
    } // StokesSpecifications

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    StokesSpecifications(const StokesSpecifications& from) : d_rho(from.d_rho), d_mu(from.d_mu), d_lambda(from.d_lambda)
    {
        // intentionally blank
        return;
    } // StokesSpecifications

    /*!
     * \brief Destructor.
     */
    virtual inline ~StokesSpecifications()
    {
        // intentionally blank
        return;
    } // ~StokesSpecifications

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StokesSpecifications& operator=(const StokesSpecifications& that)
    {
        if (&that != this)
        {
            d_rho = that.d_rho;
            d_mu = that.d_mu;
            d_lambda = that.d_lambda;
        }
        return *this;
    } // operator=

    /*!
     * \return The mass density coefficient of the fluid.
     */
    inline double getRho() const
    {
        return d_rho;
    } // getRho

    /*!
     * \brief Set the mass density coefficient of the fluid.
     */
    inline void setRho(double rho)
    {
        d_rho = rho;
        return;
    } // setRho

    /*!
     * \return The dynamic viscosity coefficient of the fluid.
     */
    inline double getMu() const
    {
        return d_mu;
    } // getMu

    /*!
     * \brief Set the dynamic viscosity coefficient of the fluid.
     */
    inline void setMu(double mu)
    {
        d_mu = mu;
        return;
    } // setMu

    /*!
     * \return The drag coefficient of the fluid.
     */
    inline double getLambda() const
    {
        return d_lambda;
    } // getLambda

    /*!
     * \brief Set the drag coefficient of the fluid.
     */
    inline void setLambda(double lambda)
    {
        d_lambda = lambda;
        return;
    } // setLambda

protected:
private:
    /*!
     * \brief The mass density (rho), dynamic viscosity (mu), and drag (lambda)
     * coefficients.
     */
    double d_rho;
    double d_mu;
    double d_lambda;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_StokesSpecifications
