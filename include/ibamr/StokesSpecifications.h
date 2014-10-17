// Filename: StokesSpecifications.h
// Created on 26 Aug 2007 by Boyce Griffith
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

#ifndef included_StokesSpecifications
#define included_StokesSpecifications

/////////////////////////////// INCLUDES /////////////////////////////////////

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
    inline ~StokesSpecifications()
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

#endif //#ifndef included_StokesSpecifications
