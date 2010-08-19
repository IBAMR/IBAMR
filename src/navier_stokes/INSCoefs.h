// Filename: INSCoefs.h
// Created on 26 Aug 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef included_INSCoefs
#define included_INSCoefs

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/DescribedClass.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSCoefs is a lightweight utility class which is used to specify
 * the physical parameters of the incompressible Navier-Stokes equations.
 */
class INSCoefs
    : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    inline
    INSCoefs(
        const double rho,
        const double mu,
        const double lambda)
        : d_rho(rho),
          d_mu(mu),
          d_lambda(lambda)
        {
            // intentionally blank
            return;
        }// INSCoefs

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    INSCoefs(
        const INSCoefs& from)
        : d_rho(from.d_rho),
          d_mu(from.d_mu),
          d_lambda(from.d_lambda)
        {
            // intentionally blank
            return;
        }// INSCoefs

    /*!
     * \brief Destructor.
     */
    inline
    ~INSCoefs()
        {
            // intentionally blank
            return;
        }// ~INSCoefs

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSCoefs& operator=(
        const INSCoefs& that)
        {
            if (&that != this)
            {
                d_rho = that.d_rho;
                d_mu = that.d_mu;
                d_lambda = that.d_lambda;
            }
            return *this;
        }// operator=

    /*!
     * \return The mass density coefficient of the fluid.
     */
    inline double
    getRho() const
        {
            return d_rho;
        }// getRho

    /*!
     * \return The dynamic viscosity coefficient of the fluid.
     */
    inline double
    getMu() const
        {
            return d_mu;
        }// getMu

    /*!
     * \return The drag coefficient of the fluid.
     */
    inline double
    getLambda() const
        {
            return d_lambda;
        }// getLambda

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    INSCoefs();

    /*!
     * \brief The mass density (rho), dynamic viscosity (mu), and drag (lambda)
     * coefficients.
     */
    double d_rho;
    double d_mu;
    double d_lambda;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "INSCoefs.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSCoefs
