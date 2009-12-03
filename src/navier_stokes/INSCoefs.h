#ifndef included_INSCoefs
#define included_INSCoefs

// Filename: INSCoefs.h
// Last modified: <10.Sep.2009 23:15:52 griffith@griffith-macbook-pro.local>
// Created on 26 Aug 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

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
     * \brief Destructor.
     */
    inline
    ~INSCoefs()
        {
            // intentionally blank
            return;
        }// ~INSCoefs


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
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    INSCoefs(
        const INSCoefs& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSCoefs& operator=(
        const INSCoefs& that);

    /*!
     * \brief The mass density (rho), dynamic viscosity (mu), and drag (lambda)
     * coefficients.
     */
    const double d_rho;
    const double d_mu;
    const double d_lambda;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "INSCoefs.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSCoefs
