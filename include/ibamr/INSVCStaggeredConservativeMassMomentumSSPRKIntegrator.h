// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_INSVCStaggeredConservativeMassMomentumSSPRKIntegrator
#define included_IBAMR_INSVCStaggeredConservativeMassMomentumSSPRKIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/INSVCStaggeredConservativeMassMomentumRKIntegrator.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSVCStaggeredConservativeMassMomentumSSPRKIntegrator is a derived class that integrates
 * the staggered density field
 *
 *  \f$ \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) = S(x,t) \f$
 *
 * and computes the conservative form of the convective operator
 * \f$ \nabla \cdot (\rho u u)\f$.
 *
 * This class implements the SSP-RK2 and SSP-RK3 time stepping schemes.
 *
 * Class INSVCStaggeredConservativeMassMomentumSSPRKIntegrator computes the convective
 * derivative of a side-centered velocity field using various bounded-limiters
 * described by Patel and Natarajan.
 *
 * A side-centered density update is provided by this class, which is used in the
 * conservative discretization of the incompressible Navier-Stokes equation.
 * There is an optional mass density source term \f$ S(x,t) \f$ that can be set to check the order
 * of accuracy via manufactured solutions.
 *
 * References
 * Patel, JK. and Natarajan, G., <A HREF="https://www.sciencedirect.com/science/article/pii/S0045793014004009">
 * A generic framework for design of interface capturing schemes for multi-fluid flows</A>
 *
 * \note This class is specialized in that it computes a conservative discretization of the form
 * \f$N = \nabla \cdot (u \rho u)\f$, where the density \f$\rho\f$ can vary in space and time.
 * This operator is to be used in conjuction with the conservative form of the variable coefficient
 * Navier-Stokes equations, which will produce better results for high density ratio flows.
 *
 * \see INSVCStaggeredHierarchyIntegrator
 */
class INSVCStaggeredConservativeMassMomentumSSPRKIntegrator : public INSVCStaggeredConservativeMassMomentumRKIntegrator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSVCStaggeredConservativeMassMomentumSSPRKIntegrator(std::string object_name,
                                                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~INSVCStaggeredConservativeMassMomentumSSPRKIntegrator() = default;

    /*!
     * \brief Integrate density and momentum field.
     */
    void integrate(double dt) override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredConservativeMassMomentumSSPRKIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredConservativeMassMomentumSSPRKIntegrator(
        const INSVCStaggeredConservativeMassMomentumSSPRKIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredConservativeMassMomentumSSPRKIntegrator&
    operator=(const INSVCStaggeredConservativeMassMomentumSSPRKIntegrator& that) = delete;

    // Number of SSP-RK steps to take. Default is set for SSP-RK3.
    int d_num_steps = 3;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSVCStaggeredConservativeMassMomentumSSPRKIntegrator
