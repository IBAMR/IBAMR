// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_IBImplicitStrategy
#define included_IBAMR_IBImplicitStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBStrategy.h"

#include "petscmat.h"
#include "petscvec.h"

#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK
namespace SAMRAI
{
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBImplicitStrategy provides a generic interface for specifying
 * the implementation details of a particular implicit version of the IB method.
 */
class IBImplicitStrategy : public IBStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBImplicitStrategy() = default;

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBImplicitStrategy() = default;

    /*!
     * Create solution and rhs data.
     */
    virtual void createSolverVecs(Vec* X_vec, Vec* F_vec) = 0;

    /*!
     * Setup solution and rhs data.
     */
    virtual void setupSolverVecs(Vec* X_vec, Vec* F_vec) = 0;

    /*!
     * Set the value of the updated position vector.
     */
    virtual void setUpdatedPosition(Vec& X_new_vec) = 0;

    /*!
     * Compute the nonlinear residual for backward Euler time stepping.
     */
    virtual void computeResidualBackwardEuler(Vec& R_vec) = 0;

    /*!
     * Compute the nonlinear residual for midpoint rule time stepping.
     */
    virtual void computeResidualMidpointRule(Vec& R_vec) = 0;

    /*!
     * Compute the nonlinear residual for trapezoidal rule time stepping.
     */
    virtual void computeResidualTrapezoidalRule(Vec& R_vec) = 0;

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitStrategy(const IBImplicitStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitStrategy& operator=(const IBImplicitStrategy& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_IBImplicitStrategy
