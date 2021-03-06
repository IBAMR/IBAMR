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
     * Set the value of the intermediate position vector used in evaluating the
     * linearized problem.
     */
    virtual void setLinearizedPosition(Vec& X_vec, double data_time) = 0;

    /*!
     * Compute the nonlinear residual.
     */
    virtual void computeResidual(Vec& R_vec) = 0;

    /*!
     * Compute the linearized residual for the given intermediate position
     * vector.
     */
    virtual void computeLinearizedResidual(Vec& X_vec, Vec& R_vec) = 0;

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval for use in evaluating the
     * residual of the linearized problem.
     */
    virtual void interpolateLinearizedVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time) = 0;

    /*!
     * Compute the Lagrangian force of the linearized problem for the specified
     * configuration of the updated position vector.
     */
    virtual void computeLinearizedLagrangianForce(Vec& X_vec, double data_time) = 0;

    /*!
     * Construct the linearized Lagrangian force Jacobian.
     */
    virtual void constructLagrangianForceJacobian(Mat& A, MatType mat_type, double data_time) = 0;

    /*!
     * Spread the Lagrangian force of the linearized problem to the Cartesian
     * grid at the specified time within the current time interval.
     */
    virtual void spreadLinearizedForce(
        int f_data_idx,
        IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
        double data_time) = 0;

    /*!
     * Construct the IB interpolation operator.
     */
    virtual void constructInterpOp(Mat& J,
                                   void (*spread_fnc)(const double, double*),
                                   int stencil_width,
                                   const std::vector<int>& num_dofs_per_proc,
                                   int dof_index_idx,
                                   double data_time) = 0;

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

#endif //#ifndef included_IBAMR_IBImplicitStrategy
