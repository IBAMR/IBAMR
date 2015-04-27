// Filename: IBImplicitStrategy.h
// Created on 21 Sep 2011 by Boyce Griffith
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

#ifndef included_IBImplicitStrategy
#define included_IBImplicitStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "ibamr/IBStrategy.h"
#include "petscvec.h"
#include "petscmat.h"

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
    IBImplicitStrategy();

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBImplicitStrategy();

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
    virtual void setLinearizedPosition(Vec& X_vec) = 0;

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
    virtual void constructLagrangianForceJacobian(Mat& A, MatType mat_type) = 0;

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
                                   const int stencil_width,
                                   const std::vector<int>& num_dofs_per_proc,
                                   const int dof_index_idx) = 0;

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitStrategy(const IBImplicitStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitStrategy& operator=(const IBImplicitStrategy& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBImplicitStrategy
