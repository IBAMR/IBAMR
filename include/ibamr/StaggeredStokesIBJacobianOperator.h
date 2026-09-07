// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#ifndef included_IBAMR_StaggeredStokesIBJacobianOperator
#define included_IBAMR_StaggeredStokesIBJacobianOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/StaggeredStokesIBOperator.h>

#include <ibtk/JacobianOperator.h>

#include <tbox/Pointer.h>

#include <petscmat.h>

#include <string>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Stokes-IB Jacobian action on Eulerian velocity-pressure increments.
 *
 * Uses the velocity-pressure formulation described by \ref StaggeredStokesIBOperator.
 * formJacobian() uses the base velocity to reconstruct the force-evaluation position.
 * apply() adds the resulting force-derivative contribution to the Stokes
 * momentum action, leaving its pressure/divergence action unchanged. The
 * strategy path holds interpolation and spreading fixed: it does not
 * differentiate their dependence on moving coupling positions. A supplied
 * coupling matrix may instead provide the already-scaled IB contribution.
 */
class StaggeredStokesIBJacobianOperator : public IBTK::JacobianOperator
{
public:
    /*!
     * \brief Constructor.
     */
    explicit StaggeredStokesIBJacobianOperator(const std::string& object_name);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesIBJacobianOperator() override;

    /*!
     * \brief Set context data required by this operator.
     *
     * \see StaggeredStokesIBOperator::Context
     * \see StaggeredStokesIBOperator::setOperatorContext
     */
    void setOperatorContext(const StaggeredStokesIBOperator::Context& ctx);

    /*!
     * \brief Set the IB velocity-coupling contribution in coupled velocity-pressure ordering.
     *
     * The assembled matrix must include the sign and time/grid factors;
     * apply() adds it directly to the Stokes action without further scaling.
     * It must use the full coupled global numbering, communicator, and local
     * ownership of the single patch level and DOF fields in the Context (see
     * StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices()).
     * Only velocity entries may be nonzero: pressure rows and columns are not
     * masked by apply(). Multilevel supplied-matrix application is unsupported.
     * Passing nullptr selects the strategy action instead. formJacobian() does
     * not rebuild or rescale this matrix.
     * The operator retains a PETSc reference until replacement or deallocation.
     */
    void setIBCouplingJacobian(Mat& SAJ_mat);

    /*!
     * \brief Form and cache Jacobian state at the specified point.
     *
     * Requires initialized state and prepared strategy data as described in
     * StaggeredStokesIBOperator::Context. Copies x and updates the shared
     * strategy's positions and linearization state. Call again when the base
     * state, time interval, or force data change before using the strategy action.
     */
    void formJacobian(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) override;

    /*!
     * \brief Get the vector used as the Jacobian base state.
     *
     * Returns the operator-owned copy, or nullptr before the first formJacobian()
     * and after deallocation. A later formJacobian() overwrites it; deallocation
     * frees its patch components even if a caller retains the vector Pointer.
     * Do not use mutation of this vector to reform the Jacobian.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> getBaseVector() const override;

    /*!
     * \brief Compute \f$y = J[x]\f$.
     *
     * Requires initialized state and, for the strategy action, a current
     * formJacobian(). Both the Stokes action and increment interpolation use
     * homogeneous boundary data, independently of the wrapper's boundary flag.
     * Uses and updates shared IB strategy state; see StaggeredStokesIBOperator::Context.
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Compute \f$z = J[x] + y\f$.
     */
    void applyAdd(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& z) override;

    /*!
     * \brief Initialize hierarchy-dependent operator state.
     *
     * Enables and updates fixed coupling on the shared IB strategy regardless
     * of StaggeredStokesIBOperator::Context::use_fixed_le_operators, and initializes the shared Stokes
     * operator. The caller must first prepare the strategy's time-step data;
     * see StaggeredStokesIBOperator::Context.
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Deallocate hierarchy-dependent operator state.
     */
    void deallocateOperatorState() override;

    /*!
     * \brief Modify right-hand side values to account for inhomogeneous
     * boundary conditions.
     * Uses this wrapper's boundary flag and times in the shared Stokes operator.
     */
    void modifyRhsForBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Impose solution boundary conditions.
     * Uses this wrapper's boundary flag and times in the shared Stokes operator.
     */
    void imposeSolBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u) override;

private:
    StaggeredStokesIBJacobianOperator() = delete;
    StaggeredStokesIBJacobianOperator(const StaggeredStokesIBJacobianOperator& from) = delete;
    StaggeredStokesIBJacobianOperator& operator=(const StaggeredStokesIBJacobianOperator& that) = delete;

    StaggeredStokesIBOperator::Context d_ctx;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_base_vector = nullptr;
    Mat d_SAJ_mat = nullptr;
    Vec d_input_vec = nullptr;
    Vec d_output_vec = nullptr;
    Vec d_solver_X = nullptr;
    Vec d_solver_X0 = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBJacobianOperator
