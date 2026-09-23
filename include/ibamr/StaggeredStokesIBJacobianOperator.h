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
 * differentiate their dependence on moving coupling positions. It uses only
 * the position derivative of the force; velocity-dependent terms such as
 * target-point damping are not differentiated. A supplied
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
     * The assembled velocity contribution must include its sign and time/grid factors.
     * It must use the full coupled global numbering, communicator, and local
     * ownership of the single patch level and DOF fields in the Context (see
     * StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices()).
     * Only velocity entries may be nonzero. Multilevel supplied-matrix application
     * is unsupported. Passing nullptr selects the strategy action instead, which
     * requires a formJacobian() with the strategy action since the state was
     * initialized; formJacobian() with a supplied matrix does not count. formJacobian()
     * does not update this matrix.
     *
     * The operator retains a PETSc reference until the matrix is replaced,
     * cleared with nullptr, or the operator is destroyed. Deallocating or
     * reinitializing the operator state does not release it, so replace the
     * matrix whenever the hierarchy or DOF numbering changes. Supplied-matrix-only
     * use requires just the Stokes operator and DOF fields in Context. With a
     * Context that also supports the strategy action, the matrix may be installed
     * or cleared while the operator is initialized.
     */
    void setIBCouplingJacobian(Mat SAJ_mat);

    /*!
     * \brief Form and cache Jacobian state at the specified point.
     *
     * Requires initialized state and copies x into the base vector. With the
     * strategy action, it also requires the prepared strategy data described in
     * StaggeredStokesIBOperator::Context, and updates the strategy's linearization
     * point (not the fixed coupling configuration, which this operator never
     * changes) to x; it must be called again when the base state, time interval, or
     * force data change. With a supplied coupling matrix, only the base vector is
     * recorded and the strategy is not used.
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
     * The strategy action uses and updates shared IB strategy state; see
     * StaggeredStokesIBOperator::Context.
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Initialize hierarchy-dependent operator state.
     *
     * When an IB strategy is supplied, this operator applies the linearization of the
     * residual with the interpolation and spreading operators held fixed at the coupling
     * configuration that the strategy holds; it does not include their dependence on the
     * structure position and does not change the configuration. It is the exact Jacobian
     * of StaggeredStokesIBOperator when that operator also uses fixed coupling operators
     * (see IBStrategy::setUseFixedLEOperators()) and the force depends only on the
     * positions. Otherwise it is an approximation that can serve, for example, in a
     * preconditioner, or next to a matrix-free finite-difference Jacobian that the
     * solver uses. The owner calls IBStrategy::updateFixedLEOperators() when the
     * configuration changes and then calls formJacobian() again; this operator caches
     * no coupling data. Reinitialization deallocates the previous state, including the
     * Jacobian base, but retains the supplied coupling matrix; see formJacobian() and
     * setIBCouplingJacobian().
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Deallocate hierarchy-dependent operator state.
     *
     * Releases the cached base and work vectors. The supplied coupling matrix is retained.
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
    /*! \brief Check Context dependencies and indices for the selected action. */
    void validateContext(bool supplied_matrix) const;

    /*! \brief Default construction is disabled. */
    StaggeredStokesIBJacobianOperator() = delete;
    /*! \brief Copy construction is disabled. */
    StaggeredStokesIBJacobianOperator(const StaggeredStokesIBJacobianOperator& from) = delete;
    /*! \brief Copy assignment is disabled. */
    StaggeredStokesIBJacobianOperator& operator=(const StaggeredStokesIBJacobianOperator& that) = delete;

    StaggeredStokesIBOperator::Context d_ctx;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_base_vector = nullptr;
    Mat d_SAJ_mat = nullptr;
    Vec d_input_vec = nullptr;
    Vec d_output_vec = nullptr;
    Vec d_solver_X = nullptr;
    Vec d_solver_X0 = nullptr;
    bool d_strategy_linearization_formed = false;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBJacobianOperator
