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

#ifndef included_IBAMR_StaggeredStokesIBOperator
#define included_IBAMR_StaggeredStokesIBOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/IBImplicitStrategy.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/GeneralOperator.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <IntVector.h>

#include <string>
#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK
namespace SAMRAI
{
namespace math
{
template <int DIM, class TYPE>
class HierarchyDataOpsReal;
} // namespace math
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
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
 * \brief Nonlinear Stokes-IB operator on Eulerian velocity-pressure vectors.
 *
 * In this velocity-pressure formulation, the nonlinear solver unknowns are
 * Eulerian velocity and pressure. Lagrangian positions are determined from the
 * velocity through the time-stepping equations, not independent solver unknowns.
 * apply() advances the IB strategy with the input velocity and selected
 * time-stepping rule, then subtracts the resulting spread force from the Stokes
 * momentum action, with weight 1 for
 * backward Euler/midpoint and 1/2 for trapezoidal stepping. The Stokes divergence
 * component is unchanged. Known current-time terms and the complete time-step
 * right-hand side are the caller's responsibility.
 */
class StaggeredStokesIBOperator : public IBTK::GeneralOperator
{
public:
    /*!
     * \brief Shared coupling state for the nonlinear and Jacobian operators.
     *
     * Operator vectors have side-centered velocity as component 0 and
     * cell-centered pressure as component 1. The vectors, data operations, and
     * schedules below must refer to the same hierarchy.
     *
     * Configure stokes_op's coefficients and boundary objects, and initialize
     * its state, before initializing either operator; the owner of the Context
     * deallocates it after both operators. The operators only apply stokes_op,
     * and set its times and boundary mode. Prepare ib_implicit_ops for the active hierarchy and time
     * step. The operators use the coupling configuration that the strategy holds and never change it:
     * with fixed coupling operators enabled (IBStrategy::setUseFixedLEOperators()), the interpolation
     * and spreading operators stay at the configuration of the last IBStrategy::updateFixedLEOperators(),
     * which the owner calls when the configuration changes, and otherwise they follow the current
     * positions.
     * Trapezoidal stepping also requires the current Lagrangian velocity stored
     * by the strategy; u_current_idx instead identifies current Eulerian velocity.
     *
     * Allocate side-centered velocity/force scratch data at u_idx/f_idx and
     * current velocity at u_current_idx with the ghost widths required by the
     * strategy. Supply the corresponding velocity synchronization/ghost-fill
     * and force-prolongation schedules. The optional u_phys_bdry_op is borrowed
     * and must remain valid during use; omit it only if the boundary setup needs
     * no such strategy, e.g. a periodic domain.
     *
     * Supplied-matrix-only Jacobian use requires just stokes_op and coupled velocity/pressure
     * DOF fields u_dof_index_idx/p_dof_index_idx on the input vector's single level; this
     * supports a Jacobian operator that takes its coupling matrix directly from its owner,
     * introduced in the next PR in this stack.
     * Copying Context copies handles and indices, not the shared objects or
     * patch data. Keep that shared state valid throughout operator use.
     */
    struct Context
    {
        SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_implicit_ops = nullptr;
        SAMRAI::tbox::Pointer<StaggeredStokesOperator> stokes_op = nullptr;
        IBTK::RobinPhysBdryPatchStrategy* u_phys_bdry_op = nullptr;
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double>> hier_velocity_data_ops = nullptr;
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM>>> u_synch_scheds;
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> u_ghost_fill_scheds;
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> f_prolongation_scheds;
        int u_idx = IBTK::invalid_index;
        int f_idx = IBTK::invalid_index;
        int u_current_idx = IBTK::invalid_index;
        int u_dof_index_idx = IBTK::invalid_index;
        int p_dof_index_idx = IBTK::invalid_index;
        //! Supported rules are BACKWARD_EULER, TRAPEZOIDAL_RULE, and MIDPOINT_RULE.
        TimeSteppingType time_stepping_type = MIDPOINT_RULE;
    };

    /*!
     * \brief Constructor.
     */
    explicit StaggeredStokesIBOperator(const std::string& object_name, bool homogeneous_bc = false);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesIBOperator() override;

    /*!
     * \brief Set context data required by this operator.
     *
     * Copies ctx; later edits to the caller's Context do not update this copy.
     * Deallocate before replacing dependencies or hierarchy-dependent data,
     * then initialize with vectors matching the new Context.
     */
    void setOperatorContext(const Context& ctx);

    /*!
     * \brief Compute \f$y = A[x]\f$.
     *
     * Requires initialized state and the prepared Context. Updates the shared
     * IB strategy's velocity, positions, and force. The Stokes action uses
     * homogeneous boundary conditions; velocity interpolation uses physical
     * boundary data, independently of this wrapper's homogeneous-boundary flag.
     */
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    /*!
     * \brief Initialize hierarchy-dependent operator state.
     *
     * The operator uses the coupling configuration that the strategy holds and does not change it.
     * With fixed coupling operators enabled, the owner decides when to call
     * IBStrategy::updateFixedLEOperators() to freeze the interpolation and spreading positions for a
     * time step or stage.
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
    /*! \brief Default construction is disabled. */
    StaggeredStokesIBOperator() = delete;
    /*! \brief Copy construction is disabled. */
    StaggeredStokesIBOperator(const StaggeredStokesIBOperator& from) = delete;
    /*! \brief Copy assignment is disabled. */
    StaggeredStokesIBOperator& operator=(const StaggeredStokesIBOperator& that) = delete;

    Context d_ctx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBOperator
