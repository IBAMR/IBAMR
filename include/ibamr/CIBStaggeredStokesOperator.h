// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_CIBStaggeredStokesOperator
#define included_IBAMR_CIBStaggeredStokesOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/CIBStrategy.h>
#include <ibamr/StaggeredStokesOperator.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/LinearOperator.h>

#include <tbox/Pointer.h>

#include <petscvec.h>

#include <string>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

namespace IBAMR
{
class CIBStrategy;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CIBStaggeredStokesOperator is a concrete IBTK::LinearOperator which
 * implements a staggered-grid (MAC) discretization of the incompressible Stokes
 * operator while maintaining the constraint of rigidity for the immersed structures.
 *
 * This class is intended to be used with an iterative (Krylov or Newton-Krylov)
 * incompressible flow solver.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class CIBStaggeredStokesOperator : public IBAMR::StaggeredStokesOperator
{
    //////////////////////////////////////////////////////////////////////////////
public:
    /*!
     * \brief Class constructor.
     */
    CIBStaggeredStokesOperator(std::string object_name,
                               SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy,
                               bool homogeneous_bc = true,
                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr);

    /*!
     * \brief Destructor.
     */
    ~CIBStaggeredStokesOperator();

    /*!
     * \brief Compute hierarchy dependent data required for computing y=Ax, including the ghost cell fill of the
     * velocity for interpolation.
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by initializeOperatorState().
     */
    void deallocateOperatorState() override;

    //\{ // Additional functionality of CIBStaggeredStokesOperator.

    /*!
     * \name Linear operator functionality.
     */
    using LinearOperator::apply;
    virtual void apply(Vec x, Vec y);

    // Set scaling factors for various operators to improve the condition number
    // of the system.

    /*!
     * \brief Set scale factor for interp operator.
     */
    void setInterpScaleFactor(const double beta);

    /*!
     * \brief Set scale factor for spread operator.
     */
    void setSpreadScaleFactor(const double gamma);

    /*!
     * \brief Set scale factor for regularizing mobility matrix.
     */
    void setRegularizeMobilityFactor(const double delta);

    /*!
     * \brief Set if the mean of the Lagrangian force is to be subtracted
     * from the Eulerian force variable.
     *
     * \note This operation is needed for certain situations like Stokes flow
     * with periodic BCs.
     */
    void setNormalizeSpreadForce(const bool normalize_force);

    /*
     * Set y := y - A*0, i.e., shift the right-hand-side vector to account for
     * inhomogeneous boundary conditions.
     */
    using LinearOperator::modifyRhsForBcs;
    virtual void modifyRhsForBcs(Vec y);

    /*!
     * \brief Impose boudary conditions in the solution vector.
     */
    using LinearOperator::imposeSolBcs;
    virtual void imposeSolBcs(Vec x);
    //\}

    //////////////////////////////////////////////////////////////////////////////
protected:
    // Pointer to a constraint based rigid IB Method.
    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;

    // Scaling factors for various operators.
    double d_scale_interp = 1.0, d_scale_spread = 1.0, d_reg_mob_factor = 1.0;
    bool d_normalize_spread_force = false;
    //////////////////////////////////////////////////////////////////////////////

private:
    /*!
     * \brief Return the transaction component of the ghost cell fill of the velocity \a U_idx for interpolation.
     */
    IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent
    getVelocityInterpComponent(int U_idx) const;

    /*!
     * \brief Return the transaction component of the ghost cell fill of the velocity \a U_idx with no Robin
     * coefficients, for the interpolation in modifyRhsForBcs().
     */
    IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent
    getVelocityRhsInterpComponent(int U_idx) const;

    /*!
     * \brief Fill the same-level and physical-boundary ghost cells of the side-centered velocity, for interpolation.
     *
     * \note The coarse-fine ghost cells are not set, because IB interpolation uses fine-side values only.
     */
    void fillVelocityGhostCells(int U_idx);

    /*!
     * \brief The ghost cell fill of the velocity for interpolation, built in initializeOperatorState().
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_interp_fill;

    /*!
     * \brief The transaction component of d_U_interp_fill from initializeOperatorState(), to which
     * fillVelocityGhostCells() resets the fill after use.
     */
    IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent d_U_interp_comp;

    /*!
     * \brief The ghost cell fill of the velocity with no Robin coefficients, for the interpolation in
     * modifyRhsForBcs(), built in initializeOperatorState().
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_rhs_interp_fill;

    /*!
     * \brief The transaction component of d_U_rhs_interp_fill from initializeOperatorState(), to which
     * modifyRhsForBcs() resets the fill after use.
     */
    IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent d_U_rhs_interp_comp;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CIBStaggeredStokesOperator(const CIBStaggeredStokesOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CIBStaggeredStokesOperator& operator=(const CIBStaggeredStokesOperator& that) = delete;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_CIBStaggeredStokesOperator
