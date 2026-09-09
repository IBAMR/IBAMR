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

#ifndef included_IBAMR_StaggeredStokesIBJacobianFACPreconditioner
#define included_IBAMR_StaggeredStokesIBJacobianFACPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/IBImplicitStrategy.h>
#include <ibamr/StaggeredStokesSolver.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/FACPreconditioner.h>

#include <tbox/Pointer.h>

#include <petscmat.h>

#include <PoissonSpecifications.h>

#include <string>
#include <vector>

namespace IBAMR
{
class StaggeredStokesIBLevelRelaxationFACOperator;
class StaggeredStokesPhysicalBoundaryHelper;
} // namespace IBAMR
namespace IBTK
{
class FACPreconditionerStrategy;
} // namespace IBTK
namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief FAC preconditioner for the Stokes-IB Jacobian.
 *
 * Uses the velocity-pressure formulation described by \ref StaggeredStokesIBOperator.
 * The supplied FAC strategy must be a StaggeredStokesIBLevelRelaxationFACOperator;
 * this wrapper forwards its Stokes coefficients, boundary conditions, nullspaces,
 * and IB configuration to that strategy. Before initializeSolverState(), configure
 * these data, set the time interval and solution time, select an IB time rule,
 * and supply the force Jacobian and interpolation matrices. See
 * StaggeredStokesIBLevelRelaxationFACOperator for matrix and hierarchy requirements.
 */
class StaggeredStokesIBJacobianFACPreconditioner : public IBTK::FACPreconditioner, public StaggeredStokesSolver
{
public:
    /*!
     * \brief Constructor.
     *
     * fac_strategy must be a StaggeredStokesIBLevelRelaxationFACOperator.
     */
    StaggeredStokesIBJacobianFACPreconditioner(const std::string& object_name,
                                               SAMRAI::tbox::Pointer<IBTK::FACPreconditionerStrategy> fac_strategy,
                                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                               const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesIBJacobianFACPreconditioner() override = default;

    /*!
     * \brief Set the velocity block Poisson coefficients.
     */
    void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& U_problem_coefs) override;

    /*!
     * \brief Set whether velocity and pressure each contain a null space.
     */
    void setComponentsHaveNullSpace(bool has_velocity_nullspace, bool has_pressure_nullspace) override;

    /*!
     * \brief Set physical boundary condition coefficient objects.
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef) override;

    /*!
     * \brief Set helper object used for physical-boundary operations.
     */
    void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper) override;

    /*!
     * \brief Set IB time stepping type used by preconditioning operators.
     *
     * Must be called before initialization; see
     * StaggeredStokesIBLevelRelaxationFACOperator::setIBTimeSteppingType().
     */
    void setIBTimeSteppingType(TimeSteppingType time_stepping_type);

    /*!
     * \brief Set the Lagrangian force Jacobian matrix.
     *
     * \see StaggeredStokesIBLevelRelaxationFACOperator::setIBForceJacobian
     */
    void setIBForceJacobian(Mat& A_mat);

    /*!
     * \brief Set the Lagrangian-Eulerian interpolation matrix.
     *
     * \see StaggeredStokesIBLevelRelaxationFACOperator::setIBInterpOp
     */
    void setIBInterpOp(Mat& J_mat);

    /*!
     * \brief Set an optional IB strategy whose fixed coupling is updated at initialization.
     *
     * Initialization enables and updates fixed Lagrangian-Eulerian operators on
     * this shared strategy. It does not construct the matrices supplied through
     * setIBForceJacobian() and setIBInterpOp(). The caller must prepare the
     * strategy's hierarchy and time-step data first; with IBMethod, enable fixed
     * coupling before IBMethod::preprocessIntegrateData(). Passing nullptr omits
     * this update, not the requirement to supply the two matrices.
     */
    void setIBImplicitStrategy(SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_implicit_ops);

    /*!
     * \brief Initialize hierarchy-dependent solver state.
     *
     * Requires configured Stokes/IB data and compatible velocity-pressure
     * vectors. Updates the optional IB strategy, then initializes FAC state.
     * See setIBImplicitStrategy() and StaggeredStokesIBLevelRelaxationFACOperator.
     */
    void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Get the concrete IB FAC strategy object.
     *
     * Returns the shared strategy, not a copy or independent initialized state.
     * Returns nullptr if the constructor's concrete-strategy requirement was
     * not met.
     */
    SAMRAI::tbox::Pointer<StaggeredStokesIBLevelRelaxationFACOperator> getIBFACPreconditionerStrategy() const;

private:
    /*! \brief Default construction is disabled. */
    StaggeredStokesIBJacobianFACPreconditioner() = delete;
    /*! \brief Copy construction is disabled. */
    StaggeredStokesIBJacobianFACPreconditioner(const StaggeredStokesIBJacobianFACPreconditioner& from) = delete;
    /*! \brief Copy assignment is disabled. */
    StaggeredStokesIBJacobianFACPreconditioner&
    operator=(const StaggeredStokesIBJacobianFACPreconditioner& that) = delete;

    SAMRAI::tbox::Pointer<IBImplicitStrategy> d_ib_implicit_ops = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBJacobianFACPreconditioner
