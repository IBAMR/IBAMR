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
 * \brief FAC preconditioner wrapper for velocity-path Stokes-IB Jacobian
 * solves.
 */
class StaggeredStokesIBJacobianFACPreconditioner : public IBTK::FACPreconditioner, public StaggeredStokesSolver
{
public:
    /*!
     * \brief Constructor.
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
    void setComponentsHaveNullSpace(const bool has_velocity_nullspace, const bool has_pressure_nullspace) override;

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
     */
    void setIBTimeSteppingType(TimeSteppingType time_stepping_type);

    /*!
     * \brief Set the Lagrangian force Jacobian matrix.
     */
    void setIBForceJacobian(Mat& A_mat);

    /*!
     * \brief Set the Lagrangian-Eulerian interpolation matrix.
     */
    void setIBInterpOp(Mat& J_mat);

    /*!
     * \brief Set the IB strategy used by this preconditioner.
     */
    void setIBImplicitStrategy(SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_implicit_ops);

    /*!
     * \brief Initialize hierarchy-dependent solver state.
     */
    void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Get the concrete IB FAC strategy object.
     */
    SAMRAI::tbox::Pointer<StaggeredStokesIBLevelRelaxationFACOperator> getIBFACPreconditionerStrategy() const;

private:
    StaggeredStokesIBJacobianFACPreconditioner() = delete;
    StaggeredStokesIBJacobianFACPreconditioner(const StaggeredStokesIBJacobianFACPreconditioner& from) = delete;
    StaggeredStokesIBJacobianFACPreconditioner&
    operator=(const StaggeredStokesIBJacobianFACPreconditioner& that) = delete;

    SAMRAI::tbox::Pointer<IBImplicitStrategy> d_ib_implicit_ops = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBJacobianFACPreconditioner
