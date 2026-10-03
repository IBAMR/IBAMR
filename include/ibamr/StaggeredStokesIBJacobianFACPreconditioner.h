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

#include <ibamr/StaggeredStokesFACPreconditioner.h>
#include <ibamr/ibamr_enums.h>

#include <tbox/Pointer.h>

#include <petscmat.h>

#include <string>

namespace IBAMR
{
class StaggeredStokesIBLevelRelaxationFACOperator;
} // namespace IBAMR
namespace SAMRAI
{
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
 * The supplied strategy must be nonnull. Subclasses replacing the protected strategy
 * must preserve its StaggeredStokesIBLevelRelaxationFACOperator type.
 * See that class for configuration and matrix requirements. The coupling matrices
 * describe one configuration of the structure, that is, interpolation and spreading
 * held fixed; the owner supplies new matrices and reinitializes the preconditioner
 * whenever the configuration changes. The coupling then is the Jacobian's when the
 * residual's coupling operators are held fixed the same way and the force depends
 * only on the positions, and otherwise an approximation. Either way the
 * preconditioner can serve beside a matrix-free finite-difference Jacobian.
 */
class StaggeredStokesIBJacobianFACPreconditioner : public StaggeredStokesFACPreconditioner
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesIBJacobianFACPreconditioner(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesIBJacobianFACPreconditioner() override = default;

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
    void setIBForceJacobian(Mat A_mat);

    /*!
     * \brief Set the Lagrangian-Eulerian interpolation matrix.
     *
     * \see StaggeredStokesIBLevelRelaxationFACOperator::setIBInterpOp
     */
    void setIBInterpOp(Mat J_mat);

    /*!
     * \brief Return the installed strategy as a StaggeredStokesIBLevelRelaxationFACOperator.
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
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBJacobianFACPreconditioner
