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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
/*! \brief Validate the strategy before the base constructor dereferences it. */
Pointer<StaggeredStokesIBLevelRelaxationFACOperator>
checked_strategy(const std::string& object_name, Pointer<StaggeredStokesIBLevelRelaxationFACOperator> strategy)
{
    if (!strategy)
    {
        TBOX_ERROR(object_name << "::StaggeredStokesIBJacobianFACPreconditioner():\n"
                               << "  fac_strategy must be nonnull.\n");
    }
    return strategy;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesIBJacobianFACPreconditioner::StaggeredStokesIBJacobianFACPreconditioner(
    const std::string& object_name,
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy,
    Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : StaggeredStokesFACPreconditioner(object_name,
                                       checked_strategy(object_name, fac_strategy),
                                       input_db,
                                       default_options_prefix)
{
    return;
} // StaggeredStokesIBJacobianFACPreconditioner

void
StaggeredStokesIBJacobianFACPreconditioner::setIBTimeSteppingType(const TimeSteppingType time_stepping_type)
{
    getIBFACPreconditionerStrategy()->setIBTimeSteppingType(time_stepping_type);
    return;
} // setIBTimeSteppingType

void
StaggeredStokesIBJacobianFACPreconditioner::setIBForceJacobian(Mat A_mat)
{
    getIBFACPreconditionerStrategy()->setIBForceJacobian(A_mat);
    return;
} // setIBForceJacobian

void
StaggeredStokesIBJacobianFACPreconditioner::setIBInterpOp(Mat J_mat)
{
    getIBFACPreconditionerStrategy()->setIBInterpOp(J_mat);
    return;
} // setIBInterpOp

Pointer<StaggeredStokesIBLevelRelaxationFACOperator>
StaggeredStokesIBJacobianFACPreconditioner::getIBFACPreconditionerStrategy() const
{
    return d_fac_strategy;
} // getIBFACPreconditionerStrategy

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
