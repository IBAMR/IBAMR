// ---------------------------------------------------------------------
//
// Copyright (c) 2025 - 2025 by the IBAMR developers
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

#include "ibtk/IdentityOperator.h"

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IdentityOperator::IdentityOperator(std::string object_name, bool homogeneous_bc)
    : LinearOperator(std::move(object_name), homogeneous_bc), LinearSolver()
{
    // intentionally blank
    return;
} // IdentityOperator()

void
IdentityOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    y.copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    return;
} // apply

bool
IdentityOperator::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    x.copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));
    return true;
} // solveSystem

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
