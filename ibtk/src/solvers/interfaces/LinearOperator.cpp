// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibtk/LinearOperator.h"

#include "Box.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LinearOperator::LinearOperator(std::string object_name, bool homogeneous_bc)
    : GeneralOperator(std::move(object_name), homogeneous_bc)
{
    // intentionally blank
    return;
} // LinearOperator()

LinearOperator::~LinearOperator()
{
    deallocateOperatorState();
    return;
} // ~LinearOperator()

void
LinearOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
{
    if (d_homogeneous_bc) return;

    // Set y := y - A*0, i.e., shift the right-hand-side vector to account for
    // inhomogeneous boundary conditions.
    Pointer<SAMRAIVectorReal<NDIM, double> > x = y.cloneVector("");
    Pointer<SAMRAIVectorReal<NDIM, double> > b = y.cloneVector("");
    x->allocateVectorData();
    b->allocateVectorData();
    x->setToScalar(0.0);
    apply(*x, *b);
    y.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false), b);
    x->freeVectorComponents();
    b->freeVectorComponents();
    return;
} // modifyRhsForBcs

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
