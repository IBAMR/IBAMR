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

#include "ibamr/ConvectiveOperator.h"

#include "SAMRAIVectorReal.h"

#include <utility>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ConvectiveOperator::ConvectiveOperator(std::string object_name, const ConvectiveDifferencingType difference_form)
    : GeneralOperator(std::move(object_name)), d_difference_form(difference_form)
{
    // intentionally blank
    return;
} // ConvectiveOperator

ConvectiveOperator::~ConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~ConvectiveOperator

void
ConvectiveOperator::setAdvectionVelocity(const int u_idx)
{
    d_u_idx = u_idx;
    return;
} // setAdvectionVelocity

int
ConvectiveOperator::getAdvectionVelocity() const
{
    return d_u_idx;
} // getAdvectionVelocity

void
ConvectiveOperator::setConvectiveDifferencingType(const ConvectiveDifferencingType difference_form)
{
    d_difference_form = difference_form;
    return;
} // setConvectiveDifferencingType

ConvectiveDifferencingType
ConvectiveOperator::getConvectiveDifferencingType() const
{
    return d_difference_form;
} // getConvectiveDifferencingType

void
ConvectiveOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    // Get the vector components.
    const int Q_idx = x.getComponentDescriptorIndex(0);
    const int N_idx = y.getComponentDescriptorIndex(0);

    // Compute the action of the operator.
    applyConvectiveOperator(Q_idx, N_idx);
    return;
} // apply

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
