// Filename: ConvectiveOperator.cpp
// Created on 21 Aug 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/ConvectiveOperator.h"
#include "SAMRAIVectorReal.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ConvectiveOperator::ConvectiveOperator(const std::string& object_name, const ConvectiveDifferencingType difference_form)
    : GeneralOperator(object_name), d_difference_form(difference_form), d_u_idx(-1)
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
