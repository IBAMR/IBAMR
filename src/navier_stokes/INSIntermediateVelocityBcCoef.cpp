// Filename: INSIntermediateVelocityBcCoef.cpp
// Created on 24 Jul 2008 by Boyce Griffith
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

#include <stddef.h>
#include <limits>
#include <ostream>
#include <vector>

#include "ArrayData.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class Variable;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSIntermediateVelocityBcCoef::INSIntermediateVelocityBcCoef(const int comp_idx,
                                                             const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                             const bool homogeneous_bc)
    : d_comp_idx(comp_idx), d_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL))
{
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // INSIntermediateVelocityBcCoef

INSIntermediateVelocityBcCoef::~INSIntermediateVelocityBcCoef()
{
    // intentionally blank
    return;
} // ~INSIntermediateVelocityBcCoef

void
INSIntermediateVelocityBcCoef::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
} // setPhysicalBcCoefs

void
INSIntermediateVelocityBcCoef::setSolutionTime(double /*solution_time*/)
{
    // intentionally blank
    return;
} // setSolutionTime

void
INSIntermediateVelocityBcCoef::setTimeInterval(double /*current_time*/, double /*new_time*/)
{
    // intentionally blank
    return;
} // setTimeInterval

void
INSIntermediateVelocityBcCoef::setTargetPatchDataIndex(int target_idx)
{
    ExtendedRobinBcCoefStrategy::setTargetPatchDataIndex(target_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ExtendedRobinBcCoefStrategy* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPatchDataIndex(target_idx);
    }
    return;
} // setTargetPatchDataIndex

void
INSIntermediateVelocityBcCoef::clearTargetPatchDataIndex()
{
    ExtendedRobinBcCoefStrategy::clearTargetPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ExtendedRobinBcCoefStrategy* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // clearTargetPatchDataIndex

void
INSIntermediateVelocityBcCoef::setHomogeneousBc(bool homogeneous_bc)
{
    ExtendedRobinBcCoefStrategy::setHomogeneousBc(homogeneous_bc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ExtendedRobinBcCoefStrategy* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    return;
} // setHomogeneousBc

void
INSIntermediateVelocityBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                          Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                          Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                          const Pointer<Variable<NDIM> >& variable,
                                          const Patch<NDIM>& patch,
                                          const BoundaryBox<NDIM>& bdry_box,
                                          double fill_time) const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
#endif
    // Set the unmodified velocity bc coefs.
    d_bc_coefs[d_comp_idx]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
    if (d_homogeneous_bc && gcoef_data) gcoef_data->fillAll(0.0);
    return;
} // setBcCoefs

IntVector<NDIM>
INSIntermediateVelocityBcCoef::numberOfExtensionsFillable() const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
#endif
    IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ret_val = IntVector<NDIM>::min(ret_val, d_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
