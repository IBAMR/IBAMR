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

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{

class BoundaryBox;

class Variable;

class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSIntermediateVelocityBcCoef::INSIntermediateVelocityBcCoef(
    const int comp_idx,
    const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& bc_coefs,
    const bool homogeneous_bc)
    : d_comp_idx(comp_idx), d_bc_coefs(NDIM)
{
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}

INSIntermediateVelocityBcCoef::~INSIntermediateVelocityBcCoef()
{
    // intentionally blank
    return;
}

void INSIntermediateVelocityBcCoef::setPhysicalBcCoefs(
    const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& bc_coefs)
{
    TBOX_ASSERT(d_bc_coefs.size() == NDIM);
    d_bc_coefs = bc_coefs;
    return;
}

void INSIntermediateVelocityBcCoef::setSolutionTime(double /*solution_time*/)
{
    // intentionally blank
    return;
}

void INSIntermediateVelocityBcCoef::setTimeInterval(double /*current_time*/, double /*new_time*/)
{
    // intentionally blank
    return;
}

void INSIntermediateVelocityBcCoef::setTargetPatchDataIndex(int target_idx)
{
    ExtendedRobinBcCoefStrategy::setTargetPatchDataIndex(target_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = boost::dynamic_pointer_cast<ExtendedRobinBcCoefStrategy>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPatchDataIndex(target_idx);
    }
    return;
}

void INSIntermediateVelocityBcCoef::clearTargetPatchDataIndex()
{
    ExtendedRobinBcCoefStrategy::clearTargetPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = boost::dynamic_pointer_cast<ExtendedRobinBcCoefStrategy>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPatchDataIndex();
    }
    return;
}

void INSIntermediateVelocityBcCoef::setHomogeneousBc(bool homogeneous_bc)
{
    ExtendedRobinBcCoefStrategy::setHomogeneousBc(homogeneous_bc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = boost::dynamic_pointer_cast<ExtendedRobinBcCoefStrategy>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    return;
}

void INSIntermediateVelocityBcCoef::setBcCoefs(const boost::shared_ptr<ArrayData<double>>& acoef_data,
                                               const boost::shared_ptr<ArrayData<double>>& bcoef_data,
                                               const boost::shared_ptr<ArrayData<double>>& gcoef_data,
                                               const boost::shared_ptr<Variable>& variable,
                                               const Patch& patch,
                                               const BoundaryBox& bdry_box,
                                               double fill_time) const
{
    // Set the unmodified velocity bc coefs.
    TBOX_ASSERT(d_bc_coefs[d_comp_idx]);
    d_bc_coefs[d_comp_idx]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
    if (d_homogeneous_bc && gcoef_data) gcoef_data->fillAll(0.0);
    return;
}

IntVector INSIntermediateVelocityBcCoef::numberOfExtensionsFillable() const
{
    IntVector ret_val(DIM, std::numeric_limits<int>::max());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
        ret_val = IntVector::min(ret_val, d_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
