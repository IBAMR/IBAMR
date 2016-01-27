// Filename: INSProjectionBcCoef.cpp
// Created on 23 Jul 2008 by Boyce Griffith
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
#include "BoundaryBox.h"
#include "Box.h"
#include "Index.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSProjectionBcCoef::INSProjectionBcCoef(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                         const bool homogeneous_bc)
    : d_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_solution_time(std::numeric_limits<double>::quiet_NaN())
{
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // INSProjectionBcCoef

INSProjectionBcCoef::~INSProjectionBcCoef()
{
    // intentionally blank
    return;
} // ~INSProjectionBcCoef

void
INSProjectionBcCoef::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
} // setPhysicalBcCoefs

void
INSProjectionBcCoef::setSolutionTime(double solution_time)
{
    d_solution_time = solution_time;
    return;
} // setSolutionTime

void
INSProjectionBcCoef::setTimeInterval(double /*current_time*/, double /*new_time*/)
{
    // intentionally blank
    return;
} // setTimeInterval

void
INSProjectionBcCoef::setTargetPatchDataIndex(int target_idx)
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
INSProjectionBcCoef::clearTargetPatchDataIndex()
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
INSProjectionBcCoef::setHomogeneousBc(bool homogeneous_bc)
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
INSProjectionBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                const Pointer<Variable<NDIM> >& variable,
                                const Patch<NDIM>& patch,
                                const BoundaryBox<NDIM>& bdry_box,
                                double /*fill_time*/) const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
    TBOX_ASSERT(acoef_data);
    TBOX_ASSERT(bcoef_data);
#endif
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(!gcoef_data || (bc_coef_box == gcoef_data->getBox()));
#endif
    // Set the unmodified velocity bc coefs.
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    d_bc_coefs[bdry_normal_axis]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, d_solution_time);

    // Ensure homogeneous boundary conditions are enforced.
    if (d_homogeneous_bc && gcoef_data) gcoef_data->fillAll(0.0);

    // Update the boundary condition coefficients.  Specifically, normal
    // velocity boundary conditions are converted into Neumann conditions for
    // the pressure, and normal traction boundary conditions are converted into
    // Dirichlet conditions for the pressure.
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const Index<NDIM>& i = it();
        double dummy_val;
        double& alpha = acoef_data ? (*acoef_data)(i, 0) : dummy_val;
        double& beta = bcoef_data ? (*bcoef_data)(i, 0) : dummy_val;
        double& gamma = gcoef_data ? (*gcoef_data)(i, 0) : dummy_val;
        const bool velocity_bc = MathUtilities<double>::equalEps(alpha, 1.0);
        const bool traction_bc = MathUtilities<double>::equalEps(beta, 1.0);
#if !defined(NDEBUG)
        TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
        if (velocity_bc)
        {
            alpha = 0.0;
            beta = 1.0;
            gamma = 0.0;
        }
        else if (traction_bc)
        {
            alpha = 1.0;
            beta = 0.0;
            gamma = -gamma;
        }
        else
        {
            TBOX_ERROR("this statement should not be reached!\n");
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
INSProjectionBcCoef::numberOfExtensionsFillable() const
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
