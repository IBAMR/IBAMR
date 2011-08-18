// Filename: INSStaggeredVelocityBcCoef.C
// Created on 22 Jul 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "INSStaggeredVelocityBcCoef.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <SideData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredVelocityBcCoef::INSStaggeredVelocityBcCoef(
    const unsigned int comp_idx,
    const INSProblemCoefs& problem_coefs,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
    const bool homogeneous_bc)
    : d_comp_idx(comp_idx),
      d_problem_coefs(problem_coefs),
      d_bc_coefs(static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_target_idx(-1),
      d_homogeneous_bc(false)
{
    setPhysicalBoundaryConditions(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSStaggeredVelocityBcCoef

INSStaggeredVelocityBcCoef::~INSStaggeredVelocityBcCoef()
{
    // intentionally blank
    return;
}// ~INSStaggeredVelocityBcCoef

void
INSStaggeredVelocityBcCoef::setPhysicalBoundaryConditions(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
    d_bc_coefs = bc_coefs;
    return;
}// setPhysicalBoundaryConditions

void
INSStaggeredVelocityBcCoef::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

void
INSStaggeredVelocityBcCoef::setTargetPatchDataIndex(
    const int target_idx)
{
    d_target_idx = target_idx;
    return;
}// setTargetPatchDataIndex

void
INSStaggeredVelocityBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredVelocityBcCoef::setBcCoefs(
    Pointer<ArrayData<NDIM,double> >& acoef_data,
    Pointer<ArrayData<NDIM,double> >& bcoef_data,
    Pointer<ArrayData<NDIM,double> >& gcoef_data,
    const Pointer<Variable<NDIM> >& variable,
    const Patch<NDIM>& patch,
    const BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d] != NULL);
    }
#endif
    // Set the unmodified velocity bc coefs.
    d_bc_coefs[d_comp_idx]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // We do not make any further modifications to the values of acoef_data and
    // bcoef_data beyond this point.
    if (gcoef_data.isNull()) return;

    // Ensure homogeneous boundary conditions are enforced.
    if (d_homogeneous_bc) gcoef_data->fillAll(0.0);

    // Modify Neumann boundary conditions to correspond to traction (stress)
    // boundary conditions.
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!acoef_data.isNull());
    TBOX_ASSERT(!bcoef_data.isNull());
    TBOX_ASSERT(!gcoef_data.isNull());
#endif
    const unsigned int location_index   = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index/2;
    const bool is_lower        = location_index%2 == 0;
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif
    Pointer<SideData<NDIM,double> > u_data =
        patch.checkAllocated(d_target_idx)
        ? patch.getPatchData(d_target_idx)
        : Pointer<PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_data.isNull());
    TBOX_ASSERT(u_data->getGhostCellWidth().max() == u_data->getGhostCellWidth().min());
#endif
    const Box<NDIM>& ghost_box = u_data->getGhostBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double mu = d_problem_coefs.getMu();
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const Index<NDIM>& i = it();
        const double& alpha = (*acoef_data)(i,0);
        const double& beta  = (*bcoef_data)(i,0);
        double& gamma = (*gcoef_data)(i,0);

        const bool velocity_bc = MathUtilities<double>::equalEps(alpha,1.0);
        const bool traction_bc = MathUtilities<double>::equalEps(beta ,1.0);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
        if (velocity_bc)
        {
            // intentionally blank
        }
        else if (traction_bc)
        {
            if (d_comp_idx == bdry_normal_axis)
            {
                Index<NDIM> i_intr0 = i;
                Index<NDIM> i_intr1 = i;

                if (is_lower)
                {
                    i_intr0(bdry_normal_axis) += 0;
                    i_intr1(bdry_normal_axis) += 1;
                }
                else
                {
                    i_intr0(bdry_normal_axis) -= 1;
                    i_intr1(bdry_normal_axis) -= 2;
                }

                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (d != bdry_normal_axis)
                    {
                        i_intr0(d) = std::max(i_intr0(d),ghost_box.lower()(d));
                        i_intr0(d) = std::min(i_intr0(d),ghost_box.upper()(d));

                        i_intr1(d) = std::max(i_intr1(d),ghost_box.lower()(d));
                        i_intr1(d) = std::min(i_intr1(d),ghost_box.upper()(d));
                    }
                }

                // Specify a Neumann boundary condition which corresponds to a
                // finite difference approximation to the divergence free
                // condition at the boundary of the domain using extrapolated
                // values of the tangential velocities.
                double du_norm_dx_norm = 0.0;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    if (axis != bdry_normal_axis)
                    {
                        const SideIndex<NDIM> i_s_intr0_upper(i_intr0, axis, SideIndex<NDIM>::Upper);
                        const SideIndex<NDIM> i_s_intr1_upper(i_intr1, axis, SideIndex<NDIM>::Upper);
                        const double u_tan_upper = 1.5*(*u_data)(i_s_intr0_upper)-0.5*(*u_data)(i_s_intr1_upper);

                        const SideIndex<NDIM> i_s_intr0_lower(i_intr0, axis, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> i_s_intr1_lower(i_intr1, axis, SideIndex<NDIM>::Lower);
                        const double u_tan_lower = 1.5*(*u_data)(i_s_intr0_lower)-0.5*(*u_data)(i_s_intr1_lower);

                        du_norm_dx_norm -= (u_tan_upper-u_tan_lower)/dx[axis];
                    }
                }
                gamma = (is_lower ? -1.0 : +1.0)*du_norm_dx_norm;
            }
            else
            {
                // Compute the tangential derivative of the normal component of
                // the velocity at the boundary.
                Index<NDIM> i_lower(i), i_upper(i);
                i_lower(d_comp_idx) = std::max(ghost_box.lower()(d_comp_idx),i(d_comp_idx)-1);
                i_upper(d_comp_idx) = std::min(ghost_box.upper()(d_comp_idx),i(d_comp_idx)  );
                const SideIndex<NDIM> i_s_lower(i_lower, bdry_normal_axis, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> i_s_upper(i_upper, bdry_normal_axis, SideIndex<NDIM>::Lower);
                const double du_norm_dx_tan = ((*u_data)(i_s_upper)-(*u_data)(i_s_lower))/dx[d_comp_idx];

                // Correct the boundary condition value.
                gamma = (is_lower ? -1.0 : +1.0)*(gamma/mu - du_norm_dx_tan);
            }
        }
        else
        {
            TBOX_ERROR("this statement should not be reached!\n");
        }
    }
    return;
}// setBcCoefs

IntVector<NDIM>
INSStaggeredVelocityBcCoef::numberOfExtensionsFillable() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d] != NULL);
    }
#endif
    IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ret_val = IntVector<NDIM>::min(ret_val, d_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
