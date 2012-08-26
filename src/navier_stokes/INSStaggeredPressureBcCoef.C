// Filename: INSStaggeredPressureBcCoef.C
// Created on 23 Jul 2008 by Boyce Griffith
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

#include "INSStaggeredPressureBcCoef.h"

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

INSStaggeredPressureBcCoef::INSStaggeredPressureBcCoef(
    const StokesSpecifications* problem_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    const bool homogeneous_bc)
    : d_problem_coefs(problem_coefs),
      d_u_current_idx(-1),
      d_u_new_idx(-1),
      d_bc_coefs(NDIM,static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_target_idx(-1),
      d_homogeneous_bc(false)
{
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSStaggeredPressureBcCoef

INSStaggeredPressureBcCoef::~INSStaggeredPressureBcCoef()
{
    // intentionally blank
    return;
}// ~INSStaggeredPressureBcCoef

void
INSStaggeredPressureBcCoef::setVelocityCurrentPatchDataIndex(
    const int u_current_idx)
{
    d_u_current_idx = u_current_idx;
    return;
}// setVelocityCurrentPatchDataIndex

void
INSStaggeredPressureBcCoef::setVelocityNewPatchDataIndex(
    const int u_new_idx)
{
    d_u_new_idx = u_new_idx;
    return;
}// setVelocityNewPatchDataIndex

void
INSStaggeredPressureBcCoef::setStokesSpecifications(
    const StokesSpecifications* problem_coefs)
{
    d_problem_coefs = problem_coefs;
    return;
}// setStokesSpecifications

void
INSStaggeredPressureBcCoef::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
}// setPhysicalBcCoefs

void
INSStaggeredPressureBcCoef::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

void
INSStaggeredPressureBcCoef::setTargetPatchDataIndex(
    const int target_idx)
{
    d_target_idx = target_idx;
    return;
}// setTargetPatchDataIndex

void
INSStaggeredPressureBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredPressureBcCoef::setBcCoefs(
    Pointer<ArrayData<NDIM,double> >& acoef_data,
    Pointer<ArrayData<NDIM,double> >& bcoef_data,
    Pointer<ArrayData<NDIM,double> >& gcoef_data,
    const Pointer<Variable<NDIM> >& variable,
    const Patch<NDIM>& patch,
    const BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
    const double half_time = 0.5*(d_current_time+d_new_time);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d] != NULL);
    }
    TBOX_ASSERT(MathUtilities<double>::equalEps(fill_time,d_new_time) ||
                MathUtilities<double>::equalEps(fill_time,half_time));
    TBOX_ASSERT(!acoef_data.isNull());
    TBOX_ASSERT(!bcoef_data.isNull());
    TBOX_ASSERT(!gcoef_data.isNull());
#else
    NULL_USE(fill_time);
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

    // Set the unmodified velocity bc coefs.
    d_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, half_time);

    // Modify the velocity boundary conditions to correspond to pressure
    // boundary conditions.
    Pointer<SideData<NDIM,double> > u_current_data =
        patch.checkAllocated(d_u_current_idx)
        ? patch.getPatchData(d_u_current_idx)
        : Pointer<PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_current_data.isNull());
    TBOX_ASSERT(u_current_data->getGhostCellWidth().max() == u_current_data->getGhostCellWidth().min());
#endif
    Pointer<SideData<NDIM,double> > u_new_data =
        patch.checkAllocated(d_u_new_idx)
        ? patch.getPatchData(d_u_new_idx)
        : Pointer<PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_new_data.isNull());
    TBOX_ASSERT(u_new_data->getGhostCellWidth().max() == u_new_data->getGhostCellWidth().min());
#endif
    const Box<NDIM> ghost_box = u_current_data->getGhostBox() * u_new_data->getGhostBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double mu = d_problem_coefs->getMu();
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const Index<NDIM>& i = it();
        double& alpha = (*acoef_data)(i,0);
        double& beta  = (*bcoef_data)(i,0);
        double& gamma = (*gcoef_data)(i,0);

        const bool velocity_bc = MathUtilities<double>::equalEps(alpha,1.0);
        const bool traction_bc = MathUtilities<double>::equalEps(beta ,1.0);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
        if (velocity_bc)
        {
            // Set the boundary condition coefficients to correspond to
            // homogeneous Neumann boundary conditions on the pressure.
            alpha = 0.0;
            beta  = 1.0;
            gamma = 0.0;
        }
        else if (traction_bc)
        {
            // Compute (d/dn)u_norm at the boundary by extrapolating the
            // divergence free condition to the boundary.
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

            double du_norm_current_dx_norm = 0.0;
            double du_norm_new_dx_norm     = 0.0;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                if (axis != bdry_normal_axis)
                {
                    const SideIndex<NDIM> i_s_intr0_upper(i_intr0, axis, SideIndex<NDIM>::Upper);
                    const SideIndex<NDIM> i_s_intr1_upper(i_intr1, axis, SideIndex<NDIM>::Upper);
                    const double u_tan_current_upper = 1.5*(*u_current_data)(i_s_intr0_upper)-0.5*(*u_current_data)(i_s_intr1_upper);
                    const double u_tan_new_upper     = 1.5*(*u_new_data    )(i_s_intr0_upper)-0.5*(*u_new_data    )(i_s_intr1_upper);

                    const SideIndex<NDIM> i_s_intr0_lower(i_intr0, axis, SideIndex<NDIM>::Lower);
                    const SideIndex<NDIM> i_s_intr1_lower(i_intr1, axis, SideIndex<NDIM>::Lower);
                    const double u_tan_current_lower = 1.5*(*u_current_data)(i_s_intr0_lower)-0.5*(*u_current_data)(i_s_intr1_lower);
                    const double u_tan_new_lower     = 1.5*(*u_new_data    )(i_s_intr0_lower)-0.5*(*u_new_data    )(i_s_intr1_lower);

                    du_norm_current_dx_norm -= (u_tan_current_upper-u_tan_current_lower)/dx[axis];
                    du_norm_new_dx_norm     -= (u_tan_new_upper    -u_tan_new_lower    )/dx[axis];
                }
            }

            // Set the boundary condition coefficients to correspond to either
            // homogeneous or inhomogeneous Dirichlet boundary conditions on the
            // pressure.
            alpha = 1.0;
            beta  = 0.0;
            gamma = (d_homogeneous_bc ? 0.0 : -gamma + mu*du_norm_current_dx_norm) + mu*du_norm_new_dx_norm;
        }
        else
        {
            TBOX_ERROR("this statement should not be reached!\n");
        }
    }
    return;
}// setBcCoefs

IntVector<NDIM>
INSStaggeredPressureBcCoef::numberOfExtensionsFillable() const
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

//////////////////////////////////////////////////////////////////////////////
