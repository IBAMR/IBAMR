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
    const StokesSpecifications* problem_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    const bool homogeneous_bc)
    : d_comp_idx(comp_idx),
      d_problem_coefs(NULL),
      d_bc_coefs(NDIM,static_cast<RobinBcCoefStrategy<NDIM>*>(NULL))
{
    setStokesSpecifications(problem_coefs);
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSStaggeredVelocityBcCoef

INSStaggeredVelocityBcCoef::~INSStaggeredVelocityBcCoef()
{
    // intentionally blank
    return;
}// ~INSStaggeredVelocityBcCoef

void
INSStaggeredVelocityBcCoef::setStokesSpecifications(
    const StokesSpecifications* problem_coefs)
{
    d_problem_coefs = problem_coefs;
    return;
}// setStokesSpecifications

void
INSStaggeredVelocityBcCoef::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
}// setPhysicalBcCoefs

void
INSStaggeredVelocityBcCoef::setSolutionTime(
    const double /*solution_time*/)
{
    // intentionally blank
    return;
}// setSolutionTime

void
INSStaggeredVelocityBcCoef::setTimeInterval(
    const double /*current_time*/,
    const double /*new_time*/)
{
    // intentionally blank
    return;
}// setTimeInterval

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
    if (!gcoef_data) return;

    // Ensure homogeneous boundary conditions are enforced.
    if (d_homogeneous_bc) gcoef_data->fillAll(0.0);

    // Where appropriate, update Neumann boundary conditions to correspond to
    // tangential traction (stress) boundary conditions.
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(acoef_data);
    TBOX_ASSERT(bcoef_data);
    TBOX_ASSERT(gcoef_data);
#endif
    const unsigned int location_index   = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index/2;
    const bool is_lower = location_index%2 == 0;
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif
    Pointer<SideData<NDIM,double> > u_data = patch.getPatchData(d_target_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_data);
    TBOX_ASSERT(u_data->getGhostCellWidth().max() == u_data->getGhostCellWidth().min());
#endif
    const Box<NDIM>& ghost_box = u_data->getGhostBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double mu = d_problem_coefs->getMu();
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
            // intentionally blank.
        }
        else if (traction_bc)
        {
            if (d_comp_idx == bdry_normal_axis)
            {
                // intentionally blank; "exact" divergence-free conditions are
                // imposed after filling all other ghost cell values.
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

//////////////////////////////////////////////////////////////////////////////
