// Filename: AdvDiffPhysicalBoundaryUtilities.C
// Created on 24 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

#include "AdvDiffPhysicalBoundaryUtilities.h"

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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

void
AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(
    Pointer<CellData<NDIM,double> > Q_data,
    Pointer<FaceData<NDIM,double> > u_ADV_data,
    Pointer<Patch<NDIM> > patch,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    const double fill_time,
    const bool inflow_boundaries_only,
    const bool homogeneous_bc)
{
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    if (!pgeom->getTouchesRegularBoundary()) return;
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    if (physical_codim1_boxes.size() == 0) return;

    // Loop over the boundary fill boxes and set boundary conditions.
    const Box<NDIM>& patch_box = patch->getBox();
    const double* const dx = pgeom->getDx();

    // Setup any extended Robin BC coef objects.
    for (int depth = 0; depth < Q_data->getDepth(); ++depth)
    {
        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[depth]);
        if (extended_bc_coef)
        {
            extended_bc_coef->clearTargetPatchDataIndex();
            extended_bc_coef->setHomogeneousBc(homogeneous_bc);
        }
    }

    // Set the boundary conditions.
    const IntVector<NDIM>& gcw = Q_data->getGhostCellWidth();
    for (int n = 0; n < physical_codim1_boxes.size(); ++n)
    {
        const BoundaryBox<NDIM>& bdry_box   = physical_codim1_boxes[n];
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index/2;
        const bool is_lower                 = location_index%2 == 0;
        static const IntVector<NDIM> gcw_to_fill = 1;
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
        Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d != bdry_normal_axis)
            {
                bc_coef_box.lower(d) = std::max(bc_coef_box.lower(d), patch_box.lower(d));
                bc_coef_box.upper(d) = std::min(bc_coef_box.upper(d), patch_box.upper(d));
            }
        }
        Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
        for (int depth = 0; depth < Q_data->getDepth(); ++depth)
        {
            bc_coefs[depth]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, fill_time);
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[depth]);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i = bc();
                const FaceIndex<NDIM> i_f(i, bdry_normal_axis, FaceIndex<NDIM>::Lower);
                bool is_inflow_bdry = (is_lower && (*u_ADV_data)(i_f) > 0.0) || (!is_lower && (*u_ADV_data)(i_f) < 0.0);
                if (!inflow_boundaries_only || is_inflow_bdry)
                {
                    const double& a = (*acoef_data)(i,0);
                    const double& b = (*bcoef_data)(i,0);
                    const double& g = (*gcoef_data)(i,0);
                    const double& h = dx[bdry_normal_axis];
                    int sgn;
                    Index<NDIM> i_intr(i);
                    if (is_lower)
                    {
                        sgn = -1;
                    }
                    else
                    {
                        sgn = +1;
                        i_intr(bdry_normal_axis) -= 1;
                    }
                    Index<NDIM> i_true(i_intr), i_ghost(i_intr);
                    for (int k = 1; k <= gcw(bdry_normal_axis); ++k)
                    {
                        i_ghost(bdry_normal_axis) = i_intr(bdry_normal_axis) + sgn*k;
                        i_true (bdry_normal_axis) = i_intr(bdry_normal_axis) - sgn*(k-1);
                        const double Q_i = (*Q_data)(i_true,depth);
                        (*Q_data)(i_ghost,depth) = -(-4.0*g*h*k-2.0*g*h+2.0*a*Q_i*h*k+a*Q_i*h-2.0*b*Q_i)/(2.0*a*h*k+a*h+2.0*b);
                    }
                }
            }
        }
    }
    return;
}// setPhysicalBoundaryConditions

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
