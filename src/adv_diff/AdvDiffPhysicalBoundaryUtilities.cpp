// Filename: AdvDiffPhysicalBoundaryUtilities.cpp
// Created on 24 Aug 2012 by Boyce Griffith
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
#include <algorithm>
#include <vector>

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h" // IWYU pragma: keep
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/FaceData.h" // IWYU pragma: keep
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/hier/Variable.h"
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

void
AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(Pointer<CellData<double> > Q_data,
                                                                Pointer<FaceData<double> > u_ADV_data,
                                                                Pointer<Patch > patch,
                                                                const std::vector<RobinBcCoefStrategy*>& bc_coefs,
                                                                const double fill_time,
                                                                const bool inflow_boundaries_only,
                                                                const bool homogeneous_bc)
{
    Pointer<CartesianPatchGeometry > pgeom = patch->getPatchGeometry();
    if (!pgeom->getTouchesRegularBoundary()) return;
    const Array<BoundaryBox > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    if (physical_codim1_boxes.size() == 0) return;

    // Loop over the boundary fill boxes and set boundary conditions.
    const Box& patch_box = patch->getBox();
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
    const IntVector& gcw = Q_data->getGhostCellWidth();
    for (int n = 0; n < physical_codim1_boxes.size(); ++n)
    {
        const BoundaryBox& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower = location_index % 2 == 0;
        static const IntVector gcw_to_fill = 1;
        const Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox trimmed_bdry_box(
            bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
        Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d != bdry_normal_axis)
            {
                bc_coef_box.lower(d) = std::max(bc_coef_box.lower(d), patch_box.lower(d));
                bc_coef_box.upper(d) = std::min(bc_coef_box.upper(d), patch_box.upper(d));
            }
        }
        Pointer<ArrayData<double> > acoef_data = new ArrayData<double>(bc_coef_box, 1);
        Pointer<ArrayData<double> > bcoef_data = new ArrayData<double>(bc_coef_box, 1);
        Pointer<ArrayData<double> > gcoef_data = new ArrayData<double>(bc_coef_box, 1);
        for (int depth = 0; depth < Q_data->getDepth(); ++depth)
        {
            bc_coefs[depth]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, fill_time);
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[depth]);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
            for (Box::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index& i = bc();
                const FaceIndex i_f(i, bdry_normal_axis, FaceIndex::Lower);
                bool is_inflow_bdry = (is_lower && (*u_ADV_data)(i_f) > 0.0) || (!is_lower && (*u_ADV_data)(i_f) < 0.0);
                if (!inflow_boundaries_only || is_inflow_bdry)
                {
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const double& g = (*gcoef_data)(i, 0);
                    const double& h = dx[bdry_normal_axis];
                    int sgn;
                    Index i_i(i), i_g(i);
                    if (is_lower)
                    {
                        sgn = -1;
                        i_g(bdry_normal_axis) = patch_box.lower(bdry_normal_axis) - 1;
                        i_i(bdry_normal_axis) = patch_box.lower(bdry_normal_axis);
                    }
                    else
                    {
                        sgn = +1;
                        i_g(bdry_normal_axis) = patch_box.upper(bdry_normal_axis) + 1;
                        i_i(bdry_normal_axis) = patch_box.upper(bdry_normal_axis);
                    }
                    for (int k = 0; k < gcw(bdry_normal_axis); ++k)
                    {
                        const double n = 1.0 + 2.0 * k;
                        const double f_i = -(a * n * h - 2.0 * b) / (a * n * h + 2.0 * b);
                        const double f_g = 2.0 * n * h / (a * n * h + 2.0 * b);
                        const double Q_i = (*Q_data)(i_i, depth);
                        (*Q_data)(i_g, depth) = f_i * Q_i + f_g * g;
                        i_g(bdry_normal_axis) += sgn;
                        i_i(bdry_normal_axis) -= sgn;
                    }
                }
            }
        }
    }
    return;
} // setPhysicalBoundaryConditions

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
