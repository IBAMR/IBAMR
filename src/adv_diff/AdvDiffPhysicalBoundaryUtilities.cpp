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

#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"

#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h" // IWYU pragma: keep
#include "FaceData.h" // IWYU pragma: keep
#include "FaceIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "RobinBcCoefStrategy.h"
#include "Variable.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

#include <algorithm>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

void
AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(Pointer<CellData<NDIM, double> > Q_data,
                                                                Pointer<FaceData<NDIM, double> > u_ADV_data,
                                                                Pointer<Patch<NDIM> > patch,
                                                                const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                                const double fill_time,
                                                                const bool inflow_boundaries_only,
                                                                const bool homogeneous_bc)
{
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    if (!pgeom->getTouchesRegularBoundary()) return;
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    if (physical_codim1_boxes.size() == 0) return;

    // Loop over the boundary fill boxes and set boundary conditions.
    const Box<NDIM>& patch_box = patch->getBox();
    const double* const dx = pgeom->getDx();

    // Setup any extended Robin BC coef objects.
    for (int depth = 0; depth < Q_data->getDepth(); ++depth)
    {
        auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[depth]);
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
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower = location_index % 2 == 0;
        static const IntVector<NDIM> gcw_to_fill = 1;
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const BoundaryBox<NDIM> trimmed_bdry_box(
            bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
        Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d != bdry_normal_axis)
            {
                bc_coef_box.lower(d) = std::max(bc_coef_box.lower(d), patch_box.lower(d));
                bc_coef_box.upper(d) = std::min(bc_coef_box.upper(d), patch_box.upper(d));
            }
        }
        Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        for (int depth = 0; depth < Q_data->getDepth(); ++depth)
        {
            if (bc_coefs[depth] == nullptr) continue;

            bc_coefs[depth]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, fill_time);
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[depth]);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const hier::Index<NDIM>& i = bc();
                const FaceIndex<NDIM> i_f(i, bdry_normal_axis, FaceIndex<NDIM>::Lower);
                bool is_inflow_bdry = (is_lower && (*u_ADV_data)(i_f) > 0.0) || (!is_lower && (*u_ADV_data)(i_f) < 0.0);
                if (!inflow_boundaries_only || is_inflow_bdry)
                {
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const double& g = (*gcoef_data)(i, 0);
                    const double& h = dx[bdry_normal_axis];
                    int sgn;
                    hier::Index<NDIM> i_i(i), i_g(i);
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
