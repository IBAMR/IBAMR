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
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIArray.h"
#include "SAMRAIArrayData.h"
#include "SAMRAIBoundaryBox.h"
#include "SAMRAIBox.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h" // IWYU pragma: keep
#include "SAMRAIFaceData.h" // IWYU pragma: keep
#include "SAMRAIFaceIndex.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAIVariable.h"

#include <algorithm>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

void
AdvDiffPhysicalBoundaryUtilities::setPhysicalBoundaryConditions(SAMRAIPointer<SAMRAICellData<double> > Q_data,
                                                                SAMRAIPointer<SAMRAIFaceData<double> > u_ADV_data,
                                                                SAMRAIPointer<SAMRAIPatch> patch,
                                                                const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
                                                                const double fill_time,
                                                                const bool inflow_boundaries_only,
                                                                const bool homogeneous_bc)
{
    SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch->getPatchGeometry();
    if (!pgeom->getTouchesRegularBoundary()) return;
    const SAMRAIArray<SAMRAIBoundaryBox> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    if (physical_codim1_boxes.size() == 0) return;

    // Loop over the boundary fill boxes and set boundary conditions.
    const SAMRAIBox& patch_box = patch->getBox();
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
    const SAMRAIIntVector& gcw = Q_data->getGhostCellWidth();
    for (int n = 0; n < physical_codim1_boxes.size(); ++n)
    {
        const SAMRAIBoundaryBox& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower = location_index % 2 == 0;
        static const SAMRAIIntVector gcw_to_fill = 1;
        const SAMRAIBox bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
        const SAMRAIBoundaryBox trimmed_bdry_box(
            bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
        SAMRAIBox bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d != bdry_normal_axis)
            {
                bc_coef_box.lower(d) = std::max(bc_coef_box.lower(d), patch_box.lower(d));
                bc_coef_box.upper(d) = std::min(bc_coef_box.upper(d), patch_box.upper(d));
            }
        }
        SAMRAIPointer<SAMRAIArrayData<double> > acoef_data = new SAMRAIArrayData<double>(bc_coef_box, 1);
        SAMRAIPointer<SAMRAIArrayData<double> > bcoef_data = new SAMRAIArrayData<double>(bc_coef_box, 1);
        SAMRAIPointer<SAMRAIArrayData<double> > gcoef_data = new SAMRAIArrayData<double>(bc_coef_box, 1);
        for (int depth = 0; depth < Q_data->getDepth(); ++depth)
        {
            if (bc_coefs[depth] == nullptr) continue;

            bc_coefs[depth]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, fill_time);
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[depth]);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
            for (SAMRAIBox::Iterator bc(bc_coef_box); bc; bc++)
            {
                const SAMRAIIndex& i = bc();
                const SAMRAIFaceIndex i_f(i, bdry_normal_axis, SAMRAIFaceIndex::Lower);
                bool is_inflow_bdry = (is_lower && (*u_ADV_data)(i_f) > 0.0) || (!is_lower && (*u_ADV_data)(i_f) < 0.0);
                if (!inflow_boundaries_only || is_inflow_bdry)
                {
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const double& g = (*gcoef_data)(i, 0);
                    const double& h = dx[bdry_normal_axis];
                    int sgn;
                    SAMRAIIndex i_i(i), i_g(i);
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
