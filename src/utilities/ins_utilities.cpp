// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2024 by the IBAMR developers
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
#include "ibamr/ins_utilities.h"

#include "ibtk/PhysicalBoundaryUtilities.h"
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
double
computeNetInflowPhysicalBoundary(Pointer<HierarchyMathOps> hier_math_ops, int u_idx, int bdry_loc_idx)
{
    Pointer<PatchHierarchy<NDIM> > patch_hier = hier_math_ops->getPatchHierarchy();
    const int wgt_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();
    const int hier_finest_ln = patch_hier->getFinestLevelNumber();

    double integral = 0.0;
    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            if (!pgeom->getTouchesRegularBoundary()) continue;

            const tbox::Array<BoundaryBox<NDIM> > physical_codim1_boxes =
                PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            if (physical_codim1_boxes.size() == 0) continue;

            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
            Pointer<SideData<NDIM, double> > wgt_data = patch->getPatchData(wgt_sc_idx);
            const double* const dx = pgeom->getDx();

            for (int n = 0; n < physical_codim1_boxes.size(); ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const int bdry_box_loc_idx = bdry_box.getLocationIndex();
                if (bdry_box_loc_idx != bdry_loc_idx) continue;

                const unsigned int axis = bdry_loc_idx / 2;
                const unsigned int side = bdry_loc_idx % 2;
                const bool is_lower = side == 0;
                const double axis_normal = is_lower ? -1.0 : 1.0;

                BoundaryBox<NDIM> trimmed_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                Box<NDIM> integration_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_box);

                // Ensure the integration box is of width one in the normal direction
                integration_box.upper(axis) = integration_box.lower(axis);

                for (Box<NDIM>::Iterator b(integration_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);

                    const double u = (*u_data)(i_s);

                    // At a physical boundary a side-centered control volume is half
                    // the size of an internal one.
                    const double ds = (*wgt_data)(i_s) / (0.5 * dx[axis]);

                    integral += -u * ds * axis_normal;
                }
            }
        }
    }

    return IBTK_MPI::sumReduction(integral);

} // computeNetInFlowPhysicalBoundary

} // namespace IBAMR
