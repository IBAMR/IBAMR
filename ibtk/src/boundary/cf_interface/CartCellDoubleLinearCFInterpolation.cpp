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

#include <ibtk/CartCellDoubleLinearCFInterpolation.h>
#include <ibtk/samrai_compatibility_names.h>

#include <BoxArray.h>
#include <SAMRAIArray.h>
#include <SAMRAIBoundaryBox.h>
#include <SAMRAIBox.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICoarseFineBoundary.h>
#include <SAMRAIGridGeometry.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAIRefineOperator.h>

#include <set>
#include <string>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_LINEAR_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(cclinearnormalinterpolation2d, CCLINEARNORMALINTERPOLATION2D)
#endif
#if (NDIM == 3)
#define CC_LINEAR_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(cclinearnormalinterpolation3d, CCLINEARNORMALINTERPOLATION3D)
#endif

// Function interfaces
extern "C"
{
    void CC_LINEAR_NORMAL_INTERPOLATION_FC(double* U,
                                           const int& U_gcw,
                                           const int& ilower0,
                                           const int& iupper0,
                                           const int& ilower1,
                                           const int& iupper1,
#if (NDIM == 3)
                                           const int& ilower2,
                                           const int& iupper2,
#endif
                                           const int& loc_index,
                                           const int& ratio,
                                           const int* blower,
                                           const int* bupper);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;
static const int GHOST_WIDTH_TO_FILL = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartCellDoubleLinearCFInterpolation::~CartCellDoubleLinearCFInterpolation()
{
    clearPatchHierarchy();
    return;
} // ~CartCellDoubleLinearCFInterpolation

void
CartCellDoubleLinearCFInterpolation::setPhysicalBoundaryConditions(SAMRAIPatch& /*patch*/,
                                                                   const double /*fill_time*/,
                                                                   const SAMRAIIntVector& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
} // setPhysicalBoundaryConditions

SAMRAIIntVector
CartCellDoubleLinearCFInterpolation::getRefineOpStencilWidth() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_refine_op->getStencilWidth().max() <= REFINE_OP_STENCIL_WIDTH);
#endif
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartCellDoubleLinearCFInterpolation::preprocessRefine(SAMRAIPatch& /*fine*/,
                                                      const SAMRAIPatch& /*coarse*/,
                                                      const SAMRAIBox& /*fine_box*/,
                                                      const SAMRAIIntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartCellDoubleLinearCFInterpolation::postprocessRefine(SAMRAIPatch& fine,
                                                       const SAMRAIPatch& coarse,
                                                       const SAMRAIBox& fine_box,
                                                       const SAMRAIIntVector& ratio)
{
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        d_refine_op->refine(fine, coarse, patch_data_index, patch_data_index, fine_box, ratio);
    }
    return;
} // postprocessRefine

void
CartCellDoubleLinearCFInterpolation::setConsistentInterpolationScheme(const bool consistent_type_2_bdry)
{
    d_consistent_type_2_bdry = consistent_type_2_bdry;
    return;
} // setConsistentInterpolationScheme

void
CartCellDoubleLinearCFInterpolation::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
CartCellDoubleLinearCFInterpolation::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
CartCellDoubleLinearCFInterpolation::setPatchDataIndices(const ComponentSelector& patch_data_indices)
{
    std::set<int> patch_data_index_set;
    for (int l = 0; l < patch_data_indices.getSize(); ++l)
    {
        if (patch_data_indices.isSet(l))
        {
            const int patch_data_index = l;
            patch_data_index_set.insert(patch_data_index);
        }
    }
    setPatchDataIndices(patch_data_index_set);
    return;
} // setPatchDataIndices

void
CartCellDoubleLinearCFInterpolation::setPatchHierarchy(SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    if (d_hierarchy) clearPatchHierarchy();
    d_hierarchy = hierarchy;
    const int finest_level_number = d_hierarchy->getFinestLevelNumber();

    d_cf_boundary.resize(finest_level_number + 1);
    const SAMRAIIntVector& max_ghost_width = getRefineOpStencilWidth();
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        d_cf_boundary[ln] = SAMRAICoarseFineBoundary(*d_hierarchy, ln, max_ghost_width);
    }

    SAMRAIPointer<SAMRAIGridGeometry> grid_geom = d_hierarchy->getGridGeometry();
    const BoxArray<NDIM>& domain_boxes = grid_geom->getPhysicalDomain();

    d_domain_boxes.resize(finest_level_number + 1);
    d_periodic_shift.resize(finest_level_number + 1);
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
        const SAMRAIIntVector& ratio = level->getRatio();
        d_domain_boxes[ln] = BoxArray<NDIM>(domain_boxes);
        d_domain_boxes[ln].refine(ratio);
        d_periodic_shift[ln] = grid_geom->getPeriodicShift(ratio);
    }
    return;
} // setPatchHierarchy

void
CartCellDoubleLinearCFInterpolation::clearPatchHierarchy()
{
    d_hierarchy.setNull();
    d_cf_boundary.clear();
    d_domain_boxes.clear();
    d_periodic_shift.clear();
    return;
} // clearPatchHierarchy

void
CartCellDoubleLinearCFInterpolation::computeNormalExtension(SAMRAIPatch& patch,
                                                            const SAMRAIIntVector& ratio,
                                                            const SAMRAIIntVector& /*ghost_width_to_fill*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
    TBOX_ASSERT(!d_consistent_type_2_bdry);
    TBOX_ASSERT(ratio.min() == ratio.max());
#endif
    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!patch.inHierarchy()) return;

    // Get the co-dimension 1 cf boundary boxes.
    const int patch_num = patch.getPatchNumber();
    const int patch_level_num = patch.getPatchLevelNumber();
#if !defined(NDEBUG)
    SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(patch_level_num);
    TBOX_ASSERT(&patch == level->getPatch(patch_num).getPointer());
#endif
    const SAMRAIArray<SAMRAIBoundaryBox>& cf_bdry_codim1_boxes =
        d_cf_boundary[patch_level_num].getBoundaries(patch_num, 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        SAMRAIPointer<SAMRAICellData<double>> data = patch.getPatchData(patch_data_index);
#if !defined(NDEBUG)
        TBOX_ASSERT(data);
#endif
        const int U_ghosts = (data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_ghosts != (data->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleLinearCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        const int data_depth = data->getDepth();
        const SAMRAIIntVector ghost_width_to_fill = GHOST_WIDTH_TO_FILL;
        SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch.getPatchGeometry();
        const SAMRAIBox& patch_box = patch.getBox();
        for (int k = 0; k < n_cf_bdry_codim1_boxes; ++k)
        {
            const SAMRAIBoundaryBox& bdry_box = cf_bdry_codim1_boxes[k];
            const SAMRAIBox bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U = data->getPointer(depth);
                const int r = ratio.min();
                CC_LINEAR_NORMAL_INTERPOLATION_FC(U,
                                                  U_ghosts,
                                                  patch_box.lower(0),
                                                  patch_box.upper(0),
                                                  patch_box.lower(1),
                                                  patch_box.upper(1),
#if (NDIM == 3)
                                                  patch_box.lower(2),
                                                  patch_box.upper(2),
#endif
                                                  location_index,
                                                  r,
                                                  bc_fill_box.lower(),
                                                  bc_fill_box.upper());
            }
        }
    }
    return;
} // computeNormalExtension

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
