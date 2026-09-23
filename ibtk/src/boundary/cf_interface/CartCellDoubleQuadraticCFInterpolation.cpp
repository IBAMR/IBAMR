// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>

#include <tbox/Array.h>

#include <BoundaryBox.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CoarseFineBoundary.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineOperator.h>

#include <set>

#include <ibtk/namespaces.h> // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(ccquadtangentialinterpolation2d, CCQUADTANGENTIALINTERPOLATION2D)
#define CC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadnormalinterpolation2d, CCQUADNORMALINTERPOLATION2D)
#endif
#if (NDIM == 3)
#define CC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(ccquadtangentialinterpolation3d, CCQUADTANGENTIALINTERPOLATION3D)
#define CC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadnormalinterpolation3d, CCQUADNORMALINTERPOLATION3D)
#endif

// Function interfaces
extern "C"
{
    void CC_QUAD_TANGENTIAL_INTERPOLATION_FC(double* U_fine,
                                             const int& U_fine_gcw,
                                             const double* U_coarse,
                                             const int& U_crse_gcw,
                                             const int& ilowerf0,
                                             const int& iupperf0,
                                             const int& ilowerf1,
                                             const int& iupperf1,
#if (NDIM == 3)
                                             const int& ilowerf2,
                                             const int& iupperf2,
#endif
                                             const int& ilowerc0,
                                             const int& iupperc0,
                                             const int& ilowerc1,
                                             const int& iupperc1,
#if (NDIM == 3)
                                             const int& ilowerc2,
                                             const int& iupperc2,
#endif
                                             const int& loc_index,
                                             const int* ratio_to_coarser,
                                             const int* blower,
                                             const int* bupper);

    void CC_QUAD_NORMAL_INTERPOLATION_FC(double* U,
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
                                         const int* ratio_to_coarser,
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

CartCellDoubleQuadraticCFInterpolation::~CartCellDoubleQuadraticCFInterpolation()
{
    clearPatchHierarchy();
    return;
} // ~CartCellDoubleQuadraticCFInterpolation

void
CartCellDoubleQuadraticCFInterpolation::setPhysicalBoundaryConditions(Patch<NDIM>& /*patch*/,
                                                                      const double /*fill_time*/,
                                                                      const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
CartCellDoubleQuadraticCFInterpolation::getRefineOpStencilWidth() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_refine_op->getStencilWidth().max() <= REFINE_OP_STENCIL_WIDTH);
#endif
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartCellDoubleQuadraticCFInterpolation::preprocessRefine(Patch<NDIM>& /*fine*/,
                                                         const Patch<NDIM>& /*coarse*/,
                                                         const Box<NDIM>& /*fine_box*/,
                                                         const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartCellDoubleQuadraticCFInterpolation::postprocessRefine(Patch<NDIM>& fine,
                                                          const Patch<NDIM>& coarse,
                                                          const Box<NDIM>& fine_box,
                                                          const IntVector<NDIM>& ratio)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!fine.inHierarchy())
    {
        for (int patch_data_index : d_patch_data_indices)
        {
            d_refine_op->refine(fine, coarse, patch_data_index, patch_data_index, fine_box, ratio);
        }
        return;
    }
#if !defined(NDEBUG)
    else
    {
        // Ensure the fine patch corresponds to the expected patch in the cached
        // patch hierarchy.
        const int patch_num = fine.getPatchNumber();
        const int fine_patch_level_num = fine.getPatchLevelNumber();
        Pointer<PatchLevel<NDIM>> fine_level = d_hierarchy->getPatchLevel(fine_patch_level_num);
        TBOX_ASSERT(&fine == fine_level->getPatch(patch_num).getPointer());
    }
#endif
    postprocessRefine_optimized(fine, coarse, ratio);
    return;
} // postprocessRefine

void
CartCellDoubleQuadraticCFInterpolation::setConsistentInterpolationScheme(const bool consistent_type_2_bdry)
{
    d_consistent_type_2_bdry = consistent_type_2_bdry;
    return;
} // setConsistentInterpolationScheme

void
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
CartCellDoubleQuadraticCFInterpolation::setPatchDataIndices(const ComponentSelector& patch_data_indices)
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
CartCellDoubleQuadraticCFInterpolation::setPatchHierarchy(Pointer<PatchHierarchy<NDIM>> hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    if (d_hierarchy) clearPatchHierarchy();
    d_hierarchy = hierarchy;
    const int finest_level_number = d_hierarchy->getFinestLevelNumber();

    d_cf_boundary.resize(finest_level_number + 1);
    const IntVector<NDIM>& max_ghost_width = getRefineOpStencilWidth();
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        d_cf_boundary[ln] = CoarseFineBoundary<NDIM>(*d_hierarchy, ln, max_ghost_width);
    }
    return;
} // setPatchHierarchy

void
CartCellDoubleQuadraticCFInterpolation::clearPatchHierarchy()
{
    d_hierarchy.setNull();
    d_cf_boundary.clear();
    return;
} // clearPatchHierarchy

void
CartCellDoubleQuadraticCFInterpolation::computeNormalExtension(Patch<NDIM>& patch,
                                                               const IntVector<NDIM>& ratio,
                                                               const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    // Ensure that the fine patch is located on the expected destination level;
    // if not, we are not guaranteed to have appropriate coarse-fine interface
    // boundary box information.
    if (!patch.inHierarchy())
    {
        return;
    }
#if !defined(NDEBUG)
    else
    {
        const int patch_num = patch.getPatchNumber();
        const int patch_level_num = patch.getPatchLevelNumber();
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(patch_level_num);
        TBOX_ASSERT(&patch == level->getPatch(patch_num).getPointer());
    }
#endif
    computeNormalExtension_optimized(patch, ratio);
    return;
} // computeNormalExtension

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CartCellDoubleQuadraticCFInterpolation::postprocessRefine_optimized(Patch<NDIM>& fine,
                                                                    const Patch<NDIM>& coarse,
                                                                    const IntVector<NDIM>& ratio)
{
    // Get the co-dimension 1 cf boundary boxes.
    const int patch_num = fine.getPatchNumber();
    const int fine_patch_level_num = fine.getPatchLevelNumber();
    const Array<BoundaryBox<NDIM>>& cf_bdry_codim1_boxes =
        d_cf_boundary[fine_patch_level_num].getBoundaries(patch_num, 1);
    if (cf_bdry_codim1_boxes.size() == 0) return;

    // Get the patch data.
    for (const auto& patch_data_index : d_patch_data_indices)
    {
        Pointer<CellData<NDIM, double>> fdata = fine.getPatchData(patch_data_index);
        Pointer<CellData<NDIM, double>> cdata = coarse.getPatchData(patch_data_index);
#if !defined(NDEBUG)
        TBOX_ASSERT(fdata);
        TBOX_ASSERT(cdata);
        TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif
        const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
        const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
        if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        const int data_depth = fdata->getDepth();
        const IntVector<NDIM> ghost_width_to_fill = GHOST_WIDTH_TO_FILL;
        Pointer<CartesianPatchGeometry<NDIM>> pgeom_fine = fine.getPatchGeometry();
        const Box<NDIM>& patch_box_fine = fine.getBox();
        const Box<NDIM>& patch_box_crse = coarse.getBox();
        for (int k = 0; k < cf_bdry_codim1_boxes.size(); ++k)
        {
            const BoundaryBox<NDIM>& bdry_box = cf_bdry_codim1_boxes[k];
            const Box<NDIM> bc_fill_box = pgeom_fine->getBoundaryFillBox(bdry_box, patch_box_fine, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U_fine = fdata->getPointer(depth);
                const double* const U_crse = cdata->getPointer(depth);
                CC_QUAD_TANGENTIAL_INTERPOLATION_FC(U_fine,
                                                    U_fine_ghosts,
                                                    U_crse,
                                                    U_crse_ghosts,
                                                    patch_box_fine.lower(0),
                                                    patch_box_fine.upper(0),
                                                    patch_box_fine.lower(1),
                                                    patch_box_fine.upper(1),
#if (NDIM == 3)
                                                    patch_box_fine.lower(2),
                                                    patch_box_fine.upper(2),
#endif
                                                    patch_box_crse.lower(0),
                                                    patch_box_crse.upper(0),
                                                    patch_box_crse.lower(1),
                                                    patch_box_crse.upper(1),
#if (NDIM == 3)
                                                    patch_box_crse.lower(2),
                                                    patch_box_crse.upper(2),
#endif
                                                    location_index,
                                                    ratio,
                                                    bc_fill_box.lower(),
                                                    bc_fill_box.upper());
            }
        }
    }
    return;
} // postprocessRefine_optimized

void
CartCellDoubleQuadraticCFInterpolation::computeNormalExtension_optimized(Patch<NDIM>& patch,
                                                                         const IntVector<NDIM>& ratio)
{
    // Get the co-dimension 1 cf boundary boxes.
    const int patch_num = patch.getPatchNumber();
    const int patch_level_num = patch.getPatchLevelNumber();
    const Array<BoundaryBox<NDIM>>& cf_bdry_codim1_boxes = d_cf_boundary[patch_level_num].getBoundaries(patch_num, 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Get the patch data.
    for (int patch_data_index : d_patch_data_indices)
    {
        Pointer<CellData<NDIM, double>> data = patch.getPatchData(patch_data_index);
#if !defined(NDEBUG)
        TBOX_ASSERT(data);
#endif
        const int U_ghosts = (data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_ghosts != (data->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartCellDoubleQuadraticCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths" << std::endl);
        }
#endif
        const int data_depth = data->getDepth();
        const IntVector<NDIM> ghost_width_to_fill = GHOST_WIDTH_TO_FILL;
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
        const Box<NDIM>& patch_box = patch.getBox();
        for (int k = 0; k < n_cf_bdry_codim1_boxes; ++k)
        {
            const BoundaryBox<NDIM>& bdry_box = cf_bdry_codim1_boxes[k];
            const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U = data->getPointer(depth);
                CC_QUAD_NORMAL_INTERPOLATION_FC(U,
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
                                                ratio,
                                                bc_fill_box.lower(),
                                                bc_fill_box.upper());
            }
        }
    }
    return;
} // computeNormalExtension_optimized

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
