// Filename: CartSideDoubleQuadraticCFInterpolation.cpp
// Created on 30 Apr 2008 by Boyce Griffith
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
#include <ostream>
#include <set>
#include <vector>

#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CartesianSideDoubleConservativeLinearRefine.h"
#include "CoarseFineBoundary.h"
#include "ComponentSelector.h"
#include "IBTK_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibtk/CartSideDoubleQuadraticCFInterpolation.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
#if (NDIM == 2)
#define SC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(scquadtangentialinterpolation2d, SCQUADTANGENTIALINTERPOLATION2D)
#define SC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(scquadnormalinterpolation2d, SCQUADNORMALINTERPOLATION2D)
#endif
#if (NDIM == 3)
#define SC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(scquadtangentialinterpolation3d, SCQUADTANGENTIALINTERPOLATION3D)
#define SC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(scquadnormalinterpolation3d, SCQUADNORMALINTERPOLATION3D)
#endif

// Function interfaces
extern "C" {
void SC_QUAD_TANGENTIAL_INTERPOLATION_FC(double* U_fine0,
                                         double* U_fine1,
#if (NDIM == 3)
                                         double* U_fine2,
#endif
                                         const int& U_fine_gcw,
                                         const double* U_crse0,
                                         const double* U_crse1,
#if (NDIM == 3)
                                         const double* U_crse2,
#endif
                                         const int& U_crse_gcw,
                                         const int* sc_indicator0,
                                         const int* sc_indicator1,
#if (NDIM == 3)
                                         const int* sc_indicator2,
#endif
                                         const int& sc_indicator_gcw,
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

void SC_QUAD_NORMAL_INTERPOLATION_FC(double* U0,
                                     double* U1,
#if (NDIM == 3)
                                     double* U2,
#endif
                                     const int& U_gcw,
                                     const double* W0,
                                     const double* W1,
#if (NDIM == 3)
                                     const double* W2,
#endif
                                     const int& W_gcw,
                                     const int* sc_indicator0,
                                     const int* sc_indicator1,
#if (NDIM == 3)
                                     const int* sc_indicator2,
#endif
                                     const int& sc_indicator_gcw,
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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleQuadraticCFInterpolation::CartSideDoubleQuadraticCFInterpolation()
    : d_patch_data_indices(),
      d_consistent_type_2_bdry(false),
      d_refine_op(new CartesianSideDoubleConservativeLinearRefine<NDIM>()),
      d_hierarchy(NULL),
      d_cf_boundary(),
      d_sc_indicator_var(new SideVariable<NDIM, int>("CartSideDoubleQuadraticCFInterpolation::sc_indicator_var"))
{
    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("CartSideDoubleQuadraticCFInterpolation::CONTEXT");
    if (var_db->checkVariableExists(d_sc_indicator_var->getName()))
    {
        d_sc_indicator_var = var_db->getVariable(d_sc_indicator_var->getName());
        d_sc_indicator_idx = var_db->mapVariableAndContextToIndex(d_sc_indicator_var, context);
    }
    else
    {
        d_sc_indicator_idx = var_db->registerVariableAndContext(d_sc_indicator_var, context, GHOST_WIDTH_TO_FILL);
    }
    return;
} // CartSideDoubleQuadraticCFInterpolation

CartSideDoubleQuadraticCFInterpolation::~CartSideDoubleQuadraticCFInterpolation()
{
    clearPatchHierarchy();
    return;
} // ~CartSideDoubleQuadraticCFInterpolation

void
CartSideDoubleQuadraticCFInterpolation::setPhysicalBoundaryConditions(Patch<NDIM>& /*patch*/,
                                                                      const double /*fill_time*/,
                                                                      const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
CartSideDoubleQuadraticCFInterpolation::getRefineOpStencilWidth() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_refine_op->getStencilWidth().max() <= REFINE_OP_STENCIL_WIDTH);
#endif
    return REFINE_OP_STENCIL_WIDTH;
} // getRefineOpStencilWidth

void
CartSideDoubleQuadraticCFInterpolation::preprocessRefine(Patch<NDIM>& /*fine*/,
                                                         const Patch<NDIM>& /*coarse*/,
                                                         const Box<NDIM>& /*fine_box*/,
                                                         const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartSideDoubleQuadraticCFInterpolation::postprocessRefine(Patch<NDIM>& fine,
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
        for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
        {
            const int& patch_data_index = *cit;
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
        Pointer<PatchLevel<NDIM> > fine_level = d_hierarchy->getPatchLevel(fine_patch_level_num);
        TBOX_ASSERT(&fine == fine_level->getPatch(patch_num).getPointer());
    }
#endif
    // Get the co-dimension 1 cf boundary boxes.
    const int patch_num = fine.getPatchNumber();
    const int fine_patch_level_num = fine.getPatchLevelNumber();
    const Array<BoundaryBox<NDIM> >& cf_bdry_codim1_boxes =
        d_cf_boundary[fine_patch_level_num]->getBoundaries(patch_num, 1);
    if (cf_bdry_codim1_boxes.size() == 0) return;

    // Get the patch data.
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int& patch_data_index = *cit;
        Pointer<SideData<NDIM, double> > fdata = fine.getPatchData(patch_data_index);
        Pointer<SideData<NDIM, double> > cdata = coarse.getPatchData(patch_data_index);
        Pointer<SideData<NDIM, int> > indicator_data = fine.getPatchData(d_sc_indicator_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(fdata);
        TBOX_ASSERT(cdata);
        TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
        TBOX_ASSERT(indicator_data);
#endif
        const int U_fine_ghosts = (fdata->getGhostCellWidth()).max();
        const int U_crse_ghosts = (cdata->getGhostCellWidth()).max();
        const int indicator_ghosts = (indicator_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_fine_ghosts != (fdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths"
                       << std::endl);
        }
        if (U_crse_ghosts != (cdata->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideDoubleQuadraticCFInterpolation::postprocessRefine():\n"
                       << "   patch data does not have uniform ghost cell widths"
                       << std::endl);
        }
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).max() == GHOST_WIDTH_TO_FILL);
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).min() == GHOST_WIDTH_TO_FILL);
#endif
        const int data_depth = fdata->getDepth();
        const IntVector<NDIM> ghost_width_to_fill = GHOST_WIDTH_TO_FILL;
        Pointer<CartesianPatchGeometry<NDIM> > pgeom_fine = fine.getPatchGeometry();
        const Box<NDIM>& patch_box_fine = fine.getBox();
        const Box<NDIM>& patch_box_crse = coarse.getBox();
        for (int k = 0; k < cf_bdry_codim1_boxes.size(); ++k)
        {
            const BoundaryBox<NDIM>& bdry_box = cf_bdry_codim1_boxes[k];
            const Box<NDIM> bc_fill_box = pgeom_fine->getBoundaryFillBox(bdry_box, patch_box_fine, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            const int* const indicator0 = indicator_data->getPointer(0);
            const int* const indicator1 = indicator_data->getPointer(1);
#if (NDIM == 3)
            const int* const indicator2 = indicator_data->getPointer(2);
#endif
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U_fine0 = fdata->getPointer(0, depth);
                double* const U_fine1 = fdata->getPointer(1, depth);
#if (NDIM == 3)
                double* const U_fine2 = fdata->getPointer(2, depth);
#endif
                const double* const U_crse0 = cdata->getPointer(0, depth);
                const double* const U_crse1 = cdata->getPointer(1, depth);
#if (NDIM == 3)
                const double* const U_crse2 = cdata->getPointer(2, depth);
#endif
                SC_QUAD_TANGENTIAL_INTERPOLATION_FC(U_fine0,
                                                    U_fine1,
#if (NDIM == 3)
                                                    U_fine2,
#endif
                                                    U_fine_ghosts,
                                                    U_crse0,
                                                    U_crse1,
#if (NDIM == 3)
                                                    U_crse2,
#endif
                                                    U_crse_ghosts,
                                                    indicator0,
                                                    indicator1,
#if (NDIM == 3)
                                                    indicator2,
#endif
                                                    indicator_ghosts,
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
} // postprocessRefine

void
CartSideDoubleQuadraticCFInterpolation::setConsistentInterpolationScheme(const bool consistent_type_2_bdry)
{
    d_consistent_type_2_bdry = consistent_type_2_bdry;
    return;
} // setConsistentInterpolationScheme

void
CartSideDoubleQuadraticCFInterpolation::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
CartSideDoubleQuadraticCFInterpolation::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
CartSideDoubleQuadraticCFInterpolation::setPatchDataIndices(const ComponentSelector& patch_data_indices)
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
CartSideDoubleQuadraticCFInterpolation::setPatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    if (d_hierarchy) clearPatchHierarchy();
    d_hierarchy = hierarchy;
    const int finest_level_number = d_hierarchy->getFinestLevelNumber();

    d_cf_boundary.resize(finest_level_number + 1);
    const IntVector<NDIM>& max_ghost_width = GHOST_WIDTH_TO_FILL;
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        d_cf_boundary[ln] = new CoarseFineBoundary<NDIM>(*d_hierarchy, ln, max_ghost_width);
    }

    Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op = NULL;
    refine_alg->registerRefine(d_sc_indicator_idx, // destination
                               d_sc_indicator_idx, // source
                               d_sc_indicator_idx, // temporary work space
                               refine_op);
    for (int ln = 0; ln <= finest_level_number; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_sc_indicator_idx))
        {
            level->allocatePatchData(d_sc_indicator_idx, 0.0);
        }
        else
        {
            level->setTime(0.0, d_sc_indicator_idx);
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, int> > sc_indicator_data = patch->getPatchData(d_sc_indicator_idx);
            sc_indicator_data->fillAll(0, sc_indicator_data->getGhostBox());
            sc_indicator_data->fillAll(1, sc_indicator_data->getBox());
        }
        refine_alg->createSchedule(d_hierarchy->getPatchLevel(ln))->fillData(0.0);
    }
    return;
} // setPatchHierarchy

void
CartSideDoubleQuadraticCFInterpolation::clearPatchHierarchy()
{
    d_hierarchy.setNull();
    for (std::vector<CoarseFineBoundary<NDIM>*>::iterator it = d_cf_boundary.begin(); it != d_cf_boundary.end(); ++it)
    {
        delete (*it);
        (*it) = NULL;
    }
    d_cf_boundary.clear();
    return;
} // clearPatchHierarchy

void
CartSideDoubleQuadraticCFInterpolation::computeNormalExtension(Patch<NDIM>& patch,
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
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(patch_level_num);
        TBOX_ASSERT(&patch == level->getPatch(patch_num).getPointer());
    }
#endif
    // Get the co-dimension 1 cf boundary boxes.
    const int patch_num = patch.getPatchNumber();
    const int patch_level_num = patch.getPatchLevelNumber();
    const Array<BoundaryBox<NDIM> >& cf_bdry_codim1_boxes = d_cf_boundary[patch_level_num]->getBoundaries(patch_num, 1);
    const int n_cf_bdry_codim1_boxes = cf_bdry_codim1_boxes.size();

    // Check to see if there are any co-dimension 1 coarse-fine boundary boxes
    // associated with the patch; if not, there is nothing to do.
    if (n_cf_bdry_codim1_boxes == 0) return;

    // Get the patch data.
    for (std::set<int>::const_iterator cit = d_patch_data_indices.begin(); cit != d_patch_data_indices.end(); ++cit)
    {
        const int& patch_data_index = *cit;
        Pointer<SideData<NDIM, double> > data = patch.getPatchData(patch_data_index);
        SideData<NDIM, double> data_copy(data->getBox(), data->getDepth(), data->getGhostCellWidth());
        data_copy.copyOnBox(*data, data->getGhostBox());
        Pointer<SideData<NDIM, int> > indicator_data = patch.getPatchData(d_sc_indicator_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(data);
        TBOX_ASSERT(indicator_data);
#endif
        const int U_ghosts = (data->getGhostCellWidth()).max();
        const int W_ghosts = (data_copy.getGhostCellWidth()).max();
        const int indicator_ghosts = (indicator_data->getGhostCellWidth()).max();
#if !defined(NDEBUG)
        if (U_ghosts != (data->getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideDoubleQuadraticCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths"
                       << std::endl);
        }
        if (W_ghosts != (data_copy.getGhostCellWidth()).min())
        {
            TBOX_ERROR("CartSideDoubleQuadraticCFInterpolation::computeNormalExtension():\n"
                       << "   patch data does not have uniform ghost cell widths"
                       << std::endl);
        }
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).max() == GHOST_WIDTH_TO_FILL);
        TBOX_ASSERT((indicator_data->getGhostCellWidth()).min() == GHOST_WIDTH_TO_FILL);
#endif
        const int data_depth = data->getDepth();
        const IntVector<NDIM> ghost_width_to_fill = GHOST_WIDTH_TO_FILL;
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
        const Box<NDIM>& patch_box = patch.getBox();
        for (int k = 0; k < n_cf_bdry_codim1_boxes; ++k)
        {
            const BoundaryBox<NDIM>& bdry_box = cf_bdry_codim1_boxes[k];
            const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);
            const unsigned int location_index = bdry_box.getLocationIndex();
            const int* const indicator0 = indicator_data->getPointer(0);
            const int* const indicator1 = indicator_data->getPointer(1);
#if (NDIM == 3)
            const int* const indicator2 = indicator_data->getPointer(2);
#endif
            for (int depth = 0; depth < data_depth; ++depth)
            {
                double* const U0 = data->getPointer(0, depth);
                double* const U1 = data->getPointer(1, depth);
#if (NDIM == 3)
                double* const U2 = data->getPointer(2, depth);
#endif
                const double* const W0 = data_copy.getPointer(0, depth);
                const double* const W1 = data_copy.getPointer(1, depth);
#if (NDIM == 3)
                const double* const W2 = data_copy.getPointer(2, depth);
#endif
                SC_QUAD_NORMAL_INTERPOLATION_FC(U0,
                                                U1,
#if (NDIM == 3)
                                                U2,
#endif
                                                U_ghosts,
                                                W0,
                                                W1,
#if (NDIM == 3)
                                                W2,
#endif
                                                W_ghosts,
                                                indicator0,
                                                indicator1,
#if (NDIM == 3)
                                                indicator2,
#endif
                                                indicator_ghosts,
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
} // computeNormalExtension

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
