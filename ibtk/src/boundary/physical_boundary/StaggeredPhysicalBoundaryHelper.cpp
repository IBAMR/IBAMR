// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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

#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "ibtk/ibtk_utilities.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "Variable.h"
#include "tbox/Array.h"
#include "tbox/MathUtilities.h"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StaggeredPhysicalBoundaryHelper::copyDataAtDirichletBoundaries(const int u_out_data_idx,
                                                               const int u_in_data_idx,
                                                               const int coarsest_ln,
                                                               const int finest_ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == invalid_level_number ? 0 : coarsest_ln);
         ln <= (finest_ln == invalid_level_number ? finest_hier_level : finest_ln);
         ++ln)
    {
        Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            Pointer<PatchNd> patch = level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                Pointer<SideDataNd<double> > u_out_data = patch->getPatchData(u_out_data_idx);
                Pointer<SideDataNd<double> > u_in_data = patch->getPatchData(u_in_data_idx);
                copyDataAtDirichletBoundaries(u_out_data, u_in_data, patch);
            }
        }
    }
    return;
} // copyDataAtDirichletBoundaries

void
StaggeredPhysicalBoundaryHelper::copyDataAtDirichletBoundaries(Pointer<SideDataNd<double> > u_out_data,
                                                               Pointer<SideDataNd<double> > u_in_data,
                                                               Pointer<PatchNd> patch) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    const Array<BoundaryBoxNd>& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayDataNd<bool> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        const ArrayDataNd<bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (BoxNd::Iterator it(bdry_locs_data.getBox()); it; it++)
        {
            const hier::IndexNd& i = it();
            if (bdry_locs_data(i, 0))
            {
                (*u_out_data)(SideIndexNd(i, bdry_normal_axis, SideIndexNd::Lower)) =
                    (*u_in_data)(SideIndexNd(i, bdry_normal_axis, SideIndexNd::Lower));
            }
        }
    }
    return;
} // copyDataAtDirichletBoundaries

void
StaggeredPhysicalBoundaryHelper::setupMaskingFunction(const int mask_data_idx,
                                                      const int coarsest_ln,
                                                      const int finest_ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == invalid_level_number ? 0 : coarsest_ln);
         ln <= (finest_ln == invalid_level_number ? finest_hier_level : finest_ln);
         ++ln)
    {
        Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            Pointer<PatchNd> patch = level->getPatch(p());
            Pointer<SideDataNd<int> > mask_data = patch->getPatchData(mask_data_idx);
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                setupMaskingFunction(mask_data, patch);
            }
            else
            {
                mask_data->fillAll(0);
            }
        }
    }
    return;
} // setupMaskingFunction

void
StaggeredPhysicalBoundaryHelper::setupMaskingFunction(Pointer<SideDataNd<int> > mask_data, Pointer<PatchNd> patch) const
{
    mask_data->fillAll(0);
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    const Array<BoundaryBoxNd>& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayDataNd<bool> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        const ArrayDataNd<bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (BoxNd::Iterator it(bdry_locs_data.getBox()); it; it++)
        {
            const hier::IndexNd& i = it();
            if (bdry_locs_data(i, 0)) (*mask_data)(SideIndexNd(i, bdry_normal_axis, SideIndexNd::Lower)) = 1;
        }
    }
    return;
} // setupMaskingFunction

bool
StaggeredPhysicalBoundaryHelper::patchTouchesDirichletBoundary(Pointer<PatchNd> patch) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return false;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        if (patchTouchesDirichletBoundaryAxis(patch, axis)) return true;
    }
    return false;
} // patchTouchesDirichletBoundary

bool
StaggeredPhysicalBoundaryHelper::patchTouchesDirichletBoundaryAxis(Pointer<PatchNd> patch,
                                                                   const unsigned int axis) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return false;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    const Array<BoundaryBoxNd>& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayDataNd<bool> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const unsigned int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        if (bdry_normal_axis == axis)
        {
            const ArrayDataNd<bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
            for (BoxNd::Iterator it(bdry_locs_data.getBox()); it; it++)
            {
                if (bdry_locs_data(it(), 0)) return true;
            }
        }
    }
    return false;
} // patchTouchesDirichletBoundaryAxis

void
StaggeredPhysicalBoundaryHelper::cacheBcCoefData(const std::vector<RobinBcCoefStrategyNd*>& u_bc_coefs,
                                                 const double fill_time,
                                                 const Pointer<PatchHierarchyNd> hierarchy)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
    TBOX_ASSERT(hierarchy);
#endif
    if (d_hierarchy) clearBcCoefData();

    // Cache boundary values.
    d_hierarchy = hierarchy;
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    d_physical_codim1_boxes.resize(finest_hier_level + 1);
    d_dirichlet_bdry_locs.resize(finest_hier_level + 1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        Pointer<PatchLevelNd> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<PatchNd> patch = level->getPatch(patch_num);
            Pointer<CartesianPatchGeometryNd> pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Array<BoundaryBoxNd>& physical_codim1_boxes = d_physical_codim1_boxes[ln][patch_num];
                physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                std::vector<Pointer<ArrayDataNd<bool> > >& dirichlet_bdry_locs = d_dirichlet_bdry_locs[ln][patch_num];
                dirichlet_bdry_locs.resize(n_physical_codim1_boxes);
                BoxNd bc_coef_box;
                BoundaryBoxNd trimmed_bdry_box;
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBoxNd& bdry_box = physical_codim1_boxes[n];
                    StaggeredPhysicalBoundaryHelper::setupBcCoefBoxes(bc_coef_box, trimmed_bdry_box, bdry_box, patch);
                    const unsigned int bdry_normal_axis = bdry_box.getLocationIndex() / 2;
                    Pointer<ArrayDataNd<double> > acoef_data = new ArrayDataNd<double>(bc_coef_box, 1);
                    Pointer<ArrayDataNd<double> > bcoef_data = new ArrayDataNd<double>(bc_coef_box, 1);
                    Pointer<ArrayDataNd<double> > gcoef_data;
                    u_bc_coefs[bdry_normal_axis]->setBcCoefs(
                        acoef_data, bcoef_data, gcoef_data, Pointer<VariableNd>(), *patch, trimmed_bdry_box, fill_time);
                    dirichlet_bdry_locs[n] = new ArrayDataNd<bool>(bc_coef_box, 1);
                    ArrayDataNd<bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
                    for (BoxNd::Iterator it(bc_coef_box); it; it++)
                    {
                        const hier::IndexNd& i = it();
                        const double& alpha = (*acoef_data)(i, 0);
                        const double& beta = (*bcoef_data)(i, 0);
#if !defined(NDEBUG)
                        TBOX_ASSERT(IBTK::rel_equal_eps(alpha + beta, 1.0));
                        TBOX_ASSERT(IBTK::rel_equal_eps(alpha, 1.0) || IBTK::rel_equal_eps(beta, 1.0));
#endif
                        bdry_locs_data(i, 0) = IBTK::rel_equal_eps(alpha, 1.0) && IBTK::abs_equal_eps(beta, 0.0);
                    }
                }
            }
            else
            {
                d_physical_codim1_boxes[ln][patch_num].resizeArray(0);
                d_dirichlet_bdry_locs[ln][patch_num].clear();
            }
        }
    }
    return;
} // cacheBcCoefData

void
StaggeredPhysicalBoundaryHelper::clearBcCoefData()
{
    d_hierarchy.setNull();
    d_physical_codim1_boxes.clear();
    d_dirichlet_bdry_locs.clear();
    return;
} // clearBcCoefData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredPhysicalBoundaryHelper::setupBcCoefBoxes(BoxNd& bc_coef_box,
                                                  BoundaryBoxNd& trimmed_bdry_box,
                                                  const BoundaryBoxNd& bdry_box,
                                                  Pointer<PatchNd> patch)
{
    Pointer<PatchGeometryNd> pgeom = patch->getPatchGeometry();
    const BoxNd& patch_box = patch->getBox();
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    BoxNd bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, /* gcw_to_fill */ IntVectorNd(1));
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d != bdry_normal_axis)
        {
            bc_fill_box.lower(d) = std::max(bc_fill_box.lower(d), patch_box.lower(d));
            bc_fill_box.upper(d) = std::min(bc_fill_box.upper(d), patch_box.upper(d));
        }
    }
    trimmed_bdry_box = BoundaryBoxNd(bdry_box.getBox() * bc_fill_box, /* codimension */ 1, location_index);
    bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
    return;
} // setupBcCoefBoxes

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
