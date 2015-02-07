// Filename: StaggeredPhysicalBoundaryHelper.cpp
// Created on 22 Jul 2008 by Boyce Griffith
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
#include <map>
#include <ostream>
#include <utility>
#include <vector>

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/hier/Variable.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredPhysicalBoundaryHelper::StaggeredPhysicalBoundaryHelper() : d_hierarchy(NULL), d_dirichlet_bdry_locs()
{
    // intentionally blank
    return;
} // StaggeredPhysicalBoundaryHelper

StaggeredPhysicalBoundaryHelper::~StaggeredPhysicalBoundaryHelper()
{
    // intentionally blank
    return;
} // ~StaggeredPhysicalBoundaryHelper

void StaggeredPhysicalBoundaryHelper::copyDataAtDirichletBoundaries(const int u_out_data_idx,
                                                                    const int u_in_data_idx,
                                                                    const int coarsest_ln,
                                                                    const int finest_ln) const
{
    TBOX_ASSERT(d_hierarchy);
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p(level); p; p++)
        {
            Pointer<Patch > patch = level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                Pointer<SideData<double> > u_out_data = patch->getPatchData(u_out_data_idx);
                Pointer<SideData<double> > u_in_data = patch->getPatchData(u_in_data_idx);
                copyDataAtDirichletBoundaries(u_out_data, u_in_data, patch);
            }
        }
    }
    return;
} // copyDataAtDirichletBoundaries

void StaggeredPhysicalBoundaryHelper::copyDataAtDirichletBoundaries(Pointer<SideData<double> > u_out_data,
                                                                    Pointer<SideData<double> > u_in_data,
                                                                    Pointer<Patch > patch) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    const Array<BoundaryBox >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayData<bool> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        const ArrayData<bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (Box::Iterator it(bdry_locs_data.getBox()); it; it++)
        {
            const Index& i = it();
            if (bdry_locs_data(i, 0))
            {
                (*u_out_data)(SideIndex(i, bdry_normal_axis, SideIndex::Lower)) =
                    (*u_in_data)(SideIndex(i, bdry_normal_axis, SideIndex::Lower));
            }
        }
    }
    return;
} // copyDataAtDirichletBoundaries

void StaggeredPhysicalBoundaryHelper::setupMaskingFunction(const int mask_data_idx,
                                                           const int coarsest_ln,
                                                           const int finest_ln) const
{
    TBOX_ASSERT(d_hierarchy);
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p(level); p; p++)
        {
            Pointer<Patch > patch = level->getPatch(p());
            Pointer<SideData<int> > mask_data = patch->getPatchData(mask_data_idx);
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

void StaggeredPhysicalBoundaryHelper::setupMaskingFunction(Pointer<SideData<int> > mask_data,
                                                           Pointer<Patch > patch) const
{
    mask_data->fillAll(0);
    if (patch->getPatchGeometry()->getTouchesRegularBoundary()) return;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    const Array<BoundaryBox >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayData<bool> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        const ArrayData<bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (Box::Iterator it(bdry_locs_data.getBox()); it; it++)
        {
            const Index& i = it();
            if (bdry_locs_data(i, 0)) (*mask_data)(SideIndex(i, bdry_normal_axis, SideIndex::Lower)) = 1;
        }
    }
    return;
} // setupMaskingFunction

bool StaggeredPhysicalBoundaryHelper::patchTouchesDirichletBoundary(Pointer<Patch > patch) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return false;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        if (patchTouchesDirichletBoundaryAxis(patch, axis)) return true;
    }
    return false;
} // patchTouchesDirichletBoundary

bool StaggeredPhysicalBoundaryHelper::patchTouchesDirichletBoundaryAxis(Pointer<Patch > patch,
                                                                        const unsigned int axis) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return false;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    const Array<BoundaryBox >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayData<bool> > >& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const unsigned int bdry_normal_axis = physical_codim1_boxes[n].getLocationIndex() / 2;
        if (bdry_normal_axis == axis)
        {
            const ArrayData<bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
            for (Box::Iterator it(bdry_locs_data.getBox()); it; it++)
            {
                if (bdry_locs_data(it(), 0)) return true;
            }
        }
    }
    return false;
} // patchTouchesDirichletBoundaryAxis

void StaggeredPhysicalBoundaryHelper::cacheBcCoefData(const std::vector<RobinBcCoefStrategy*>& u_bc_coefs,
                                                      const double fill_time,
                                                      const Pointer<PatchHierarchy > hierarchy)
{
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
    TBOX_ASSERT(hierarchy);
    if (d_hierarchy) clearBcCoefData();

    // Cache boundary values.
    d_hierarchy = hierarchy;
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    d_physical_codim1_boxes.resize(finest_hier_level + 1);
    d_dirichlet_bdry_locs.resize(finest_hier_level + 1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        Pointer<PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch > patch = level->getPatch(patch_num);
            Pointer<CartesianPatchGeometry > pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Array<BoundaryBox >& physical_codim1_boxes = d_physical_codim1_boxes[ln][patch_num];
                physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                std::vector<Pointer<ArrayData<bool> > >& dirichlet_bdry_locs =
                    d_dirichlet_bdry_locs[ln][patch_num];
                dirichlet_bdry_locs.resize(n_physical_codim1_boxes);
                Box bc_coef_box;
                BoundaryBox trimmed_bdry_box;
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox& bdry_box = physical_codim1_boxes[n];
                    StaggeredPhysicalBoundaryHelper::setupBcCoefBoxes(bc_coef_box, trimmed_bdry_box, bdry_box, patch);
                    const unsigned int bdry_normal_axis = bdry_box.getLocationIndex() / 2;
                    Pointer<ArrayData<double> > acoef_data = new ArrayData<double>(bc_coef_box, 1);
                    Pointer<ArrayData<double> > bcoef_data = new ArrayData<double>(bc_coef_box, 1);
                    Pointer<ArrayData<double> > gcoef_data;
                    u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data,
                                                             bcoef_data,
                                                             gcoef_data,
                                                             Pointer<Variable >(),
                                                             *patch,
                                                             trimmed_bdry_box,
                                                             fill_time);
                    dirichlet_bdry_locs[n] = new ArrayData<bool>(bc_coef_box, 1);
                    ArrayData<bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
                    for (Box::Iterator it(bc_coef_box); it; it++)
                    {
                        const Index& i = it();
                        const double& alpha = (*acoef_data)(i, 0);
                        const double& beta = (*bcoef_data)(i, 0);
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha + beta, 1.0));
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha, 1.0) ||
                                    MathUtilities<double>::equalEps(beta, 1.0));
                        bdry_locs_data(i, 0) = MathUtilities<double>::equalEps(alpha, 1.0) &&
                                               (beta == 0.0 || MathUtilities<double>::equalEps(beta, 0.0));
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

void StaggeredPhysicalBoundaryHelper::clearBcCoefData()
{
    d_hierarchy.setNull();
    d_physical_codim1_boxes.clear();
    d_dirichlet_bdry_locs.clear();
    return;
} // clearBcCoefData

/////////////////////////////// PROTECTED ////////////////////////////////////

void StaggeredPhysicalBoundaryHelper::setupBcCoefBoxes(Box& bc_coef_box,
                                                       BoundaryBox& trimmed_bdry_box,
                                                       const BoundaryBox& bdry_box,
                                                       Pointer<Patch > patch)
{
    Pointer<PatchGeometry > pgeom = patch->getPatchGeometry();
    const Box& patch_box = patch->getBox();
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    Box bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, /* gcw_to_fill */ IntVector::getOne(DIM));
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d != bdry_normal_axis)
        {
            bc_fill_box.lower(d) = std::max(bc_fill_box.lower(d), patch_box.lower(d));
            bc_fill_box.upper(d) = std::min(bc_fill_box.upper(d), patch_box.upper(d));
        }
    }
    trimmed_bdry_box = BoundaryBox(bdry_box.getBox() * bc_fill_box, /* codimension */ 1, location_index);
    bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
    return;
} // setupBcCoefBoxes

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
