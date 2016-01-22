// Filename: StaggeredStokesPhysicalBoundaryHelper.cpp
// Created on 28 Aug 2012 by Boyce Griffith
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
#include <map>
#include <ostream>
#include <utility>
#include <vector>

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideIndex.h"
#include "Variable.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StokesBcCoefStrategy.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "tbox/Array.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesPhysicalBoundaryHelper::StaggeredStokesPhysicalBoundaryHelper()
{
    // intentionally blank
    return;
} // StaggeredStokesPhysicalBoundaryHelper

StaggeredStokesPhysicalBoundaryHelper::~StaggeredStokesPhysicalBoundaryHelper()
{
    // intentionally blank
    return;
} // ~StaggeredStokesPhysicalBoundaryHelper

void
StaggeredStokesPhysicalBoundaryHelper::enforceNormalVelocityBoundaryConditions(
    const int u_data_idx,
    const int p_data_idx,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const double fill_time,
    const bool homogeneous_bc,
    const int coarsest_ln,
    const int finest_ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
    TBOX_ASSERT(d_hierarchy);
#endif
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        u_bc_coefs, /*p_bc_coef*/ NULL, u_data_idx, p_data_idx, homogeneous_bc);
    std::vector<int> target_data_idxs(2);
    target_data_idxs[0] = u_data_idx;
    target_data_idxs[1] = p_data_idx;
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_data_idx);
                Box<NDIM> bc_coef_box;
                BoundaryBox<NDIM> trimmed_bdry_box;
                const Array<BoundaryBox<NDIM> >& physical_codim1_boxes =
                    d_physical_codim1_boxes[ln].find(patch_num)->second;
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    StaggeredPhysicalBoundaryHelper::setupBcCoefBoxes(bc_coef_box, trimmed_bdry_box, bdry_box, patch);
                    const unsigned int bdry_normal_axis = bdry_box.getLocationIndex() / 2;
                    Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data,
                                                             bcoef_data,
                                                             gcoef_data,
                                                             Pointer<Variable<NDIM> >(),
                                                             *patch,
                                                             trimmed_bdry_box,
                                                             fill_time);
                    ExtendedRobinBcCoefStrategy* const extended_bc_coef =
                        dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
                    if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
                    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                    {
                        const Index<NDIM>& i = it();
                        const double& alpha = (*acoef_data)(i, 0);
                        const double gamma = homogeneous_bc && !extended_bc_coef ? 0.0 : (*gcoef_data)(i, 0);
#if !defined(NDEBUG)
                        const double& beta = (*bcoef_data)(i, 0);
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha + beta, 1.0));
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha, 1.0) ||
                                    MathUtilities<double>::equalEps(beta, 1.0));
#endif
                        if (MathUtilities<double>::equalEps(alpha, 1.0))
                            (*u_data)(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower)) = gamma;
                    }
                }
            }
        }
    }
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(u_bc_coefs, /*p_bc_coef*/ NULL);
    return;
} // enforceNormalVelocityBoundaryConditions
#if 0
void
StaggeredStokesPhysicalBoundaryHelper::enforceDivergenceFreeConditionAtBoundary(
    const int u_data_idx,
    const int coarsest_ln,
    const int finest_ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM,double> > u_data = patch->getPatchData(u_data_idx);
                enforceDivergenceFreeConditionAtBoundary(u_data, patch);
            }
        }
    }
    return;
}// enforceDivergenceFreeConditionAtBoundary

void
StaggeredStokesPhysicalBoundaryHelper::enforceDivergenceFreeConditionAtBoundary(
    Pointer<SideData<NDIM,double> > u_data,
    Pointer<Patch<NDIM> > patch) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayData<NDIM,bool> > >& dirichlet_bdry_locs = d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box   = physical_codim1_boxes[n];
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower                 = location_index % 2 == 0;
        const Box<NDIM>& bc_coef_box        = dirichlet_bdry_locs[n]->getBox();
        const ArrayData<NDIM,bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
        {
            const Index<NDIM>& i = it();
            if (!bdry_locs_data(i,0))
            {
                // Place i_g in the ghost cell abutting the boundary.
                Index<NDIM> i_g = i;
                if (is_lower)
                {
                    i_g(bdry_normal_axis) -= 1;
                }
                else
                {
                    // intentionally blank
                }

                // Work out from the physical boundary to fill the ghost cell
                // values so that the velocity field satisfies the discrete
                // divergence-free condition.
                for (int k = 0; k < u_data->getGhostCellWidth()(bdry_normal_axis); ++k, i_g(bdry_normal_axis) += (is_lower ? -1 : +1))
                {
                    // Determine the ghost cell value so that the divergence of
                    // the velocity field is zero in the ghost cell.
                    SideIndex<NDIM> i_g_s(i_g, bdry_normal_axis, is_lower ? SideIndex<NDIM>::Lower : SideIndex<NDIM>::Upper);
                    (*u_data)(i_g_s) = 0.0;
                    double div_u_g = 0.0;
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        const SideIndex<NDIM> i_g_s_upper(i_g,axis,SideIndex<NDIM>::Upper);
                        const SideIndex<NDIM> i_g_s_lower(i_g,axis,SideIndex<NDIM>::Lower);
                        div_u_g += ((*u_data)(i_g_s_upper)-(*u_data)(i_g_s_lower))*dx[bdry_normal_axis]/dx[axis];
                    }
                    (*u_data)(i_g_s) = (is_lower ? +1.0 : -1.0)*div_u_g;
                }
            }
        }
    }
    return;
}// enforceDivergenceFreeConditionAtBoundary
#endif
void
StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                                          RobinBcCoefStrategy<NDIM>* p_bc_coef,
                                                          int u_target_data_idx,
                                                          int p_target_data_idx,
                                                          bool homogeneous_bc)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
#endif
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ExtendedRobinBcCoefStrategy* extended_u_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[d]);
        if (extended_u_bc_coef)
        {
            extended_u_bc_coef->clearTargetPatchDataIndex();
            extended_u_bc_coef->setHomogeneousBc(homogeneous_bc);
        }
        StokesBcCoefStrategy* stokes_u_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(u_bc_coefs[d]);
        if (stokes_u_bc_coef)
        {
            stokes_u_bc_coef->setTargetVelocityPatchDataIndex(u_target_data_idx);
            stokes_u_bc_coef->setTargetPressurePatchDataIndex(p_target_data_idx);
        }
    }
    ExtendedRobinBcCoefStrategy* extended_p_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(p_bc_coef);
    if (extended_p_bc_coef)
    {
        extended_p_bc_coef->clearTargetPatchDataIndex();
        extended_p_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    StokesBcCoefStrategy* stokes_p_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(p_bc_coef);
    if (stokes_p_bc_coef)
    {
        stokes_p_bc_coef->setTargetVelocityPatchDataIndex(u_target_data_idx);
        stokes_p_bc_coef->setTargetPressurePatchDataIndex(p_target_data_idx);
    }
    return;
} // setupBcCoefObjects

void
StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                                          RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
#endif
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        StokesBcCoefStrategy* stokes_u_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(u_bc_coefs[d]);
        if (stokes_u_bc_coef)
        {
            stokes_u_bc_coef->clearTargetVelocityPatchDataIndex();
            stokes_u_bc_coef->clearTargetPressurePatchDataIndex();
        }
    }
    StokesBcCoefStrategy* stokes_p_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(p_bc_coef);
    if (stokes_p_bc_coef)
    {
        stokes_p_bc_coef->clearTargetVelocityPatchDataIndex();
        stokes_p_bc_coef->clearTargetPressurePatchDataIndex();
    }
    return;
} // resetBcCoefObjects

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
