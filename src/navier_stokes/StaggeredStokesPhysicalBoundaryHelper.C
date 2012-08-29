// Filename: StaggeredStokesPhysicalBoundaryHelper.C
// Created on 28 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "StaggeredStokesPhysicalBoundaryHelper.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesPhysicalBoundaryHelper::StaggeredStokesPhysicalBoundaryHelper()
    : d_neumann_bdry_locs(),
      d_neumann_bdry_vals(),
      d_problem_coefs(NULL),
      d_u_current_data_idx(-1)
{
    // intentionally blank
    return;
}// StaggeredStokesPhysicalBoundaryHelper

StaggeredStokesPhysicalBoundaryHelper::~StaggeredStokesPhysicalBoundaryHelper()
{
    // intentionally blank
    return;
}// ~StaggeredStokesPhysicalBoundaryHelper

void
StaggeredStokesPhysicalBoundaryHelper::enforceDivergenceFreeConditionAtBoundary(
    const int u_data_idx,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM,double> > u_data = patch->getPatchData(u_data_idx);
            enforceDivergenceFreeConditionAtBoundary(u_data, patch);
        }
    }
    return;
}// enforceDivergenceFreeConditionAtBoundary

void
StaggeredStokesPhysicalBoundaryHelper::enforceDivergenceFreeConditionAtBoundary(
    Pointer<SideData<NDIM,double> > u_data,
    Pointer<Patch<NDIM> > patch) const
{
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayData<NDIM,bool> > >& neumann_bdry_locs = d_neumann_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box   = physical_codim1_boxes[n];
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower                 = location_index % 2 == 0;
        const Box<NDIM>& bc_coef_box        = neumann_bdry_locs[n]->getBox();
        const ArrayData<NDIM,bool>& bdry_locs_data = *neumann_bdry_locs[n];
        for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
        {
            const Index<NDIM>& i = it();
            if (bdry_locs_data(i,0))
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

                // Determine the ghost cell value so that the divergence of the
                // velocity field is zero in the ghost cell.
                double div_u_g = 0.0;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    if (axis == bdry_normal_axis)
                    {
                        if (is_lower)
                        {
                            div_u_g += +(*u_data)(SideIndex<NDIM>(i_g,axis,SideIndex<NDIM>::Upper))/dx[axis];
                        }
                        else
                        {
                            div_u_g += -(*u_data)(SideIndex<NDIM>(i_g,axis,SideIndex<NDIM>::Lower))/dx[axis];
                        }
                    }
                    else
                    {
                        div_u_g += ((*u_data)(SideIndex<NDIM>(i_g,axis,SideIndex<NDIM>::Upper)) - (*u_data)(SideIndex<NDIM>(i_g,axis,SideIndex<NDIM>::Lower)))/dx[axis];
                    }
                }
                if (is_lower)
                {
                    (*u_data)(SideIndex<NDIM>(i_g,bdry_normal_axis,SideIndex<NDIM>::Lower)) = +div_u_g*dx[bdry_normal_axis];
                }
                else
                {
                    (*u_data)(SideIndex<NDIM>(i_g,bdry_normal_axis,SideIndex<NDIM>::Upper)) = -div_u_g*dx[bdry_normal_axis];
                }
            }
        }
    }
    return;
}// enforceDivergenceFreeConditionAtBoundary

void
StaggeredStokesPhysicalBoundaryHelper::enforceNormalTractionBoundaryConditions(
    const int p_data_idx,
    const int u_new_data_idx,
    const bool homogeneous_bcs,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM,double> > p_data = patch->getPatchData(p_data_idx);
            Pointer<SideData<NDIM,double> > u_new_data = patch->getPatchData(u_new_data_idx);
            enforceNormalTractionBoundaryConditions(p_data, u_new_data, homogeneous_bcs, patch);
        }
    }
    return;
}// enforceNormalTractionBoundaryConditions

void
StaggeredStokesPhysicalBoundaryHelper::enforceNormalTractionBoundaryConditions(
    Pointer<CellData<NDIM,double> > p_data,
    Pointer<SideData<NDIM,double> > u_new_data,
    const bool homogeneous_bcs,
    Pointer<Patch<NDIM> > patch) const
{
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    const double mu = d_problem_coefs->getMu();
    Pointer<SideData<NDIM,double> > u_current_data = patch->getPatchData(d_u_current_data_idx);
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayData<NDIM,bool  > > >& neumann_bdry_locs = d_neumann_bdry_locs[ln].find(patch_num)->second;
    const std::vector<Pointer<ArrayData<NDIM,double> > >& neumann_bdry_vals = d_neumann_bdry_vals[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box   = physical_codim1_boxes[n];
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower                 = location_index % 2 == 0;
        const Box<NDIM>& bc_coef_box        = neumann_bdry_locs[n]->getBox();
        const ArrayData<NDIM,bool  >& bdry_locs_data = *neumann_bdry_locs[n];
        const ArrayData<NDIM,double>& bdry_vals_data = *neumann_bdry_vals[n];
        for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
        {
            const Index<NDIM>& i = it();
            if (bdry_locs_data(i,0))
            {
                // Place i_i in the interior cell abutting the boundary, and
                // place i_g in the ghost    cell abutting the boundary.
                Index<NDIM> i_i(i), i_g(i);
                if (is_lower)
                {
                    i_g(bdry_normal_axis) -= 1;
                }
                else
                {
                    i_i(bdry_normal_axis) -= 1;
                }

                // The boundary condition is -p + 2*mu*du_n/dx_n = g.
                const SideIndex<NDIM> i_s_upper(is_lower ? i_i : i_g, bdry_normal_axis, SideIndex<NDIM>::Upper);
                const SideIndex<NDIM> i_s_lower(is_lower ? i_g : i_i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                const double du_norm_current_dx_norm = ((*u_current_data)(i_s_upper)-(*u_current_data)(i_s_lower))/(2.0*dx[bdry_normal_axis]);
                const double du_norm_new_dx_norm     = ((*u_new_data    )(i_s_upper)-(*u_new_data    )(i_s_lower))/(2.0*dx[bdry_normal_axis]);
                const double gamma = (homogeneous_bcs ? 0.0 : -bdry_vals_data(i,0) + mu*du_norm_current_dx_norm) + mu*du_norm_new_dx_norm;
                (*p_data)(i_g) = 2.0*gamma-(*p_data)(i_i);
            }
        }
    }
    return;
}// enforceNormalTractionBoundaryConditions

void
StaggeredStokesPhysicalBoundaryHelper::setStokesSpecifications(
    const StokesSpecifications* problem_coefs)
{
    d_problem_coefs = problem_coefs;
    return;
}// setStokesSpecifications

void
StaggeredStokesPhysicalBoundaryHelper::setCurrentVelocityDataIndex(
    const int u_current_data_idx)
{
    d_u_current_data_idx = u_current_data_idx;
    return;
}// setCurrentVelocityDataIndex

void
StaggeredStokesPhysicalBoundaryHelper::cacheBcCoefData(
    const int u_data_idx,
    const Pointer<Variable<NDIM> > u_var,
    std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const double velocity_bc_fill_time,
    const double traction_bc_fill_time,
    const IntVector<NDIM>& gcw_to_fill,
    const Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    StaggeredPhysicalBoundaryHelper::cacheBcCoefData(u_data_idx, u_var, u_bc_coefs, velocity_bc_fill_time, gcw_to_fill, hierarchy);
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    d_neumann_bdry_locs.resize(finest_hier_level+1);
    d_neumann_bdry_vals.resize(finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

            // Look up the boundary fill boxes.
            const Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln][patch_num];
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();

            // Compute the locations of the Neumann boundaries.
            std::vector<Pointer<ArrayData<NDIM,bool  > > >& neumann_bdry_locs = d_neumann_bdry_locs[ln][patch_num];
            std::vector<Pointer<ArrayData<NDIM,double> > >& neumann_bdry_vals = d_neumann_bdry_vals[ln][patch_num];
            neumann_bdry_locs.resize(n_physical_codim1_boxes,Pointer<ArrayData<NDIM,bool  > >(NULL));
            neumann_bdry_vals.resize(n_physical_codim1_boxes,Pointer<ArrayData<NDIM,double> >(NULL));
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index   = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (d != bdry_normal_axis)
                    {
                        bc_fill_box.lower(d) = std::max(bc_fill_box.lower(d), patch_box.lower(d));
                        bc_fill_box.upper(d) = std::min(bc_fill_box.upper(d), patch_box.upper(d));
                    }
                }
                const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());

                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
                ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
                ArrayData<NDIM,double> gcoef_data(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data, false);
                Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data, false);
                Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(&gcoef_data, false);
                u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, u_var, *patch, trimmed_bdry_box, traction_bc_fill_time);

                neumann_bdry_locs[n] = new ArrayData<NDIM,bool  >(bc_coef_box, 1);
                neumann_bdry_vals[n] = new ArrayData<NDIM,double>(bc_coef_box, 1);
                ArrayData<NDIM,bool  >& neumann_bdry_locs_data = *neumann_bdry_locs[n];
                ArrayData<NDIM,double>& neumann_bdry_vals_data = *neumann_bdry_vals[n];
                for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const Index<NDIM>& i = it();
                    const double& alpha = acoef_data(i,0);
                    const double& beta  = bcoef_data(i,0);
                    const double& gamma = gcoef_data(i,0);
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(MathUtilities<double>::equalEps(alpha+beta,1.0));
                    TBOX_ASSERT(MathUtilities<double>::equalEps(alpha,1.0) || MathUtilities<double>::equalEps(beta,1.0));
#endif
                    neumann_bdry_locs_data(i,0) = MathUtilities<double>::equalEps(beta,1.0) && (alpha == 0.0 || MathUtilities<double>::equalEps(alpha,0.0));
                    neumann_bdry_vals_data(i,0) = neumann_bdry_locs_data(i,0) ? gamma : std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }
    return;
}// cacheBcCoefData

void
StaggeredStokesPhysicalBoundaryHelper::clearBcCoefData()
{
    StaggeredPhysicalBoundaryHelper::clearBcCoefData();
    d_neumann_bdry_locs.clear();
    d_neumann_bdry_vals.clear();
    return;
}// clearBcCoefData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
