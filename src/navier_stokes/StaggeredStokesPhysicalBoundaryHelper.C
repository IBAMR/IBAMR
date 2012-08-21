// Filename: StaggeredStokesPhysicalBoundaryHelper.C
// Created on 22 Jul 2008 by Boyce Griffith
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
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <SideIndex.h>
#include <SideData.h>
#include <tbox/MathUtilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesPhysicalBoundaryHelper::StaggeredStokesPhysicalBoundaryHelper()
    : d_hierarchy(NULL),
      d_dirichlet_bdry_locs(),
      d_dirichlet_bdry_vals()
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
StaggeredStokesPhysicalBoundaryHelper::enforceDirichletBcs(
    const int u_data_idx,
    const bool homogeneous_bcs,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            Pointer<SideData<NDIM,double> > u_data = patch->getPatchData(u_data_idx);

            // Lookup the boundary fill boxes.
            const Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();

            // Enforce homogeneous or inhomogeneous Dirichlet boundary
            // conditions.
            const std::vector<Pointer<ArrayData<NDIM,bool  > > >& dirichlet_bdry_locs = d_dirichlet_bdry_locs[ln].find(patch_num)->second;
            const std::vector<Pointer<ArrayData<NDIM,double> > >& dirichlet_bdry_vals = d_dirichlet_bdry_vals[ln].find(patch_num)->second;
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index   = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const Box<NDIM>& bc_coef_box = dirichlet_bdry_locs[n]->getBox();
                ArrayData<NDIM,bool  >& bdry_locs_data = *dirichlet_bdry_locs[n];
                ArrayData<NDIM,double>& bdry_vals_data = *dirichlet_bdry_vals[n];
                for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const Index<NDIM>& i = it();
                    if (bdry_locs_data(i,0))
                    {
                        (*u_data)(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower)) = homogeneous_bcs ? 0.0 : bdry_vals_data(i,0);
                    }
                }
            }
        }
    }
    return;
}// enforceDirichletBcs

void
StaggeredStokesPhysicalBoundaryHelper::copyDataAtDirichletBoundaries(
    const int u_out_data_idx,
    const int u_in_data_idx,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln); ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            Pointer<SideData<NDIM,double> > u_out_data = patch->getPatchData(u_out_data_idx);
            Pointer<SideData<NDIM,double> >  u_in_data = patch->getPatchData( u_in_data_idx);

            // Lookup the boundary fill boxes.
            const Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();

            // Copy data at Dirichlet boundaries.
            const std::vector<Pointer<ArrayData<NDIM,bool> > >& dirichlet_bdry_locs = d_dirichlet_bdry_locs[ln].find(patch_num)->second;
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index   = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const Box<NDIM>& bc_coef_box = dirichlet_bdry_locs[n]->getBox();
                ArrayData<NDIM,bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
                for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const Index<NDIM>& i = it();
                    if (bdry_locs_data(i,0))
                    {
                        (*u_out_data)(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower)) = (*u_in_data)(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower));
                    }
                }
            }
        }
    }
    return;
}// copyDataAtDirichletBoundaries

void
StaggeredStokesPhysicalBoundaryHelper::cacheBcCoefData(
    const int u_data_idx,
    const Pointer<Variable<NDIM> > u_var,
    std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const double fill_time,
    const IntVector<NDIM>& gcw_to_fill,
    const Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    // Indicate whether we are employing homogeneous or inhomogeneous boundary
    // conditions for all extended Robin BC coef strategy objects employed by
    // this object.
    for (std::vector<RobinBcCoefStrategy<NDIM>*>::iterator it = u_bc_coefs.begin(); it != u_bc_coefs.end(); ++it)
    {
        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(*it);
        if (extended_bc_coef != NULL)
        {
            extended_bc_coef->setTargetPatchDataIndex(u_data_idx);
            extended_bc_coef->setHomogeneousBc(false);
        }
    }

    if (!d_hierarchy.isNull()) clearBcCoefData();

    d_hierarchy = hierarchy;
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    d_physical_codim1_boxes.resize(finest_hier_level+1);
    d_dirichlet_bdry_locs  .resize(finest_hier_level+1);
    d_dirichlet_bdry_vals  .resize(finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

            // Compute the boundary fill boxes.
            Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln][patch_num];
            physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();

            // Compute the locations of the Dirichlet boundary.
            std::vector<Pointer<ArrayData<NDIM,bool> > >& dirichlet_bdry_locs = d_dirichlet_bdry_locs[ln][patch_num];
            dirichlet_bdry_locs.resize(n_physical_codim1_boxes,Pointer<ArrayData<NDIM,bool> >(NULL));
            std::vector<Pointer<ArrayData<NDIM,double> > >& dirichlet_bdry_vals = d_dirichlet_bdry_vals[ln][patch_num];
            dirichlet_bdry_vals.resize(n_physical_codim1_boxes,Pointer<ArrayData<NDIM,double> >(NULL));
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                const BoundaryBox<NDIM> trimmed_bdry_box(bdry_box.getBox()*bc_fill_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
                const unsigned int location_index   = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;

                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
                ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
                ArrayData<NDIM,double> gcoef_data(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data, false);
                Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data, false);
                Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(&gcoef_data, false);
                u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, u_var, *patch, trimmed_bdry_box, fill_time);

                dirichlet_bdry_locs[n] = new ArrayData<NDIM,bool  >(bc_coef_box, 1);
                dirichlet_bdry_vals[n] = new ArrayData<NDIM,double>(bc_coef_box, 1);
                ArrayData<NDIM,bool  >& bdry_locs_data = *dirichlet_bdry_locs[n];
                ArrayData<NDIM,double>& bdry_vals_data = *dirichlet_bdry_vals[n];
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
                    bdry_locs_data(i,0) = MathUtilities<double>::equalEps(alpha,1.0) || !MathUtilities<double>::equalEps(beta,1.0);
                    bdry_vals_data(i,0) = gamma;
                }
            }
        }
    }
    return;
}// cacheBcCoefData

void
StaggeredStokesPhysicalBoundaryHelper::clearBcCoefData()
{
    d_hierarchy.setNull();
    d_physical_codim1_boxes.clear();
    d_dirichlet_bdry_locs  .clear();
    d_dirichlet_bdry_vals  .clear();
    return;
}// clearBcCoefData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
