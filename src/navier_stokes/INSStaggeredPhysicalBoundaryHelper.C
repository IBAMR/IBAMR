// Filename: INSStaggeredPhysicalBoundaryHelper.C
// Created on 22 Jul 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "INSStaggeredPhysicalBoundaryHelper.h"

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

INSStaggeredPhysicalBoundaryHelper::INSStaggeredPhysicalBoundaryHelper()
    : d_hierarchy(NULL),
      d_dirichlet_bdry_locs(),
      d_dirichlet_bdry_vals()
{
    // intentionally blank
    return;
}// INSStaggeredPhysicalBoundaryHelper

INSStaggeredPhysicalBoundaryHelper::~INSStaggeredPhysicalBoundaryHelper()
{
    // intentionally blank
    return;
}// ~INSStaggeredPhysicalBoundaryHelper

void
INSStaggeredPhysicalBoundaryHelper::zeroValuesAtDirichletBoundaries(
    const int patch_data_idx,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln);
         ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            Pointer<SideData<NDIM,double> > patch_data = patch->getPatchData(patch_data_idx);

            // Compute the boundary fill boxes.
            const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();

            // Compute the locations of the Dirichlet boundary.
            const std::vector<Pointer<ArrayData<NDIM,bool> > >& dirichlet_bdry_locs = (*d_dirichlet_bdry_locs[ln].find(patch_num)).second;
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const int location_index   = bdry_box.getLocationIndex();
                const int bdry_normal_axis = location_index / 2;
                const Box<NDIM>& bc_coef_box = dirichlet_bdry_locs[n]->getBox();
                ArrayData<NDIM,bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
                for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const Index<NDIM>& i = it();
                    if (bdry_locs_data(i,0))
                    {
                        const SideIndex<NDIM> i_s(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                        (*patch_data)(i_s) = 0.0;
                    }
                }
            }
        }
    }
    return;
}// zeroValuesAtDirichletBoundaries

void
INSStaggeredPhysicalBoundaryHelper::resetValuesAtDirichletBoundaries(
    const int patch_data_idx,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == -1 ? 0 : coarsest_ln);
         ln <= (finest_ln == -1 ? finest_hier_level : finest_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch<NDIM> > patch = level->getPatch(patch_num);
            Pointer<SideData<NDIM,double> > patch_data = patch->getPatchData(patch_data_idx);

            // Compute the boundary fill boxes.
            const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();

            // Compute the locations of the Dirichlet boundary.
            const std::vector<Pointer<ArrayData<NDIM,bool> > >& dirichlet_bdry_locs = (*d_dirichlet_bdry_locs[ln].find(patch_num)).second;
            const std::vector<Pointer<ArrayData<NDIM,double> > >& dirichlet_bdry_vals = (*d_dirichlet_bdry_vals[ln].find(patch_num)).second;
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const int location_index   = bdry_box.getLocationIndex();
                const int bdry_normal_axis = location_index / 2;
                const Box<NDIM>& bc_coef_box = dirichlet_bdry_locs[n]->getBox();
                ArrayData<NDIM,bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
                ArrayData<NDIM,double>& bdry_vals_data = *dirichlet_bdry_vals[n];
                for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const Index<NDIM>& i = it();
                    if (bdry_locs_data(i,0))
                    {
                        const SideIndex<NDIM> i_s(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                        (*patch_data)(i_s) = bdry_vals_data(i,0);
                    }
                }
            }
        }
    }
    return;
}// resetValuesAtDirichletBoundaries

void
INSStaggeredPhysicalBoundaryHelper::cacheBcCoefData(
    const int u_idx,
    const Pointer<Variable<NDIM> >& u_var,
    std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const double fill_time,
    const IntVector<NDIM>& gcw_to_fill,
    const Pointer<PatchHierarchy<NDIM> >& hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    // Indicate whether we are employing homogeneous or inhomogeneous boundary
    // conditions for all extended Robin BC coef strategy objects employed by
    // this object.
    for (std::vector<RobinBcCoefStrategy<NDIM>*>::iterator it = u_bc_coefs.begin();
         it != u_bc_coefs.end(); ++it)
    {
        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(*it);
        if (extended_bc_coef != NULL)
        {
            extended_bc_coef->setTargetPatchDataIndex(u_idx);
            extended_bc_coef->setHomogeneousBc(false);
        }
    }

    if (!d_hierarchy.isNull()) clearBcCoefData();
    d_hierarchy = hierarchy;
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    d_dirichlet_bdry_locs.resize(finest_hier_level+1);
    d_dirichlet_bdry_vals.resize(finest_hier_level+1);
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
            const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
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
                const int location_index   = bdry_box.getLocationIndex();
                const int bdry_normal_axis = location_index / 2;

                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, u_var, *patch, trimmed_bdry_box, fill_time);

                dirichlet_bdry_locs[n] = new ArrayData<NDIM,bool>(bc_coef_box, 1);
                ArrayData<NDIM,bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
                dirichlet_bdry_vals[n] = new ArrayData<NDIM,double>(bc_coef_box, 1);
                ArrayData<NDIM,double>& bdry_vals_data = *dirichlet_bdry_vals[n];
                for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const Index<NDIM>& i = it();
                    const double& alpha = (*acoef_data)(i,0);
                    const double& beta  = (*bcoef_data)(i,0);
                    const double& gamma = (*gcoef_data)(i,0);
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
INSStaggeredPhysicalBoundaryHelper::clearBcCoefData()
{
    d_hierarchy.setNull();
    d_dirichlet_bdry_locs.clear();
    d_dirichlet_bdry_vals.clear();
    return;
}// clearBcCoefData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::INSStaggeredPhysicalBoundaryHelper>;

//////////////////////////////////////////////////////////////////////////////
