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
#include <ibamr/StokesBcCoefStrategy.h>
#include <ibamr/namespaces.h>

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
    : d_problem_coefs(NULL),
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
StaggeredStokesPhysicalBoundaryHelper::enforceNormalVelocityBoundaryConditions(
    const int u_data_idx,
    const int p_data_idx,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const double fill_time,
    const bool homogeneous_bc,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
    TBOX_ASSERT(d_hierarchy);
#endif
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
                Pointer<SideData<NDIM,double> > u_data = patch->getPatchData(u_data_idx);
                Box<NDIM> bc_coef_box;
                BoundaryBox<NDIM> trimmed_bdry_box;
                const Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    setupBcCoefBoxes(bc_coef_box, trimmed_bdry_box, bdry_box, patch);
                    const unsigned int bdry_normal_axis = bdry_box.getLocationIndex() / 2;
                    Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
                    if (extended_bc_coef)
                    {
                        extended_bc_coef->clearTargetPatchDataIndex();
                        extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                    }
                    StokesBcCoefStrategy* stokes_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
                    if (stokes_bc_coef)
                    {
                        stokes_bc_coef->setTargetVelocityPatchDataIndex(u_data_idx);
                        stokes_bc_coef->setTargetPressurePatchDataIndex(p_data_idx);
                    }
                    u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM> >(), *patch, trimmed_bdry_box, fill_time);
                    if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
                    if (stokes_bc_coef)
                    {
                        stokes_bc_coef->clearTargetVelocityPatchDataIndex();
                        stokes_bc_coef->clearTargetPressurePatchDataIndex();
                    }
                    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                    {
                        const Index<NDIM>& i = it();
                        const double& alpha = (*acoef_data)(i,0);
                        const double  gamma = homogeneous_bc && !extended_bc_coef ? 0.0 : (*gcoef_data)(i,0);
#ifdef DEBUG_CHECK_ASSERTIONS
                        const double& beta  = (*bcoef_data)(i,0);
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha+beta,1.0));
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha,1.0) || MathUtilities<double>::equalEps(beta,1.0));
#endif
                        if (MathUtilities<double>::equalEps(alpha,1.0)) (*u_data)(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower)) = gamma;
                    }
                }
            }
        }
    }
    return;
}// enforceNormalVelocityBoundaryConditions
void
StaggeredStokesPhysicalBoundaryHelper::enforceNormalTractionBoundaryConditions(
    const int u_data_idx,
    const int p_data_idx,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const double fill_time,
    const bool homogeneous_bc,
    const int coarsest_ln,
    const int finest_ln) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
    TBOX_ASSERT(d_hierarchy);
#endif
    const double mu = d_problem_coefs->getMu();
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
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            if (pgeom->getTouchesRegularBoundary())
            {
                Pointer<CellData<NDIM,double> > p_half_data    = patch->getPatchData(p_data_idx);
                Pointer<SideData<NDIM,double> > u_new_data     = patch->getPatchData(u_data_idx);
                Pointer<SideData<NDIM,double> > u_current_data = patch->getPatchData(d_u_current_data_idx);
                Box<NDIM> bc_coef_box;
                BoundaryBox<NDIM> trimmed_bdry_box;
                const Array<BoundaryBox<NDIM> >& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    setupBcCoefBoxes(bc_coef_box, trimmed_bdry_box, bdry_box, patch);
                    const unsigned int location_index   = bdry_box.getLocationIndex();
                    const unsigned int bdry_normal_axis = location_index / 2;
                    const bool is_lower                 = location_index % 2 == 0;
                    Pointer<ArrayData<NDIM,double> > acoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > bcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > gcoef_data = new ArrayData<NDIM,double>(bc_coef_box, 1);
                    ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
                    if (extended_bc_coef)
                    {
                        extended_bc_coef->clearTargetPatchDataIndex();
                        extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                    }
                    StokesBcCoefStrategy* stokes_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
                    if (stokes_bc_coef)
                    {
                        stokes_bc_coef->setTargetVelocityPatchDataIndex(u_data_idx);
                        stokes_bc_coef->setTargetPressurePatchDataIndex(p_data_idx);
                    }
                    u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM> >(), *patch, trimmed_bdry_box, fill_time);
                    if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
                    if (stokes_bc_coef)
                    {
                        stokes_bc_coef->clearTargetVelocityPatchDataIndex();
                        stokes_bc_coef->clearTargetPressurePatchDataIndex();
                    }
                    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                    {
                        const Index<NDIM>& i = it();
                        const double& beta  = (*bcoef_data)(i,0);
                        const double  gamma = homogeneous_bc && !extended_bc_coef ? 0.0 : (*gcoef_data)(i,0);
#ifdef DEBUG_CHECK_ASSERTIONS
                        const double& alpha = (*acoef_data)(i,0);
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha+beta,1.0));
                        TBOX_ASSERT(MathUtilities<double>::equalEps(alpha,1.0) || MathUtilities<double>::equalEps(beta,1.0));
#endif
                        if (MathUtilities<double>::equalEps(beta,1.0))
                        {
                            // Place i_i in the interior cell abutting the
                            // boundary, and place i_g in the ghost cell
                            // abutting the boundary.
                            Index<NDIM> i_i(i), i_g(i);
                            if (is_lower)
                            {
                                i_g(bdry_normal_axis) -= 1;
                            }
                            else
                            {
                                i_i(bdry_normal_axis) -= 1;
                            }

                            // The boundary condition is -p + 2*mu*du_n/dx_n =
                            // g.
                            static const int NVALS = 3;
                            double u_current[NVALS], u_new[NVALS];
                            SideIndex<NDIM> i_s(i_i, bdry_normal_axis, is_lower ? SideIndex<NDIM>::Lower : SideIndex<NDIM>::Upper);
                            for (int k = 0; k < NVALS; ++k, i_s(bdry_normal_axis) += (is_lower ? 1 : -1))
                            {
                                u_current[k] = (*u_current_data)(i_s);
                                u_new    [k] = (*u_new_data    )(i_s);
                            }
                            const double h = dx[bdry_normal_axis];
                            const double du_norm_current_dx_norm = (is_lower ? +1.0 : -1.0)*(2.0*u_current[1]-1.5*u_current[0]-0.5*u_current[2])/h;
                            const double du_norm_new_dx_norm     = (is_lower ? +1.0 : -1.0)*(2.0*u_new    [1]-1.5*u_new    [0]-0.5*u_new    [2])/h;
                            const double g = -gamma + (homogeneous_bc ? 0.0 : mu*du_norm_current_dx_norm) + mu*du_norm_new_dx_norm;
                            (*p_half_data)(i_g) = 2.0*g-(*p_half_data)(i_i);
                        }
                    }
                }
            }
        }
    }
    return;
}// enforceNormalTractionBoundaryConditions

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

                // Determine the ghost cell value so that the divergence of the
                // velocity field is zero in the ghost cell.
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
    return;
}// enforceDivergenceFreeConditionAtBoundary

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

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
