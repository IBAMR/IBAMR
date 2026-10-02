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

#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/StokesBcCoefStrategy.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/StaggeredPhysicalBoundaryHelper.h>

#include <tbox/Array.h>
#include <tbox/MathUtilities.h>
#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>
#include <Variable.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

const short int StaggeredStokesPhysicalBoundaryHelper::NORMAL_TRACTION_BDRY = 0x100;
const short int StaggeredStokesPhysicalBoundaryHelper::NORMAL_VELOCITY_BDRY = 0x200;
const short int StaggeredStokesPhysicalBoundaryHelper::ALL_BDRY = 0x100 | 0x200;

/////////////////////////////// PUBLIC ///////////////////////////////////////

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
        u_bc_coefs, /*p_bc_coef*/ nullptr, u_data_idx, p_data_idx, homogeneous_bc);
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == IBTK::invalid_level_number ? 0 : coarsest_ln);
         ln <= (finest_ln == IBTK::invalid_level_number ? finest_hier_level : finest_ln);
         ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const int patch_num = p();
            Pointer<Patch<NDIM>> patch = level->getPatch(patch_num);
            Pointer<PatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_data_idx);
                Box<NDIM> bc_coef_box;
                BoundaryBox<NDIM> trimmed_bdry_box;
                const Array<BoundaryBox<NDIM>>& physical_codim1_boxes =
                    d_physical_codim1_boxes[ln].find(patch_num)->second;
                const int n_physical_codim1_boxes = physical_codim1_boxes.size();
                for (int n = 0; n < n_physical_codim1_boxes; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    StaggeredPhysicalBoundaryHelper::setupBcCoefBoxes(bc_coef_box, trimmed_bdry_box, bdry_box, patch);
                    const unsigned int bdry_normal_axis = bdry_box.getLocationIndex() / 2;
                    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double>> gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data,
                                                             bcoef_data,
                                                             gcoef_data,
                                                             Pointer<Variable<NDIM>>(),
                                                             *patch,
                                                             trimmed_bdry_box,
                                                             fill_time);
                    auto const extended_bc_coef =
                        dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
                    if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);
                    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                    {
                        const hier::Index<NDIM>& i = it();
                        const double& alpha = (*acoef_data)(i, 0);
                        const double gamma = homogeneous_bc && !extended_bc_coef ? 0.0 : (*gcoef_data)(i, 0);
#if !defined(NDEBUG)
                        const double& beta = (*bcoef_data)(i, 0);
                        TBOX_ASSERT(IBTK::rel_equal_eps(alpha + beta, 1.0));
                        TBOX_ASSERT(IBTK::rel_equal_eps(alpha, 1.0) || IBTK::rel_equal_eps(beta, 1.0));
#endif
                        if (IBTK::rel_equal_eps(alpha, 1.0))
                            (*u_data)(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower)) = gamma;
                    }
                }
            }
        }
    }
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(u_bc_coefs, /*p_bc_coef*/ nullptr);
    return;
} // enforceNormalVelocityBoundaryConditions

void
StaggeredStokesPhysicalBoundaryHelper::enforceDivergenceFreeConditionAtBoundary(const int u_data_idx,
                                                                                const int coarsest_ln,
                                                                                const int finest_ln,
                                                                                const short int bdry_tag) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == IBTK::invalid_level_number ? 0 : coarsest_ln);
         ln <= (finest_ln == IBTK::invalid_level_number ? finest_hier_level : finest_ln);
         ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_data_idx);
                enforceDivergenceFreeConditionAtBoundary(u_data, patch, bdry_tag);
            }
        }
    }
    return;
} // enforceDivergenceFreeConditionAtBoundary

void
StaggeredStokesPhysicalBoundaryHelper::enforceDivergenceFreeConditionAtBoundary(Pointer<SideData<NDIM, double>> u_data,
                                                                                Pointer<Patch<NDIM>> patch,
                                                                                const short int bdry_tag) const
{
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary()) return;
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Array<BoundaryBox<NDIM>>& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayData<NDIM, bool>>>& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower = location_index % 2 == 0;
        const Box<NDIM>& bc_coef_box = dirichlet_bdry_locs[n]->getBox();
        const ArrayData<NDIM, bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
        {
            const hier::Index<NDIM>& i = it();
            if ((bdry_locs_data(i, 0) == 0.0 && (bdry_tag & NORMAL_TRACTION_BDRY)) ||
                (bdry_locs_data(i, 0) == 1.0 && (bdry_tag & NORMAL_VELOCITY_BDRY)))
            {
                // Place i_g in the ghost cell abutting the boundary.
                hier::Index<NDIM> i_g = i;
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
                for (int k = 0; k < u_data->getGhostCellWidth()(bdry_normal_axis);
                     ++k, i_g(bdry_normal_axis) += (is_lower ? -1 : +1))
                {
                    // Determine the ghost cell value so that the divergence of
                    // the velocity field is zero in the ghost cell.
                    SideIndex<NDIM> i_g_s(
                        i_g, bdry_normal_axis, is_lower ? SideIndex<NDIM>::Lower : SideIndex<NDIM>::Upper);
                    (*u_data)(i_g_s) = divergenceFreeNormalGhostValue(*u_data, i_g, bdry_normal_axis, is_lower, dx);
                }
            }
        }
    }
    return;
} // enforceDivergenceFreeConditionAtBoundary

void
StaggeredStokesPhysicalBoundaryHelper::setNormalTractionGhostValues(
    const int u_data_idx,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const int coarsest_ln,
    const int finest_ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy);
#endif
    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    for (int ln = (coarsest_ln == IBTK::invalid_level_number ? 0 : coarsest_ln);
         ln <= (finest_ln == IBTK::invalid_level_number ? finest_hier_level : finest_ln);
         ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_data_idx);
                setNormalTractionGhostValues(u_data, patch, u_bc_coefs);
            }
        }
    }
    return;
} // setNormalTractionGhostValues

void
StaggeredStokesPhysicalBoundaryHelper::setNormalTractionGhostValues(
    Pointer<SideData<NDIM, double>> u_data,
    Pointer<Patch<NDIM>> patch,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs) const
{
    // TRACTION conditions prescribe -p + 2*mu*du_n/dx_n = g at an open boundary face.  The pressure boundary value is
    // p_b = -g, and the normal velocity ghost value set here supplies the viscous term.
    //
    // Notation, for a lower boundary: h is the grid spacing normal to the boundary; u_B, u_I, and u_G are the normal
    // velocities on the boundary face, on the next face inside the domain, and on the next face outside the domain;
    // p_I is the pressure in the cell abutting the boundary; and the pressure ghost value is p_G = 2*p_b - p_I.  For an
    // upper boundary, read h as -h.
    //
    // 1. A divergence-free normal velocity ghost value is not needed.
    //
    // The normal momentum equation at the boundary face involves ghost values only through the terms
    //
    //     mu*(u_I - 2*u_B + u_G)/h^2 - (p_I - p_G)/h,
    //
    // and so it depends on u_G and p_b only through p_b + mu*u_G/(2*h).  The direct discretization of the boundary
    // condition uses u_G = u_div, the value that makes the discrete divergence vanish in the ghost cell, together
    // with
    //
    //     p_b = 2*mu*D - g,    D = (u_I - u_div)/(2*h).
    //
    // The same equation results from p_b = -g and
    //
    //     u_G = u_div + 4*h*D = 2*u_I - u_div,
    //
    // which is the value set here.  It also results from the reflected value u_G = u_I and p_b = mu*D - g.  The
    // viscosity cancels from u_G, so each viscous term that is evaluated with this ghost value imposes the condition
    // with its own coefficient and at its own time level.  For PSEUDO_TRACTION conditions, -p = g, the reflected
    // value u_G = u_I and p_b = -g impose the condition and nothing is set here.
    //
    // 2. The result is a mimetic discretization.
    //
    // The equation at the boundary face is a momentum balance over the half cell between the boundary and the centers
    // of the cells abutting it.  With tangential velocity ghost values that satisfy mu*(du_t/dx_n + du_n/dx_t) = g_t
    // on the boundary, the terms above together with the tangential viscous terms equal
    //
    //     (2/h)*(-p_I + 2*mu*(u_I - u_B)/h - g) + D_t g_t - (mu/h)*(div u)_I,
    //
    // in which D_t g_t is the centered difference of the tangential traction data along the boundary, summed over the
    // tangential directions, and (div u)_I is the discrete divergence in the cell abutting the boundary, which is zero
    // for a discretely divergence-free velocity.  The normal stress -p + 2*mu*du_n/dx_n is evaluated at the cell
    // center and is the data g on the boundary; the shear stress on the boundary is the data g_t, at the nodes (edges
    // in three dimensions) at which the discrete shear stress is defined.  No velocity difference is taken across the
    // boundary.  With PSEUDO_TRACTION conditions the terms equal
    //
    //     (2/h)*(-p_I + mu*(u_I - u_B)/h - g) + mu*(D_t D_t u_n)_B.
    if (!patch->getPatchGeometry()->getTouchesRegularBoundary())
    {
        return;
    }
    const int ln = patch->getPatchLevelNumber();
    const int patch_num = patch->getPatchNumber();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Array<BoundaryBox<NDIM>>& physical_codim1_boxes = d_physical_codim1_boxes[ln].find(patch_num)->second;
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    const std::vector<Pointer<ArrayData<NDIM, bool>>>& dirichlet_bdry_locs =
        d_dirichlet_bdry_locs[ln].find(patch_num)->second;
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const unsigned int location_index = physical_codim1_boxes[n].getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower = location_index % 2 == 0;
        auto stokes_u_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
        if (!stokes_u_bc_coef || stokes_u_bc_coef->getTractionBcType() != TRACTION)
        {
            continue;
        }
        const Box<NDIM>& bc_coef_box = dirichlet_bdry_locs[n]->getBox();
        const ArrayData<NDIM, bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
        {
            const hier::Index<NDIM>& i = it();
            if (bdry_locs_data(i, 0) == 0.0)
            {
                // Place i_g in the ghost cell abutting the boundary, i_s_ghost on the ghost face of that cell, and
                // i_s_inside on the next face inside the domain.
                hier::Index<NDIM> i_g = i;
                if (is_lower)
                {
                    i_g(bdry_normal_axis) -= 1;
                }
                SideIndex<NDIM> i_s_ghost(
                    i_g, bdry_normal_axis, is_lower ? SideIndex<NDIM>::Lower : SideIndex<NDIM>::Upper);
                SideIndex<NDIM> i_s_inside(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                i_s_inside(bdry_normal_axis) += is_lower ? 1 : -1;
                const double u_inside = (*u_data)(i_s_inside);
                const double u_div = divergenceFreeNormalGhostValue(*u_data, i_g, bdry_normal_axis, is_lower, dx);
                (*u_data)(i_s_ghost) = 2.0 * u_inside - u_div;
            }
        }
    }
    return;
} // setNormalTractionGhostValues

double
StaggeredStokesPhysicalBoundaryHelper::divergenceFreeNormalGhostValue(const SideData<NDIM, double>& u_data,
                                                                      const hier::Index<NDIM>& i_g,
                                                                      const unsigned int normal_axis,
                                                                      const bool is_lower,
                                                                      const double* const dx)
{
    double div_u_g = 0.0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const SideIndex<NDIM> i_g_s_upper(i_g, axis, SideIndex<NDIM>::Upper);
        const SideIndex<NDIM> i_g_s_lower(i_g, axis, SideIndex<NDIM>::Lower);
        double u_upper = u_data(i_g_s_upper);
        double u_lower = u_data(i_g_s_lower);
        if (axis == normal_axis)
        {
            (is_lower ? u_lower : u_upper) = 0.0;
        }
        div_u_g += (u_upper - u_lower) * dx[normal_axis] / dx[axis];
    }
    return (is_lower ? +1.0 : -1.0) * div_u_g;
} // divergenceFreeNormalGhostValue

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
        auto extended_u_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[d]);
        if (extended_u_bc_coef)
        {
            extended_u_bc_coef->clearTargetPatchDataIndex();
            extended_u_bc_coef->setHomogeneousBc(homogeneous_bc);
        }
        auto stokes_u_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(u_bc_coefs[d]);
        if (stokes_u_bc_coef)
        {
            stokes_u_bc_coef->setTargetVelocityPatchDataIndex(u_target_data_idx);
            stokes_u_bc_coef->setTargetPressurePatchDataIndex(p_target_data_idx);
        }
    }
    auto extended_p_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(p_bc_coef);
    if (extended_p_bc_coef)
    {
        extended_p_bc_coef->clearTargetPatchDataIndex();
        extended_p_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    auto stokes_p_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(p_bc_coef);
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
        auto stokes_u_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(u_bc_coefs[d]);
        if (stokes_u_bc_coef)
        {
            stokes_u_bc_coef->clearTargetVelocityPatchDataIndex();
            stokes_u_bc_coef->clearTargetPressurePatchDataIndex();
        }
    }
    auto stokes_p_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(p_bc_coef);
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
