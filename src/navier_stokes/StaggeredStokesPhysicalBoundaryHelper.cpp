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
#include <ibtk/ibtk_enums.h>

#include <tbox/Array.h>
#include <tbox/MathUtilities.h>
#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <EdgeData.h>
#include <Index.h>
#include <IntVector.h>
#include <NodeData.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>
#include <Variable.h>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

// FORTRAN ROUTINES
#define H_AVG2_FC IBAMR_FC_FUNC_(h_avg2, H_AVG2)
#define H_AVG4_FC IBAMR_FC_FUNC_(h_avg4, H_AVG4)
#define H_AVG12_FC IBAMR_FC_FUNC_(h_avg12, H_AVG12)

extern "C"
{
    double H_AVG2_FC(const double&, const double&);

    double H_AVG4_FC(const double&, const double&, const double&, const double&);

    double H_AVG12_FC(const double&,
                      const double&,
                      const double&,
                      const double&,
                      const double&,
                      const double&,
                      const double&,
                      const double&,
                      const double&,
                      const double&,
                      const double&,
                      const double&);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace
{
// Return the normal velocity on the ghost face of the ghost cell i_g that makes the discrete divergence of u_data
// vanish in that cell.  The ghost face is the face of i_g normal to normal_axis that is farther from the physical
// boundary: the lower face if is_lower is true and the upper face otherwise.  The value of u_data on the ghost face is
// not read; the other faces of i_g are.
double
divergence_free_normal_ghost_value(const SideData<NDIM, double>& u_data,
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
}

// The data centering of the coefficient of the variable-coefficient viscous term, and the numbers of its values on a
// cell and on a face.
#if (NDIM == 2)
using ViscousCoefData = NodeData<NDIM, double>;
constexpr int NUM_CELL_COEFS = 4;
constexpr int NUM_FACE_COEFS = 2;
#endif
#if (NDIM == 3)
using ViscousCoefData = EdgeData<NDIM, double>;
constexpr int NUM_CELL_COEFS = 12;
constexpr int NUM_FACE_COEFS = 4;
#endif

// Return the average of the first n entries of v as the variable-coefficient viscous kernels compute it: the
// arithmetic mean, or the harmonic mean given by the functions that those kernels call, so that the two agree in how
// they treat zero values.  The number n is that of the coefficients on a face or on a cell.
double
average_viscous_coefs(const std::array<double, NUM_CELL_COEFS>& v, const int n, const VCInterpType interp_type)
{
    if (interp_type == VC_HARMONIC_INTERP)
    {
        switch (n)
        {
        case 2:
            return H_AVG2_FC(v[0], v[1]);
        case 4:
            return H_AVG4_FC(v[0], v[1], v[2], v[3]);
#if (NDIM == 3)
        case 12:
            return H_AVG12_FC(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11]);
#endif
        default:
            TBOX_ERROR("average_viscous_coefs(): unsupported number of values: " << n << "\n");
        }
    }
    double sum = 0.0;
    for (int k = 0; k < n; ++k)
    {
        sum += v[k];
    }
    return sum / n;
}

// Return the cell-centered viscous coefficient in the cell i_c as the variable-coefficient viscous kernels compute
// it: the average of the node-centered (two dimensions) or edge-centered (three dimensions) coefficients on the cell.
double
cell_viscous_coef(const ViscousCoefData& coef_data, const hier::Index<NDIM>& i_c, const VCInterpType interp_type)
{
    std::array<double, NUM_CELL_COEFS> values;
    int n = 0;
#if (NDIM == 2)
    const ArrayData<NDIM, double>& coef_array = coef_data.getArrayData();
    for (int j = 0; j <= 1; ++j)
    {
        for (int i = 0; i <= 1; ++i)
        {
            values[n++] = coef_array(i_c + hier::Index<NDIM>(i, j), 0);
        }
    }
#endif
#if (NDIM == 3)
    for (unsigned int edge_axis = 0; edge_axis < NDIM; ++edge_axis)
    {
        const ArrayData<NDIM, double>& coef_array = coef_data.getArrayData(edge_axis);
        const unsigned int axis_1 = (edge_axis + 1) % NDIM;
        const unsigned int axis_2 = (edge_axis + 2) % NDIM;
        for (int j = 0; j <= 1; ++j)
        {
            for (int i = 0; i <= 1; ++i)
            {
                hier::Index<NDIM> i_e = i_c;
                i_e(axis_1) += i;
                i_e(axis_2) += j;
                values[n++] = coef_array(i_e, 0);
            }
        }
    }
#endif
    return average_viscous_coefs(values, n, interp_type);
}

// Return the average of the node-centered (two dimensions) or edge-centered (three dimensions) viscous coefficients
// on the face with index i_f that is normal to normal_axis.
double
face_viscous_coef(const ViscousCoefData& coef_data,
                  const hier::Index<NDIM>& i_f,
                  const unsigned int normal_axis,
                  const VCInterpType interp_type)
{
    std::array<double, NUM_CELL_COEFS> values;
    int n = 0;
#if (NDIM == 2)
    const ArrayData<NDIM, double>& coef_array = coef_data.getArrayData();
    const unsigned int tangential_axis = (normal_axis + 1) % NDIM;
    for (int i = 0; i <= 1; ++i)
    {
        hier::Index<NDIM> i_n = i_f;
        i_n(tangential_axis) += i;
        values[n++] = coef_array(i_n, 0);
    }
#endif
#if (NDIM == 3)
    // The edges of the face that are parallel to one tangential axis are offset along the other.
    for (unsigned int k = 1; k < NDIM; ++k)
    {
        const unsigned int edge_axis = (normal_axis + k) % NDIM;
        const unsigned int offset_axis = (normal_axis + NDIM - k) % NDIM;
        const ArrayData<NDIM, double>& coef_array = coef_data.getArrayData(edge_axis);
        for (int i = 0; i <= 1; ++i)
        {
            hier::Index<NDIM> i_e = i_f;
            i_e(offset_axis) += i;
            values[n++] = coef_array(i_e, 0);
        }
    }
#endif
#if !defined(NDEBUG)
    TBOX_ASSERT(n == NUM_FACE_COEFS);
#endif
    return average_viscous_coefs(values, n, interp_type);
}
} // namespace

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
                    (*u_data)(i_g_s) = divergence_free_normal_ghost_value(*u_data, i_g, bdry_normal_axis, is_lower, dx);
                }
            }
        }
    }
    return;
} // enforceDivergenceFreeConditionAtBoundary

void
StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm(
    const int f_data_idx,
    const int u_data_idx,
    const double viscous_coef,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool linear_pressure_extrapolation,
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
                Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_data_idx);
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_data_idx);
                addNormalTractionViscousTerm(
                    f_data, u_data, patch, viscous_coef, u_bc_coefs, linear_pressure_extrapolation);
            }
        }
    }
    return;
} // addNormalTractionViscousTerm

void
StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm(
    Pointer<SideData<NDIM, double>> f_data,
    Pointer<SideData<NDIM, double>> u_data,
    Pointer<Patch<NDIM>> patch,
    const double viscous_coef,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool linear_pressure_extrapolation) const
{
    // TRACTION conditions prescribe -p + 2*mu*du_n/dx_n = g at a boundary face where the normal velocity is not
    // prescribed.  The pressure boundary value is p_b = -g, and the velocity boundary conditions set the reflected
    // normal velocity ghost value.  This function adds the part of the viscous term at the boundary face that those
    // leave out.
    //
    // Notation, for a lower boundary: h is the grid spacing normal to the boundary; u_B, u_I, and u_G are the normal
    // velocities on the boundary face, on the next face inside the domain, and on the next face outside the domain; p_I
    // is the pressure in the cell abutting the boundary; and the pressure ghost value is p_G = 2*p_b - p_I.  For an
    // upper boundary, read h as -h.
    //
    // 1. A divergence-free normal velocity ghost value is not needed.
    //
    // The normal momentum equation at the boundary face involves ghost values only through the terms
    //
    //     mu*(u_I - 2*u_B + u_G)/h^2 - (p_I - p_G)/h,
    //
    // and so it depends on u_G and p_b only through p_b + mu*u_G/(2*h).  The direct discretization of the boundary
    // condition uses u_G = u_div, the value that makes the discrete divergence vanish in the ghost cell, together with
    //
    //     p_b = 2*mu*E - g,    E = (u_I - u_div)/(2*h).
    //
    // The same equation results from p_b = -g and u_G = 2*u_I - u_div, that is, from p_b = -g, the reflected value
    // u_G = u_I, and the additional term
    //
    //     mu*(u_I - u_div)/h^2
    //
    // at the boundary face, which is what is added here, with the coefficient of the viscous term in place of mu.  Each
    // viscous term to which this function is applied therefore imposes the condition with its own coefficient and at
    // its own time level.  PSEUDO_TRACTION conditions prescribe -p + mu*du_n/dx_n = g.  Their direct discretization
    // uses u_G = u_div and p_b = mu*E - g, which is the same equation as the reflected value with p_b = -g, so nothing
    // is added.
    //
    // 2. The result is a mimetic discretization.
    //
    // The equation at the boundary face is a momentum balance over the half cell between the boundary and the centers
    // of the cells abutting it.  With tangential velocity ghost values that satisfy mu*(du_t/dx_n + du_n/dx_t) = g_t on
    // the boundary, the terms above together with the tangential viscous terms equal
    //
    //     (2/h)*(-p_I + 2*mu*(u_I - u_B)/h - g) + D_t g_t - (mu/h)*(div u)_I,
    //
    // in which D_t g_t is the centered difference of the tangential traction data along the boundary, summed over the
    // tangential directions, and (div u)_I is the discrete divergence in the cell abutting the boundary, which is zero
    // for a discretely divergence-free velocity.  The normal stress -p + 2*mu*du_n/dx_n is evaluated at the cell center
    // and is the data g on the boundary; the shear stress on the boundary is the data g_t, at the nodes (edges in three
    // dimensions) at which the discrete shear stress is defined.  No velocity difference is taken across the boundary.
    // With PSEUDO_TRACTION conditions the terms equal
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
        const double h = dx[bdry_normal_axis];
        const Box<NDIM>& bc_coef_box = dirichlet_bdry_locs[n]->getBox();
        const ArrayData<NDIM, bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
        for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
        {
            const hier::Index<NDIM>& i = it();
            if (bdry_locs_data(i, 0) == 0.0)
            {
                // The added term assumes that the pressure ghost value is p_G = 2*p_b - p_I.
                if (!linear_pressure_extrapolation && viscous_coef != 0.0)
                {
                    TBOX_ERROR("StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm():\n"
                               << "  TRACTION conditions at a boundary at which the normal velocity is not prescribed\n"
                               << "  require the boundary interpolation type \"LINEAR\", set by the input parameter\n"
                               << "  bdry_interp_type of StaggeredStokesOperator and by U_P_bdry_interp_type of\n"
                               << "  INSStaggeredHierarchyIntegrator.  The type \"QUADRATIC\" is not supported.\n");
                }
                // Place i_g in the ghost cell abutting the boundary, i_s_boundary on the boundary face, and
                // i_s_inside on the next face inside the domain.
                hier::Index<NDIM> i_g = i;
                if (is_lower)
                {
                    i_g(bdry_normal_axis) -= 1;
                }
                const SideIndex<NDIM> i_s_boundary(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                SideIndex<NDIM> i_s_inside = i_s_boundary;
                i_s_inside(bdry_normal_axis) += is_lower ? 1 : -1;
                const double u_inside = (*u_data)(i_s_inside);
                const double u_div = divergence_free_normal_ghost_value(*u_data, i_g, bdry_normal_axis, is_lower, dx);
                (*f_data)(i_s_boundary) += viscous_coef * (u_inside - u_div) / (h * h);
            }
        }
    }
    return;
} // addNormalTractionViscousTerm

void
StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm(
    const int f_data_idx,
    const int u_data_idx,
    const int viscous_coef_data_idx,
    const VCInterpType viscous_coef_interp_type,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool linear_pressure_extrapolation,
    const int coarsest_ln,
    const int finest_ln) const
{
    // The variable-coefficient viscous term is the divergence of the viscous stress D*(grad u + grad u^T).  With the
    // notation of the version of this function for a constant coefficient, and with D_I and D_G the cell-centered
    // coefficients in the cell abutting the boundary and in the ghost cell, the normal momentum equation at a boundary
    // face involves ghost values only through the terms
    //
    //     (2/h^2)*(D_I*(u_I - u_B) - D_G*(u_B - u_G)) - (p_I - p_G)/h + D_t tau,
    //
    // in which tau = D*(du_n/dx_t + du_t/dx_n) is the shear stress at the nodes (edges in three dimensions) on the
    // boundary and D_t tau is its centered difference along the boundary, summed over the tangential directions.
    //
    // The momentum balance over the half cell between the boundary and the centers of the cells abutting it, with the
    // normal stress -p + 2*mu*du_n/dx_n = g on the boundary, is
    //
    //     (2/h)*(-p_I + 2*D_I*(u_I - u_B)/h - g) + D_t tau.
    //
    // With p_b = -g, the two agree if (2/h^2)*(D_I*(u_I - u_B) + D_G*(u_B - u_G)) is added at the boundary face, which
    // replaces the normal viscous stress in the ghost cell by the reflection of that in the cell abutting the
    // boundary.  The result does not depend on u_G or on D_G, and TRACTION velocity conditions make tau the
    // tangential traction data.
    //
    // PSEUDO_TRACTION conditions prescribe -p + mu*du_n/dx_n = g, so the normal stress on the boundary is
    // g + mu*du_n/dx_n, in which the normal derivative on the boundary is (u_I - u_div)/(2*h).  The term
    // D_B*(u_I - u_div)/h^2 is then subtracted, in which D_B is the average of the coefficient over the boundary face.
    //
    // For a constant coefficient and either type of condition, the result equals the constant-coefficient viscous
    // term D*(Laplacian u) with the term added by the version of this function for a constant coefficient, plus
    // (D/h)*(div u)_I, in which (div u)_I is the discrete divergence in the cell abutting the boundary.
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
            if (!patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                continue;
            }
            Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_data_idx);
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_data_idx);
            Pointer<ViscousCoefData> coef_data = patch->getPatchData(viscous_coef_data_idx);
#if !defined(NDEBUG)
            TBOX_ASSERT(f_data);
            TBOX_ASSERT(u_data);
            TBOX_ASSERT(coef_data);
#endif
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
                if (!stokes_u_bc_coef)
                {
                    continue;
                }
                const bool pseudo_traction = stokes_u_bc_coef->getTractionBcType() == PSEUDO_TRACTION;
                const double h = dx[bdry_normal_axis];
                const Box<NDIM>& bc_coef_box = dirichlet_bdry_locs[n]->getBox();
                const ArrayData<NDIM, bool>& bdry_locs_data = *dirichlet_bdry_locs[n];
                for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                {
                    const hier::Index<NDIM>& i = it();
                    if (bdry_locs_data(i, 0) != 0.0)
                    {
                        continue;
                    }
                    // Place i_g in the ghost cell abutting the boundary, i_c in the cell abutting the boundary inside
                    // the domain, i_s_boundary on the boundary face, and i_s_inside and i_s_outside on the next faces
                    // inside and outside the domain.
                    hier::Index<NDIM> i_g = i, i_c = i;
                    (is_lower ? i_g : i_c)(bdry_normal_axis) -= 1;
                    const SideIndex<NDIM> i_s_boundary(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                    SideIndex<NDIM> i_s_inside = i_s_boundary, i_s_outside = i_s_boundary;
                    i_s_inside(bdry_normal_axis) += is_lower ? 1 : -1;
                    i_s_outside(bdry_normal_axis) += is_lower ? -1 : 1;
                    const double u_boundary = (*u_data)(i_s_boundary);
                    const double u_inside = (*u_data)(i_s_inside);
                    const double u_outside = (*u_data)(i_s_outside);
                    const double coef_inside = cell_viscous_coef(*coef_data, i_c, viscous_coef_interp_type);
                    const double coef_outside = cell_viscous_coef(*coef_data, i_g, viscous_coef_interp_type);
                    // The added term assumes that the pressure ghost value is p_G = 2*p_b - p_I.
                    if (!linear_pressure_extrapolation)
                    {
                        TBOX_ERROR(
                            "StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm():\n"
                            << "  TRACTION and PSEUDO_TRACTION conditions at a boundary at which the normal velocity\n"
                            << "  is not prescribed require the boundary interpolation type \"LINEAR\", set by the\n"
                            << "  input parameter bdry_interp_type of VCStaggeredStokesOperator.  The type\n"
                            << "  \"QUADRATIC\" is not supported.\n");
                    }
                    double term =
                        2.0 * (coef_inside * (u_inside - u_boundary) + coef_outside * (u_boundary - u_outside));
                    if (pseudo_traction)
                    {
                        const double coef_boundary =
                            face_viscous_coef(*coef_data, i, bdry_normal_axis, viscous_coef_interp_type);
                        const double u_div =
                            divergence_free_normal_ghost_value(*u_data, i_g, bdry_normal_axis, is_lower, dx);
                        term -= coef_boundary * (u_inside - u_div);
                    }
                    (*f_data)(i_s_boundary) += term / (h * h);
                }
            }
        }
    }
    return;
} // addNormalTractionViscousTerm

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
