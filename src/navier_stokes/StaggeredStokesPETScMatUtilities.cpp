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

#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StokesBcCoefStrategy.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/compiler_hints.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Array.h>
#include <tbox/MathUtilities.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <petsclog.h>
#include <petscmat.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellGeometry.h>
#include <CellIndex.h>
#include <Index.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <set>
#include <utility>
#include <vector>

#include "./ins_staggered_traction_stencil.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
} // compute_tangential_extension

// The key of a row of the velocity part of the matrix: the component and the index of the side.
using RowKey = std::array<int, NDIM + 1>;

RowKey
get_row_key(const SideIndex<NDIM>& i_s)
{
    RowKey key;
    key[0] = static_cast<int>(i_s.getAxis());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d + 1] = i_s(d);
    }
    return key;
} // get_row_key

// The terms of the rows of the velocity part of the matrix that couple a tangential velocity to the normal velocities
// on a physical boundary, as pairs of the global DOF index of the column and the coefficient.
using TractionCouplings = std::map<RowKey, std::vector<std::pair<int, double>>>;

// Add value to the coefficient of column in row, which is appended to the row if the row does not contain the column.
void
add_to_row(std::vector<std::pair<int, double>>& row, const int column, const double value)
{
    auto it = std::find_if(row.begin(), row.end(), [column](const auto& c) { return c.first == column; });
    if (it == row.end())
    {
        row.emplace_back(column, value);
    }
    else
    {
        it->second += value;
    }
} // add_to_row

// Return the couplings of the rows of the tangential velocity components next to the physical boundaries of the patch
// at which the boundary condition object of the component imposes TRACTION conditions.
//
// With a TRACTION condition, the ghost value of the tangential velocity next to the boundary is u_G = u_I + h*gamma,
// in which h is the grid spacing normal to the boundary and gamma, the inhomogeneous Robin coefficient, includes the
// tangential difference of the normal velocity on the boundary. The row of u_I contains the ghost coefficient D/h^2
// times u_G, so the derivative of gamma with respect to a normal velocity multiplied by D/h is the coefficient of that
// velocity in the row. The part of gamma that does not depend on the solution is data and belongs to the right-hand
// side. A normal velocity that is not a DOF of the level, which is a velocity in a coarse-fine ghost cell, has no
// column; the right-hand side accounts for it.
TractionCouplings
compute_traction_couplings(Patch<NDIM>& patch,
                           const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                           const double data_time,
                           const double D,
                           const SideData<NDIM, int>& u_dof_index_data)
{
    TractionCouplings couplings;
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Box<NDIM> ghost_box = Box<NDIM>::grow(patch.getBox(), IntVector<NDIM>(1));
    const Array<BoundaryBox<NDIM>> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        auto u_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(u_bc_coefs[axis]);
        if (!u_bc_coef)
        {
            continue;
        }
        std::array<double, NDIM> shift;
        shift.fill(0.0);
        shift[axis] = -0.5 * dx[axis];
        const traction_stencil::ShiftedPatchGeometry shifted_geometry(patch, shift);
        for (int n = 0; n < physical_codim1_boxes.size(); ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;
            if (bdry_normal_axis == axis)
            {
                continue;
            }
            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, patch);
            const Box<NDIM> bc_coef_box = compute_tangential_extension(
                PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);
            Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double>> gcoef_data;
            static const bool homogeneous_bc = true;
            u_bc_coef->clearTargetPatchDataIndex();
            u_bc_coef->setHomogeneousBc(homogeneous_bc);
            const traction_stencil::ShiftedPatchGeometry::Scope scope(patch, shifted_geometry);
            u_bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, nullptr, patch, trimmed_bdry_box, data_time);
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const hier::Index<NDIM>& i = bc();
                if (!((*acoef_data)(i, 0) == 0.0 && (*bcoef_data)(i, 0) == 1.0))
                {
                    continue;
                }
                const auto stencil =
                    u_bc_coef->getNormalVelocityStencil(i, patch, trimmed_bdry_box, ghost_box, data_time);
                if (stencil.empty())
                {
                    continue;
                }
                hier::Index<NDIM> i_intr = i;
                if (!is_lower)
                {
                    i_intr(bdry_normal_axis) -= 1;
                }
                std::vector<std::pair<int, double>>& row =
                    couplings[get_row_key(SideIndex<NDIM>(i_intr, axis, SideIndex<NDIM>::Lower))];
                const double h = dx[bdry_normal_axis];
                for (const auto& entry : stencil)
                {
                    const int column = u_dof_index_data(entry.first);
                    if (column < 0)
                    {
                        continue;
                    }
                    add_to_row(row, column, D / h * entry.second);
                }
            }
        }
    }
    return couplings;
} // compute_traction_couplings

// Add the terms of the rows of the normal velocity on the physical boundaries of the patch at which TRACTION
// conditions are imposed and the normal velocity is not prescribed, as pairs of the global DOF index of the column and
// the coefficient.
//
// At such a boundary face, StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm() adds D/h^2 times
// u_I - u_div to the viscous term, in which h is the grid spacing normal to the boundary, u_I is the normal velocity on
// the next face inside the domain, and u_div is the normal velocity that makes the discrete divergence vanish in the
// ghost cell abutting the boundary. The stencil of u_I - u_div, which traction_stencil::get_normal_stress_stencil()
// provides, reads the normal velocity on the boundary face and the tangential velocities on the faces of the ghost
// cell, which are ghost values. A ghost value of a tangential velocity is f_i*u + f_g*gamma, in which u is the velocity
// on the face inside the domain, f_i and f_g are the coefficients of the linear extrapolation that defines the ghost
// value, and gamma is the inhomogeneous Robin coefficient. The derivative of gamma with respect to the normal
// velocities on the boundary, which the boundary condition object of the tangential component provides as for the
// couplings of the tangential rows, gives the remaining terms. The part of gamma that does not depend on the solution
// is data and belongs to the right-hand side, as are velocities that are not DOFs of the level.
void
add_normal_stress_couplings(TractionCouplings& couplings,
                            Patch<NDIM>& patch,
                            const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                            const double data_time,
                            const double D,
                            const SideData<NDIM, int>& u_dof_index_data)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Box<NDIM> ghost_box = Box<NDIM>::grow(patch.getBox(), IntVector<NDIM>(1));
    const Array<BoundaryBox<NDIM>> physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(patch);
    for (int n = 0; n < physical_codim1_boxes.size(); ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const bool is_lower = location_index % 2 == 0;
        const auto normal_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(u_bc_coefs[bdry_normal_axis]);
        if (!normal_bc_coef || normal_bc_coef->getTractionBcType() != TRACTION)
        {
            continue;
        }
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, patch);
        const Box<NDIM> normal_bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        // Find the faces at which the normal velocity is not prescribed.
        Pointer<ArrayData<NDIM, double>> normal_acoef_data = new ArrayData<NDIM, double>(normal_bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double>> normal_bcoef_data = new ArrayData<NDIM, double>(normal_bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double>> normal_gcoef_data;
        normal_bc_coef->clearTargetPatchDataIndex();
        normal_bc_coef->setHomogeneousBc(true);
        normal_bc_coef->setBcCoefs(
            normal_acoef_data, normal_bcoef_data, normal_gcoef_data, nullptr, patch, trimmed_bdry_box, data_time);
        std::vector<hier::Index<NDIM>> open_faces;
        for (Box<NDIM>::Iterator it(normal_bc_coef_box); it; it++)
        {
            if ((*normal_acoef_data)(it(), 0) == 0.0 && (*normal_bcoef_data)(it(), 0) == 1.0)
            {
                open_faces.push_back(it());
            }
        }
        if (open_faces.empty())
        {
            continue;
        }

        // Determine the boundary conditions of the tangential components at the locations of their ghost values, and
        // the dependence of the inhomogeneous coefficient on the normal velocities.
        struct TangentialBc
        {
            Pointer<ArrayData<NDIM, double>> acoef_data, bcoef_data;
            std::map<RowKey, std::vector<std::pair<SideIndex<NDIM>, double>>> stencils;
        };
        std::array<TangentialBc, NDIM> tangential_bcs;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            if (axis == bdry_normal_axis)
            {
                continue;
            }
            const Box<NDIM> bc_coef_box = compute_tangential_extension(normal_bc_coef_box, axis);
            TangentialBc& bc = tangential_bcs[axis];
            bc.acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            bc.bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double>> gcoef_data;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(true);
            }
            std::array<double, NDIM> shift;
            shift.fill(0.0);
            shift[axis] = -0.5 * dx[axis];
            const traction_stencil::ShiftedPatchGeometry shifted_geometry(patch, shift);
            const traction_stencil::ShiftedPatchGeometry::Scope scope(patch, shifted_geometry);
            u_bc_coefs[axis]->setBcCoefs(
                bc.acoef_data, bc.bcoef_data, gcoef_data, nullptr, patch, trimmed_bdry_box, data_time);
            auto ins_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(u_bc_coefs[axis]);
            if (!ins_bc_coef)
            {
                continue;
            }
            for (Box<NDIM>::Iterator bc_it(bc_coef_box); bc_it; bc_it++)
            {
                const hier::Index<NDIM>& i = bc_it();
                if ((*bc.acoef_data)(i, 0) == 0.0 && (*bc.bcoef_data)(i, 0) == 1.0)
                {
                    bc.stencils[get_row_key(SideIndex<NDIM>(i, axis, SideIndex<NDIM>::Lower))] =
                        ins_bc_coef->getNormalVelocityStencil(i, patch, trimmed_bdry_box, ghost_box, data_time);
                }
            }
        }

        // Add the terms of the rows.
        const double h = dx[bdry_normal_axis];
        const double coef = D / (h * h);
        for (const hier::Index<NDIM>& i : open_faces)
        {
            std::vector<std::pair<int, double>>& row =
                couplings[get_row_key(SideIndex<NDIM>(i, bdry_normal_axis, SideIndex<NDIM>::Lower))];
            const auto add_term = [&row, &u_dof_index_data](const SideIndex<NDIM>& i_s, const double value)
            {
                const int column = u_dof_index_data(i_s);
                if (column >= 0)
                {
                    add_to_row(row, column, value);
                }
            };
            for (const auto& entry : traction_stencil::get_normal_stress_stencil(i, bdry_normal_axis, is_lower, dx))
            {
                const SideIndex<NDIM>& i_s = entry.first;
                const unsigned int axis = i_s.getAxis();
                if (axis == bdry_normal_axis)
                {
                    add_term(i_s, coef * entry.second);
                    continue;
                }

                // The face is a ghost face of the tangential velocity. Its location on the boundary has the index of
                // the boundary face along the axis normal to the boundary, and the face inside the domain from which
                // the ghost value is extrapolated is the one next to it along that axis.
                hier::Index<NDIM> i_bc = i_s;
                SideIndex<NDIM> i_s_inside = i_s;
                if (is_lower)
                {
                    i_bc(bdry_normal_axis) += 1;
                    i_s_inside(bdry_normal_axis) += 1;
                }
                else
                {
                    i_s_inside(bdry_normal_axis) -= 1;
                }
                const TangentialBc& bc = tangential_bcs[axis];
                const double a = (*bc.acoef_data)(i_bc, 0);
                const double b = (*bc.bcoef_data)(i_bc, 0);
                if (!((a == 1.0 && b == 0.0) || (a == 0.0 && b == 1.0)))
                {
                    TBOX_ERROR("StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp():\n"
                               << "  unsupported boundary condition coefficients (a, b) = (" << a << ", " << b
                               << ") for the tangential velocity.\n"
                               << "  Only a prescribed velocity, (a, b) = (1, 0), or a prescribed traction, (a, b) = "
                                  "(0, 1),\n"
                               << "  is supported.\n");
                }
                const double f_i = -(a * h - 2.0 * b) / (a * h + 2.0 * b);
                const double f_g = 2.0 * h / (a * h + 2.0 * b);
                add_term(i_s_inside, coef * entry.second * f_i);
                const auto stencil = bc.stencils.find(get_row_key(SideIndex<NDIM>(i_bc, axis, SideIndex<NDIM>::Lower)));
                if (stencil != bc.stencils.end())
                {
                    for (const auto& dgamma : stencil->second)
                    {
                        add_term(dgamma.first, coef * entry.second * f_g * dgamma.second);
                    }
                }
            }
        }
    }
    return;
} // add_normal_stress_couplings
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
    Mat& mat,
    const PoissonSpecifications& u_problem_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    double data_time,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    Pointer<PatchLevel<NDIM>> patch_level,
    RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Setup the finite difference stencils.
    static const int uu_stencil_sz = 2 * NDIM + 1;
    std::array<hier::Index<NDIM>, uu_stencil_sz> uu_stencil(
        array_constant<hier::Index<NDIM>, uu_stencil_sz>(hier::Index<NDIM>(0)));
    for (unsigned int axis = 0, uu_stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
        {
            uu_stencil[uu_stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }
    static const int up_stencil_sz = 2;
    std::array<std::array<hier::Index<NDIM>, up_stencil_sz>, NDIM> up_stencil(
        array_constant<std::array<hier::Index<NDIM>, up_stencil_sz>, NDIM>(
            array_constant<hier::Index<NDIM>, up_stencil_sz>(hier::Index<NDIM>(0))));
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side)
        {
            up_stencil[axis][side](axis) = (side == 0 ? -1 : 0);
        }
    }
    static const int pu_stencil_sz = 2 * NDIM;
    std::array<hier::Index<NDIM>, pu_stencil_sz> pu_stencil(
        array_constant<hier::Index<NDIM>, pu_stencil_sz>(hier::Index<NDIM>(0)));
    for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
        {
            pu_stencil[pu_stencil_index](axis) = (side == 0 ? 0 : +1);
        }
    }

    // Determine the index ranges.
    const int mpi_rank = IBTK_MPI::getRank();
    const int nlocal = num_dofs_per_proc[mpi_rank];
    const int ilower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int iupper = ilower + nlocal;
    const int ntotal = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the couplings of tangential velocities to the normal velocities on the physical boundary.
    const double C = (u_problem_coefs.cIsZero() ? 0.0 : u_problem_coefs.getCConstant());
    const double D = u_problem_coefs.getDConstant();
    std::map<int, TractionCouplings> traction_couplings;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        if (!patch->getPatchGeometry()->getTouchesRegularBoundary())
        {
            continue;
        }
        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        TractionCouplings couplings = compute_traction_couplings(*patch, u_bc_coefs, data_time, D, *u_dof_index_data);
        add_normal_stress_couplings(couplings, *patch, u_bc_coefs, data_time, D, *u_dof_index_data);
        traction_couplings[patch->getPatchNumber()] = std::move(couplings);
    }

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(nlocal, 0), o_nnz(nlocal, 0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        const auto patch_couplings = traction_couplings.find(patch->getPatchNumber());
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& ic = b();
                const SideIndex<NDIM> is(ic, axis, SideIndex<NDIM>::Lower);
                const int u_dof_index = (*u_dof_index_data)(is);
                if (UNLIKELY(ilower > u_dof_index || u_dof_index >= iupper)) continue;
                const int u_local_idx = u_dof_index - ilower;
                d_nnz[u_local_idx] += 1;
                for (unsigned int d = 0, uu_stencil_index = 1; d < NDIM; ++d)
                {
                    for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
                    {
                        const int uu_dof_index = (*u_dof_index_data)(is + uu_stencil[uu_stencil_index]);
                        if (LIKELY(uu_dof_index >= ilower && uu_dof_index < iupper))
                        {
                            d_nnz[u_local_idx] += 1;
                        }
                        else
                        {
                            o_nnz[u_local_idx] += 1;
                        }
                    }
                }
                for (int side = 0, up_stencil_index = 0; side <= 1; ++side, ++up_stencil_index)
                {
                    const int up_dof_index = (*p_dof_index_data)(ic + up_stencil[axis][up_stencil_index]);
                    if (LIKELY(up_dof_index >= ilower && up_dof_index < iupper))
                    {
                        d_nnz[u_local_idx] += 1;
                    }
                    else
                    {
                        o_nnz[u_local_idx] += 1;
                    }
                }
                if (patch_couplings != traction_couplings.end())
                {
                    const auto row = patch_couplings->second.find(get_row_key(is));
                    if (row != patch_couplings->second.end())
                    {
                        for (const auto& coupling : row->second)
                        {
                            if (coupling.first >= ilower && coupling.first < iupper)
                            {
                                d_nnz[u_local_idx] += 1;
                            }
                            else
                            {
                                o_nnz[u_local_idx] += 1;
                            }
                        }
                    }
                }
                d_nnz[u_local_idx] = std::min(nlocal, d_nnz[u_local_idx]);
                o_nnz[u_local_idx] = std::min(ntotal - nlocal, o_nnz[u_local_idx]);
            }
        }
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            const int p_dof_index = (*p_dof_index_data)(ic);
            if (UNLIKELY(ilower > p_dof_index || p_dof_index >= iupper)) continue;
            const int p_local_idx = p_dof_index - ilower;
            d_nnz[p_local_idx] += 1;
            for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
            {
                for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
                {
                    const int pu_dof_index = (*u_dof_index_data)(
                        SideIndex<NDIM>(ic + pu_stencil[pu_stencil_index], axis, SideIndex<NDIM>::Lower));
                    if (LIKELY(pu_dof_index >= ilower && pu_dof_index < iupper))
                    {
                        d_nnz[p_local_idx] += 1;
                    }
                    else
                    {
                        o_nnz[p_local_idx] += 1;
                    }
                }
            }
            d_nnz[p_local_idx] = std::min(nlocal, d_nnz[p_local_idx]);
            o_nnz[p_local_idx] = std::min(ntotal - nlocal, o_nnz[p_local_idx]);
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        nlocal,
                        nlocal,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        0,
                        get_data_or_null(d_nnz),
                        0,
                        get_data_or_null(o_nnz),
                        &mat);
    IBTK_CHKERRQ(ierr);

// Set some general matrix options.
#if !defined(NDEBUG)
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        TractionCouplings empty_couplings;
        const auto patch_couplings = traction_couplings.find(patch->getPatchNumber());
        TractionCouplings& couplings =
            patch_couplings != traction_couplings.end() ? patch_couplings->second : empty_couplings;
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        const IntVector<NDIM> no_ghosts(0);
        SideData<NDIM, double> uu_matrix_coefs(patch_box, uu_stencil_sz, no_ghosts);
        SideData<NDIM, double> up_matrix_coefs(patch_box, up_stencil_sz, no_ghosts);
        CellData<NDIM, double> pu_matrix_coefs(patch_box, pu_stencil_sz, no_ghosts);

        // Compute all matrix coefficients, including those on the physical
        // boundary; however, do not yet take physical boundary conditions into
        // account.  Boundary conditions are handled subsequently.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            std::vector<double> uu_mat_vals(uu_stencil_sz, 0.0);
            uu_mat_vals[0] = C; // diagonal
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double dx_sq = dx[d] * dx[d];
                uu_mat_vals[0] -= 2 * D / dx_sq;    // diagonal
                uu_mat_vals[2 * d + 1] = D / dx_sq; // lower off-diagonal
                uu_mat_vals[2 * d + 2] = D / dx_sq; // upper off-diagonal
            }
            for (int uu_stencil_index = 0; uu_stencil_index < uu_stencil_sz; ++uu_stencil_index)
            {
                uu_matrix_coefs.fill(uu_mat_vals[uu_stencil_index], uu_stencil_index);
            }

            // grad p
            for (int d = 0; d < NDIM; ++d)
            {
                up_matrix_coefs.getArrayData(d).fill(-1.0 / dx[d], 0);
                up_matrix_coefs.getArrayData(d).fill(+1.0 / dx[d], 1);
            }

            // -div u
            std::vector<double> pu_mat_vals(pu_stencil_sz, 0.0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                pu_matrix_coefs.fill(+1.0 / dx[d], 2 * d);
                pu_matrix_coefs.fill(-1.0 / dx[d], 2 * d + 1);
            }
        }

        // Data structures required to set physical boundary conditions.
        const Array<BoundaryBox<NDIM>> physical_codim1_boxes =
            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();
        const double* const patch_x_lower = pgeom->getXLower();
        const double* const patch_x_upper = pgeom->getXUpper();
        const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
        Array<Array<bool>> touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            touches_regular_bdry[axis].resizeArray(2);
            touches_periodic_bdry[axis].resizeArray(2);
            for (int upperlower = 0; upperlower < 2; ++upperlower)
            {
                touches_regular_bdry[axis][upperlower] = pgeom->getTouchesRegularBoundary(axis, upperlower);
                touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis, upperlower);
            }
        }

        // Modify matrix coefficients to account for physical boundary
        // conditions along boundaries which ARE NOT aligned with the data axis.
        //
        // NOTE: It important to set these values first to avoid problems at
        // corners in the physical domain.  In particular, since Dirichlet
        // boundary conditions for values located on the physical boundary
        // override all other boundary conditions, we set those values last.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const bool is_lower = location_index % 2 == 0;

                if (bdry_normal_axis == axis) continue;

                const Box<NDIM> bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
                const BoundaryBox<NDIM> trimmed_bdry_box =
                    PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box<NDIM> bc_coef_box = compute_tangential_extension(
                    PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

                Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef_data;

                // Temporarily reset the patch geometry object associated with
                // the patch so that boundary conditions are set at the correct
                // spatial locations.
                std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    shifted_patch_x_lower[d] = patch_x_lower[d];
                    shifted_patch_x_upper[d] = patch_x_upper[d];
                }
                shifted_patch_x_lower[axis] -= 0.5 * dx[axis];
                shifted_patch_x_upper[axis] -= 0.5 * dx[axis];
                patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                         touches_regular_bdry,
                                                                         touches_periodic_bdry,
                                                                         dx,
                                                                         shifted_patch_x_lower.data(),
                                                                         shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Restore the original patch geometry object.
                patch->setPatchGeometry(pgeom);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const hier::Index<NDIM>& i = bc();
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const bool velocity_bc = (a == 1.0 && b == 0.0);
                    const bool traction_bc = (a == 0.0 && b == 1.0);
                    hier::Index<NDIM> i_intr = i;
                    if (is_lower)
                    {
                        i_intr(bdry_normal_axis) += 0;
                    }
                    else
                    {
                        i_intr(bdry_normal_axis) -= 1;
                    }
                    const SideIndex<NDIM> i_s(i_intr, axis, SideIndex<NDIM>::Lower);

                    if (velocity_bc)
                    {
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 0) -= uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 0) -= uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }
                    }
                    else if (traction_bc)
                    {
                        // A traction condition prescribes the normal derivative of a tangential velocity
                        // component, so the ghost value is the interior value plus boundary data. Add the ghost
                        // coefficient to the diagonal.
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 0) += uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 0) += uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }
                    }
                    else
                    {
                        TBOX_ERROR(
                            "StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp():\n"
                            << "  unsupported boundary condition coefficients (a, b) = (" << a << ", " << b
                            << ") for the tangential velocity.\n"
                            << "  Only a prescribed velocity, (a, b) = (1, 0), or a prescribed traction, (a, b) = "
                               "(0, 1),\n"
                            << "  is supported.\n");
                    }
                }
            }
        }

        // Modify matrix coefficients to account for physical boundary
        // conditions along boundaries which ARE aligned with the data axis.
        //
        // NOTE: It important to set these values last to avoid problems at corners
        // in the physical domain.  In particular, since Dirichlet boundary
        // conditions for values located on the physical boundary override all other
        // boundary conditions, we set those values last.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const bool is_lower = location_index % 2 == 0;

                if (bdry_normal_axis != axis) continue;

                const Box<NDIM> bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
                const BoundaryBox<NDIM> trimmed_bdry_box =
                    PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

                Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef_data;

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Set the pressure boundary condition coefficients.
                Pointer<ArrayData<NDIM, double>> p_acoef_data, p_bcoef_data;
                if (p_bc_coef)
                {
                    p_acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    p_bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double>> p_gcoef_data;
                    auto extended_p_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(p_bc_coef);
                    if (extended_p_bc_coef)
                    {
                        extended_p_bc_coef->clearTargetPatchDataIndex();
                        extended_p_bc_coef->setHomogeneousBc(homogeneous_bc);
                    }
                    p_bc_coef->setBcCoefs(
                        p_acoef_data, p_bcoef_data, p_gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                }

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const hier::Index<NDIM>& i = bc();
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const bool velocity_bc = (a == 1.0 && b == 0.0);
                    const bool traction_bc = (a == 0.0 && b == 1.0);
                    if (velocity_bc)
                    {
                        // The row prescribes the velocity, so the couplings of the row are dropped.
                        couplings.erase(get_row_key(i_s));
                        uu_matrix_coefs(i_s, 0) = 1.0;
                        for (int k = 1; k < uu_stencil_sz; ++k)
                        {
                            uu_matrix_coefs(i_s, k) = 0.0;
                        }
                        for (int k = 0; k < up_stencil_sz; ++k)
                        {
                            up_matrix_coefs(i_s, k) = 0.0;
                        }
                    }
                    else if (traction_bc)
                    {
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) +=
                                uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) +=
                                uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }

                        // The pressure in the ghost cell is p_G = f_i*p_I + f_g*g, in which p_I is the pressure in
                        // the interior cell abutting the boundary. Add the ghost coefficient times f_i to the
                        // coefficient of p_I. The term f_g*g is moved to the right-hand side.
                        if (p_bc_coef)
                        {
                            const double h = dx[bdry_normal_axis];
                            const double p_a = (*p_acoef_data)(i, 0);
                            const double p_b = (*p_bcoef_data)(i, 0);
                            const double f_i = -(p_a * h - 2.0 * p_b) / (p_a * h + 2.0 * p_b);
                            const int ghost_index = is_lower ? 0 : 1;
                            const int interior_index = is_lower ? 1 : 0;
                            up_matrix_coefs(i_s, interior_index) += f_i * up_matrix_coefs(i_s, ghost_index);
                            up_matrix_coefs(i_s, ghost_index) = 0.0;
                        }
                    }
                    else
                    {
                        TBOX_ERROR(
                            "StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp():\n"
                            << "  unsupported boundary condition coefficients (a, b) = (" << a << ", " << b
                            << ") for the normal velocity.\n"
                            << "  Only a prescribed velocity, (a, b) = (1, 0), or a prescribed traction, (a, b) = "
                               "(0, 1),\n"
                            << "  is supported.\n");
                    }
                }
            }
        }

        // Set matrix coefficients.
        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& ic = b();
                const SideIndex<NDIM> is(ic, axis, SideIndex<NDIM>::Lower);
                const int u_dof_index = (*u_dof_index_data)(is);
                if (UNLIKELY(ilower > u_dof_index || u_dof_index >= iupper)) continue;

                const int u_stencil_sz = uu_stencil_sz + up_stencil_sz;
                std::vector<double> u_mat_vals(u_stencil_sz);
                std::vector<int> u_mat_cols(u_stencil_sz);

                u_mat_vals[0] = uu_matrix_coefs(is, 0);
                u_mat_cols[0] = u_dof_index;
                for (unsigned int d = 0, uu_stencil_index = 1; d < NDIM; ++d)
                {
                    for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
                    {
                        u_mat_vals[uu_stencil_index] = uu_matrix_coefs(is, uu_stencil_index);
                        u_mat_cols[uu_stencil_index] = (*u_dof_index_data)(is + uu_stencil[uu_stencil_index]);
                    }
                }
                for (int side = 0, up_stencil_index = 0; side <= 1; ++side, ++up_stencil_index)
                {
                    u_mat_vals[uu_stencil_sz + side] = up_matrix_coefs(is, up_stencil_index);
                    u_mat_cols[uu_stencil_sz + side] = (*p_dof_index_data)(ic + up_stencil[axis][up_stencil_index]);
                }
                const auto row_couplings = couplings.find(get_row_key(is));
                if (row_couplings != couplings.end())
                {
                    // A coupling to a velocity next to the row is a coupling to a column that the row has already.
                    const int num_stencil_cols = static_cast<int>(u_mat_cols.size());
                    for (const auto& coupling : row_couplings->second)
                    {
                        const auto stencil_end = u_mat_cols.begin() + num_stencil_cols;
                        const auto column = std::find(u_mat_cols.begin(), stencil_end, coupling.first);
                        if (column == stencil_end)
                        {
                            u_mat_cols.push_back(coupling.first);
                            u_mat_vals.push_back(coupling.second);
                        }
                        else
                        {
                            u_mat_vals[column - u_mat_cols.begin()] += coupling.second;
                        }
                    }
                }

                const int u_num_vals = static_cast<int>(u_mat_vals.size());
                ierr =
                    MatSetValues(mat, 1, &u_dof_index, u_num_vals, u_mat_cols.data(), u_mat_vals.data(), INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }

        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            const int p_dof_index = (*p_dof_index_data)(ic);
            if (UNLIKELY(ilower > p_dof_index || p_dof_index >= iupper)) continue;

            const int p_stencil_sz = pu_stencil_sz + 1;
            std::vector<double> p_mat_vals(p_stencil_sz);
            std::vector<int> p_mat_cols(p_stencil_sz);

            for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
            {
                for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
                {
                    p_mat_vals[pu_stencil_index] = pu_matrix_coefs(ic, pu_stencil_index);
                    p_mat_cols[pu_stencil_index] = (*u_dof_index_data)(
                        SideIndex<NDIM>(ic + pu_stencil[pu_stencil_index], axis, SideIndex<NDIM>::Lower));
                }
            }
            p_mat_vals[pu_stencil_sz] = 0.0;
            p_mat_cols[pu_stencil_sz] = p_dof_index;

            ierr =
                MatSetValues(mat, 1, &p_dof_index, p_stencil_sz, p_mat_cols.data(), p_mat_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;
} // constructPatchLevelMACStokesOp

void
StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains(std::vector<std::set<int>>& is_overlap,
                                                                   std::vector<std::set<int>>& is_nonoverlap,
                                                                   const IntVector<NDIM>& box_size,
                                                                   const IntVector<NDIM>& overlap_size,
                                                                   const std::vector<int>& /*num_dofs_per_proc*/,
                                                                   int u_dof_index_idx,
                                                                   int p_dof_index_idx,
                                                                   Pointer<PatchLevel<NDIM>> patch_level,
                                                                   Pointer<CoarseFineBoundary<NDIM>> /*cf_boundary*/)
{
    // Clear previously stored index sets.
    for (auto& k : is_overlap)
    {
        k.clear();
    }
    is_overlap.clear();
    for (auto& k : is_nonoverlap)
    {
        k.clear();
    }
    is_nonoverlap.clear();

    // Create variables to keep track of whether a particular velocity location
    // is the "master" location.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, int>> patch_num_var = new SideVariable<NDIM, int>(
        "StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains()::"
        "patch_num_var");
    static const int patch_num_idx = var_db->registerPatchDataIndex(patch_num_var);
    patch_level->allocatePatchData(patch_num_idx);
    Pointer<SideVariable<NDIM, bool>> u_mastr_loc_var = new SideVariable<NDIM, bool>(
        "StaggeredStokesPETScMatUtilities::"
        "constructPatchLevelASMSubdomains()::u_"
        "mastr_loc_var");
    static const int u_mastr_loc_idx = var_db->registerPatchDataIndex(u_mastr_loc_var);
    patch_level->allocatePatchData(u_mastr_loc_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        Pointer<SideData<NDIM, int>> patch_num_data = patch->getPatchData(patch_num_idx);
        Pointer<SideData<NDIM, bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        patch_num_data->fillAll(patch_num);
        u_mastr_loc_data->fillAll(false);
    }

    // Synchronize the patch number at patch boundaries to determine which patch
    // owns a given DOF along patch boundaries.
    RefineAlgorithm<NDIM> bdry_synch_alg;
    bdry_synch_alg.registerRefine(patch_num_idx, patch_num_idx, patch_num_idx, nullptr, new SideSynchCopyFillPattern());
    bdry_synch_alg.createSchedule(patch_level)->fillData(0.0);

    // For a single patch in a periodic domain, the far side DOFs are not master.
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = patch_level->getGridGeometry();
    IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift(patch_level->getRatio());
    const BoxArray<NDIM>& domain_boxes = patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(domain_boxes.size() == 1);
#endif
    const hier::Index<NDIM>& domain_upper = domain_boxes[0].upper();

    // Determine the number of local DOFs.
    int local_dof_count = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const Box<NDIM>& patch_box = patch->getBox();
        const IntVector<NDIM> patch_size = patch_box.numberCells();

        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<SideData<NDIM, int>> patch_num_data = patch->getPatchData(patch_num_idx);
        Pointer<SideData<NDIM, bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const int upper_domain_side_idx = domain_upper(component_axis) + 1;
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> is(i, component_axis, SideIndex<NDIM>::Lower);
                bool fully_periodic_patch_in_axis =
                    periodic_shift(component_axis) && (patch_size(component_axis) == periodic_shift(component_axis));
                bool periodic_image = fully_periodic_patch_in_axis && (i(component_axis) == upper_domain_side_idx);
                if ((*patch_num_data)(is) == patch_num && !periodic_image)
                {
                    (*u_mastr_loc_data)(is) = true;
                    ++local_dof_count;
                }
            }
        }
        local_dof_count += CellGeometry<NDIM>::toCellBox(patch_box).size();
    }

    // Determine the subdomains associated with this processor.
    const int n_local_patches = patch_level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::vector<Box<NDIM>>> overlap_boxes(n_local_patches), nonoverlap_boxes(n_local_patches);
    int patch_counter = 0, subdomain_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        IndexUtilities::partitionPatchBox(
            overlap_boxes[patch_counter], nonoverlap_boxes[patch_counter], patch_box, box_size, overlap_size);
        subdomain_counter += overlap_boxes[patch_counter].size();
    }
    is_overlap.resize(subdomain_counter);
    is_nonoverlap.resize(subdomain_counter);

    // Fill in the IS'es.
    int nonoverlap_dof_counter = 0;
    subdomain_counter = 0, patch_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_patch_box[NDIM];
        for (int axis = 0; axis < NDIM; ++axis)
        {
            side_patch_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        }
        Pointer<SideData<NDIM, bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
#if !defined(NDEBUG)
        {
            const int u_data_depth = u_dof_data->getDepth();
            const int p_data_depth = p_dof_data->getDepth();
            TBOX_ASSERT(u_data_depth == 1);
            TBOX_ASSERT(p_data_depth == 1);
            TBOX_ASSERT(u_dof_data->getGhostCellWidth().min() >= overlap_size.max());
            TBOX_ASSERT(p_dof_data->getGhostCellWidth().min() >= overlap_size.max());
        }
#endif
        int n_patch_subdomains = static_cast<int>(nonoverlap_boxes[patch_counter].size());
        for (int k = 0; k < n_patch_subdomains; ++k, ++subdomain_counter)
        {
            // The nonoverlapping subdomains.
            const Box<NDIM>& sub_box = nonoverlap_boxes[patch_counter][k];
            Box<NDIM> side_sub_box[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_sub_box[axis] = SideGeometry<NDIM>::toSideBox(sub_box, axis);
            }

            // Get the local DOFs.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_sub_box[axis]); b; b++)
                {
                    const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                    const bool at_upper_subdomain_bdry = (i_s(axis) == side_sub_box[axis].upper(axis));
                    const bool at_upper_patch_bdry = (i_s(axis) == side_patch_box[axis].upper(axis));
                    if (!at_upper_subdomain_bdry || (at_upper_patch_bdry && (*u_mastr_loc_data)(i_s)))
                    {
                        const int dof_idx = (*u_dof_data)(i_s);
                        if (dof_idx >= 0) is_nonoverlap[subdomain_counter].insert(dof_idx);
                    }
                }
            }
            for (Box<NDIM>::Iterator b(sub_box); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const int dof_idx = (*p_dof_data)(i);
                if (dof_idx >= 0) is_nonoverlap[subdomain_counter].insert(dof_idx);
            }
            const int n_nonoverlap = static_cast<int>(is_nonoverlap[subdomain_counter].size());
            nonoverlap_dof_counter += n_nonoverlap;

            // The overlapping subdomains.
            const Box<NDIM>& overlap_sub_box = overlap_boxes[patch_counter][k];
            Box<NDIM> side_overlap_sub_box[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_overlap_sub_box[axis] = SideGeometry<NDIM>::toSideBox(overlap_sub_box, axis);
            }

            // Get the overlap DOFs.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_overlap_sub_box[axis]); b; b++)
                {
                    const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                    const int dof_idx = (*u_dof_data)(i_s);
                    if (dof_idx >= 0) is_overlap[subdomain_counter].insert(dof_idx);
                }
            }
            for (Box<NDIM>::Iterator b(overlap_sub_box); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const int dof_idx = (*p_dof_data)(i);
                if (dof_idx >= 0) is_overlap[subdomain_counter].insert(dof_idx);
            }
        }
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(local_dof_count == nonoverlap_dof_counter);
#else
    NULL_USE(local_dof_count);
    NULL_USE(nonoverlap_dof_counter);
#endif

    // Deallocate patch_num variable data.
    patch_level->deallocatePatchData(patch_num_idx);
    patch_level->deallocatePatchData(u_mastr_loc_idx);
    return;
} // constructPatchLevelASMSubdomains

void
StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
    std::vector<std::set<int>>& is_field,
    std::vector<std::string>& is_field_name,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    // Destroy the previously stored IS'es
    for (auto& k : is_field)
    {
        k.clear();
    }
    is_field.clear();
    is_field_name.clear();

    // Resize vectors
    is_field.resize(2);
    is_field_name.resize(2);

    // Name of the fields.
    static const int U_FIELD_IDX = 0;
    static const int P_FIELD_IDX = 1;
    is_field_name[U_FIELD_IDX] = "velocity";
    is_field_name[P_FIELD_IDX] = "pressure";

    // DOFs on this processor.
    const int mpi_rank = IBTK_MPI::getRank();
    const int n_local_dofs = num_dofs_per_proc[mpi_rank];

    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int last_local_dof = first_local_dof + n_local_dofs;

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_patch_box[NDIM];
        for (int axis = 0; axis < NDIM; ++axis)
        {
            side_patch_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        }

        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
#if !defined(NDEBUG)
        const int u_data_depth = u_dof_data->getDepth();
        const int p_data_depth = p_dof_data->getDepth();
        TBOX_ASSERT(u_data_depth == 1);
        TBOX_ASSERT(p_data_depth == 1);
#endif

        // Get the local velocity DOFs.
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(side_patch_box[axis]); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                const int dof_idx = (*u_dof_data)(i_s);
                if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
                {
                    is_field[0].insert(dof_idx);
                }
            }
        }

        // Get the local pressure DOFs.
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            const int dof_idx = (*p_dof_data)(i);
            if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
            {
                is_field[1].insert(dof_idx);
            }
        }
    }

    return;
} // constructPatchLevelFields

void
StaggeredStokesPETScMatUtilities::constructProlongationOp(Mat& mat,
                                                          const std::string& u_op_type,
                                                          const std::string& p_op_type,
                                                          int u_dof_index_idx,
                                                          int p_dof_index_idx,
                                                          const std::vector<int>& num_fine_dofs_per_proc,
                                                          const std::vector<int>& num_coarse_dofs_per_proc,
                                                          Pointer<PatchLevel<NDIM>> fine_patch_level,
                                                          Pointer<PatchLevel<NDIM>> coarse_patch_level,
                                                          const AO& coarse_level_ao,
                                                          const int u_coarse_ao_offset,
                                                          const int p_coarse_ao_offset)
{
    int ierr;
    Mat p_prolong_mat = nullptr;
    PETScMatUtilities::constructProlongationOp(mat,
                                               u_op_type,
                                               u_dof_index_idx,
                                               num_fine_dofs_per_proc,
                                               num_coarse_dofs_per_proc,
                                               fine_patch_level,
                                               coarse_patch_level,
                                               coarse_level_ao,
                                               u_coarse_ao_offset);

    PETScMatUtilities::constructProlongationOp(p_prolong_mat,
                                               p_op_type,
                                               p_dof_index_idx,
                                               num_fine_dofs_per_proc,
                                               num_coarse_dofs_per_proc,
                                               fine_patch_level,
                                               coarse_patch_level,
                                               coarse_level_ao,
                                               p_coarse_ao_offset);

    // P{u,p} = (P_u + P_p){u,p}
    ierr = MatAXPY(mat, 1.0, p_prolong_mat, DIFFERENT_NONZERO_PATTERN);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&p_prolong_mat);
    IBTK_CHKERRQ(ierr);

} // constructPatchLevelProlongationOp

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
