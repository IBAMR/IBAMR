// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibamr/AcousticStreamingPETScMatUtilities.h"

#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/PETScMatUtilities.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/ibtk_utilities.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellGeometry.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchLevel.h"
#include "ProcessorMapping.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "tbox/Array.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "petscmat.h"
#include <petsclog.h>

#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{

static const int REAL = 0;
static const int IMAG = 1;
static const int CENTER = 0;
static const int PLUS = 1;
static const int MINUS = 2;

inline hier::Index<NDIM>
get_shift(int dir, int shift)
{
    SAMRAI::hier::Index<NDIM> iv(0);
    iv(dir) = shift;
    return iv;
} // get_shift

inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
} // compute_tangential_extension

void
compute_grad_p_sc_matrix_coefficients(SideData<NDIM, double>& matrix_coefs,
                                      Pointer<Patch<NDIM> > patch,
                                      const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      double data_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
#endif
    const int stencil_sz = 2;

#if !defined(NDEBUG)
    TBOX_ASSERT(matrix_coefs.getDepth() == stencil_sz);
#endif
    matrix_coefs.fillAll(0.0);

    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Compute all matrix coefficients, including those on the physical
    // boundary; however, do not yet take physical boundary conditions into
    // account.  Boundary conditions are handled subsequently.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        // grad p
        matrix_coefs.getArrayData(axis).fill(-1.0 / dx[axis], /*side*/ 0);
        matrix_coefs.getArrayData(axis).fill(+1.0 / dx[axis], /*side*/ 1);
    }

    // Data structures required to set physical boundary conditions.
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
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
    // conditions along boundaries which ARE aligned with the data axis.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;

            if (bdry_normal_axis != axis) continue;

            const Box<NDIM> bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data;

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

            // Modify the matrix coefficients to account for homogeneous
            // boundary conditions.
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const hier::Index<NDIM>& i = bc();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                const double& a = (*acoef_data)(i, 0);
                const bool velocity_bc = (a == 1.0 || IBTK::rel_equal_eps(a, 1.0));
#if !defined(NDEBUG)
                TBOX_ASSERT(velocity_bc);
#endif
                for (int k = 0; k < stencil_sz; ++k)
                {
                    matrix_coefs(i_s, k) = 0.0;
                }
            }
        }
    }

    return;

} // compute_grad_p_sc_matrix_coefficients

// grad/div slot
inline int
gd_slot(const int axis, const int part, const int k)
{
    // k = 0 center, 1 plus, 2 minus
    return (part * NDIM + axis) * 3 + k;
} // gd_slot

inline double
robin_factor(const double a, const double b, const double h)
{
    return b == 0.0 ? -1.0 : -(a * h - 2.0 * b) / (a * h + 2.0 * b);
} // robin_factor

template <class GetBcCoef>
void
set_cc_centered_derivative_coefs(CellData<NDIM, double>& coefs,
                                 Pointer<Patch<NDIM> > patch,
                                 GetBcCoef get_bc_coef,
                                 const double data_time,
                                 const double* const dx)
{
    const Box<NDIM>& patch_box = patch->getBox();

    for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
    {
        const CellIndex<NDIM>& i = b();

        for (int part = 0; part < 2; ++part)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                coefs(i, gd_slot(axis, part, CENTER)) = 0.0;
                coefs(i, gd_slot(axis, part, PLUS)) = +1.0 / (2.0 * dx[axis]);
                coefs(i, gd_slot(axis, part, MINUS)) = -1.0 / (2.0 * dx[axis]);
            }
        }
    }

    const auto physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();

    for (int part = 0; part < 2; ++part)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const hier::BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const int axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            RobinBcCoefStrategy<NDIM>* bc_coef = get_bc_coef(part, axis);
            if (!bc_coef) continue;

            const auto trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const auto bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            static const bool homogeneous_bc = true;

            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }

            bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);

            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            for (hier::Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const CellIndex<NDIM>& i_bdry = bc();

                CellIndex<NDIM> i_intr = i_bdry;
                if (!is_lower) i_intr(axis) -= 1;

                const double a = (*acoef_data)(i_bdry, 0);
                const double b = (*bcoef_data)(i_bdry, 0);
                const double r = robin_factor(a, b, dx[axis]);

                const int c = gd_slot(axis, part, CENTER);
                const int p = gd_slot(axis, part, PLUS);
                const int m = gd_slot(axis, part, MINUS);

                if (is_lower)
                {
                    coefs(i_intr, c) += r * coefs(i_intr, m);
                    coefs(i_intr, m) = 0.0;
                }
                else
                {
                    coefs(i_intr, c) += r * coefs(i_intr, p);
                    coefs(i_intr, p) = 0.0;
                }
            }
        }
    }
}

} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
AcousticStreamingPETScMatUtilities::constructPatchLevelFOAcousticStreamingOp(
    Mat& mat,
    double omega,
    double sound_speed,
    int rho_idx,
    int mu_idx,
    int lambda_idx,
    int chi_idx,
    const std::array<std::vector<RobinBcCoefStrategy<NDIM>*>, 2>& u_bc_coefs,
    const std::array<RobinBcCoefStrategy<NDIM>*, 2>& p_bc_coefs,
    double data_time,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    Pointer<PatchLevel<NDIM> > patch_level,
    IBTK::VCInterpType mu_interp_type)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > u_dof_index_var;
    var_db->mapIndexToVariable(u_dof_index_idx, u_dof_index_var);
    Pointer<SideVariable<NDIM, int> > u_dof_index_sc_var = u_dof_index_var;
    Pointer<CellVariable<NDIM, int> > u_dof_index_cc_var = u_dof_index_var;
    Pointer<Variable<NDIM> > p_dof_index_var;
    var_db->mapIndexToVariable(p_dof_index_idx, p_dof_index_var);
    Pointer<CellVariable<NDIM, int> > p_dof_index_cc_var = p_dof_index_var;

    if (u_dof_index_sc_var && p_dof_index_cc_var)
    {
        constructPatchLevelFOAcousticStreamingOp_MAC(mat,
                                                     omega,
                                                     sound_speed,
                                                     rho_idx,
                                                     mu_idx,
                                                     lambda_idx,
                                                     chi_idx,
                                                     u_bc_coefs,
                                                     data_time,
                                                     num_dofs_per_proc,
                                                     u_dof_index_idx,
                                                     p_dof_index_idx,
                                                     patch_level,
                                                     mu_interp_type);
    }
    else if (u_dof_index_cc_var && p_dof_index_cc_var)
    {
        constructPatchLevelFOAcousticStreamingOp_cell(mat,
                                                      omega,
                                                      sound_speed,
                                                      rho_idx,
                                                      mu_idx,
                                                      lambda_idx,
                                                      chi_idx,
                                                      u_bc_coefs,
                                                      p_bc_coefs,
                                                      data_time,
                                                      num_dofs_per_proc,
                                                      u_dof_index_idx,
                                                      p_dof_index_idx,
                                                      patch_level,
                                                      mu_interp_type);
    }
    else
    {
        TBOX_ERROR("AcousticStreamingPETScMatUtilities::constructPatchLevelFOAcousticStreamingOp():\n"
                   << "  unsupported data centering types for variables " << u_dof_index_var->getName() << " and "
                   << p_dof_index_var->getName() << "\n");
    }

    return;
} // constructPatchLevelFOAcousticStreamingOp

/////////////////////////////// PRIVATE ///////////////////////////////////////

void
AcousticStreamingPETScMatUtilities::constructPatchLevelFOAcousticStreamingOp_MAC(
    Mat& mat,
    double omega,
    double sound_speed,
    int rho_idx,
    int mu_idx,
    int lambda_idx,
    int chi_idx,
    const std::array<std::vector<RobinBcCoefStrategy<NDIM>*>, 2>& u_bc_coefs,
    double data_time,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    Pointer<PatchLevel<NDIM> > patch_level,
    IBTK::VCInterpType mu_interp_type)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Determine the index ranges.
    const int mpi_rank = IBTK_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int proc_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int proc_upper = proc_lower + n_local;
    const int n_total = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(n_local, 0), o_nnz(n_local, 0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        const int u_dof_index_depth = u_dof_index_data->getDepth();
        const int p_dof_index_depth = p_dof_index_data->getDepth();
#if !defined(NDEBUG)
        TBOX_ASSERT(u_dof_index_depth == 2);
        TBOX_ASSERT(p_dof_index_depth == 2);
#endif
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const hier::Index<NDIM>& cc = b();
                const SideIndex<NDIM> i(cc, axis, SideIndex<NDIM>::Lower);

                for (int comp = 0; comp < 2; ++comp)
                {
                    const int i_dof_index = (*u_dof_index_data)(i, comp);

                    // Non-zero structure for u
                    if (proc_lower <= i_dof_index && i_dof_index < proc_upper)
                    {
                        const int u_local_idx = i_dof_index - proc_lower;
                        d_nnz[u_local_idx] += 1;

                        const int other_comp = (comp == REAL ? IMAG : REAL);

                        // Stencil for divergence of viscous and dilatational stress tensor
                        // acting on the other component of velocity.
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (d == axis)
                            {
                                const int i_dof = (*u_dof_index_data)(i, other_comp);
                                if (i_dof >= proc_lower && i_dof < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }

                                hier::Index<NDIM> shift_axis = get_shift(axis, 1);

                                const int i_dof_hi = (*u_dof_index_data)(i + shift_axis, other_comp);
                                if (i_dof_hi >= proc_lower && i_dof_hi < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }
                                const int i_dof_lo = (*u_dof_index_data)(i - shift_axis, other_comp);
                                if (i_dof_lo >= proc_lower && i_dof_lo < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }
                            }
                            else
                            {
                                hier::Index<NDIM> shift_d_plus = get_shift(d, 1);
                                hier::Index<NDIM> shift_d_minus = get_shift(d, -1);
                                hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);

                                const int i_dof_hi = (*u_dof_index_data)(i + shift_d_plus, other_comp);
                                if (i_dof_hi >= proc_lower && i_dof_hi < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }
                                const int i_dof_lo = (*u_dof_index_data)(i + shift_d_minus, other_comp);
                                if (i_dof_lo >= proc_lower && i_dof_lo < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }

                                const SideIndex<NDIM> j_se(cc, d, SideIndex<NDIM>::Lower);
                                const int j_se_dof_index = (*u_dof_index_data)(j_se, other_comp);
                                if (j_se_dof_index >= proc_lower && j_se_dof_index < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }

                                const SideIndex<NDIM> j_sw(cc + shift_axis_minus, d, SideIndex<NDIM>::Lower);
                                const int j_sw_dof_index = (*u_dof_index_data)(j_sw, other_comp);
                                if (j_sw_dof_index >= proc_lower && j_sw_dof_index < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }

                                const SideIndex<NDIM> j_ne(cc, d, SideIndex<NDIM>::Upper);
                                const int j_ne_dof_index = (*u_dof_index_data)(j_ne, other_comp);
                                if (j_ne_dof_index >= proc_lower && j_ne_dof_index < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }

                                const SideIndex<NDIM> j_nw(cc + shift_axis_minus, d, SideIndex<NDIM>::Upper);
                                const int j_nw_dof_index = (*u_dof_index_data)(j_nw, other_comp);
                                if (j_nw_dof_index >= proc_lower && j_nw_dof_index < proc_upper)
                                {
                                    d_nnz[u_local_idx] += 1;
                                }
                                else
                                {
                                    o_nnz[u_local_idx] += 1;
                                }
                            }
                        }

                        // Stencil for pressure gradient acting on the other component of pressure.
                        for (int side = 0; side <= 1; ++side)
                        {
                            hier::Index<NDIM> shift_axis = get_shift(axis, side - 1);
                            const int p_dof_index = (*p_dof_index_data)(cc + shift_axis, other_comp);
                            if (LIKELY(p_dof_index >= proc_lower && p_dof_index < proc_upper))
                            {
                                d_nnz[u_local_idx] += 1;
                            }
                            else
                            {
                                o_nnz[u_local_idx] += 1;
                            }
                        }

                        d_nnz[u_local_idx] = std::min(n_local, d_nnz[u_local_idx]);
                        o_nnz[u_local_idx] = std::min(n_total - n_local, o_nnz[u_local_idx]);
                    }
                }
            }
        }

        // Stencil for divergence of velocity acting on the other component.
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            for (int comp = 0; comp < 2; ++comp)
            {
                const int p_dof_index = (*p_dof_index_data)(ic, comp);
                if (LIKELY(p_dof_index >= proc_lower && p_dof_index < proc_upper))
                {
                    const int p_local_idx = p_dof_index - proc_lower;
                    d_nnz[p_local_idx] += 1;

                    const int other_comp = (comp == REAL ? IMAG : REAL);

                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (int side = 0; side <= 1; ++side)
                        {
                            hier::Index<NDIM> shift_axis = get_shift(axis, side);
                            const int u_dof_index = (*u_dof_index_data)(
                                SideIndex<NDIM>(ic + shift_axis, axis, SideIndex<NDIM>::Lower), other_comp);
                            if (LIKELY(u_dof_index >= proc_lower && u_dof_index < proc_upper))
                            {
                                d_nnz[p_local_idx] += 1;
                            }
                            else
                            {
                                o_nnz[p_local_idx] += 1;
                            }
                        }
                    }
                    d_nnz[p_local_idx] = std::min(n_local, d_nnz[p_local_idx]);
                    o_nnz[p_local_idx] = std::min(n_total - n_local, o_nnz[p_local_idx]);
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        n_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        0,
                        n_local ? &d_nnz[0] : nullptr,
                        0,
                        n_local ? &o_nnz[0] : nullptr,
                        &mat);
    IBTK_CHKERRQ(ierr);

// Set some general matrix options.
#if !defined(NDEBUG)
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    using StencilMapType = std::map<hier::Index<NDIM>, int, IndexFortranOrder>;
    static std::vector<StencilMapType> stencil_map_vec;
    static const int uu_stencil_sz = (2 * NDIM + 1) + 4 * (NDIM - 1); // 9 for 2D and 15 for 3D
    static const int up_stencil_sz = 2;
    static const int pu_stencil_sz = 2 * NDIM;
    static const hier::Index<NDIM> ORIGIN(0);

#if (NDIM == 2)
    // Create stencil dictionary.
    enum DIRECTIONS
    {
        CENTER = 0,
        EAST = 1,
        WEST = 2,
        NORTH = 3,
        SOUTH = 4,
        NORTHEAST = 5,
        NORTHWEST = 6,
        SOUTHEAST = 7,
        SOUTHWEST = 8,
        X = 0,
        Y = 1
    };
    IBTK_DO_ONCE(static StencilMapType sm; sm[ORIGIN] = CENTER; sm[get_shift(X, 1)] = EAST; sm[get_shift(X, -1)] = WEST;
                 sm[get_shift(Y, 1)] = NORTH;
                 sm[get_shift(Y, -1)] = SOUTH;
                 sm[get_shift(Y, 1) + get_shift(X, 1)] = NORTHEAST;
                 sm[get_shift(Y, 1) + get_shift(X, -1)] = NORTHWEST;
                 sm[get_shift(Y, -1) + get_shift(X, 1)] = SOUTHEAST;
                 sm[get_shift(Y, -1) + get_shift(X, -1)] = SOUTHWEST;
                 stencil_map_vec.push_back(sm););

#elif (NDIM == 3)
    // In 3D, the shifted directions depend on the axis under consideration
    enum COMMONDIRECTIONS
    {
        CENTER = 0,
        EAST = 1,
        WEST = 2,
        NORTH = 3,
        SOUTH = 4,
        TOP = 5,
        BOTTOM = 6,
        X = 0,
        Y = 1,
        Z = 2
    };
    IBTK_DO_ONCE(for (int axis = 0; axis < NDIM; ++axis) {
        static StencilMapType sm;
        // Common to all axes
        sm[ORIGIN] = CENTER;
        sm[get_shift(X, 1)] = EAST;
        sm[get_shift(X, -1)] = WEST;
        sm[get_shift(Y, 1)] = NORTH;
        sm[get_shift(Y, -1)] = SOUTH;
        sm[get_shift(Z, 1)] = TOP;
        sm[get_shift(Z, -1)] = BOTTOM;

        // Specific to certain axes
        int idx = BOTTOM;
        for (int d = 0; d < NDIM; ++d)
        {
            if (d == axis) continue;
            idx += 1;
            sm[get_shift(axis, 1) + get_shift(d, 1)] = idx;
            idx += 1;
            sm[get_shift(axis, -1) + get_shift(d, 1)] = idx;
            idx += 1;
            sm[get_shift(axis, 1) + get_shift(d, -1)] = idx;
            idx += 1;
            sm[get_shift(axis, -1) + get_shift(d, -1)] = idx;
        }
        stencil_map_vec.push_back(sm);
    });
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to first order acoustic streaming system.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

        Pointer<SideData<NDIM, double> > chi_data = nullptr;
        if (chi_idx != IBTK::invalid_index) chi_data = patch->getPatchData(chi_idx);

        const IntVector<NDIM> no_ghosts(0);
        SideData<NDIM, double> uur_matrix_coefs(patch_box, uu_stencil_sz, no_ghosts);
        SideData<NDIM, double> uui_matrix_coefs(patch_box, uu_stencil_sz, no_ghosts);
        SideData<NDIM, double> upr_matrix_coefs(patch_box, up_stencil_sz, no_ghosts);
        SideData<NDIM, double> upi_matrix_coefs(patch_box, up_stencil_sz, no_ghosts);
        CellData<NDIM, double> pu_matrix_coefs(patch_box, pu_stencil_sz, no_ghosts);

        // Compute matrix coefficients corresponding to viscous and dilatational stress for both velocity components.
        // Do not consider Dirichlet BCs in the matrix coefficients yet, as we will consider them after matrix
        // assembly.
        PoissonUtilities::computeVCSCViscousDilatationalOpMatrixCoefficients(uur_matrix_coefs,
                                                                             patch,
                                                                             stencil_map_vec,
                                                                             u_bc_coefs[REAL],
                                                                             data_time,
                                                                             mu_idx,
                                                                             lambda_idx,
                                                                             mu_interp_type,
                                                                             /*consider_dirichlet_bcs*/ false);
        PoissonUtilities::computeVCSCViscousDilatationalOpMatrixCoefficients(uui_matrix_coefs,
                                                                             patch,
                                                                             stencil_map_vec,
                                                                             u_bc_coefs[IMAG],
                                                                             data_time,
                                                                             mu_idx,
                                                                             lambda_idx,
                                                                             mu_interp_type,
                                                                             /*consider_dirichlet_bcs*/ false);

        // Compute matrix coefficients corresponding to the gradient of both components of pressure.
        // Note: here we give BCs of the other component of velocity to modify the pressure gradient.
        compute_grad_p_sc_matrix_coefficients(upr_matrix_coefs, patch, u_bc_coefs[IMAG], data_time);
        compute_grad_p_sc_matrix_coefficients(upi_matrix_coefs, patch, u_bc_coefs[REAL], data_time);

// Set matrix coefficients for velocity DOFs.
#if (NDIM == 2)
        StencilMapType& stencil_map = stencil_map_vec[0];
#endif
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
#if (NDIM == 3)
            StencilMapType& stencil_map = stencil_map_vec[axis];
#endif
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& ic = b();
                const SideIndex<NDIM> is(ic, axis, SideIndex<NDIM>::Lower);
                double chi = chi_data.isNull() ? 0.0 : (*chi_data)(is, 0);

                for (int comp = 0; comp < 2; ++comp)
                {
                    const int other_comp = (comp == REAL ? IMAG : REAL);
                    const double freq = (comp == REAL ? omega : -omega);
                    SideData<NDIM, double>& uu_matrix_coefs = (comp == REAL ? uui_matrix_coefs : uur_matrix_coefs);
                    SideData<NDIM, double>& up_matrix_coefs = (comp == REAL ? upi_matrix_coefs : upr_matrix_coefs);

                    const int u_dof_index = (*u_dof_index_data)(is, comp);
                    if (UNLIKELY(proc_lower > u_dof_index || u_dof_index >= proc_upper)) continue;

                    const int u_stencil_sz = 1 + uu_stencil_sz + up_stencil_sz; // 1 for vel comp + uu_stencil_size for
                                                                                // other vel comps + 2 pressure gradient
                    std::vector<double> u_mat_vals(u_stencil_sz);
                    std::vector<int> u_mat_cols(u_stencil_sz);

                    int idx = 0;
                    u_mat_vals[idx] = freq * (*rho_data)(is, 0);
                    u_mat_cols[idx] = u_dof_index;

                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[ORIGIN]) + chi;
                            u_mat_cols[idx] = (*u_dof_index_data)(is, other_comp);

                            const hier::Index<NDIM> shift_axis_plus = get_shift(axis, 1);

                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[shift_axis_plus]);
                            u_mat_cols[idx] = (*u_dof_index_data)(is + shift_axis_plus, other_comp);

                            const hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);

                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[shift_axis_minus]);
                            u_mat_cols[idx] = (*u_dof_index_data)(is + shift_axis_minus, other_comp);
                        }
                        else
                        {
                            const hier::Index<NDIM> shift_d_plus = get_shift(d, 1);
                            const hier::Index<NDIM> shift_d_minus = get_shift(d, -1);
                            const hier::Index<NDIM> shift_axis_plus = get_shift(axis, 1);
                            const hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);

                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[shift_d_plus]);
                            u_mat_cols[idx] = (*u_dof_index_data)(is + shift_d_plus, other_comp);

                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[shift_d_minus]);
                            u_mat_cols[idx] = (*u_dof_index_data)(is + shift_d_minus, other_comp);

                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[shift_d_plus + shift_axis_plus]);
                            const SideIndex<NDIM> ne(ic, d, SideIndex<NDIM>::Upper);
                            u_mat_cols[idx] = (*u_dof_index_data)(ne, other_comp);

                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[shift_d_plus + shift_axis_minus]);
                            const SideIndex<NDIM> nw(ic + shift_axis_minus, d, SideIndex<NDIM>::Upper);
                            u_mat_cols[idx] = (*u_dof_index_data)(nw, other_comp);

                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[shift_d_minus + shift_axis_plus]);
                            const SideIndex<NDIM> se(ic, d, SideIndex<NDIM>::Lower);
                            u_mat_cols[idx] = (*u_dof_index_data)(se, other_comp);

                            idx += 1;
                            u_mat_vals[idx] = uu_matrix_coefs(is, stencil_map[shift_d_minus + shift_axis_minus]);
                            const SideIndex<NDIM> sw(ic + shift_axis_minus, d, SideIndex<NDIM>::Lower);
                            u_mat_cols[idx] = (*u_dof_index_data)(sw, other_comp);
                        }
                    }

                    // Stencil for pressure gradient acting on the other component of pressure.
                    for (int side = 0; side <= 1; ++side)
                    {
                        hier::Index<NDIM> shift_axis = get_shift(axis, side - 1);

                        idx += 1;
                        u_mat_vals[idx] = up_matrix_coefs(is, side);
                        u_mat_cols[idx] = (*p_dof_index_data)(ic + shift_axis, other_comp);
                    }

#if !defined(NDEBUG)
                    TBOX_ASSERT(idx == (u_stencil_sz - 1));
#endif
                    ierr =
                        MatSetValues(mat, 1, &u_dof_index, u_stencil_sz, &u_mat_cols[0], &u_mat_vals[0], INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
        }

        // Set matrix coefficients for pressure DOFs.
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();

            for (int comp = 0; comp < 2; ++comp)
            {
                const int p_dof_index = (*p_dof_index_data)(ic, comp);
                if (LIKELY(p_dof_index >= proc_lower && p_dof_index < proc_upper))
                {
                    const int p_stencil_sz = pu_stencil_sz + 1;
                    std::vector<double> p_mat_vals(p_stencil_sz);
                    std::vector<int> p_mat_cols(p_stencil_sz);

                    const int other_comp = (comp == REAL ? IMAG : REAL);
                    const double p_coef =
                        (comp == REAL ? omega / (sound_speed * sound_speed) : -omega / (sound_speed * sound_speed));

                    int idx = 0;
                    p_mat_vals[idx] = p_coef;
                    p_mat_cols[idx] = p_dof_index;

                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (int side = 0; side <= 1; ++side)
                        {
                            hier::Index<NDIM> shift_axis = get_shift(axis, side);
                            SideIndex<NDIM> is(ic + shift_axis, axis, SideIndex<NDIM>::Lower);
                            double sgn = (side == 0 ? -1.0 : 1.0);

                            idx += 1;
                            p_mat_vals[idx] = sgn * (*rho_data)(is, 0) / dx[axis];
                            p_mat_cols[idx] = (*u_dof_index_data)(is, other_comp);
                        }
                    }

                    ierr =
                        MatSetValues(mat, 1, &p_dof_index, p_stencil_sz, &p_mat_cols[0], &p_mat_vals[0], INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    // Modify the matrix coefficients for Dirichlet BCs for the normal component of the velocity
    std::vector<int> dirichlet_rows;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();

        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;

                if (bdry_normal_axis != axis) continue;

                const Box<NDIM> bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
                const BoundaryBox<NDIM> trimmed_bdry_box =
                    PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

                Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                for (int comp = 0; comp < 2; ++comp)
                {
                    auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[comp][axis]);
                    if (extended_bc_coef)
                    {
                        extended_bc_coef->clearTargetPatchDataIndex();
                        extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                    }
                    u_bc_coefs[comp][axis]->setBcCoefs(
                        acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                    if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                    // Modify the matrix coefficients to account for homogeneous
                    // boundary conditions.
                    //
                    // If b == 0 we enforce Dirichlet condition directly.
                    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                    {
                        const hier::Index<NDIM>& i = bc();
                        const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);

                        const int u_dof_index = (*u_dof_index_data)(i_s, comp);
                        if (UNLIKELY(proc_lower > u_dof_index || u_dof_index >= proc_upper)) continue;

                        const double& b = (*bcoef_data)(i, 0);
                        if (IBTK::abs_equal_eps(b, 0.0))
                        {
#if !defined(NDEBUG)
                            const double& a = (*acoef_data)(i, 0);
                            TBOX_ASSERT(IBTK::rel_equal_eps(a, 1.0));
#endif
                            dirichlet_rows.push_back(u_dof_index);
                        }
                    }
                }
            }
        }
    }

    // Send-in the (local) Dirichlet row information to PETSc Mat
    const int num_dirichlet_rows = dirichlet_rows.size();
    ierr = MatZeroRows(mat, num_dirichlet_rows, &dirichlet_rows[0], /*diag_value*/ 1.0, NULL, NULL);
    IBTK_CHKERRQ(ierr);
    return;
} // constructPatchLevelFOAcousticStreamingOp_MAC

void
AcousticStreamingPETScMatUtilities::constructPatchLevelFOAcousticStreamingOp_cell(
    Mat& mat,
    const double omega,
    const double sound_speed,
    const int rho_idx,
    const int mu_idx,
    const int lambda_idx,
    const int chi_idx,
    const std::array<std::vector<RobinBcCoefStrategy<NDIM>*>, 2>& u_bc_coefs,
    const std::array<RobinBcCoefStrategy<NDIM>*, 2>& p_bc_coefs,
    const double data_time,
    const std::vector<int>& num_dofs_per_proc,
    const int u_dof_index_idx,
    const int p_dof_index_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
    IBTK::VCInterpType mu_interp_type)
{
    auto other_part = [](const int part) { return part == REAL ? IMAG : REAL; };
    auto u_depth = [](const int axis, const int part) { return 2 * axis + part; };

    using StencilMapType = IBTK::PoissonUtilities::CCVCStencilMap;
    const auto make_key = IBTK::PoissonUtilities::makeCCVCStencilKey;

    static const hier::Index<NDIM> ORIGIN(0);
    static const int uu_stencil_sz = 1 + 2 * NDIM + 9 * (NDIM - 1);
    static const int up_stencil_sz = 3;        // p(i), p(i+e_axis), p(i-e_axis)
    static const int pu_stencil_sz = 3 * NDIM; // u_axis(i), u_axis(i+e_axis), u_axis(i-e_axis)

    static std::vector<StencilMapType> stencil_map_vec;
    if (stencil_map_vec.empty())
    {
        stencil_map_vec.resize(NDIM);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            StencilMapType& sm = stencil_map_vec[axis];
            int idx = 0;

            sm[make_key(ORIGIN, axis)] = idx++;

            for (int d = 0; d < NDIM; ++d)
            {
                sm[make_key(get_shift(d, +1), axis)] = idx++;
                sm[make_key(get_shift(d, -1), axis)] = idx++;
            }

            for (int comp = 0; comp < NDIM; ++comp)
            {
                if (comp == axis) continue;

                const auto ea_p = get_shift(axis, +1);
                const auto ea_m = get_shift(axis, -1);
                const auto ec_p = get_shift(comp, +1);
                const auto ec_m = get_shift(comp, -1);

                sm[make_key(ea_p + ec_p, comp)] = idx++;
                sm[make_key(ea_m + ec_p, comp)] = idx++;
                sm[make_key(ea_p, comp)] = idx++;
                sm[make_key(ea_m, comp)] = idx++;
                sm[make_key(ea_p + ec_m, comp)] = idx++;
                sm[make_key(ea_m + ec_m, comp)] = idx++;

                sm[make_key(ec_p, comp)] = idx++;
                sm[make_key(ORIGIN, comp)] = idx++;
                sm[make_key(ec_m, comp)] = idx++;
            }

            TBOX_ASSERT(idx == uu_stencil_sz);
        }
    }

    int ierr = 0;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    const int mpi_rank = IBTK_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int proc_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int proc_upper = proc_lower + n_local;
    const int n_total = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    std::vector<int> d_nnz(n_local, 0), o_nnz(n_local, 0);

    auto add_nnz = [&](const int row_local, const int col)
    {
        if (proc_lower <= col && col < proc_upper)
            d_nnz[row_local] += 1;
        else
            o_nnz[row_local] += 1;
    };

    for (hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        auto patch = patch_level->getPatch(p());
        const auto& patch_box = patch->getBox();

        Pointer<CellData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);

        for (hier::Box<NDIM>::Iterator b(pdat::CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const pdat::CellIndex<NDIM>& i = b();

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (int part = 0; part < 2; ++part)
                {
                    const int row = (*u_dof_index_data)(i, u_depth(axis, part));
                    if (!(proc_lower <= row && row < proc_upper)) continue;

                    const int row_local = row - proc_lower;
                    const int opart = other_part(part);

                    add_nnz(row_local, row);
                    add_nnz(row_local, (*u_dof_index_data)(i, u_depth(axis, opart)));

                    for (int d = 0; d < NDIM; ++d)
                    {
                        add_nnz(row_local, (*u_dof_index_data)(i + get_shift(d, +1), u_depth(axis, opart)));
                        add_nnz(row_local, (*u_dof_index_data)(i + get_shift(d, -1), u_depth(axis, opart)));
                    }

                    for (int comp = 0; comp < NDIM; ++comp)
                    {
                        if (comp == axis) continue;

                        const auto ea_p = get_shift(axis, +1);
                        const auto ea_m = get_shift(axis, -1);
                        const auto ec_p = get_shift(comp, +1);
                        const auto ec_m = get_shift(comp, -1);

                        const hier::Index<NDIM> shifts[9] = { ea_p + ec_p, ea_m + ec_p, ea_p,   ea_m, ea_p + ec_m,
                                                              ea_m + ec_m, ec_p,        ORIGIN, ec_m };

                        for (int k = 0; k < 9; ++k)
                            add_nnz(row_local, (*u_dof_index_data)(i + shifts[k], u_depth(comp, opart)));
                    }

                    add_nnz(row_local, (*p_dof_index_data)(i, opart));
                    add_nnz(row_local, (*p_dof_index_data)(i + get_shift(axis, +1), opart));
                    add_nnz(row_local, (*p_dof_index_data)(i + get_shift(axis, -1), opart));

                    d_nnz[row_local] = std::min(n_local, d_nnz[row_local]);
                    o_nnz[row_local] = std::min(n_total - n_local, o_nnz[row_local]);
                }
            }

            for (int part = 0; part < 2; ++part)
            {
                const int row = (*p_dof_index_data)(i, part);
                if (!(proc_lower <= row && row < proc_upper)) continue;

                const int row_local = row - proc_lower;
                const int opart = other_part(part);

                add_nnz(row_local, row);

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    add_nnz(row_local, (*u_dof_index_data)(i, u_depth(axis, opart)));
                    add_nnz(row_local, (*u_dof_index_data)(i + get_shift(axis, +1), u_depth(axis, opart)));
                    add_nnz(row_local, (*u_dof_index_data)(i + get_shift(axis, -1), u_depth(axis, opart)));
                }

                d_nnz[row_local] = std::min(n_local, d_nnz[row_local]);
                o_nnz[row_local] = std::min(n_total - n_local, o_nnz[row_local]);
            }
        }
    }

    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        n_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        0,
                        n_local ? &d_nnz[0] : nullptr,
                        0,
                        n_local ? &o_nnz[0] : nullptr,
                        &mat);
    IBTK_CHKERRQ(ierr);

    ierr = MatSetBlockSize(mat, 2 * NDIM + 2);
    IBTK_CHKERRQ(ierr);

#if !defined(NDEBUG)
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    for (hier::PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        auto patch = patch_level->getPatch(p());
        const auto& patch_box = patch->getBox();

        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        Pointer<CellData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        Pointer<CellData<NDIM, int> > rho_data = patch->getPatchData(rho_idx);

        Pointer<pdat::CellData<NDIM, double> > chi_data = nullptr;
        if (chi_idx != IBTK::invalid_index) chi_data = patch->getPatchData(chi_idx);

        const hier::IntVector<NDIM> no_ghosts(0);

        pdat::CellData<NDIM, double> grad_p_coefs(patch_box, 3 * NDIM * 2, no_ghosts);
        pdat::CellData<NDIM, double> div_u_coefs(patch_box, 3 * NDIM * 2, no_ghosts);

        auto get_p_bc_coef = [&p_bc_coefs](const int part,
                                           const int /*axis*/) -> SAMRAI::solv::RobinBcCoefStrategy<NDIM>*
        { return p_bc_coefs[part]; };

        auto get_u_bc_coef = [&u_bc_coefs](const int part, const int axis) -> SAMRAI::solv::RobinBcCoefStrategy<NDIM>*
        { return u_bc_coefs[part][axis]; };

        set_cc_centered_derivative_coefs(grad_p_coefs, patch, get_p_bc_coef, data_time, dx);
        set_cc_centered_derivative_coefs(div_u_coefs, patch, get_u_bc_coef, data_time, dx);

        pdat::CellData<NDIM, double> uur_matrix_coefs(patch_box, uu_stencil_sz * NDIM, no_ghosts);
        pdat::CellData<NDIM, double> uui_matrix_coefs(patch_box, uu_stencil_sz * NDIM, no_ghosts);

        IBTK::VCViscousDilatationalOpSpec op_spec;
        op_spec.d_C_is_const = true;
        op_spec.d_C_const = 0.0;
        op_spec.d_C_idx = IBTK::invalid_index;
        op_spec.d_D_is_const = false;
        op_spec.d_D_idx = mu_idx;
        op_spec.d_L_const = false;
        op_spec.d_L_idx = lambda_idx;

        IBTK::PoissonUtilities::computeVCCCViscousDilatationalOpMatrixCoefficients(
            uur_matrix_coefs, patch, stencil_map_vec, &op_spec, u_bc_coefs[REAL], data_time, mu_interp_type);

        IBTK::PoissonUtilities::computeVCCCViscousDilatationalOpMatrixCoefficients(
            uui_matrix_coefs, patch, stencil_map_vec, &op_spec, u_bc_coefs[IMAG], data_time, mu_interp_type);

        for (hier::Box<NDIM>::Iterator b(pdat::CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const pdat::CellIndex<NDIM>& i = b();

            for (int axis = 0; axis < NDIM; ++axis)
            {
                const StencilMapType& sm = stencil_map_vec[axis];
                const int uu_offset = axis * uu_stencil_sz;

                for (int part = 0; part < 2; ++part)
                {
                    const int row = (*u_dof_index_data)(i, u_depth(axis, part));
                    if (!(proc_lower <= row && row < proc_upper)) continue;

                    const int opart = other_part(part);
                    const double freq = (part == REAL ? omega : -omega);

                    auto& uu_coefs = (opart == REAL ? uur_matrix_coefs : uui_matrix_coefs);

                    const int row_stencil_sz = 1 + uu_stencil_sz + up_stencil_sz;
                    std::vector<double> mat_vals(row_stencil_sz);
                    std::vector<int> mat_cols(row_stencil_sz);

                    int idx = 0;
                    mat_vals[idx] = freq * (*rho_data)(i, 0);
                    mat_cols[idx] = row;

                    idx += 1;
                    mat_vals[idx] = uu_coefs(i, uu_offset + sm.at(make_key(ORIGIN, axis))) +
                                    (chi_data.isNull() ? 0.0 : (*chi_data)(i, 0));
                    mat_cols[idx] = (*u_dof_index_data)(i, u_depth(axis, opart));

                    for (int d = 0; d < NDIM; ++d)
                    {
                        idx += 1;
                        mat_vals[idx] = uu_coefs(i, uu_offset + sm.at(make_key(get_shift(d, +1), axis)));
                        mat_cols[idx] = (*u_dof_index_data)(i + get_shift(d, +1), u_depth(axis, opart));

                        idx += 1;
                        mat_vals[idx] = uu_coefs(i, uu_offset + sm.at(make_key(get_shift(d, -1), axis)));
                        mat_cols[idx] = (*u_dof_index_data)(i + get_shift(d, -1), u_depth(axis, opart));
                    }

                    for (int comp = 0; comp < NDIM; ++comp)
                    {
                        if (comp == axis) continue;

                        const auto ea_p = get_shift(axis, +1);
                        const auto ea_m = get_shift(axis, -1);
                        const auto ec_p = get_shift(comp, +1);
                        const auto ec_m = get_shift(comp, -1);

                        const hier::Index<NDIM> shifts[9] = { ea_p + ec_p, ea_m + ec_p, ea_p,   ea_m, ea_p + ec_m,
                                                              ea_m + ec_m, ec_p,        ORIGIN, ec_m };

                        for (int k = 0; k < 9; ++k)
                        {
                            idx += 1;
                            mat_vals[idx] = uu_coefs(i, uu_offset + sm.at(make_key(shifts[k], comp)));
                            mat_cols[idx] = (*u_dof_index_data)(i + shifts[k], u_depth(comp, opart));
                        }
                    }

                    const int c = gd_slot(axis, opart, CENTER);
                    const int pp = gd_slot(axis, opart, PLUS);
                    const int mm = gd_slot(axis, opart, MINUS);

                    idx += 1;
                    mat_vals[idx] = grad_p_coefs(i, c);
                    mat_cols[idx] = (*p_dof_index_data)(i, opart);

                    idx += 1;
                    mat_vals[idx] = grad_p_coefs(i, pp);
                    mat_cols[idx] = (*p_dof_index_data)(i + get_shift(axis, +1), opart);

                    idx += 1;
                    mat_vals[idx] = grad_p_coefs(i, mm);
                    mat_cols[idx] = (*p_dof_index_data)(i + get_shift(axis, -1), opart);

                    TBOX_ASSERT(idx == row_stencil_sz - 1);

                    ierr = MatSetValues(mat, 1, &row, row_stencil_sz, &mat_cols[0], &mat_vals[0], INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }

            for (int part = 0; part < 2; ++part)
            {
                const int row = (*p_dof_index_data)(i, part);
                if (!(proc_lower <= row && row < proc_upper)) continue;

                const int opart = other_part(part);
                const double p_coef =
                    (part == REAL ? omega / (sound_speed * sound_speed) : -omega / (sound_speed * sound_speed));

                const int row_stencil_sz = 1 + pu_stencil_sz;
                std::vector<double> mat_vals(row_stencil_sz);
                std::vector<int> mat_cols(row_stencil_sz);

                int idx = 0;
                mat_vals[idx] = p_coef;
                mat_cols[idx] = row;

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    const auto ip = i + get_shift(axis, +1);
                    const auto im = i + get_shift(axis, -1);

                    const int c = gd_slot(axis, opart, CENTER);
                    const int pp = gd_slot(axis, opart, PLUS);
                    const int mm = gd_slot(axis, opart, MINUS);

                    idx += 1;
                    mat_vals[idx] = (*rho_data)(i, 0) * div_u_coefs(i, c);
                    mat_cols[idx] = (*u_dof_index_data)(i, u_depth(axis, opart));

                    idx += 1;
                    mat_vals[idx] = (*rho_data)(ip, 0) * div_u_coefs(i, pp);
                    mat_cols[idx] = (*u_dof_index_data)(ip, u_depth(axis, opart));

                    idx += 1;
                    mat_vals[idx] = (*rho_data)(im, 0) * div_u_coefs(i, mm);
                    mat_cols[idx] = (*u_dof_index_data)(im, u_depth(axis, opart));
                }

                TBOX_ASSERT(idx == row_stencil_sz - 1);

                ierr = MatSetValues(mat, 1, &row, row_stencil_sz, &mat_cols[0], &mat_vals[0], INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
