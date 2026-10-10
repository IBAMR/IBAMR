// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2026 by the IBAMR developers
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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/solver_utilities.h>

#include <BoxList.h>
#include <GridGeometry.h>
#include <PatchLevel.h>
#include <SideGeometry.h>

namespace IBTK
{
namespace
{
// Offsets by which the periodic images of a box of the level are shifted, including no shift.
std::vector<SAMRAI::hier::IntVector<NDIM>>
periodic_image_offsets(const SAMRAI::hier::PatchLevel<NDIM>& level)
{
    const SAMRAI::hier::IntVector<NDIM> periodic_shift = level.getGridGeometry()->getPeriodicShift(level.getRatio());
    std::vector<SAMRAI::hier::IntVector<NDIM>> offsets(1, SAMRAI::hier::IntVector<NDIM>(0));
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift(d) == 0)
        {
            continue;
        }
        const auto n = offsets.size();
        for (int image = -1; image <= 1; image += 2)
        {
            for (std::size_t k = 0; k < n; ++k)
            {
                SAMRAI::hier::IntVector<NDIM> offset = offsets[k];
                offset(d) = image * periodic_shift(d);
                offsets.push_back(offset);
            }
        }
    }
    return offsets;
}

// Call zero_entry(i, k) for each index i of the box of degrees of freedom of a patch and each entry k of the stencil
// such that the index i + stencil[k] is in the physical domain but is not a degree of freedom of the level. The degrees
// of freedom are the cells if side_axis is negative, and otherwise the sides normal to the axis side_axis.
template <class ZeroEntry>
void
for_each_off_level_entry(const SAMRAI::hier::PatchLevel<NDIM>& level,
                         const SAMRAI::hier::Box<NDIM>& patch_box,
                         const int side_axis,
                         const std::vector<SAMRAI::hier::Index<NDIM>>& stencil,
                         ZeroEntry zero_entry)
{
    using namespace SAMRAI;
    if (level.getLevelNumber() == 0)
    {
        return;
    }
    const auto to_dof_box = [side_axis](const hier::Box<NDIM>& box)
    { return side_axis < 0 ? box : pdat::SideGeometry<NDIM>::toSideBox(box, side_axis); };
    const hier::Box<NDIM> dof_box = to_dof_box(patch_box);
    const hier::Box<NDIM> ghost_box = hier::Box<NDIM>::grow(dof_box, hier::IntVector<NDIM>(1));
    const std::vector<hier::IntVector<NDIM>> offsets = periodic_image_offsets(level);

    // The degrees of freedom in the ghost region of the patch that are in the physical domain, including its periodic
    // images, and that are not degrees of freedom of the level.
    hier::BoxList<NDIM> ghost_dofs(ghost_box);
    ghost_dofs.removeIntersections(dof_box);
    hier::BoxList<NDIM> off_level_dofs;
    for (hier::BoxList<NDIM>::Iterator g(ghost_dofs); g; g++)
    {
        for (int n = 0; n < level.getPhysicalDomain().getNumberOfBoxes(); ++n)
        {
            for (const hier::IntVector<NDIM>& offset : offsets)
            {
                const hier::Box<NDIM> domain_box =
                    hier::Box<NDIM>::shift(to_dof_box(level.getPhysicalDomain()[n]), offset);
                const hier::Box<NDIM> intersection = g() * domain_box;
                if (!intersection.empty())
                {
                    off_level_dofs.appendItem(intersection);
                }
            }
        }
    }
    for (int n = 0; n < level.getBoxes().getNumberOfBoxes(); ++n)
    {
        for (const hier::IntVector<NDIM>& offset : offsets)
        {
            const hier::Box<NDIM> level_box = hier::Box<NDIM>::shift(to_dof_box(level.getBoxes()[n]), offset);
            if (level_box.intersects(ghost_box))
            {
                off_level_dofs.removeIntersections(level_box);
            }
        }
    }

    for (hier::BoxList<NDIM>::Iterator b(off_level_dofs); b; b++)
    {
        for (hier::Box<NDIM>::Iterator g(b()); g; g++)
        {
            for (std::size_t k = 0; k < stencil.size(); ++k)
            {
                hier::Index<NDIM> i = g();
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    i(d) -= stencil[k](d);
                }
                if (dof_box.contains(i))
                {
                    zero_entry(i, static_cast<int>(k));
                }
            }
        }
    }
}
} // namespace

void
reportPETScKSPConvergedReason(const std::string& object_name, const KSPConvergedReason& reason, std::ostream& os)
{
    switch (static_cast<int>(reason))
    {
    case KSP_CONVERGED_RTOL:
        os << object_name
           << ": converged: |Ax-b| <= rtol*|b| --- residual norm is less than specified relative tolerance.\n";
        break;
    case KSP_CONVERGED_ATOL:
        os << object_name
           << ": converged: |Ax-b| <= atol --- residual norm is less than specified absolute tolerance.\n";
        break;
    case KSP_CONVERGED_ITS:
        os << object_name << ": converged: maximum number of iterations reached.\n";
        break;
    case KSP_CONVERGED_STEP_LENGTH:
        os << object_name << ": converged: step size less than specified tolerance.\n";
        break;
    case KSP_DIVERGED_NULL:
        os << object_name << ": diverged: null.\n";
        break;
    case KSP_DIVERGED_ITS:
        os << object_name
           << ": diverged: reached maximum number of iterations before any convergence criteria were satisfied.\n";
        break;
    case KSP_DIVERGED_DTOL:
        os << object_name
           << ": diverged: |Ax-b| >= dtol*|b| --- residual is greater than specified divergence tolerance.\n";
        break;
    case KSP_DIVERGED_BREAKDOWN:
        os << object_name << ": diverged: breakdown in the Krylov method.\n";
        break;
    case KSP_DIVERGED_BREAKDOWN_BICG:
        os << object_name << ": diverged: breakdown in the bi-congugate gradient method.\n";
        break;
    case KSP_DIVERGED_NONSYMMETRIC:
        os << object_name
           << ": diverged: it appears the operator or preconditioner is not symmetric, but this Krylov method (KSPCG, "
              "KSPMINRES, KSPCR) requires symmetry\n";
        break;
    case KSP_DIVERGED_INDEFINITE_PC:
        os << object_name
           << ": diverged: it appears the preconditioner is indefinite (has both positive and negative eigenvalues), "
              "but this Krylov method (KSPCG) requires it to be positive definite.\n";
        break;
    case KSP_CONVERGED_ITERATING:
        os << object_name << ": iterating: KSPSolve() is still running.\n";
        break;
    default:
        os << object_name << ": unknown completion code " << static_cast<int>(reason) << " reported.\n";
        break;
    }
    return;
} // reportPETScKSPConvergedReason

void
reportPETScSNESConvergedReason(const std::string& object_name, const SNESConvergedReason& reason, std::ostream& os)
{
    switch (static_cast<int>(reason))
    {
    case SNES_CONVERGED_FNORM_ABS:
        os << object_name << ": converged: |F| less than specified absolute tolerance.\n";
        break;
    case SNES_CONVERGED_FNORM_RELATIVE:
        os << object_name << ": converged: |F| less than specified relative tolerance.\n";
        break;
    case SNES_CONVERGED_SNORM_RELATIVE:
        os << object_name << ": converged: step size less than specified relative tolerance.\n";
        break;
    case SNES_CONVERGED_ITS:
        os << object_name << ": converged: maximum number of iterations reached.\n";
        break;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
        os << object_name << ": diverged: new x location passed to the function is not in the function domain.\n";
        break;
    case SNES_DIVERGED_FUNCTION_COUNT:
        os << object_name << ": diverged: exceeded maximum number of function evaluations.\n";
        break;
    case SNES_DIVERGED_LINEAR_SOLVE:
        os << object_name << ": diverged: the linear solve failed.\n";
        break;
#if PETSC_VERSION_LT(3, 25, 0)
    case SNES_DIVERGED_FNORM_NAN:
        os << object_name << ": diverged: |F| is NaN.\n";
        break;
#else
    case SNES_DIVERGED_FUNCTION_NANORINF:
        os << object_name << ": diverged: |F| is NaN or infinity.\n";
        break;
#endif
    case SNES_DIVERGED_MAX_IT:
        os << object_name << ": diverged: exceeded maximum number of iterations.\n";
        break;
    case SNES_DIVERGED_LINE_SEARCH:
        os << object_name << ": diverged: line-search failure.\n";
        break;
    case SNES_DIVERGED_INNER:
        os << object_name << ": diverged: inner solve failed.\n";
        break;
    case SNES_DIVERGED_LOCAL_MIN:
        os << object_name << ": diverged: attained non-zero local minimum.\n";
        break;
    case SNES_DIVERGED_TR_DELTA:
        os << object_name << ": diverged: trust-region delta.\n";
        break;
    case SNES_CONVERGED_ITERATING:
        os << object_name << ": iterating.\n";
        break;
    default:
        os << object_name << ": unknown completion code " << static_cast<int>(reason) << " reported.\n";
        break;
    }
} // reportPETScSNESConvergedReason

// hypre defines HYPRE_RELEASE_NUMBER in version 2.21 and newer.
#if HYPRE_RELEASE_NUMBER >= 22100
namespace
{
#if HYPRE_RELEASE_NUMBER < 22900
// These versions of hypre cannot report whether hypre is initialized.
bool s_hypre_initialized = false;
#endif

PetscErrorCode
finalize_hypre()
{
    PetscFunctionBeginUser;
    HYPRE_Finalize();
#if HYPRE_RELEASE_NUMBER < 22900
    s_hypre_initialized = false;
#endif
    PetscFunctionReturn(0);
}
} // namespace
#endif

void
initialize_hypre()
{
#if HYPRE_RELEASE_NUMBER >= 22900
    if (!HYPRE_Initialized())
    {
        HYPRE_Initialize();
        const int ierr = PetscRegisterFinalize(finalize_hypre);
        IBTK_CHKERRQ(ierr);
    }
#elif HYPRE_RELEASE_NUMBER >= 22100
    if (!s_hypre_initialized)
    {
        HYPRE_Init();
        const int ierr = PetscRegisterFinalize(finalize_hypre);
        IBTK_CHKERRQ(ierr);
        s_hypre_initialized = true;
    }
#endif
    return;
} // initialize_hypre

std::array<HYPRE_Int, NDIM>
hypre_array(const SAMRAI::hier::Index<NDIM>& index)
{
    std::array<HYPRE_Int, NDIM> result;
    for (unsigned int d = 0; d < NDIM; ++d) result[d] = index[d];
    return result;
}

void
copyFromHypre(SAMRAI::pdat::CellData<NDIM, double>& dst_data,
              const std::vector<HYPRE_StructVector>& vectors,
              const SAMRAI::hier::Box<NDIM>& box)
{
    const unsigned int depth = dst_data.getDepth();
#ifndef NDEBUG
    TBOX_ASSERT(depth == vectors.size());
#endif
    const SAMRAI::hier::Box<NDIM> transfer_box = box * dst_data.getGhostBox();
    if (transfer_box.empty())
    {
        return;
    }
    std::array<HYPRE_Int, NDIM> lower = hypre_array(transfer_box.lower());
    std::array<HYPRE_Int, NDIM> upper = hypre_array(transfer_box.upper());
    std::array<HYPRE_Int, NDIM> value_lower = hypre_array(dst_data.getGhostBox().lower());
    std::array<HYPRE_Int, NDIM> value_upper = hypre_array(dst_data.getGhostBox().upper());
    for (unsigned int k = 0; k < depth; ++k)
    {
        HYPRE_StructVectorGetBoxValues2(
            vectors[k], lower.data(), upper.data(), value_lower.data(), value_upper.data(), dst_data.getPointer(k));
    }
    return;
} // copyFromHypre

void
copyFromHypre(SAMRAI::pdat::SideData<NDIM, double>& dst_data,
              HYPRE_SStructVector vector,
              const SAMRAI::hier::Box<NDIM>& box)
{
    // SAMRAI gives a side the index of the cell above it, and hypre gives it
    // the index of the cell below it, so subtract one from each index in the
    // direction normal to the side.
    for (int var = 0; var < NDIM; ++var)
    {
        const unsigned int axis = var;
        // Intersect the boxes in side index space, because two boxes of cells
        // that do not overlap can still share sides.
        const SAMRAI::hier::Box<NDIM>& value_box = dst_data.getArrayData(axis).getBox();
        const SAMRAI::hier::Box<NDIM> transfer_box = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(box, axis) * value_box;
        if (transfer_box.empty())
        {
            continue;
        }
        std::array<HYPRE_Int, NDIM> lower = hypre_array(transfer_box.lower());
        std::array<HYPRE_Int, NDIM> upper = hypre_array(transfer_box.upper());
        std::array<HYPRE_Int, NDIM> value_lower = hypre_array(value_box.lower());
        std::array<HYPRE_Int, NDIM> value_upper = hypre_array(value_box.upper());
        --lower[axis];
        --upper[axis];
        --value_lower[axis];
        --value_upper[axis];
        HYPRE_SStructVectorGetBoxValues2(vector,
                                         0,
                                         lower.data(),
                                         upper.data(),
                                         var,
                                         value_lower.data(),
                                         value_upper.data(),
                                         dst_data.getPointer(axis));
    }
    return;
} // copyFromHypre

void
copyToHypre(const std::vector<HYPRE_StructVector>& vectors,
            SAMRAI::pdat::CellData<NDIM, double>& src_data,
            const SAMRAI::hier::Box<NDIM>& box)
{
    const unsigned int depth = src_data.getDepth();
#ifndef NDEBUG
    TBOX_ASSERT(depth == vectors.size());
#endif
    const SAMRAI::hier::Box<NDIM> transfer_box = box * src_data.getGhostBox();
    if (transfer_box.empty())
    {
        return;
    }
    std::array<HYPRE_Int, NDIM> lower = hypre_array(transfer_box.lower());
    std::array<HYPRE_Int, NDIM> upper = hypre_array(transfer_box.upper());
    std::array<HYPRE_Int, NDIM> value_lower = hypre_array(src_data.getGhostBox().lower());
    std::array<HYPRE_Int, NDIM> value_upper = hypre_array(src_data.getGhostBox().upper());
    for (unsigned int k = 0; k < depth; ++k)
    {
        HYPRE_StructVectorSetBoxValues2(
            vectors[k], lower.data(), upper.data(), value_lower.data(), value_upper.data(), src_data.getPointer(k));
    }
    return;
} // copyToHypre

void
copyToHypre(HYPRE_SStructVector& vector,
            SAMRAI::pdat::SideData<NDIM, double>& src_data,
            const SAMRAI::hier::Box<NDIM>& box)
{
    // SAMRAI gives a side the index of the cell above it, and hypre gives it
    // the index of the cell below it, so subtract one from each index in the
    // direction normal to the side.
    for (int var = 0; var < NDIM; ++var)
    {
        const unsigned int axis = var;
        // Intersect the boxes in side index space, because two boxes of cells
        // that do not overlap can still share sides.
        const SAMRAI::hier::Box<NDIM>& value_box = src_data.getArrayData(axis).getBox();
        const SAMRAI::hier::Box<NDIM> transfer_box = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(box, axis) * value_box;
        if (transfer_box.empty())
        {
            continue;
        }
        std::array<HYPRE_Int, NDIM> lower = hypre_array(transfer_box.lower());
        std::array<HYPRE_Int, NDIM> upper = hypre_array(transfer_box.upper());
        std::array<HYPRE_Int, NDIM> value_lower = hypre_array(value_box.lower());
        std::array<HYPRE_Int, NDIM> value_upper = hypre_array(value_box.upper());
        --lower[axis];
        --upper[axis];
        --value_lower[axis];
        --value_upper[axis];
        HYPRE_SStructVectorSetBoxValues2(vector,
                                         0,
                                         lower.data(),
                                         upper.data(),
                                         var,
                                         value_lower.data(),
                                         value_upper.data(),
                                         src_data.getPointer(axis));
    }
    return;
} // copyToHypre

void
clearOffLevelMatrixEntries(SAMRAI::pdat::CellData<NDIM, double>& matrix_coefficients,
                           const SAMRAI::hier::PatchLevel<NDIM>& level,
                           const std::vector<SAMRAI::hier::Index<NDIM>>& stencil)
{
    for_each_off_level_entry(level,
                             matrix_coefficients.getBox(),
                             /*side_axis*/ -1,
                             stencil,
                             [&matrix_coefficients](const SAMRAI::hier::Index<NDIM>& i, const int k)
                             { matrix_coefficients(SAMRAI::pdat::CellIndex<NDIM>(i), k) = 0.0; });
}

void
clearOffLevelMatrixEntries(SAMRAI::pdat::SideData<NDIM, double>& matrix_coefficients,
                           const SAMRAI::hier::PatchLevel<NDIM>& level,
                           const std::vector<SAMRAI::hier::Index<NDIM>>& stencil)
{
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for_each_off_level_entry(
            level,
            matrix_coefficients.getBox(),
            static_cast<int>(axis),
            stencil,
            [&matrix_coefficients, axis](const SAMRAI::hier::Index<NDIM>& i, const int k) {
                matrix_coefficients(SAMRAI::pdat::SideIndex<NDIM>(i, axis, SAMRAI::pdat::SideIndex<NDIM>::Lower), k) =
                    0.0;
            });
    }
}
} // namespace IBTK
