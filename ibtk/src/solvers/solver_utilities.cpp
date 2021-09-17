// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
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

#include "ibtk/solver_utilities.h"

#include <memory>

namespace IBTK
{
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
#if PETSC_VERSION_LT(3, 12, 0)
    case SNES_CONVERGED_TR_DELTA:
        os << object_name << ": converged: trust-region delta.\n";
        break;
#endif
    case SNES_DIVERGED_FUNCTION_DOMAIN:
        os << object_name << ": diverged: new x location passed to the function is not in the function domain.\n";
        break;
    case SNES_DIVERGED_FUNCTION_COUNT:
        os << object_name << ": diverged: exceeded maximum number of function evaluations.\n";
        break;
    case SNES_DIVERGED_LINEAR_SOLVE:
        os << object_name << ": diverged: the linear solve failed.\n";
        break;
    case SNES_DIVERGED_FNORM_NAN:
        os << object_name << ": diverged: |F| is NaN.\n";
        break;
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
#if PETSC_VERSION_GE(3, 12, 0)
    case SNES_DIVERGED_TR_DELTA:
        os << object_name << ": diverged: trust-region delta.\n";
        break;
#endif
    case SNES_CONVERGED_ITERATING:
        os << object_name << ": iterating.\n";
        break;
    default:
        os << object_name << ": unknown completion code " << static_cast<int>(reason) << " reported.\n";
        break;
    }
} // reportPETScSNESConvergedReason

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
    const bool copy_data = dst_data.getGhostBox() != box;
    std::unique_ptr<SAMRAI::pdat::CellData<NDIM, double> > dst_data_box(
        copy_data ? new SAMRAI::pdat::CellData<NDIM, double>(box, dst_data.getDepth(), 0) : nullptr);
    SAMRAI::pdat::CellData<NDIM, double>& hypre_data = copy_data ? *dst_data_box : dst_data;
    unsigned int depth = dst_data.getDepth();
#ifndef NDEBUG
    TBOX_ASSERT(depth == vectors.size());
#endif
    auto lower = hypre_array(box.lower());
    auto upper = hypre_array(box.upper());
    for (unsigned int k = 0; k < depth; ++k)
    {
        HYPRE_StructVectorGetBoxValues(vectors[k], lower.data(), upper.data(), hypre_data.getPointer(k));
    }
    if (copy_data)
    {
        dst_data.copyOnBox(hypre_data, box);
    }
    return;
} // copyFromHypre

void
copyFromHypre(SAMRAI::pdat::SideData<NDIM, double>& dst_data,
              HYPRE_SStructVector vector,
              const SAMRAI::hier::Box<NDIM>& box)
{
    const bool copy_data = dst_data.getGhostBox() != box;
    std::unique_ptr<SAMRAI::pdat::SideData<NDIM, double> > dst_data_box(
        copy_data ? new SAMRAI::pdat::SideData<NDIM, double>(box, 1, 0) : nullptr);
    SAMRAI::pdat::SideData<NDIM, double>& hypre_data = copy_data ? *dst_data_box : dst_data;
    for (int var = 0; var < NDIM; ++var)
    {
        const unsigned int axis = var;
        auto lower = hypre_array(box.lower());
        lower[axis] -= 1;
        auto upper = hypre_array(box.upper());
        HYPRE_SStructVectorGetBoxValues(vector, 0, lower.data(), upper.data(), var, hypre_data.getPointer(axis));
    }
    if (copy_data)
    {
        dst_data.copyOnBox(hypre_data, box);
    }
    return;
} // copyFromHypre

void
copyToHypre(const std::vector<HYPRE_StructVector>& vectors,
            SAMRAI::pdat::CellData<NDIM, double>& src_data,
            const SAMRAI::hier::Box<NDIM>& box)
{
    const bool copy_data = src_data.getGhostBox() != box;
    std::unique_ptr<SAMRAI::pdat::CellData<NDIM, double> > src_data_box(
        copy_data ? new SAMRAI::pdat::CellData<NDIM, double>(box, src_data.getDepth(), 0) : nullptr);
    SAMRAI::pdat::CellData<NDIM, double>& hypre_data = copy_data ? *src_data_box : src_data;
    if (copy_data) hypre_data.copyOnBox(src_data, box);
    unsigned int depth = src_data.getDepth();
#ifndef NDEBUG
    TBOX_ASSERT(depth == vectors.size());
#endif
    auto lower = hypre_array(box.lower());
    auto upper = hypre_array(box.upper());
    for (unsigned int k = 0; k < depth; ++k)
    {
        HYPRE_StructVectorSetBoxValues(vectors[k], lower.data(), upper.data(), hypre_data.getPointer(k));
    }
    return;
} // copyToHypre

void
copyToHypre(HYPRE_SStructVector& vector,
            SAMRAI::pdat::SideData<NDIM, double>& src_data,
            const SAMRAI::hier::Box<NDIM>& box)
{
    const bool copy_data = src_data.getGhostBox() != box;
    std::unique_ptr<SAMRAI::pdat::SideData<NDIM, double> > src_data_box(
        copy_data ? new SAMRAI::pdat::SideData<NDIM, double>(box, 1, 0) : nullptr);
    SAMRAI::pdat::SideData<NDIM, double>& hypre_data = copy_data ? *src_data_box : src_data;
    if (copy_data) hypre_data.copyOnBox(src_data, box);
    for (int var = 0; var < NDIM; ++var)
    {
        const unsigned int axis = var;
        auto lower = hypre_array(box.lower());
        lower[axis] -= 1;
        auto upper = hypre_array(box.upper());
        HYPRE_SStructVectorSetBoxValues(vector, 0, lower.data(), upper.data(), var, hypre_data.getPointer(axis));
    }
    return;
} // copyToHypre
} // namespace IBTK
